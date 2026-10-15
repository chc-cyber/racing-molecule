"""Conservation, periodicity and full-molecule normal tests, without LAMMPS."""

import json
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from molecular_car.periodic_au import (
    CACHE_Z_MAX, CACHE_Z_MIN, EV_TO_KCAL_MOL, ROOT, PeriodicAuPotential,
    _normal_geometry, create_periodic_surface,
)


def benzene(angle=0., height=3.4):
    phi = np.arange(6)*np.pi/3
    xy = np.stack([np.cos(phi), np.sin(phi)], axis=1)
    values = np.column_stack([np.concatenate([1.4*xy, 2.48*xy]), np.zeros(12)])
    radians = np.deg2rad(angle)
    rotation = np.array([[1, 0, 0], [0, np.cos(radians), -np.sin(radians)],
                         [0, np.sin(radians), np.cos(radians)]])
    values = values@rotation.T
    values += [.379, .713, height]
    return values


class PeriodicAuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.elements = ["C"]*6+["H"]*6
        cls.direct = PeriodicAuPotential(cls.elements, benzene(), use_cache=False, accelerator="numpy")
        # A modest table suffices to test the interpolator's own derivative;
        # production accuracy is checked against the native oracle separately.
        directory = ROOT/"outputs/periodic-au-check/unit-cache"
        directory.mkdir(parents=True, exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(prefix="molecular-periodic-au-", dir=directory)
        cls.cache = PeriodicAuPotential(cls.elements, benzene(), use_cache=True,
            cache_dir=cls.temporary.name, grid_nx=8, grid_ny=12, grid_dz=.25, accelerator="numpy")

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def assert_force_gradient(self, surface, values, axes=None, step=1e-5, atol=3e-7):
        _, force = surface.energy_and_forces_eV(values)
        gradient = np.empty_like(force)
        selected = [(i, j) for i in range(len(values)) for j in range(3)] if axes is None else axes
        for atom, axis in selected:
            plus, minus = values.copy(), values.copy()
            plus[atom, axis] += step
            minus[atom, axis] -= step
            gradient[atom, axis] = (surface.energy_and_forces_eV(plus)[0]-surface.energy_and_forces_eV(minus)[0])/(2*step)
            self.assertAlmostEqual(force[atom, axis], -gradient[atom, axis], delta=atol,
                msg=f"atom {atom}, axis {axis}")

    def test_direct_analytic_force_for_tilt_and_local_twist(self):
        for angle in (0, 47, 85):
            coordinates = benzene(angle, 5.1)
            coordinates[1, 2] += .12
            coordinates[9, 2] -= .07
            self.assert_force_gradient(self.direct, coordinates)

    def test_cached_analytic_force_including_normal_chain(self):
        coordinates = benzene(42, 5.2)
        coordinates[1, 2] += .12
        coordinates[8, 1] += .06
        self.assert_force_gradient(self.cache, coordinates)

    def test_conservative_fallback_and_blending_edges(self):
        for height in (1.9, CACHE_Z_MIN, 2.15, 2.3, 11.7, 11.85, CACHE_Z_MAX, 12.1):
            coordinates = benzene(height=height)
            # Test the molecular translation derivative so the normal-chain
            # contributions cancel, then individual atom derivatives as well.
            _, force = self.cache.energy_and_forces_eV(coordinates)
            step = 1e-5
            shift = np.array([0., 0., step])
            fd = (self.cache.energy_and_forces_eV(coordinates+shift)[0]
                 -self.cache.energy_and_forces_eV(coordinates-shift)[0])/(2*step)
            self.assertAlmostEqual(np.sum(force[:, 2]), -fd, delta=2e-6)
            self.assert_force_gradient(self.cache, coordinates, axes=[(0, 2), (7, 2)], atol=2e-6)
        for height in (1.9, 12.1, 17.5):
            direct = self.direct.energy_and_forces_eV(benzene(height=height))
            cached = self.cache.energy_and_forces_eV(benzene(height=height))
            self.assertAlmostEqual(direct[0], cached[0], delta=1e-12)
            np.testing.assert_allclose(direct[1], cached[1], atol=1e-12, rtol=0)

    def test_periodic_translation_across_many_cells(self):
        coordinates = benzene(33, 5)
        for surface in (self.direct, self.cache):
            initial = surface.energy_and_forces_eV(coordinates)
            translated = coordinates+[997*surface.period[0], -613*surface.period[1], 0]
            result = surface.energy_and_forces_eV(translated)
            self.assertAlmostEqual(result[0], initial[0], delta=1e-11)
            np.testing.assert_allclose(result[1], initial[1], atol=3e-11, rtol=0)

    def test_query_cell_seams_have_consistent_gradient(self):
        coordinates = benzene(18, 5)
        for surface in (self.direct, self.cache):
            coordinates[:, 0] -= coordinates[0, 0]
            coordinates[:, 1] -= coordinates[0, 1]
            self.assert_force_gradient(surface, coordinates, axes=[(0, 0), (0, 1)])

    def test_car_is_never_folded_for_normal_neighbours(self):
        first = benzene(21, 4.8)
        second = first+[5*self.direct.period[0], 4*self.direct.period[1], 0]
        whole = np.concatenate([first, second])
        surface = PeriodicAuPotential(self.elements*2, whole, use_cache=False)
        energy, force = surface.energy_and_forces_eV(whole)
        first_energy, first_force = self.direct.energy_and_forces_eV(first)
        self.assertAlmostEqual(energy, 2*first_energy, delta=1e-12)
        np.testing.assert_allclose(force[:12], first_force, atol=1e-12, rtol=0)
        np.testing.assert_allclose(force[12:], first_force, atol=1e-12, rtol=0)
        normals = _normal_geometry(whole, surface.normal_cutoffs)
        self.assertTrue(np.all(normals[2][:6] == 3))
        self.assertTrue(np.all(normals[2][6:12] == 1))

    def test_nontrivial_hydrogen_normals_are_not_fixed(self):
        # Two nearby C neighbours give H a real local plane. Carbon also has
        # two neighbours; this covers the less common two-neighbour chain.
        elements = ["C", "C", "H"]
        coordinates = np.array([[.17, .34, 3.4], [1.57, .34, 3.513], [.86, 1.24, 3.82]])
        direct = PeriodicAuPotential(elements, coordinates, use_cache=False)
        cached = PeriodicAuPotential(elements, coordinates, use_cache=True,
            cache_dir=self.temporary.name, grid_nx=8, grid_ny=12, grid_dz=.25)
        self.assertTrue(np.all(_normal_geometry(coordinates, direct.normal_cutoffs)[2] == 2))
        self.assert_force_gradient(direct, coordinates)
        self.assert_force_gradient(cached, coordinates)

    def test_cache_nodes_match_exact_sampled_potential(self):
        # Carbon's angular term remains exact, and hydrogen is exactly fixed-z
        # here. Catmull-Rom passes through each stored scalar potential node.
        coordinates = benzene(height=5)
        coordinates[:, 0] = np.round(coordinates[:, 0]/(self.cache.period[0]/8))*(self.cache.period[0]/8)
        coordinates[:, 1] = np.round(coordinates[:, 1]/(self.cache.period[1]/12))*(self.cache.period[1]/12)
        direct_energy = self.direct.energy_and_forces_eV(coordinates)[0]
        cached_energy = self.cache.energy_and_forces_eV(coordinates)[0]
        self.assertAlmostEqual(direct_energy, cached_energy, delta=5e-12)

    def test_cache_reload_and_metadata_are_reviewable(self):
        again = PeriodicAuPotential(self.elements, benzene(), use_cache=True,
            cache_dir=self.temporary.name, grid_nx=8, grid_ny=12, grid_dz=.25)
        self.assertTrue(again.cache_loaded)
        self.assertEqual(again.cache_build_seconds, 0)
        json.dumps(again.metadata(), allow_nan=False)
        np.testing.assert_array_equal(again._grid, self.cache._grid)
        self.assertTrue(Path(again.metadata()["cache_file"]).is_file())

    def test_factory_units_and_early_validation(self):
        atoms = [{"element":element} for element in self.elements]
        surface = create_periodic_surface(atoms, benzene(), use_cache=False)
        energy, force = surface.energy_and_forces_eV(benzene())
        kcal_energy, kcal_force = surface.energy_and_forces(benzene())
        self.assertEqual(kcal_energy, energy*EV_TO_KCAL_MOL)
        np.testing.assert_array_equal(kcal_force, force*EV_TO_KCAL_MOL)
        with self.assertRaisesRegex(ValueError, "C/H"):
            create_periodic_surface(["O"], np.zeros((1, 3)))
        with self.assertRaisesRegex(ValueError, "three normal"):
            create_periodic_surface(["C"]+["H"]*4,
                np.array([[0, 0, 4], [1, 0, 4], [-1, 0, 4], [0, 1, 4], [0, -1, 4]]), use_cache=False)
        with self.assertRaisesRegex(ValueError, "zero normal"):
            create_periodic_surface(["C"]*3, np.array([[0, 0, 4], [.8, 0, 4], [1.6, 0, 4]]), use_cache=False)

    @unittest.skipUnless(importlib.util.find_spec("numba") is not None, "Optional Numba not installed")
    def test_compiled_accelerator_matches_numpy_complete_forces(self):
        for cached in (False, True):
            options = dict(use_cache=cached, cache_dir=self.temporary.name,
                           grid_nx=8, grid_ny=12, grid_dz=.25)
            numpy_surface = PeriodicAuPotential(self.elements, benzene(), accelerator="numpy", **options)
            compiled = PeriodicAuPotential(self.elements, benzene(), accelerator="numba", **options)
            self.assertEqual(compiled.metadata()["accelerator"], "numba")
            self.assertIsNotNone(compiled.metadata()["numba_version"])
            for angle in (0, 47, 85):
                coordinates = benzene(angle, 5.1)
                coordinates[1, 2] += .12
                coordinates[9, 2] -= .07
                self.assert_surfaces_equal(numpy_surface, compiled, coordinates)
                self.assert_force_gradient(compiled, coordinates)
            for height in (1.9, 2., 2.15, 2.3, 11.7, 11.85, 12., 12.1):
                self.assert_surfaces_equal(numpy_surface, compiled, benzene(height=height))
            values = benzene(33, 5)+[997*compiled.period[0], -613*compiled.period[1], 0]
            self.assert_surfaces_equal(numpy_surface, compiled, values)
            # The contiguous input ABI is preserved for a strided view.
            signatures = len(compiled._compiled_kernel.signatures)
            self.assert_surfaces_equal(numpy_surface, compiled, np.asfortranarray(values))
            self.assertEqual(len(compiled._compiled_kernel.signatures), signatures)
        elements = ["C", "C", "H"]
        coordinates = np.array([[.17, .34, 3.4], [1.57, .34, 3.513], [.86, 1.24, 3.82]])
        numpy_surface = PeriodicAuPotential(elements, coordinates, accelerator="numpy", **options)
        compiled = PeriodicAuPotential(elements, coordinates, accelerator="numba", **options)
        self.assert_surfaces_equal(numpy_surface, compiled, coordinates)
        self.assert_force_gradient(compiled, coordinates)
        first = benzene(21, 4.8)
        whole = np.concatenate([first, first+[5*self.direct.period[0], 4*self.direct.period[1], 0]])
        numpy_surface = PeriodicAuPotential(self.elements*2, whole, accelerator="numpy", **options)
        compiled = PeriodicAuPotential(self.elements*2, whole, accelerator="numba", **options)
        self.assert_surfaces_equal(numpy_surface, compiled, whole)

    def assert_surfaces_equal(self, numpy_surface, compiled, coordinates):
        expected = numpy_surface.energy_and_forces_eV(coordinates)
        actual = compiled.energy_and_forces_eV(coordinates)
        self.assertAlmostEqual(actual[0], expected[0], delta=3e-11)
        np.testing.assert_allclose(actual[1], expected[1], atol=3e-11, rtol=0)

    def test_missing_optional_accelerator_is_explicit(self):
        with mock.patch("molecular_car.periodic_au.importlib.util.find_spec", return_value=None):
            automatic = PeriodicAuPotential(self.elements, benzene(), use_cache=False, accelerator="auto")
            self.assertEqual(automatic.metadata()["accelerator"], "numpy")
            self.assert_surfaces_equal(self.direct, automatic, benzene())
            with self.assertRaisesRegex(ImportError, "Numba"):
                PeriodicAuPotential(self.elements, benzene(), use_cache=False, accelerator="numba")


if __name__ == "__main__":
    unittest.main()
