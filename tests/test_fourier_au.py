"""Independent conservative-proxy tests using deliberately coarse small tables."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from molecular_car.fourier_au import FourierAuPotential, _box_weight
from molecular_car.periodic_au import PeriodicAuPotential, ROOT, _normal_geometry


def benzene(angle=0.,height=3.5):
    phi=np.arange(6)*np.pi/3
    xy=np.stack([np.cos(phi),np.sin(phi)],axis=1)
    values=np.column_stack([np.concatenate([1.4*xy,2.48*xy]),np.zeros(12)])
    a=np.deg2rad(angle)
    values=values@np.array([[1,0,0],[0,np.cos(a),-np.sin(a)],[0,np.sin(a),np.cos(a)]]).T
    return values+[.379,.713,height]


class FourierAuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.elements=['C']*6+['H']*6
        directory=ROOT/'outputs/fourier-au-check/unit-cache'
        directory.mkdir(parents=True,exist_ok=True)
        cls.temporary=tempfile.TemporaryDirectory(prefix='fourier-conservative-',dir=directory)
        cls.options=dict(angular_step=.5,chart_bound=2.,chart_blend=.4,
            height_min=2.8,height_max=4.6,height_blend=.3,height_step=.15,
            grid_nx=8,grid_ny=12,grid_dz=.25,cache_dir=cls.temporary.name,
            fourier_cache_dir=cls.temporary.name)
        cls.proxy=FourierAuPotential(cls.elements,benzene(),accelerator='numpy',**cls.options)
        cls.compiled=FourierAuPotential(cls.elements,benzene(),accelerator='numba',**cls.options)
        cls.baseline=PeriodicAuPotential(cls.elements,benzene(),accelerator='numpy',
            cache_dir=cls.temporary.name,grid_nx=8,grid_ny=12,grid_dz=.25)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def assert_fd(self,surface,positions,axes=None,atol=2e-6):
        energy,force=surface.energy_and_forces_eV(positions)
        axes=axes or [(i,j) for i in range(len(positions)) for j in range(3)]
        for atom,axis in axes:
            plus,minus=positions.copy(),positions.copy()
            plus[atom,axis]+=1e-5;minus[atom,axis]-=1e-5
            fd=(surface.energy_and_forces_eV(plus)[0]-surface.energy_and_forces_eV(minus)[0])/2e-5
            self.assertAlmostEqual(force[atom,axis],-fd,delta=atol,msg=f'{atom}/{axis}')

    def test_entire_force_is_its_potential_derivative(self):
        values=benzene(35,4.)
        values[1,2]+=.12;values[9,1]-=.04
        self.assert_fd(self.proxy,values)
        self.assert_fd(self.compiled,values)

    def test_conservative_height_blend_and_boundaries(self):
        for z in (2.79,2.8,2.95,3.1,4.3,4.45,4.6,4.61):
            self.assert_fd(self.compiled,benzene(height=z),[(0,2),(7,2)])

    def test_conservative_chart_blend_and_equator(self):
        for ratio in (1.59,1.6,1.8,1.99,2.,2.01):
            # Move the whole tilted molecule high enough to avoid hard contact.
            values=benzene(np.rad2deg(np.arctan(ratio)),5.)
            self.assert_fd(self.compiled,values,[(0,0),(0,1),(0,2),(1,2)])
        self.assert_fd(self.compiled,benzene(90,5.),[(0,1),(0,2)])

    def test_outside_support_matches_validated_cache(self):
        for coordinates in (benzene(85,6.),benzene(height=8.),benzene(height=2.7)):
            a=self.compiled.energy_and_forces_eV(coordinates)
            b=self.baseline.energy_and_forces_eV(coordinates)
            self.assertAlmostEqual(a[0],b[0],delta=1e-12)
            np.testing.assert_allclose(a[1],b[1],atol=1e-11,rtol=0)
            self.assertEqual(self.compiled.last_evaluation_counts['exact_direction_fallback'],6)

    def test_python_and_compiled_agree(self):
        for angle in (0,25,61,85):
            coordinates=benzene(angle,4.)
            coordinates[1,2]+=.06
            a=self.proxy.energy_and_forces_eV(coordinates)
            b=self.compiled.energy_and_forces_eV(coordinates)
            self.assertAlmostEqual(a[0],b[0],delta=2e-12)
            np.testing.assert_allclose(a[1],b[1],atol=2e-11,rtol=0)

    def test_multicell_unwrapped_translation_and_seam(self):
        coordinates=benzene(20,4.)
        coordinates[:,:2]-=coordinates[0,:2]
        self.assert_fd(self.compiled,coordinates,[(0,0),(0,1)])
        a=self.compiled.energy_and_forces_eV(coordinates)
        b=self.compiled.energy_and_forces_eV(coordinates+[997*self.proxy.period[0],-613*self.proxy.period[1],0])
        self.assertAlmostEqual(a[0],b[0],delta=2e-11)
        np.testing.assert_allclose(a[1],b[1],atol=2e-10,rtol=0)

    def test_opposite_normal_chart_has_identical_potential(self):
        point=np.array([.8,.6,3.7]);n=np.array([.3,.8,1.]);n/=np.linalg.norm(n)
        a,b=n[:2]/n[2]
        one=self.proxy._proxy(0,point,n,a,b);two=self.proxy._proxy(0,point,-n,a,b)
        self.assertEqual(one[0],two[0]);np.testing.assert_array_equal(one[1],two[1])
        np.testing.assert_array_equal(one[2],-two[2])

    def test_metadata_cache_identity_and_hydrogen_direction(self):
        self.assertTrue(self.compiled.fourier_cache_loaded)
        self.assertEqual(self.compiled.fourier_cache_path,self.proxy.fourier_cache_path)
        metadata=self.compiled.metadata()
        self.assertEqual(metadata['accelerator'],'numba')
        self.assertEqual(metadata['proxy_height_support_angstrom'],[2.8,4.6])
        self.assertEqual(metadata['fourier_mode_count'],12)
        # H with two C neighbours exercises a real normal rather than fixed z.
        points=np.array([[0.,0.,3.65],[1.1,0.,3.7],[0.,1.1,3.8]])
        hydrogen=FourierAuPotential(['H','C','C'],points,accelerator='numba',**self.options)
        self.assert_fd(hydrogen,points)

    def test_invalid_and_degenerate_input_rejected(self):
        with self.assertRaises(ValueError):
            self.compiled.energy_and_forces_eV(np.zeros((12,2)))
        with self.assertRaises(ValueError):
            FourierAuPotential(['O'],accelerator='numpy',**self.options)
        values=benzene();values[1]=values[0]+[.8,0,0];values[5]=values[0]+[1.4,0,0]
        with self.assertRaises(ValueError):
            self.compiled.energy_and_forces_eV(values)


if __name__=='__main__':
    unittest.main()
