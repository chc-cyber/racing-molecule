import unittest
import numpy as np
from molecular_car import fields
from molecular_car.surface_racing import (gradient_energy_forces,continuous_surface_energy_forces,
    paired_velocities,decode_images,KB_EV_K,MVV_TO_EV)


class SurfaceRacePhysicsTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[.3,.8,3.1],[1.7,-.2,3.3],[-.4,1.1,2.7]])
        self.q=np.array([.07,-.1,.03])
        self.a=np.array([1,.2,.9])

    def test_driver_matches_original_and_has_no_y_force(self):
        params=({"x":0,"y":0,"z":7,"strength":0},.08)
        e,f=fields.energy_and_forces(self.p,self.q,self.a,*params,drive={"mode":"gradient","strength":.75})
        e0,f0=fields.energy_and_forces(self.p,self.q,self.a,*params,drive={"mode":"gradient","strength":0})
        actual,force=gradient_energy_forces(self.p,self.q,self.a,.75,0)
        self.assertAlmostEqual(actual,e-e0,places=11)
        np.testing.assert_allclose(force,f-f0,atol=1e-12)
        np.testing.assert_array_equal(force[:,1],0)

    def test_surface_removes_only_original_boundaries(self):
        e,f=fields.energy_and_forces(self.p,self.q,self.a,{},.08,drive={"mode":"gradient","strength":0})
        for axis,lower,upper in ((0,fields.TRACK_X_MIN,fields.TRACK_X_MAX),(1,-fields.LANE_HALF_WIDTH,fields.LANE_HALF_WIDTH)):
            eb,fb=fields._soft_boundary(self.p[:,axis],lower,upper)
            e-=eb.sum();f[:,axis]-=fb
        actual,force=continuous_surface_energy_forces(self.p)
        self.assertAlmostEqual(actual,e,places=12)
        np.testing.assert_allclose(force,f,atol=1e-12)

    def test_origin_translation_preserves_field_and_gradient(self):
        e,f=gradient_energy_forces(self.p,self.q,self.a,.25,0)
        moved,force=gradient_energy_forces(self.p+[500,0,0],self.q,self.a,.25,500)
        self.assertAlmostEqual(e,moved,places=11)
        np.testing.assert_allclose(f,force,atol=1e-13)
        delta=np.array([[.11,.07,-.02],[.06,-.02,.01],[-.08,.09,.04]])
        eps=1e-5
        ep=gradient_energy_forces(self.p+eps*delta,self.q,self.a,.25,0)[0]
        em=gradient_energy_forces(self.p-eps*delta,self.q,self.a,.25,0)[0]
        self.assertAlmostEqual((ep-em)/(2*eps),-np.sum(f*delta),places=7)

    def test_matched_velocity_kinetic_units_and_zero_COM(self):
        masses=np.array([12.011]*6+[1.008]*6)
        v=paired_velocities(masses,600,1234567)
        np.testing.assert_allclose(np.average(v,axis=0,weights=masses),0,atol=1e-14)
        np.testing.assert_array_equal(v,paired_velocities(masses,600,1234567))
        actual=.5*MVV_TO_EV*np.sum(masses[:,None]*v**2)
        self.assertAlmostEqual(actual,.5*(3*12-3)*KB_EV_K*600,places=12)

    def test_image_unwrap_counts_at_seams(self):
        counts=np.array([[0,0,0],[3,-2,0],[-1,1,0]])
        packed=((counts[:,0]+512)|((counts[:,1]+512)<<10)|((counts[:,2]+512)<<20))
        np.testing.assert_array_equal(decode_images(packed,1023,512,10,20),counts)


if __name__=="__main__":unittest.main()
