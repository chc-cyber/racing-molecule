"""Energy-gradient and phase-average checks of the two game controls."""
import unittest
import numpy as np
from molecular_car.fields import energy_and_forces, gradient_coefficients
from molecular_car.external_fields import ExternalFieldEvaluator
from molecular_car.lattice_corrugation import LatticeCorrugation


class AnalyticPeriodicProbe:
    period=np.array([2.,3.])
    def energy_and_forces_eV(self,r,center_weights=None,center_dz=None):
        x,y,z=r.T;kx,ky=2*np.pi/self.period
        a=np.exp(-z);c=np.cos(kx*x)*np.cos(ky*y);h=.03*np.cos(2*kx*x+6*ky*y)
        e=(z-3)**2+a*c+h
        f=np.array([a*kx*np.sin(kx*x)*np.cos(ky*y)+.03*2*kx*np.sin(2*kx*x+6*ky*y),
                    a*ky*np.cos(kx*x)*np.sin(ky*y)+.03*6*ky*np.sin(2*kx*x+6*ky*y),
                    -2*(z-3)+a*c]).T
        if center_weights is not None:
            f*=center_weights[:,None];f[:,2]-=center_dz*e;e*=center_weights
        return float(e.sum()),f
    def metadata(self):return {}


class CrystalContrastTests(unittest.TestCase):
    def test_phase_average_binding_and_full_gradient(self):
        base=AnalyticPeriodicProbe();r=np.array([[.31,-.21,3.3],[11.2,7.4,3.1]])
        for scale in (0.,1.,3.,12.):
            potential=LatticeCorrugation(base,scale);e,f=potential.energy_and_forces_eV(r)
            for atom,axis in np.ndindex(r.shape):
                delta=np.zeros_like(r);delta[atom,axis]=1e-5
                numeric=-(potential.energy_and_forces_eV(r+delta)[0]-potential.energy_and_forces_eV(r-delta)[0])/2e-5
                self.assertAlmostEqual(numeric,f[atom,axis],places=6)
            cell=[]
            for i in range(12):
                for j in range(18):cell.append(potential.energy_and_forces_eV(r+[i*2/12,j*3/18,0])[0])
            self.assertAlmostEqual(np.mean(cell),np.sum((r[:,2]-3)**2),places=10)
        original=base.energy_and_forces_eV(r)
        np.testing.assert_array_equal(LatticeCorrugation(base,1.).energy_and_forces_eV(r)[1],original[1])

    def test_field_growth_same_initial_push_and_cartesian_gradient(self):
        r=np.array([[.31,-.21,3.3],[201.2,7.4,3.1]])
        q=np.array([.3,-.3]);alpha=np.array([1.,.85]);tip={'x':6,'y':0,'z':7,'strength':.6}
        compiled=ExternalFieldEvaluator(q,alpha)
        for growth in (.1,.25,1.,2.):
            b,g=gradient_coefficients(.6,growth)
            self.assertAlmostEqual(b*g,.5*.0015*.6**2,places=15)
            drive={'mode':'gradient','strength':.6,'growth':growth}
            e,f=energy_and_forces(r,q,alpha,tip,.08,drive)
            ce,cf=compiled.energy_and_forces(r,1,.6,6,0,7,.08,growth)
            np.testing.assert_allclose(ce,e,atol=1e-10);np.testing.assert_allclose(cf,f,atol=1e-10)
            for atom,axis in np.ndindex(r.shape):
                d=np.zeros_like(r);d[atom,axis]=1e-5
                numeric=-(energy_and_forces(r+d,q,alpha,tip,.08,drive)[0]-energy_and_forces(r-d,q,alpha,tip,.08,drive)[0])/2e-5
                self.assertAlmostEqual(numeric,f[atom,axis],places=6)

    def test_invalid_controls_rejected(self):
        for x in (True,float('nan'),-.1,13.):
            with self.assertRaises(ValueError):LatticeCorrugation(AnalyticPeriodicProbe(),x)
        for x in (True,float('nan'),0.,3.):
            with self.assertRaises(ValueError):gradient_coefficients(.6,x)

    def test_repulsive_core_gate_gradient(self):
        base=AnalyticPeriodicProbe()
        potential=LatticeCorrugation(base,12.)
        for height in (1.99,2.01,2.21,2.41,2.59,2.61):
            r=np.array([[.31,-.21,height],[1.2,.4,2.3]])
            e,f=potential.energy_and_forces_eV(r)
            for atom,axis in np.ndindex(r.shape):
                d=np.zeros_like(r);d[atom,axis]=1e-6
                numeric=-(potential.energy_and_forces_eV(r+d)[0]-potential.energy_and_forces_eV(r-d)[0])/2e-6
                self.assertAlmostEqual(numeric,f[atom,axis],places=6)
            if height<=2.:
                np.testing.assert_allclose(f[0],base.energy_and_forces_eV(r)[1][0],atol=1e-14)


if __name__=='__main__':unittest.main()
