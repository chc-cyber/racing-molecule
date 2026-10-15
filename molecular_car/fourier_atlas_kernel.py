"""Compiled conservative v2 atlas with fused six-phase crystal contrast.

Normal geometry and its chain derivative are evaluated once. The six lateral
translations act analytically on Fourier modes; radial tables and exact
fallbacks retain the same translated sums, including all height derivatives.
"""
import math
import numpy as np
from .periodic_au_kernel import (numba, NUMBA_VERSION, _taper, _radial,
    _anisotropic, _weights, _interpolate, _cache_weight, _normal_geometry)
from .fourier_au_kernel import _chain, _box_weight, _hermite

@numba.njit(cache=True,fastmath=False)
def _chart_proxy(coeff,wavevectors,multiplicity,code,x,y,z,nz,a,b,
           astep,bound,zstep,zstart,period_x,period_y,phase_mean,contrast,contrast_dz):
    az,bz=(a+bound)/astep,(b+bound)/astep
    sz=(z-zstart)/zstep
    ia,ib,kz=math.floor(az),math.floor(bz),math.floor(sz)
    ha,da=_hermite(az-ia,astep); hb,db=_hermite(bz-ib,astep)
    wz,dz=_weights(sz-kz)
    values=np.zeros((4,len(wavevectors)),dtype=np.complex128)
    # Traverse the contiguous mode axis once per interpolation node. Compute
    # scalar tensor weights once, rather than twelve times per complex mode.
    for k in range(4):
        for i in range(4):
            for j in range(4):
                component=i//2+2*(j//2)
                w0=ha[i]*hb[j]*wz[k]
                wa=da[i]*hb[j]*wz[k]
                wb=ha[i]*db[j]*wz[k]
                wzz=ha[i]*hb[j]*dz[k]/zstep
                for mode in range(len(wavevectors)):
                    c=coeff[code,kz+k-1,ia+i%2,ib+j%2,component,mode]
                    values[0,mode]+=c*w0
                    values[1,mode]+=c*wa
                    values[2,mode]+=c*wb
                    values[3,mode]+=c*wzz
    e,gx,gy,gz,ga,gb=0.,0.,0.,0.,0.,0.
    x,y=x%period_x,y%period_y
    for mode in range(len(wavevectors)):
        phase=wavevectors[mode,0]*x+wavevectors[mode,1]*y
        phasevalue=complex(math.cos(phase),math.sin(phase))*multiplicity[mode]
        factor=phase_mean[mode]+contrast*(1-phase_mean[mode])
        v=phasevalue*values[0,mode]*factor
        e+=v.real
        gx-=wavevectors[mode,0]*v.imag
        gy-=wavevectors[mode,1]*v.imag
        gz+=(phasevalue*(values[3,mode]*factor+values[0,mode]*contrast_dz*(1-phase_mean[mode]))).real
        ga+=(phasevalue*values[1,mode]*factor).real
        gb+=(phasevalue*values[2,mode]*factor).real
    return e,-gx,-gy,-gz,ga/nz,gb/nz,-(a*ga+b*gb)/nz


@numba.njit(cache=True,fastmath=False)
def _atlas_proxy(coeff,wavevectors,multiplicity,code,x,y,z,normal,
                 astep,bound,blend,zstep,zstart,px,py,phase_mean,contrast,contrast_dz):
    # Each chart rotates the normal coordinates, never the substrate x/y.
    axes=((0,1,2),(1,2,0),(2,0,1))
    energy=0.; force=np.zeros(3); gradient=np.zeros(3)
    weight_sum=0.; weight_gradient=np.zeros(3)
    for chart in range(3):
        ax,ay,az=axes[chart]
        nz=normal[az]
        if abs(nz)<1e-12:
            continue
        a,b=normal[ax]/nz,normal[ay]/nz
        wa,da=_box_weight(a,-bound,bound,blend)
        wb,db=_box_weight(b,-bound,bound,blend)
        w=wa*wb
        if w==0:
            continue
        v=_chart_proxy(coeff[:,chart],wavevectors,multiplicity,code,x,y,z,nz,a,b,
            astep,bound,zstep,zstart,px,py,phase_mean,contrast,contrast_dz)
        dw=np.zeros(3)
        dw[ax]=da*wb/nz; dw[ay]=wa*db/nz; dw[az]=-(a*da*wb+b*wa*db)/nz
        energy+=w*v[0]; weight_sum+=w
        for axis in range(3):
            force[axis]+=w*v[axis+1]
            weight_gradient[axis]+=dw[axis]
            gradient[axis]+=dw[axis]*v[0]
        gradient[ax]+=w*v[4]; gradient[ay]+=w*v[5]; gradient[az]+=w*v[6]
    if weight_sum<=0:
        raise ValueError('No v2 chart covers this normal.')
    energy/=weight_sum
    return energy,force/weight_sum,(gradient-weight_gradient*energy)/weight_sum


@numba.njit(cache=True,fastmath=False)
def energy_and_forces(positions,codes,normal_cutoffs,tile,parameters,grid,
                      px,py,grid_dz,z_start,has_cache,coefficients,
                      wavevectors,multiplicity,angular_step,chart_bound,chart_blend,
                      height_step,fourier_z_start,height_min,height_max,height_blend,
                      topology_neighbours,topology_counts,normal_epsilon,
                      phase_mean,scale,center_weights,center_dz):
    normals,neighbours,counts,first,second,length=_normal_geometry(
        positions,normal_cutoffs,topology_neighbours,topology_counts,normal_epsilon)
    count=len(positions)
    force=np.zeros((count,3)); normal_gradient=np.zeros((count,3))
    proxy_counts=np.zeros(3,dtype=np.int64)
    total=0.
    for i in range(count):
        code=codes[i]
        if code<0:
            continue
        p,gold=parameters[code,0],parameters[code,1]
        x,y,z=positions[i,0],positions[i,1],positions[i,2]
        nx,ny,nz=normals[i,0],normals[i,1],normals[i,2]
        t=min(1.,max(0.,(z-2.)/.6))
        contrast=1.+(scale-1.)*t**3*(10+t*(-15+6*t))
        contrast_dz=(scale-1.)*30*t*t*(t-1)**2/.6
        rw,rd=_cache_weight(z) if has_cache else (0.,0.)
        direction=code==0 or counts[i]>=2
        w,dw=1.,0.
        gn=np.zeros(3); pe=0.; pf=np.zeros(3); pg=np.zeros(3)
        if direction:
            wz,dz=_box_weight(z,height_min,height_max,height_blend)
            t=min(1.,max(0.,(nx*nx+ny*ny+nz*nz-.999)/.0009))
            wn=t**3*(10+t*(-15+6*t)); dn=30*t*t*(1-t)**2/.0009
            w=wz*wn; dw=dz*wn
            gn=wz*dn*2*normals[i]
            proxy_counts[0 if w==1 else 1 if w>0 else 2]+=1
            if w>0:
                pe,pf,pg=_atlas_proxy(coefficients,wavevectors,multiplicity,code,x,y,z,normals[i],
                    angular_step,chart_bound,chart_blend,height_step,fourier_z_start,
                    px,py,phase_mean,contrast,contrast_dz)
        ce=0.; cf=np.zeros(3); be=0.; bf=np.zeros(3)
        ee=0.; ef=np.zeros(3); eg=np.zeros(3)
        for phase in range(1 if scale==1. else 6):
            weight=(1+5*contrast)/6 if phase==0 else (1-contrast)/6
            derivative=5*contrast_dz/6 if phase==0 else -contrast_dz/6
            xq=(x+(phase//3)*px/2)%px
            yq=(y+(phase%3)*py/3)%py
            if rw>0:
                e,fx,fy,fz=_interpolate(grid,code,xq,yq,z,px,py,grid_dz,z_start)
                ce+=weight*e
                cf[0]+=weight*fx; cf[1]+=weight*fy; cf[2]+=weight*fz-derivative*e
            if rw<1 or (direction and w<1):
                for site in range(len(tile)):
                    dx,dy,dzz=xq-tile[site,0],yq-tile[site,1],z-tile[site,2]
                    square=dx*dx+dy*dy+dzz*dzz
                    if square>=256.:
                        continue
                    if square<1e-16:
                        raise ValueError('A molecular atom coincides with a gold lattice site.')
                    radius=math.sqrt(square); taper,td=_taper(radius)
                    fe,ffx,ffy,ffz=0.,0.,0.,0.
                    if code==1:
                        fixed=_anisotropic(dx,dy,dzz,radius,0.,0.,1.,p,taper,td)
                        fe,ffx,ffy,ffz=fixed[0],fixed[1],fixed[2],fixed[3]
                    if rw<1:
                        radial,rf=_radial(radius,p,gold,taper,td)
                        be+=weight*(radial+fe)
                        bf[0]+=weight*(rf*dx+ffx); bf[1]+=weight*(rf*dy+ffy)
                        bf[2]+=weight*(rf*dzz+ffz)-derivative*(radial+fe)
                    if direction and w<1:
                        actual=_anisotropic(dx,dy,dzz,radius,nx,ny,nz,p,taper,td)
                        ee+=weight*(actual[0]-fe)
                        ef[0]+=weight*(actual[1]-ffx); ef[1]+=weight*(actual[2]-ffy)
                        ef[2]+=weight*(actual[3]-ffz)-derivative*(actual[0]-fe)
                        for axis in range(3):
                            eg[axis]+=weight*actual[axis+4]
        atomic_energy=(1-rw)*be+rw*ce+w*pe+(1-w)*ee
        total+=center_weights[i]*atomic_energy
        for axis in range(3):
            force[i,axis]=center_weights[i]*((1-rw)*bf[axis]+rw*cf[axis]+w*pf[axis]+(1-w)*ef[axis])
            normal_gradient[i,axis]=center_weights[i]*(w*pg[axis]+(1-w)*eg[axis]+gn[axis]*(pe-ee))
        force[i,2]-=center_weights[i]*(rd*(ce-be)+dw*(pe-ee))+center_dz[i]*atomic_energy
    return total,_chain(force,normal_gradient,normals,neighbours,counts,first,second,length),proxy_counts


