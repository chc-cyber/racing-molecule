"""Single-thread conservative bounded Fourier/Hermite SAIP proxy kernel."""
from __future__ import annotations
import math
import numpy as np
from .periodic_au_kernel import numba,NUMBA_VERSION,_taper,_radial,_anisotropic,_weights,_interpolate,_cache_weight,_normal_geometry

@numba.njit(cache=True,fastmath=False)
def _chain(force,normal_gradient,normals,neighbours,neighbour_count,first,second,normal_length):
    count=len(normals)
    # Reverse-mode derivative of the normalized local plane. The molecule
    # is never folded into an Au cell while finding or applying this chain.
    for i in range(count):
        if neighbour_count[i] < 2:
            continue
        nx, ny, nz = normals[i, 0], normals[i, 1], normals[i, 2]
        dot = normal_gradient[i, 0]*nx+normal_gradient[i, 1]*ny+normal_gradient[i, 2]*nz
        gx = (normal_gradient[i, 0]-dot*nx)/normal_length[i]
        gy = (normal_gradient[i, 1]-dot*ny)/normal_length[i]
        gz = (normal_gradient[i, 2]-dot*nz)/normal_length[i]
        ax, ay, az = first[i, 0], first[i, 1], first[i, 2]
        bx, by, bz = second[i, 0], second[i, 1], second[i, 2]
        gax, gay, gaz = by*gz-bz*gy, bz*gx-bx*gz, bx*gy-by*gx
        gbx, gby, gbz = gy*az-gz*ay, gz*ax-gx*az, gx*ay-gy*ax
        n0, n1 = neighbours[i, 0], neighbours[i, 1]
        pivot = neighbours[i, 2] if neighbour_count[i] == 3 else i
        force[n0, 0] -= gax
        force[n0, 1] -= gay
        force[n0, 2] -= gaz
        force[n1, 0] -= gbx
        force[n1, 1] -= gby
        force[n1, 2] -= gbz
        force[pivot, 0] += gax+gbx
        force[pivot, 1] += gay+gby
        force[pivot, 2] += gaz+gbz
    return force


@numba.njit(cache=True,fastmath=False,inline="always")
def _box_weight(value,low,high,width):
    tl=min(1.,max(0.,(value-low)/width))
    th=min(1.,max(0.,(high-value)/width))
    wl=tl**3*(10+tl*(-15+6*tl)); wh=th**3*(10+th*(-15+6*th))
    dl=30*tl*tl*(1-tl)**2; dh=30*th*th*(1-th)**2
    return wl*wh,(dl*wh-wl*dh)/width


@numba.njit(cache=True,fastmath=False,inline="always")
def _hermite(t,spacing):
    return ((2*t**3-3*t*t+1,-2*t**3+3*t*t,
             spacing*(t**3-2*t*t+t),spacing*(t**3-t*t)),
            ((6*t*t-6*t)/spacing,(-6*t*t+6*t)/spacing,
             3*t*t-4*t+1,3*t*t-2*t))


@numba.njit(cache=True,fastmath=False)
def _proxy_weight(z,nx,ny,nz,zmin,zmax,zblend,bound,blend):
    if abs(nz)<1e-12:
        return 0.,0.,0.,0.,0.,0.,0.
    a,b=nx/nz,ny/nz
    wz,dz=_box_weight(z,zmin,zmax,zblend)
    wa,da=_box_weight(a,-bound,bound,blend)
    wb,db=_box_weight(b,-bound,bound,blend)
    # The angle proxy is fitted on the unit sphere. A smooth conservative
    # weight returns regularized, non-unit planes to the exact lattice sum.
    t=min(1.,max(0.,(nx*nx+ny*ny+nz*nz-.999)/.0009))
    wn=t**3*(10+t*(-15+6*t)); dn=30*t*t*(1-t)**2/.0009
    w=wz*wa*wb
    return (w*wn,dz*wa*wb*wn,wz*da*wb/nz*wn+w*dn*2*nx,
            wz*wa*db/nz*wn+w*dn*2*ny,
            -wz*(da*wb*a+wa*db*b)/nz*wn+w*dn*2*nz,a,b)


@numba.njit(cache=True,fastmath=False)
def _proxy(coeff,wavevectors,multiplicity,code,x,y,z,nz,a,b,
           astep,bound,zstep,zstart,period_x,period_y):
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
        v=phasevalue*values[0,mode]
        e+=v.real
        gx-=wavevectors[mode,0]*v.imag
        gy-=wavevectors[mode,1]*v.imag
        gz+=(phasevalue*values[3,mode]).real
        ga+=(phasevalue*values[1,mode]).real
        gb+=(phasevalue*values[2,mode]).real
    return e,-gx,-gy,-gz,ga/nz,gb/nz,-(a*ga+b*gb)/nz


@numba.njit(cache=True,fastmath=False)
def energy_and_forces(positions,codes,normal_cutoffs,tile,parameters,grid,
                      period_x,period_y,grid_dz,z_start,has_cache,coefficients,
                      wavevectors,multiplicity,angular_step,chart_bound,chart_blend,
                      height_step,fourier_z_start,height_min,height_max,height_blend,
                      topology_neighbours,topology_counts,normal_epsilon,center_weights=None,center_dz=None):
    normals,neighbours,counts,first,second,length=_normal_geometry(positions,normal_cutoffs,topology_neighbours,topology_counts,normal_epsilon)
    count=len(positions)
    force=np.zeros((count,3)); normal_gradient=np.zeros((count,3))
    proxy_counts=np.zeros(3,dtype=np.int64)
    total=0.
    for i in range(count):
        code=codes[i]
        if code < 0:
            continue
        p,gold=parameters[code,0],parameters[code,1]
        x,y,z=positions[i,0],positions[i,1],positions[i,2]
        nx,ny,nz=normals[i,0],normals[i,1],normals[i,2]
        rw,rd=_cache_weight(z) if has_cache else (0.,0.)
        ce,cfx,cfy,cfz=0.,0.,0.,0.
        if rw>0:
            ce,cfx,cfy,cfz=_interpolate(grid,code,x,y,z,period_x,period_y,grid_dz,z_start)
        direction=code==0 or counts[i]>=2
        w,dw,gnx,gny,gnz,a,b=1.,0.,0.,0.,0.,0.,0.
        pe,pfx,pfy,pfz,pgx,pgy,pgz=0.,0.,0.,0.,0.,0.,0.
        if direction:
            w,dw,gnx,gny,gnz,a,b=_proxy_weight(z,nx,ny,nz,height_min,height_max,
                height_blend,chart_bound,chart_blend)
            proxy_counts[0 if w==1 else 1 if w>0 else 2]+=1
            if w>0:
                pe,pfx,pfy,pfz,pgx,pgy,pgz=_proxy(coefficients,wavevectors,multiplicity,
                    code,x,y,z,nz,a,b,angular_step,chart_bound,height_step,fourier_z_start,period_x,period_y)
        be,bfx,bfy,bfz,ee,efx,efy,efz,egx,egy,egz=0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
        if rw<1 or (direction and w<1):
            xq,yq=x%period_x,y%period_y
            for site in range(len(tile)):
                dx,dy,dzz=xq-tile[site,0],yq-tile[site,1],z-tile[site,2]
                square=dx*dx+dy*dy+dzz*dzz
                if square>=256.:
                    continue
                if square<1e-16:
                    raise ValueError("A molecular atom coincides with a gold lattice site.")
                radius=math.sqrt(square); taper,td=_taper(radius)
                fe,ffx,ffy,ffz=0.,0.,0.,0.
                if code==1:
                    fixed=_anisotropic(dx,dy,dzz,radius,0.,0.,1.,p,taper,td)
                    fe,ffx,ffy,ffz=fixed[0],fixed[1],fixed[2],fixed[3]
                if rw<1:
                    radial,rf=_radial(radius,p,gold,taper,td)
                    be+=radial+fe; bfx+=rf*dx+ffx; bfy+=rf*dy+ffy; bfz+=rf*dzz+ffz
                if direction and w<1:
                    actual=_anisotropic(dx,dy,dzz,radius,nx,ny,nz,p,taper,td)
                    ee+=actual[0]-fe; efx+=actual[1]-ffx; efy+=actual[2]-ffy; efz+=actual[3]-ffz
                    egx+=actual[4]; egy+=actual[5]; egz+=actual[6]
        total+=(1-rw)*be+rw*ce+w*pe+(1-w)*ee
        force[i,0]+=(1-rw)*bfx+rw*cfx+w*pfx+(1-w)*efx
        force[i,1]+=(1-rw)*bfy+rw*cfy+w*pfy+(1-w)*efy
        force[i,2]+=(1-rw)*bfz+rw*cfz+w*pfz+(1-w)*efz-rd*(ce-be)-dw*(pe-ee)
        normal_gradient[i,0]=w*pgx+(1-w)*egx+gnx*(pe-ee)
        normal_gradient[i,1]=w*pgy+(1-w)*egy+gny*(pe-ee)
        normal_gradient[i,2]=w*pgz+(1-w)*egz+gnz*(pe-ee)
        if center_weights is not None:
            atomic_energy=(1-rw)*be+rw*ce+w*pe+(1-w)*ee
            # Replace the unweighted center contribution before the chain rule.
            total+=(center_weights[i]-1)*atomic_energy
            for axis in range(3):
                force[i,axis]*=center_weights[i]
                normal_gradient[i,axis]*=center_weights[i]
            force[i,2]-=center_dz[i]*atomic_energy
    return total,_chain(force,normal_gradient,normals,neighbours,counts,first,second,length),proxy_counts
