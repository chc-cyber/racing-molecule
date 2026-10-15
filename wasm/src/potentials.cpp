// Float64 translations of the existing conservative Numba kernels.
// No timestep, thermostat, cutoff, normal policy, or response model is changed.
#include "potentials.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace molecular_car {
namespace {
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
Vec3 subtract(const Vec3& a, const Vec3& b) { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
double modulo(double x, double period) { double r=std::fmod(x,period); return r<0?r+period:r; }
int modulo(int x,int period) { int r=x%period;return r<0?r+period:r; }
std::array<double,2> taper(double r) {
  double x=r/16.0, square=x*x, cube=square*x, offset=x-1;
  return {1+square*square*(-35+x*(84+x*(-70+20*x))),140*cube*offset*offset*offset/16.0};
}
std::array<double,2> cache_weight(double z) {
  double l=std::clamp((z-2.0)/.3,0.0,1.0),h=std::clamp((12.0-z)/.3,0.0,1.0);
  double lw=l*l*l*(10+l*(-15+6*l)),hw=h*h*h*(10+h*(-15+6*h));
  double ld=30*l*l*(l-1)*(l-1)/.3,hd=30*h*h*(h-1)*(h-1)/.3;
  return {lw*hw,ld*hw-lw*hd};
}
std::array<std::array<double,4>,2> weights(double t) {
  double s=t*t,c=s*t;
  return {{{-.5*t+s-.5*c,1-2.5*s+1.5*c,.5*t+2*s-1.5*c,-.5*s+.5*c},
           {-.5+2*t-1.5*s,-5*t+4.5*s,.5+4*t-4.5*s,-t+1.5*s}}};
}
std::array<double,2> softplus(double t) {
  if(t>=0) { double e=std::exp(-t);return {t+std::log1p(e),1/(1+e)}; }
  double e=std::exp(t);return {std::log1p(e),e/(1+e)};
}
std::array<double,2> wall(double x,double lower,double upper,const Fields& p) {
  auto a=softplus((lower-x)/p.wall_width),b=softplus((x-upper)/p.wall_width);
  return {.5*p.wall_stiffness*p.wall_width*p.wall_width*(a[0]*a[0]+b[0]*b[0]),
          p.wall_stiffness*p.wall_width*(a[0]*a[1]-b[0]*b[1])};
}
std::array<double,2> radial(double r,const std::array<double,8>& p,
                          const std::array<double,8>& gold,double t,double dt) {
  double organic=.5*p[3]*std::exp(-p[1]*(r-p[0]));
  double rep=(.5*gold[3]+gold[4])*std::exp(-gold[1]*(r-gold[0]));
  double argument=p[5]*(r/p[6]-1);
  double damp;
  if(argument>=0) damp=1/(1+std::exp(-argument));
  else { double e=std::exp(argument);damp=e/(1+e); }
  double raw=organic+rep-p[7]*std::pow(r,-6)*damp;
  double derivative=-p[1]*organic-gold[1]*rep+p[7]*std::pow(r,-7)*damp*(6-r*p[5]/p[6]*(1-damp));
  return {t*raw,-(dt*raw+t*derivative)/r};
}
std::array<double,7> anisotropic(const Vec3& d,double r,const Vec3& n,
                               const std::array<double,8>& p,double t,double dt) {
  double projection=dot(d,n),rho=dot(d,d)-projection*projection;
  double exp=p[4]*std::exp(-p[1]*(r-p[0])-rho/p[2]),energy=t*exp;
  double rf=-(dt-p[1]*t)*exp/r,lf=2*energy/p[2],nf=lf*projection;
  return {energy,rf*d[0]+lf*(d[0]-projection*n[0]),rf*d[1]+lf*(d[1]-projection*n[1]),
          rf*d[2]+lf*(d[2]-projection*n[2]),nf*d[0],nf*d[1],nf*d[2]};
}
void finite(double e,const std::vector<Vec3>& f) {
  if(!std::isfinite(e)) throw std::runtime_error("Nonfinite potential energy");
  for(const auto& v:f) for(double c:v) if(!std::isfinite(c)) throw std::runtime_error("Nonfinite force");
}
}  // namespace

std::array<double,4> Grid::interpolate(int code,const Vec3& p) const {
  double dx=px/nx,dy=py/ny,sx=modulo(p[0],px)/dx,sy=modulo(p[1],py)/dy,sz=(p[2]-zstart)/dz;
  int bx=std::floor(sx),by=std::floor(sy),bz=std::floor(sz);
  auto wx=weights(sx-bx),wy=weights(sy-by),wz=weights(sz-bz);
  std::array<double,4> result{};
  for(int k=0;k<4;k++) {
    int zi=bz+k-1;
    if(zi<0||zi>=nz) throw std::runtime_error("Potential interpolation outside padded grid");
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) {
      int xi=modulo(bx+i-1,nx),yi=modulo(by+j-1,ny);
      double value=values[((code*nz+zi)*nx+xi)*ny+yi];
      result[0]+=value*wz[0][k]*wx[0][i]*wy[0][j];
      result[1]+=value*wz[0][k]*wx[1][i]*wy[0][j];
      result[2]+=value*wz[0][k]*wx[0][i]*wy[1][j];
      result[3]+=value*wz[1][k]*wx[0][i]*wy[0][j];
    }
  }
  result[1]/=-dx;result[2]/=-dy;result[3]/=-dz;
  return result;
}

double fields(const std::vector<Vec3>& x,const std::vector<Response>& atoms,const Fields& p,
              int mode,double strength,const Vec3& tip,double corrugation,bool continuum,
              std::vector<Vec3>& force, double growth) {
  if(mode!=0&&mode!=1) throw std::runtime_error("Unknown drive mode");
  if(!std::isfinite(growth)||growth<.1||growth>2)throw std::runtime_error("Field growth must be in [0.1, 2]");
  double root=std::sqrt(growth);
  double total=0,amplitude=p.scale*strength,base=p.base*strength/root,slope=p.slope*strength*root;
  double kx=2*std::acos(-1.0)/p.spacing,ky=kx/std::sqrt(3.0);
  force.assign(x.size(),{0,0,0});
  for(size_t i=0;i<x.size();i++) {
    double px=x[i][0],py=x[i][1],z=x[i][2],q=atoms[i].charge,alpha=atoms[i].alpha,e=0;
    auto& f=force[i];
    if(continuum) {
      double h=z-p.height,ex=std::exp(-p.decay*h),square=ex*ex;
      e+=p.depth*(square-2*ex);f[2]+=2*p.depth*p.decay*(square-ex);
      if(corrugation!=0) {
        double cx=std::cos(kx*px),sx=std::sin(kx*px),cy=std::cos(ky*py),sy=std::sin(ky*py);
        double pattern=(2*cx*cy+std::cos(2*ky*py))/3;
        double envelope=corrugation*std::exp(-p.corr_decay*h);
        e+=envelope*pattern;
        f[0]+=envelope*(2*kx/3)*sx*cy;
        f[1]+=envelope*(2*ky/3)*(cx*sy+std::sin(2*ky*py));
        f[2]+=p.corr_decay*envelope*pattern;
      }
    }
    auto bx=wall(px,p.xmin,p.xmax,p),by=wall(py,-p.halfwidth,p.halfwidth,p);
    e+=bx[0]+by[0];f[0]+=bx[1];f[1]+=by[1];
    if(mode==1) {
      double ex=slope*z,ez=base+slope*px,response=p.gain*alpha;
      e+=-q*ez*z-.5*response*(ex*ex+ez*ez);
      f[0]+=q*ex+response*slope*ez;f[2]+=q*ez+response*slope*slope*z;
    } else if(amplitude!=0) {
      Vec3 d=subtract(x[i],tip);double r2=dot(d,d),s2=p.softening*p.softening;
      double inv=1/std::sqrt(r2+s2),inv2=inv*inv,ef=amplitude*inv*inv2;
      e+=q*amplitude*inv-.5*alpha*ef*ef*r2;
      double factor=q*ef+alpha*amplitude*amplitude*(s2-2*r2)*std::pow(inv2,4);
      for(int a=0;a<3;a++) f[a]+=factor*d[a];
    }
    total+=e;
  }
  finite(total,force);return total;
}

double raw_substrate(const std::vector<Vec3>& x,const std::vector<Response>& atoms,const Substrate& s,
                 const NormalTopology& topology,std::vector<Vec3>& force,
                 const std::vector<double>& center_weights={},const std::vector<double>& center_dz={}) {
  size_t count=x.size();
  std::vector<Vec3> normals(count,{0,0,1}),first(count),second(count),ng(count);
  std::vector<double> lengths(count,1);
  force.assign(count,{0,0,0});double total=0;
  for(size_t i=0;i<count;i++) if(topology.counts[i]>=2) {
    auto ids=topology.neighbours[i];int pivot=topology.counts[i]==3?ids[2]:i;
    first[i]=subtract(x[ids[0]],x[pivot]);second[i]=subtract(x[ids[1]],x[pivot]);
    auto n=cross(first[i],second[i]);lengths[i]=std::sqrt(dot(n,n)+s.epsilon*s.epsilon);
    for(int a=0;a<3;a++) normals[i][a]=n[a]/lengths[i];
  }
  for(size_t i=0;i<count;i++) {
    int code=atoms[i].code;auto wd=cache_weight(x[i][2]);double w=wd[0],dw=wd[1];
    double base_e=0,anis_e=0,cached_e=0;Vec3 base{},anis{},cached{};
    if(w>0) {
      if(code>=0) {
        auto v=s.saip.interpolate(code,x[i]);cached_e=v[0];for(int a=0;a<3;a++) cached[a]=v[a+1];
      } else {
        for(int c=0;c<2;c++) {
          auto v=s.hetero.interpolate(c,x[i]);double scale=c==0?atoms[i].hetero_a:atoms[i].hetero_b;
          cached_e+=scale*v[0];for(int a=0;a<3;a++) cached[a]+=scale*v[a+1];
        }
      }
    }
    bool direction=code==0||(code==1&&topology.counts[i]>=2);
    if(w<1||direction) {
      Vec3 local={modulo(x[i][0],s.saip.px),modulo(x[i][1],s.saip.py),x[i][2]};
      for(const auto& site:s.tile) {
        Vec3 d=subtract(local,site);double square=dot(d,d);
        if(square>=256) continue;
        if(square<1e-16) throw std::runtime_error("Molecular atom coincides with Au lattice site");
        double r=std::sqrt(square);auto td=taper(r);double t=td[0],dt=td[1];
        if(code<0) {
          double r6=std::pow(square,-3),a=atoms[i].hetero_a,b=atoms[i].hetero_b;
          double raw=a*r6*r6+b*r6,derivative=(-12*a*r6*r6-6*b*r6)/r;
          double factor=-(dt*raw+t*derivative)/r;
          base_e+=t*raw;for(int axis=0;axis<3;axis++) base[axis]+=factor*d[axis];
        } else {
          const auto& p=s.parameters[code][0];const auto& gold=s.parameters[code][1];
          std::array<double,7> fixed{};
          if(code==1) fixed=anisotropic(d,r,{0,0,1},p,t,dt);
          if(w<1) {
            auto radial_value=radial(r,p,gold,t,dt);base_e+=radial_value[0]+fixed[0];
            for(int axis=0;axis<3;axis++) base[axis]+=radial_value[1]*d[axis]+fixed[axis+1];
          }
          if(direction) {
            auto actual=anisotropic(d,r,normals[i],p,t,dt);anis_e+=actual[0]-fixed[0];
            for(int axis=0;axis<3;axis++) { anis[axis]+=actual[axis+1]-fixed[axis+1];ng[i][axis]+=actual[axis+4]; }
          }
        }
      }
    }
    double atomic_energy=(1-w)*base_e+w*cached_e+anis_e;
    total+=center_weights.empty()?atomic_energy:center_weights[i]*atomic_energy;
    for(int a=0;a<3;a++) force[i][a]+=(1-w)*base[a]+w*cached[a]+anis[a];
    force[i][2]-=dw*(cached_e-base_e);
    if(!center_weights.empty()) {
      for(int a=0;a<3;a++){force[i][a]*=center_weights[i];ng[i][a]*=center_weights[i];}
      force[i][2]-=center_dz[i]*atomic_energy;
    }
  }
  for(size_t i=0;i<count;i++) if(topology.counts[i]>=2) {
    double projection=dot(ng[i],normals[i]);Vec3 g;
    for(int a=0;a<3;a++) g[a]=(ng[i][a]-projection*normals[i][a])/lengths[i];
    auto ga=cross(second[i],g),gb=cross(g,first[i]);
    auto ids=topology.neighbours[i];int pivot=topology.counts[i]==3?ids[2]:i;
    for(int a=0;a<3;a++) {force[ids[0]][a]-=ga[a];force[ids[1]][a]-=gb[a];force[pivot][a]+=ga[a]+gb[a];}
  }
  constexpr double ev_to_kcal_mol=23.0605478306;
  total*=ev_to_kcal_mol;for(auto& f:force) for(auto& c:f) c*=ev_to_kcal_mol;
  finite(total,force);return total;
}

double reference_substrate(const std::vector<Vec3>& x,const std::vector<Response>& atoms,const Substrate& s,
                 const NormalTopology& topology,std::vector<Vec3>& force,double scale) {
  if(!std::isfinite(scale)||scale<0||scale>12)throw std::runtime_error("Lattice scale must be in [0, 12]");
  if(scale==1.)return raw_substrate(x,atoms,s,topology,force);
  std::vector<double> factors(x.size()),derivatives(x.size()),weights(x.size()),dz(x.size());
  for(size_t k=0;k<x.size();k++) {
    double t=std::clamp((x[k][2]-2.)/.6,0.,1.);
    factors[k]=1+(scale-1)*t*t*t*(10+t*(-15+6*t));
    derivatives[k]=(scale-1)*30*t*t*(t-1)*(t-1)/.6;
    weights[k]=(1+5*factors[k])/6;dz[k]=5*derivatives[k]/6;
  }
  double total=raw_substrate(x,atoms,s,topology,force,weights,dz);
  for(size_t k=0;k<x.size();k++){weights[k]=(1-factors[k])/6;dz[k]=-derivatives[k]/6;}
  for(int i=0;i<2;i++)for(int j=0;j<3;j++)if(i!=0||j!=0) {
    auto shifted=x;for(auto& r:shifted){r[0]+=i*s.saip.px/2;r[1]+=j*s.saip.py/3;}
    std::vector<Vec3> extra;total+=raw_substrate(shifted,atoms,s,topology,extra,weights,dz);
    for(size_t k=0;k<x.size();k++)for(int a=0;a<3;a++)force[k][a]+=extra[k][a];
  }
  return total;
}

namespace {
std::array<double,2> smooth(double t) {
  t=std::clamp(t,0.,1.);return {t*t*t*(10+t*(-15+6*t)),30*t*t*(t-1)*(t-1)};
}
std::array<double,2> box_weight(double v,double lo,double hi,double width) {
  auto l=smooth((v-lo)/width),h=smooth((hi-v)/width);
  return {l[0]*h[0],(l[1]*h[0]-l[0]*h[1])/width};
}
std::array<std::array<double,4>,2> hermite(double t,double step) {
  double t2=t*t,t3=t2*t;
  return {{{2*t3-3*t2+1,-2*t3+3*t2,step*(t3-2*t2+t),step*(t3-t2)},
    {(6*t2-6*t)/step,(-6*t2+6*t)/step,3*t2-4*t+1,3*t2-2*t}}};
}
// [weight,dweight/dz,dweight/dnx,dweight/dny,dweight/dnz,a,b]
std::array<double,7> proxy_weight(const FourierProxy& p,double z,const Vec3& n) {
  if(p.charts==3) {
    auto wz=box_weight(z,p.zmin,p.zmax,p.zblend),wn=smooth((dot(n,n)-.999)/.0009);
    double dn=wz[0]*wn[1]/.0009;
    return {wz[0]*wn[0],wz[1]*wn[0],2*n[0]*dn,2*n[1]*dn,2*n[2]*dn,0,0};
  }
  if(std::abs(n[2])<1e-12)return {};
  double a=n[0]/n[2],b=n[1]/n[2];
  auto wz=box_weight(z,p.zmin,p.zmax,p.zblend),wa=box_weight(a,-p.bound,p.bound,p.blend),wb=box_weight(b,-p.bound,p.bound,p.blend);
  auto wn=smooth((dot(n,n)-.999)/.0009);double dn=wn[1]/.0009,w=wz[0]*wa[0]*wb[0];
  return {w*wn[0],wz[1]*wa[0]*wb[0]*wn[0],
    wz[0]*wa[1]*wb[0]/n[2]*wn[0]+w*dn*2*n[0],
    wz[0]*wa[0]*wb[1]/n[2]*wn[0]+w*dn*2*n[1],
    -wz[0]*(wa[1]*wb[0]*a+wa[0]*wb[1]*b)/n[2]*wn[0]+w*dn*2*n[2],a,b};
}
// One Fourier interpolation covers all six translations. Their average retains
// modes with even m AND n divisible by three, not just the constant mode.
std::array<double,7> proxy_value(const Substrate& s,int code,const Vec3& x,const Vec3& n,
    double a,double b,double factor,double factor_dz,int chart=0) {
  const auto& p=s.fourier;double ax=(a+p.bound)/p.astep,by=(b+p.bound)/p.astep,zz=(x[2]-p.zstart)/p.zstep;
  int ia=std::floor(ax),ib=std::floor(by),iz=std::floor(zz);
  if(ia<0||ib<0||ia+1>=p.na||ib+1>=p.na||iz<1||iz+2>=p.nz)throw std::runtime_error("Fourier query outside padded support");
  auto ha=hermite(ax-ia,p.astep),hb=hermite(by-ib,p.astep),wz=weights(zz-iz);
  // Fixed upper bound is checked at asset installation; no per-atom allocations.
  double value[64][8]{};
  for(int k=0;k<4;k++)for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
    int component=i/2+2*(j/2);
    double tensor[4]={ha[0][i]*hb[0][j]*wz[0][k],ha[1][i]*hb[0][j]*wz[0][k],
      ha[0][i]*hb[1][j]*wz[0][k],ha[0][i]*hb[0][j]*wz[1][k]/p.zstep};
    size_t base=((((((code*p.charts+chart)*p.nz+iz+k-1)*p.na+ia+i%2)*p.na+ib+j%2)*4+component)*p.modes)*2;
    for(int m=0;m<p.modes;m++)for(int d=0;d<4;d++) {
      value[m][2*d]+=p.coefficients[base+2*m]*tensor[d];
      value[m][2*d+1]+=p.coefficients[base+2*m+1]*tensor[d];
    }
  }
  double e=0,gx=0,gy=0,gz=0,ga=0,gb=0;
  double xx=modulo(x[0],s.saip.px),yy=modulo(x[1],s.saip.py),pi=std::acos(-1.);
  for(int m=0;m<p.modes;m++) {
    auto index=p.indices[m];double kx=2*pi*index[0]/s.saip.px,ky=2*pi*index[1]/s.saip.py;
    double phase=kx*xx+ky*yy,multiplicity=index[0]||index[1]?2:1;
    double re=std::cos(phase)*multiplicity,im=std::sin(phase)*multiplicity;
    bool retained=index[0]%2==0&&index[1]%3==0;
    double amplitude=retained?1:factor,da=retained?0:factor_dz;
    auto& v=value[m];double energy=re*v[0]-im*v[1],imaginary=re*v[1]+im*v[0];
    e+=amplitude*energy;gx-=amplitude*kx*imaginary;gy-=amplitude*ky*imaginary;
    ga+=amplitude*(re*v[2]-im*v[3]);gb+=amplitude*(re*v[4]-im*v[5]);
    gz+=amplitude*(re*v[6]-im*v[7])+da*energy;
  }
  return {e,-gx,-gy,-gz,ga/n[2],gb/n[2],-(a*ga+b*gb)/n[2]};
}
// A partition of unity across three projective charts has no equatorial hole.
// The largest |normal component| supplies at least one full chart. The weights
// and their normal derivatives blend energies, never forces independently.
std::array<double,7> atlas_value(const Substrate& s,int code,const Vec3& x,const Vec3& n,
    double factor,double factor_dz) {
  const auto& p=s.fourier;
  constexpr int axes[3][3]={{0,1,2},{1,2,0},{2,0,1}};
  double total=0;Vec3 weight_gradient{};std::array<double,7> result{};
  for(int chart=0;chart<3;chart++) {
    const auto& axis=axes[chart];Vec3 q={n[axis[0]],n[axis[1]],n[axis[2]]};
    if(std::abs(q[2])<1e-12)continue;
    double a=q[0]/q[2],b=q[1]/q[2];
    auto wa=box_weight(a,-p.bound,p.bound,p.blend),wb=box_weight(b,-p.bound,p.bound,p.blend);
    double w=wa[0]*wb[0];if(w==0)continue;
    Vec3 local_gradient={wa[1]*wb[0]/q[2],wa[0]*wb[1]/q[2],-(a*wa[1]*wb[0]+b*wa[0]*wb[1])/q[2]};
    auto v=proxy_value(s,code,x,q,a,b,factor,factor_dz,chart);
    total+=w;
    for(int k=0;k<4;k++)result[k]+=w*v[k];
    for(int k=0;k<3;k++) {
      int global=axis[k];weight_gradient[global]+=local_gradient[k];
      result[4+global]+=w*v[4+k]+local_gradient[k]*v[0];
    }
  }
  if(total<=0)throw std::runtime_error("Fourier atlas has no chart for supported normal");
  for(auto& v:result)v/=total;
  for(int k=0;k<3;k++)result[k+4]-=result[0]*weight_gradient[k]/total;
  return result;
}
// The usual carbon contact path has no radial fallback or fixed-H subtraction.
// Keep its tight loop separate from the phase accumulators and specialize the
// index path, without changing site order, cutoffs, or anisotropic arithmetic.
template<bool Filtered>
[[gnu::noinline]] std::array<double,7> carbon_sites(const Substrate& s,
    const Vec3& local,const Vec3& normal,size_t begin,size_t end) {
  std::array<double,7> result{};
  const auto& p=s.parameters[0][0];
  for(size_t k=begin;k<end;k++) {
    const auto& site=s.tile[Filtered?s.stencil.sites[k]:k];
    Vec3 d=subtract(local,site);double square=dot(d,d);
    if(square>=256)continue;
    if(square<1e-16)throw std::runtime_error("Molecular atom coincides with Au lattice site");
    double r=std::sqrt(square);auto td=taper(r);
    auto value=anisotropic(d,r,normal,p,td[0],td[1]);
    for(int a=0;a<7;a++)result[a]+=value[a];
  }
  return result;
}
Grid phase_average(const Grid& grid) {
  if(grid.nx%2 || grid.ny%3)throw std::runtime_error("Six-phase cache requires aligned grid dimensions");
  Grid mean=grid;
  const size_t plane=static_cast<size_t>(grid.nx)*grid.ny;
  for(size_t base=0;base<grid.values.size();base+=plane)
    for(int x=0;x<grid.nx;x++)for(int y=0;y<grid.ny;y++) {
      double value=0;
      for(int i=0;i<2;i++)for(int j=0;j<3;j++)
        value+=grid.values[base+modulo(x+i*grid.nx/2,grid.nx)*grid.ny+modulo(y+j*grid.ny/3,grid.ny)]/6.;
      mean.values[base+x*grid.ny+y]=value;
    }
  return mean;
}
void prepare_stencil(const Substrate& s) {
  auto& bins=s.stencil;
  double low=0,high=0;
  for(const auto& site:s.tile){low=std::min(low,site[2]);high=std::max(high,site[2]);}
  bins.zmin=std::floor(low-16.);
  bins.nz=static_cast<int>(std::ceil(high+16.-bins.zmin));
  bins.offsets.clear();bins.sites.clear();
  const double dx=s.saip.px/bins.nx,dy=s.saip.py/bins.ny;
  // Original site order is retained. No COM-based wrapping or minimum-image
  // assumption is made about the molecule's size or its covalent neighbours.
  for(int z=0;z<bins.nz;z++)for(int x=0;x<bins.nx;x++)for(int y=0;y<bins.ny;y++) {
    bins.offsets.push_back(bins.sites.size());
    Vec3 lower={x*dx,y*dy,bins.zmin+z*bins.dz},upper={(x+1)*dx,(y+1)*dy,lower[2]+bins.dz};
    for(size_t k=0;k<s.tile.size();k++) {
      double square=0;
      for(int a=0;a<3;a++) {double gap=std::max({lower[a]-s.tile[k][a],s.tile[k][a]-upper[a],0.});square+=gap*gap;}
      if(square<=256.+1e-10)bins.sites.push_back(static_cast<int>(k));
    }
  }
  bins.offsets.push_back(bins.sites.size());s.stencil_ready=true;
}
std::array<double,4> cached_value(const Substrate& s,const Response& atom,const Vec3& p,bool mean) {
  if(atom.code>=0)return (mean?s.saip_mean:s.saip).interpolate(atom.code,p);
  const auto& grid=mean?s.hetero_mean:s.hetero;
  auto a=grid.interpolate(0,p),b=grid.interpolate(1,p);
  for(int k=0;k<4;k++)a[k]=atom.hetero_a*a[k]+atom.hetero_b*b[k];
  return a;
}
}

double substrate(const std::vector<Vec3>& x,const std::vector<Response>& atoms,const Substrate& s,
                 const NormalTopology& topology,std::vector<Vec3>& force,double scale) {
  if(!std::isfinite(scale)||scale<0||scale>12)throw std::runtime_error("Lattice scale must be in [0, 12]");
  // Preserve the original default exactly, including operation order.
  s.fourier.counts={};
  if(!s.fourier.enabled && (scale==1. || s.optimization==0))return reference_substrate(x,atoms,s,topology,force,scale);
  if(s.optimization>=2 && !s.means_ready) {
    s.saip_mean=phase_average(s.saip);s.hetero_mean=phase_average(s.hetero);s.means_ready=true;
  }
  if(s.optimization>=3 && !s.stencil_ready)prepare_stencil(s);
  const size_t count=x.size();auto& work=s.workspace;
  work.normals.assign(count,{0,0,1});work.first.resize(count);work.second.resize(count);
  work.lengths.assign(count,1.);work.normal_gradient.assign(count,{0,0,0});
  auto& normals=work.normals;auto& ng=work.normal_gradient;
  force.assign(count,{0,0,0});double total=0;
  for(size_t i=0;i<count;i++)if(topology.counts[i]>=2) {
    auto ids=topology.neighbours[i];int pivot=topology.counts[i]==3?ids[2]:i;
    work.first[i]=subtract(x[ids[0]],x[pivot]);work.second[i]=subtract(x[ids[1]],x[pivot]);
    auto n=cross(work.first[i],work.second[i]);work.lengths[i]=std::sqrt(dot(n,n)+s.epsilon*s.epsilon);
    for(int a=0;a<3;a++)normals[i][a]=n[a]/work.lengths[i];
  }
  for(size_t i=0;i<count;i++) {
    const auto& atom=atoms[i];int code=atom.code;
    double q=std::clamp((x[i][2]-2.)/.6,0.,1.);
    double factor=1+(scale-1)*q*q*q*(10+q*(-15+6*q));
    double derivative=(scale-1)*30*q*q*(q-1)*(q-1)/.6;
    const auto wd=cache_weight(x[i][2]);double w=wd[0],dw=wd[1];
    double ce=0,be=0,ae=0,cdz=0,bdz=0,adz=0;Vec3 cf{},bf{},af{};
    // Interpolation is linear, and shifts are exactly 16/32 and 20/60 cells.
    // Average table VALUES, then differentiate the same interpolant. The core
    // gate is applied after lookup, never resampled into a low-order grid.
    if(w>0 && s.optimization>=2) {
      auto original=cached_value(s,atom,x[i],false),mean=cached_value(s,atom,x[i],true);
      ce=factor*original[0]+(1-factor)*mean[0];cdz=derivative*(original[0]-mean[0]);
      for(int a=0;a<3;a++)cf[a]=factor*original[a+1]+(1-factor)*mean[a+1];
    }
    bool direction=code==0||(code==1&&topology.counts[i]>=2);
    std::array<double,7> pw{},pv{};
    if(s.fourier.enabled&&direction) {
      if(s.fourier.loaded&(1<<code))pw=proxy_weight(s.fourier,x[i][2],normals[i]);
      s.fourier.counts[pw[0]==1?0:pw[0]>0?1:2]++;
      if(pw[0]>0)pv=s.fourier.charts==3 ? atlas_value(s,code,x[i],normals[i],factor,derivative) :
        proxy_value(s,code,x[i],normals[i],pw[5],pw[6],factor,derivative);
    }
    bool exact_direction=direction&&pw[0]<1;
    for(int pi=0;pi<2;pi++)for(int pj=0;pj<3;pj++) {
      bool origin=pi==0&&pj==0;
      double weight=origin?(1+5*factor)/6:(1-factor)/6;
      double dz=origin?5*derivative/6:-derivative/6;
      if(weight==0 && dz==0)continue;
      Vec3 query={x[i][0]+pi*s.saip.px/2,x[i][1]+pj*s.saip.py/3,x[i][2]};
      if(w>0 && s.optimization==1) {
        auto v=cached_value(s,atom,query,false);ce+=weight*v[0];cdz+=dz*v[0];
        for(int a=0;a<3;a++)cf[a]+=weight*v[a+1];
      }
      if(w==1 && !exact_direction)continue;
      Vec3 local={modulo(query[0],s.saip.px),modulo(query[1],s.saip.py),query[2]};
      size_t begin=0,end=s.tile.size();bool filtered=false;
      const auto& bins=s.stencil;
      if(s.optimization>=3 && local[2]>=bins.zmin && local[2]<bins.zmin+bins.nz*bins.dz) {
        int bx=std::min(static_cast<int>(local[0]/s.saip.px*bins.nx),bins.nx-1);
        int by=std::min(static_cast<int>(local[1]/s.saip.py*bins.ny),bins.ny-1);
        int bz=static_cast<int>(std::floor((local[2]-bins.zmin)/bins.dz));
        size_t cell=(bz*bins.nx+bx)*bins.ny+by;
        begin=bins.offsets[cell];end=bins.offsets[cell+1];filtered=true;
      }
      double phase_be=0,phase_ae=0;Vec3 phase_bf{},phase_af{},phase_ng{};
      if(code==0 && w==1) {
        auto v=filtered?carbon_sites<true>(s,local,normals[i],begin,end):carbon_sites<false>(s,local,normals[i],begin,end);
        phase_ae=v[0];
        for(int a=0;a<3;a++){phase_af[a]=v[a+1];phase_ng[a]=v[a+4];}
      } else for(size_t k=begin;k<end;k++) {
        const auto& site=s.tile[filtered?bins.sites[k]:k];
        Vec3 d=subtract(local,site);double square=dot(d,d);
        if(square>=256)continue;
        if(square<1e-16)throw std::runtime_error("Molecular atom coincides with Au lattice site");
        double r=std::sqrt(square);auto td=taper(r);double t=td[0],dt=td[1];
        if(code<0) {
          double r6=std::pow(square,-3),a=atom.hetero_a,b=atom.hetero_b;
          double raw=a*r6*r6+b*r6,radial_derivative=(-12*a*r6*r6-6*b*r6)/r;
          double f=-(dt*raw+t*radial_derivative)/r;
          phase_be+=t*raw;for(int axis=0;axis<3;axis++)phase_bf[axis]+=f*d[axis];
        } else {
          const auto& p=s.parameters[code][0];const auto& gold=s.parameters[code][1];
          std::array<double,7> fixed{};
          if(code==1)fixed=anisotropic(d,r,{0,0,1},p,t,dt);
          if(w<1) {
            auto radial_value=radial(r,p,gold,t,dt);phase_be+=radial_value[0]+fixed[0];
            for(int axis=0;axis<3;axis++)phase_bf[axis]+=radial_value[1]*d[axis]+fixed[axis+1];
          }
          if(exact_direction) {
            auto actual=anisotropic(d,r,normals[i],p,t,dt);phase_ae+=actual[0]-fixed[0];
            for(int axis=0;axis<3;axis++) {
              phase_af[axis]+=actual[axis+1]-fixed[axis+1];phase_ng[axis]+=actual[axis+4];
            }
          }
        }
      }
      be+=weight*phase_be;ae+=weight*phase_ae;bdz+=dz*phase_be;adz+=dz*phase_ae;
      for(int a=0;a<3;a++){bf[a]+=weight*phase_bf[a];af[a]+=weight*phase_af[a];ng[i][a]+=weight*phase_ng[a];}
    }
    if(pw[0]>0) {
      double difference=pv[0]-ae;
      ae=(1-pw[0])*ae+pw[0]*pv[0];adz*=1-pw[0];
      for(int a=0;a<3;a++) {
        af[a]=(1-pw[0])*af[a]+pw[0]*pv[a+1];
        ng[i][a]=(1-pw[0])*ng[i][a]+pw[0]*pv[a+4]+pw[a+2]*difference;
      }
      af[2]-=pw[1]*difference;
    }
    total+=(1-w)*be+w*ce+ae;
    for(int a=0;a<3;a++)force[i][a]+=(1-w)*bf[a]+w*cf[a]+af[a];
    force[i][2]-=dw*(ce-be)+(1-w)*bdz+w*cdz+adz;
  }
  // All phase normal gradients are accumulated BEFORE one chain-rule scatter.
  for(size_t i=0;i<count;i++)if(topology.counts[i]>=2) {
    double projection=dot(ng[i],normals[i]);Vec3 g;
    for(int a=0;a<3;a++)g[a]=(ng[i][a]-projection*normals[i][a])/work.lengths[i];
    auto ga=cross(work.second[i],g),gb=cross(g,work.first[i]);
    auto ids=topology.neighbours[i];int pivot=topology.counts[i]==3?ids[2]:i;
    for(int a=0;a<3;a++){force[ids[0]][a]-=ga[a];force[ids[1]][a]-=gb[a];force[pivot][a]+=ga[a]+gb[a];}
  }
  constexpr double ev_to_kcal_mol=23.0605478306;
  total*=ev_to_kcal_mol;for(auto& f:force)for(auto& c:f)c*=ev_to_kcal_mol;
  finite(total,force);return total;
}
}  // namespace molecular_car
