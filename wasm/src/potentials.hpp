#pragma once
#include <array>
#include <vector>

namespace molecular_car {
using Vec3 = std::array<double, 3>;
struct Response { double charge, alpha; int code; double hetero_a, hetero_b; };
struct Fields {
  double depth, height, decay, spacing, corr_decay, softening, scale;
  double base, slope, gain, wall_width, wall_stiffness, xmin, xmax, halfwidth;
};
struct Grid {
  int nz = 0, nx = 0, ny = 0;
  double px = 0, py = 0, dz = 0, zstart = 0;
  std::vector<double> values;
  std::array<double, 4> interpolate(int code, const Vec3& position) const;
};
// Geometry-only conservative candidate lists. Each cell contains every site
// whose 16 A sphere intersects the CLOSED cell; exact cutoff stays in the kernel.
struct SiteStencil {
  int nx=4, ny=6, nz=0;
  double zmin=0, dz=1;
  std::vector<size_t> offsets;
  std::vector<int> sites;
};
struct SurfaceWorkspace {
  std::vector<Vec3> normals, first, second, normal_gradient;
  std::vector<double> lengths;
};
struct FourierProxy {
  bool enabled=false;
  int loaded=0;
  int nz=0, na=0, modes=0, charts=1;
  double astep=.2, bound=2, blend=.4, zstep=.05, zstart=2.65;
  double zmin=2.8, zmax=8, zblend=.3;
  std::vector<double> coefficients; // [CH,chart,z,a,b,value/da/db/dadb,mode,re/im]
  std::vector<std::array<int,2>> indices;
  mutable std::array<int,3> counts{}; // full / blended / exact direction queries
};
struct Substrate {
  std::vector<Vec3> tile;
  std::array<std::array<std::array<double, 8>, 2>, 2> parameters;
  Grid saip, hetero;
  double epsilon = 1e-6;
  // 0: archived algorithm, 1: fused phases, 2: averaged grids, 3: site stencil.
  // Diagnostic switch only; production always starts at level 3.
  int optimization=3;
  mutable Grid saip_mean, hetero_mean;
  mutable bool means_ready=false, stencil_ready=false;
  mutable SiteStencil stencil;
  mutable SurfaceWorkspace workspace;
  FourierProxy fourier;
};
struct NormalTopology {
  std::vector<std::array<int, 3>> neighbours;
  std::vector<int> counts;
};
double fields(const std::vector<Vec3>& x, const std::vector<Response>& atoms,
              const Fields& p, int mode, double strength, const Vec3& tip,
              double corrugation, bool continuum, std::vector<Vec3>& force, double growth=1.);
double substrate(const std::vector<Vec3>& x, const std::vector<Response>& atoms,
                 const Substrate& surface, const NormalTopology& topology,
                 std::vector<Vec3>& force, double lattice_scale=1.);
}  // namespace molecular_car
