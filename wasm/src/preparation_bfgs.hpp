// Preparation-only counterpart of molecular_car/preparation_bfgs.py.
// All trial forces use the unchanged joint conservative objective. A return
// requires a separate final-coordinate force gate in the owning engine.
#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace molecular_car {
struct BfgsOptions {
  int maximum_iterations=400, maximum_line_search_trials=30, maximum_cholesky_attempts=12;
  double hessian_step=1e-4, maximum_atom_displacement=.1, armijo=1e-4;
  double curvature_relative_threshold=1e-12, diagonal_shift_relative_start=1e-8;
};
struct PreparationBudgetStop : std::runtime_error { using std::runtime_error::runtime_error; };
struct BfgsResult {
  std::vector<double> positions;
  int iterations=0, force_evaluations=0;
  double diagonal_shift=0;
  std::string reason="iteration-limit";
};
using PreparationEvaluation=std::pair<double,std::vector<double>>;

inline double preparation_dot(const std::vector<double>& a,const std::vector<double>& b) {
  double sum=0;for(size_t i=0;i<a.size();i++)sum+=a[i]*b[i];return sum;
}

inline BfgsResult dense_bfgs(const std::vector<double>& initial,
    const std::function<PreparationEvaluation(const std::vector<double>&)>& evaluate,
    double tolerance,const BfgsOptions& options) {
  BfgsResult result;result.positions=initial;
  const size_t dimension=initial.size();
  auto value=[&](const std::vector<double>& positions) {
    auto evaluated=evaluate(positions);result.force_evaluations++;
    if(!std::isfinite(evaluated.first)||evaluated.second.size()!=dimension)
      throw std::runtime_error("Preparation BFGS received invalid energy/forces");
    for(double f:evaluated.second)if(!std::isfinite(f))
      throw std::runtime_error("Preparation BFGS received nonfinite forces");
    return evaluated;
  };
  try {
    auto state=value(result.positions);
    double energy=state.first;auto force=state.second;
    if(std::sqrt(preparation_dot(force,force))<=tolerance) {result.reason="force-tolerance";return result;}
    std::vector<double> hessian(dimension*dimension),lower(dimension*dimension),inverse(dimension*dimension);
    for(size_t column=0;column<dimension;column++) {
      auto plus=result.positions,minus=result.positions;
      plus[column]+=options.hessian_step;minus[column]-=options.hessian_step;
      auto fp=value(plus).second,fm=value(minus).second;
      for(size_t row=0;row<dimension;row++)hessian[row*dimension+column]=-(fp[row]-fm[row])/(2*options.hessian_step);
    }
    double scale=1;
    for(size_t i=0;i<dimension;i++) {
      scale=std::max(scale,std::abs(hessian[i*dimension+i]));
      for(size_t j=0;j<i;j++)hessian[i*dimension+j]=hessian[j*dimension+i]=
        .5*(hessian[i*dimension+j]+hessian[j*dimension+i]);
    }
    bool factored=false;
    for(int attempt=0;attempt<options.maximum_cholesky_attempts;attempt++) {
      std::fill(lower.begin(),lower.end(),0);factored=true;
      for(size_t i=0;i<dimension&&factored;i++)for(size_t j=0;j<=i;j++) {
        double entry=hessian[i*dimension+j]+(i==j?result.diagonal_shift:0);
        for(size_t k=0;k<j;k++)entry-=lower[i*dimension+k]*lower[j*dimension+k];
        if(i==j) {
          if(!std::isfinite(entry)||entry<=1e-12*scale){factored=false;break;}
          lower[i*dimension+j]=std::sqrt(entry);
        } else lower[i*dimension+j]=entry/lower[j*dimension+j];
      }
      if(factored)break;
      result.diagonal_shift=result.diagonal_shift==0?options.diagonal_shift_relative_start*scale:result.diagonal_shift*10;
    }
    if(!factored){result.reason="curvature-factorization-failed";return result;}
    // Invert the triangular factor, then form L^-T L^-1.
    std::fill(hessian.begin(),hessian.end(),0);
    for(size_t column=0;column<dimension;column++)for(size_t row=column;row<dimension;row++) {
      double entry=row==column?1:0;
      for(size_t k=column;k<row;k++)entry-=lower[row*dimension+k]*hessian[k*dimension+column];
      hessian[row*dimension+column]=entry/lower[row*dimension+row];
    }
    for(size_t i=0;i<dimension;i++)for(size_t j=0;j<=i;j++) {
      double entry=0;for(size_t k=std::max(i,j);k<dimension;k++)entry+=hessian[k*dimension+i]*hessian[k*dimension+j];
      inverse[i*dimension+j]=inverse[j*dimension+i]=entry;
    }
    std::vector<double> direction(dimension),trial(dimension),displacement(dimension),gradient_change(dimension),product(dimension);
    for(int iteration=0;iteration<options.maximum_iterations;iteration++) {
      if(std::sqrt(preparation_dot(force,force))<=tolerance){result.reason="force-tolerance";break;}
      for(size_t i=0;i<dimension;i++) {
        direction[i]=0;for(size_t j=0;j<dimension;j++)direction[i]+=inverse[i*dimension+j]*force[j];
      }
      if(std::any_of(direction.begin(),direction.end(),[](double v){return !std::isfinite(v);})) {
        result.reason="nonfinite-search-direction";break;
      }
      double slope=-preparation_dot(force,direction),size=0;
      if(slope>=0){result.reason="non-descent-search-direction";break;}
      for(size_t i=0;i<dimension;i+=3) {
        double square=0;for(int axis=0;axis<3;axis++)square+=direction[i+axis]*direction[i+axis];
        size=std::max(size,std::sqrt(square));
      }
      if(size==0){result.reason="zero-search-direction";break;}
      double factor=std::min(1.,options.maximum_atom_displacement/size);
      for(auto& v:direction)v*=factor;slope*=factor;
      double alpha=1;bool accepted=false;PreparationEvaluation candidate;
      for(int search=0;search<options.maximum_line_search_trials;search++) {
        for(size_t i=0;i<dimension;i++)trial[i]=result.positions[i]+alpha*direction[i];
        candidate=value(trial);
        if(candidate.first<=energy+options.armijo*alpha*slope){accepted=true;break;}
        alpha*=.5;
      }
      if(!accepted){result.reason="line-search-stalled";break;}
      for(size_t i=0;i<dimension;i++) {
        displacement[i]=trial[i]-result.positions[i];gradient_change[i]=force[i]-candidate.second[i];
      }
      double curvature=preparation_dot(displacement,gradient_change);
      if(curvature>options.curvature_relative_threshold*std::sqrt(preparation_dot(displacement,displacement))*
                    std::sqrt(preparation_dot(gradient_change,gradient_change))) {
        for(size_t i=0;i<dimension;i++) {
          product[i]=0;for(size_t j=0;j<dimension;j++)product[i]+=inverse[i*dimension+j]*gradient_change[j];
        }
        double coefficient=(1+preparation_dot(gradient_change,product)/curvature)/curvature;
        for(size_t i=0;i<dimension;i++)for(size_t j=0;j<dimension;j++)
          inverse[i*dimension+j]+=coefficient*displacement[i]*displacement[j]-
            (product[i]*displacement[j]+displacement[i]*product[j])/curvature;
      }
      result.positions=trial;energy=candidate.first;force=candidate.second;result.iterations++;
    }
    if(std::sqrt(preparation_dot(force,force))<=tolerance)result.reason="force-tolerance";
  } catch(const PreparationBudgetStop& error) {result.reason=error.what();}
  return result;
}
} // namespace molecular_car
