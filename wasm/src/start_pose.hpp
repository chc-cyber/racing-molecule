// Counterpart of molecular_car/start_pose.py: rigid, preparation-only placement.
#pragma once
#include "potentials.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>
namespace molecular_car {
using Matrix3=std::array<Vec3,3>;
struct StartPose {std::string mode="preserve";double yaw=0,pitch=0,roll=0;std::vector<int> axis_atoms,validation_ids;};
struct PosePlacement {
  StartPose requested;std::vector<Vec3> coordinates;Matrix3 rotation{};Vec3 variances{};
  bool heading=false,plane=false;std::vector<int> heading_ids;std::string strategy="input-preserved";
};
inline Matrix3 pose_identity(){return {{{1,0,0},{0,1,0},{0,0,1}}};}
inline double pose_dot(const Vec3& a,const Vec3& b){double v=0;for(int k=0;k<3;k++)v+=a[k]*b[k];return v;}
inline Vec3 pose_cross(const Vec3& a,const Vec3& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline double pose_norm(const Vec3& a){return std::sqrt(pose_dot(a,a));}
inline Vec3 pose_unit(Vec3 a){double n=pose_norm(a);for(auto& v:a)v/=n;return a;}
inline void pose_principal(const std::vector<Vec3>& x,const std::vector<int>& reference,Vec3& variances,Matrix3& axes) {
  Vec3 center{};for(int i:reference)for(int k=0;k<3;k++)center[k]+=x[i][k];for(auto& v:center)v/=reference.size();
  Matrix3 covariance{},vectors=pose_identity();
  for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
    for(int atom:reference)covariance[i][j]+=(x[atom][i]-center[i])*(x[atom][j]-center[j]);
    covariance[i][j]/=reference.size();
  }
  for(int iteration=0;iteration<32;iteration++) {
    int p=0,q=1;double largest=std::abs(covariance[0][1]);
    for(auto pair:std::array<std::array<int,2>,2>{{{0,2},{1,2}}})if(std::abs(covariance[pair[0]][pair[1]])>largest) {
      p=pair[0];q=pair[1];largest=std::abs(covariance[p][q]);
    }
    double scale=0;for(int i=0;i<3;i++)scale+=std::abs(covariance[i][i]);scale=std::max(1.,scale);
    if(largest<=1e-14*scale)break;
    double angle=.5*std::atan2(2*covariance[p][q],covariance[q][q]-covariance[p][p]);
    double c=std::cos(angle),s=std::sin(angle),pp=covariance[p][p],qq=covariance[q][q],pq=covariance[p][q];
    for(int k=0;k<3;k++) {
      if(k!=p&&k!=q) {
        double kp=covariance[k][p],kq=covariance[k][q];
        covariance[k][p]=covariance[p][k]=c*kp-s*kq;covariance[k][q]=covariance[q][k]=s*kp+c*kq;
      }
      double vp=vectors[k][p],vq=vectors[k][q];vectors[k][p]=c*vp-s*vq;vectors[k][q]=s*vp+c*vq;
    }
    covariance[p][p]=c*c*pp-2*s*c*pq+s*s*qq;covariance[q][q]=s*s*pp+2*s*c*pq+c*c*qq;
    covariance[p][q]=covariance[q][p]=0;
  }
  std::array<int,3> order{0,1,2};std::stable_sort(order.begin(),order.end(),[&](int a,int b){return covariance[a][a]>covariance[b][b];});
  for(int i=0;i<3;i++){variances[i]=std::max(0.,covariance[order[i]][order[i]]);for(int k=0;k<3;k++)axes[i][k]=vectors[k][order[i]];}
}
inline Matrix3 pose_minimum_rotation(Vec3 a,Vec3 b) {
  a=pose_unit(a);b=pose_unit(b);double cosine=pose_dot(a,b);Matrix3 result=pose_identity();
  if(cosine>1-1e-12)return result;
  if(cosine<-1+1e-12) {
    int index=0;for(int k=1;k<3;k++)if(std::abs(a[k])<std::abs(a[index]))index=k;
    Vec3 seed{};seed[index]=1;Vec3 axis=pose_unit(pose_cross(a,seed));
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)result[i][j]=2*axis[i]*axis[j]-(i==j?1:0);
    return result;
  }
  Vec3 axis=pose_cross(a,b);Matrix3 skew{{{0,-axis[2],axis[1]},{axis[2],0,-axis[0]},{-axis[1],axis[0],0}}};
  for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
    double square=0;for(int k=0;k<3;k++)square+=skew[i][k]*skew[k][j];
    result[i][j]+=skew[i][j]+square/(1+cosine);
  }
  return result;
}
inline PosePlacement place_start_pose(const std::vector<Vec3>& x,const std::vector<double>& masses,
                                     std::vector<int> reference,const StartPose& config) {
  PosePlacement result;result.requested=config;result.rotation=pose_identity();result.coordinates=x;
  for(int id:config.validation_ids)if(id<1||id>static_cast<int>(x.size()))throw std::runtime_error("Heading reference atom is not in this molecule");
  for(int id:config.axis_atoms)if(id<1||id>static_cast<int>(x.size()))throw std::runtime_error("Heading reference atom is not in this molecule");
  if(reference.empty())for(size_t i=0;i<x.size();i++)reference.push_back(i);
  Matrix3 axes{};pose_principal(x,reference,result.variances,axes);
  result.heading=result.variances[0]>1e-12&&result.variances[0]>=1.3*result.variances[1];
  result.plane=result.variances[1]>1e-12&&result.variances[2]<=.25*result.variances[1];
  if(config.mode=="manual") {
    constexpr double radians=3.14159265358979323846/180.;
    double cy=std::cos(config.yaw*radians),sy=std::sin(config.yaw*radians),cp=std::cos(config.pitch*radians),sp=std::sin(config.pitch*radians),
           cr=std::cos(config.roll*radians),sr=std::sin(config.roll*radians);
    result.rotation={{{cy*cp,cy*sp*sr-sy*cr,cy*sp*cr+sy*sr},{sy*cp,sy*sp*sr+cy*cr,sy*sp*cr-cy*sr},{-sp,cp*sr,cp*cr}}};
    result.strategy="manual-file-axes-roll-pitch-yaw";
  } else if(config.mode=="auto") {
    Vec3 forward=axes[0],normal=axes[2];double sign=1;
    for(double v:normal)if(std::abs(v)>1e-12){sign=v;break;}
    if(normal[2]<-1e-12||(std::abs(normal[2])<=1e-12&&sign<0))for(auto& v:normal)v=-v;
    if(!config.axis_atoms.empty()) {
      for(int k=0;k<3;k++)forward[k]=x[config.axis_atoms[1]-1][k]-x[config.axis_atoms[0]-1][k];
      if(pose_norm(forward)<=1e-8)throw std::runtime_error("Heading reference atoms are too close to define a direction");
      forward=pose_unit(forward);result.heading=true;result.heading_ids=config.axis_atoms;
    } else if(result.heading) {
      Vec3 mean{};for(int atom:reference)for(int k=0;k<3;k++)mean[k]+=x[atom][k];for(auto& v:mean)v/=reference.size();
      std::vector<double> projections;
      for(int atom:reference){Vec3 delta{};for(int k=0;k<3;k++)delta[k]=x[atom][k]-mean[k];projections.push_back(pose_dot(delta,forward));}
      double forward_sign=1;for(double v:forward)if(std::abs(v)>1e-12){forward_sign=v;break;}
      if(forward_sign<0){for(auto& v:forward)v=-v;for(auto& v:projections)v=-v;}
      double low=*std::min_element(projections.begin(),projections.end()),high=*std::max_element(projections.begin(),projections.end());
      double tolerance=1e-8*std::max(1.,high-low);size_t rear=0,front=0;
      while(projections[rear]>low+tolerance)rear++;while(projections[front]<high-tolerance)front++;
      result.heading_ids={reference[rear]+1,reference[front]+1};
    }
    if(result.heading&&result.plane&&pose_norm(pose_cross(normal,forward))>1e-8) {
      Vec3 right=pose_unit(pose_cross(normal,forward)),up=pose_cross(forward,right);
      result.rotation={forward,right,up};result.strategy="automatic-plane-and-heading";
    } else if(result.heading){result.rotation=pose_minimum_rotation(forward,{1,0,0});result.strategy="automatic-heading-preserve-roll";}
    else if(result.plane){result.rotation=pose_minimum_rotation(normal,{0,0,1});result.strategy="automatic-plane-preserve-yaw";}
    else result.strategy="ambiguous-geometry-input-preserved";
  }
  Vec3 center{};double total_mass=0;
  for(size_t i=0;i<x.size();i++){total_mass+=masses[i];for(int k=0;k<3;k++)center[k]+=masses[i]*x[i][k];}for(auto& v:center)v/=total_mass;
  if(result.rotation!=pose_identity())for(size_t i=0;i<x.size();i++)for(int axis=0;axis<3;axis++) {
    double value=0;for(int k=0;k<3;k++)value+=(x[i][k]-center[k])*result.rotation[axis][k];
    result.coordinates[i][axis]=value+center[axis];
  }
  double min_z=result.coordinates[0][2];for(const auto& p:result.coordinates)min_z=std::min(min_z,p[2]);
  for(auto& p:result.coordinates){p[0]-=center[0];p[1]-=center[1];p[2]+=2.9-min_z;}
  return result;
}
} // namespace molecular_car
