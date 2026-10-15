// One compiled MD module. No Python or JavaScript is called during force loops.
#include "potentials.hpp"
#include "preparation_bfgs.hpp"
#include "start_pose.hpp"
#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <library.h>
#include <ForceField/ForceField.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include "editor.hpp"

using emscripten::val;
using namespace molecular_car;
namespace {
template<typename T> val array(const std::vector<T>& v) {
  return val::global("Array").call<val>("from",val(emscripten::typed_memory_view(v.size(),v.data())));
}
val vector3(const Vec3& v) { return array(std::vector<double>(v.begin(),v.end())); }
std::vector<Vec3> positions(val v) {
  auto flat=emscripten::convertJSArrayToNumberVector<double>(v);
  if(flat.size()%3) throw std::runtime_error("Coordinates must contain complete xyz triples");
  std::vector<Vec3> x(flat.size()/3);
  for(size_t i=0;i<x.size();i++) for(int a=0;a<3;a++) {
    if(!std::isfinite(flat[3*i+a])) throw std::runtime_error("Coordinates must be finite");
    x[i][a]=flat[3*i+a];
  }
  return x;
}
std::vector<double> flatten(const std::vector<Vec3>& x) {
  std::vector<double> v;v.reserve(3*x.size());for(const auto& p:x) for(double c:p) v.push_back(c);return v;
}
double number(val v,const char* name) { return v[name].as<double>(); }
std::string format(double v) { std::ostringstream s;s<<std::setprecision(17)<<v;return s.str(); }
double norm(const std::vector<Vec3>& f) {double n=0;for(const auto& v:f) for(double c:v)n+=c*c;return std::sqrt(n);}
}

class BrowserEngine {
  val model;
  Fields field;
  Substrate gold;
  NormalTopology topology;
  std::unique_ptr<RDKit::RWMol> molecule;
  std::unique_ptr<ForceFields::ForceField> uff;
  std::vector<Response> atoms;
  std::vector<double> masses;
  std::vector<Vec3> x, total_force, internal_force, surface_force, field_force;
  void* lmp=nullptr;
  std::string filename,name,surface_mode="periodic-au",callback_error;
  double temperature=600,damping=350,corrugation=.08,strength=.75,growth=1.,lattice_scale=1.;
  double force_tolerance=1e-4,wall_limit=30;
  double timestep_fs=.5;
  double race_length=200.;
  int evaluation_limit=25000;
  int mode=1,seed=98245165,sequence=0;
  bool drive_enabled=false,thermostat=true,preparing=false,verifying=false;
  bool isotope_ignored=false;
  int evaluations=0,verification_evaluations=0;
  std::chrono::steady_clock::time_point prepare_started;
  std::string budget_stop;
  Vec3 initial_center{};
  val preparation=val::object();
  PosePlacement pose_placement;

  StartPose parse_start_pose(val value) const {
    StartPose config;config.mode="auto";
    if(value.isNull()||value.isUndefined())return config;
    if(value.typeOf().as<std::string>()!="object"||val::global("Array").call<bool>("isArray",value))
      throw std::runtime_error("Invalid start placement configuration");
    auto keys=val::global("Object").call<val>("keys",value);
    for(unsigned i=0;i<keys["length"].as<unsigned>();i++) {
      auto key=keys[i].as<std::string>();
      if(key!="mode"&&key!="yaw_deg"&&key!="pitch_deg"&&key!="roll_deg"&&key!="axis_atoms")throw std::runtime_error("Unknown start placement parameter");
    }
    if(!value["mode"].isUndefined())config.mode=value["mode"].as<std::string>();
    if(config.mode!="auto"&&config.mode!="manual"&&config.mode!="preserve")throw std::runtime_error("Unknown start placement mode");
    auto angle=[&](const char* key) {
      if(value[key].isUndefined())return 0.;
      if(value[key].typeOf().as<std::string>()!="number")throw std::runtime_error("Start angles must be numeric");
      double result=value[key].as<double>();
      if(!std::isfinite(result)||result<-180||result>180)throw std::runtime_error("Start angles must be between -180 and 180 degrees");return result;
    };
    config.yaw=angle("yaw_deg");config.pitch=angle("pitch_deg");config.roll=angle("roll_deg");
    if(!value["axis_atoms"].isUndefined()) {
      auto ids=value["axis_atoms"];
      if(!val::global("Array").call<bool>("isArray",ids))throw std::runtime_error("Heading reference must be an array");
      unsigned count=ids["length"].as<unsigned>();if(count!=0&&count!=2)throw std::runtime_error("Heading reference needs two atoms");
      for(unsigned i=0;i<count;i++) {
        if(ids[i].typeOf().as<std::string>()!="number")throw std::runtime_error("Heading atoms need integer IDs");
        double id=ids[i].as<double>();if(!std::isfinite(id)||std::floor(id)!=id||id<1||id>256)throw std::runtime_error("Invalid heading atom ID");
        config.axis_atoms.push_back(static_cast<int>(id));
      }
      if(count==2&&config.axis_atoms[0]==config.axis_atoms[1])throw std::runtime_error("Heading reference atoms must differ");
    }
    config.validation_ids=config.axis_atoms;
    if(config.mode!="manual")config.yaw=config.pitch=config.roll=0;
    if(config.mode!="auto")config.axis_atoms.clear();return config;
  }

  val start_pose_info() const {
    const auto& p=pose_placement;val report=val::object(),requested=val::object();
    requested.set("mode",p.requested.mode);requested.set("yaw_deg",p.requested.yaw);
    requested.set("pitch_deg",p.requested.pitch);requested.set("roll_deg",p.requested.roll);requested.set("axis_atoms",array(p.requested.axis_atoms));
    report.set("protocol_version",std::string("rigid-start-pose-v1"));report.set("requested",requested);report.set("strategy",p.strategy);
    val matrix=val::array();for(const auto& row:p.rotation)matrix.call<void>("push",vector3(row));
    report.set("rotation_matrix",matrix);report.set("principal_variances_A2",vector3(p.variances));
    report.set("heading_confident",p.heading);report.set("plane_confident",p.plane);report.set("heading_atom_ids",array(p.heading_ids));
    report.set("qualification",std::string("Rigid pre-relaxation placement only; unconstrained surface relaxation and subsequent dynamics may change orientation."));
    val after=val::object();after.set("heading_defined",!p.heading_ids.empty());
    if(!p.heading_ids.empty()) {
      Vec3 direction{};for(int k=0;k<3;k++)direction[k]=x[p.heading_ids[1]-1][k]-x[p.heading_ids[0]-1][k];
      constexpr double degrees=180./3.14159265358979323846;
      after.set("heading_atom_ids",array(p.heading_ids));after.set("yaw_deg",degrees*std::atan2(direction[1],direction[0]));
      after.set("tilt_deg",degrees*std::atan2(direction[2],std::hypot(direction[0],direction[1])));
    }
    if(lmp)report.set("after_relaxation",after);return report;
  }

  void command(const std::string& s) {
    lammps_command(lmp,s.c_str());
    if(lammps_has_error(lmp)) {
      char message[2048];lammps_get_last_error_message(lmp,message,sizeof(message));
      throw std::runtime_error(message);
    }
    if(!callback_error.empty()) throw std::runtime_error(callback_error);
  }
  void gather(const char* property,std::vector<Vec3>& into) const {
    into.resize(atoms.size());lammps_gather_atoms(lmp,property,1,3,into.data());
  }
  Vec3 center_of(const std::vector<Vec3>& values) const {
    Vec3 c{};double mass=0;
    for(size_t i=0;i<masses.size();i++) {mass+=masses[i];for(int a=0;a<3;a++)c[a]+=masses[i]*values[i][a];}
    for(auto& a:c)a/=mass;return c;
  }
  double evaluate(std::vector<Vec3>& coords,int field_mode,double field_strength,bool continuum,
                  double& internal_energy,double& surface_energy,double& field_energy) {
    auto flat=flatten(coords);internal_energy=uff->calcEnergy(flat.data());
    std::vector<double> gradient(flat.size(),0);uff->calcGrad(flat.data(),gradient.data());
    internal_force.assign(coords.size(),{0,0,0});
    for(size_t i=0;i<coords.size();i++) for(int a=0;a<3;a++)internal_force[i][a]=-gradient[3*i+a];
    field_energy=molecular_car::fields(coords,atoms,field,field_mode,field_strength,Vec3{0,0,0},corrugation,continuum,field_force,growth);
    surface_force.assign(coords.size(),{0,0,0});surface_energy=0;
    if(!continuum)surface_energy=molecular_car::substrate(coords,atoms,gold,topology,surface_force,lattice_scale);
    total_force=internal_force;
    for(size_t i=0;i<coords.size();i++) for(int a=0;a<3;a++)total_force[i][a]+=surface_force[i][a]+field_force[i][a];
    double energy=internal_energy+surface_energy+field_energy;
    if(!std::isfinite(energy)||!std::isfinite(norm(total_force)))throw std::runtime_error("Nonfinite molecular energy or forces");
    return energy;
  }
  static void force_callback(void* context,int64_t step,int nlocal,int* ids,double** local_x,double** f) {
    auto& e=*static_cast<BrowserEngine*>(context);
    try {
      if(nlocal!=static_cast<int>(e.atoms.size()))throw std::runtime_error("Browser MD requires one serial LAMMPS instance");
      std::vector<Vec3> coords(e.atoms.size());std::vector<bool> seen(e.atoms.size(),false);
      for(int i=0;i<nlocal;i++) {
        int index=ids[i]-1;
        if(index<0||index>=nlocal||seen[index])throw std::runtime_error("LAMMPS atom identities are invalid");
        seen[index]=true;for(int a=0;a<3;a++)coords[index][a]=local_x[i][a];
      }
      double internal,surface,field,strength=e.drive_enabled?e.strength:0;
      double energy=e.evaluate(coords,e.mode,strength,e.surface_mode=="continuum",internal,surface,field);
      for(int i=0;i<nlocal;i++)for(int a=0;a<3;a++)f[i][a]=e.total_force[ids[i]-1][a];
      lammps_fix_external_set_energy_global(e.lmp,"molecular",energy);
      if(e.preparing) {
        if(e.verifying)e.verification_evaluations++;
        else if(++e.evaluations>=e.evaluation_limit&&e.budget_stop.empty()) {
          e.budget_stop="force-evaluation-budget";lammps_force_timeout(e.lmp);
        }
      }
    } catch(const std::exception& error) {
      e.callback_error=error.what();for(int i=0;i<nlocal;i++)for(int a=0;a<3;a++)f[i][a]=0;
      lammps_force_timeout(e.lmp);
    }
  }
  void normal_topology() {
    size_t n=atoms.size();topology.neighbours.assign(n,{-1,-1,-1});topology.counts.assign(n,0);
    std::vector<std::vector<std::pair<int,double>>> graph(n);
    for(const auto* b:molecule->bonds()) {
      int a=b->getBeginAtomIdx(),c=b->getEndAtomIdx();double order=b->getBondTypeAsDouble();
      graph[a].push_back({c,order});graph[c].push_back({a,order});
    }
    for(size_t i=0;i<n;i++) {
      auto& entries=graph[i];std::sort(entries.begin(),entries.end());int degree=entries.size();
      if(atoms[i].code<0)continue;
      if((atoms[i].code==1&&(degree>1||(degree==1&&entries[0].second!=1)))||degree>4)
        throw std::runtime_error("Invalid C/H normal topology");
      if(degree<=1||degree==4)continue;
      double orders=0;bool triple=false;
      for(const auto& b:entries){orders+=b.second;triple|=b.second==3;}
      if(triple||(degree==2&&orders>=4))continue;
      int pivot=degree==3?entries[2].first:i;Vec3 a,b,c;
      for(int axis=0;axis<3;axis++){a[axis]=x[entries[0].first][axis]-x[pivot][axis];b[axis]=x[entries[1].first][axis]-x[pivot][axis];}
      c={a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
      double aa=0,bb=0,cc=0;for(int axis=0;axis<3;axis++){aa+=a[axis]*a[axis];bb+=b[axis]*b[axis];cc+=c[axis]*c[axis];}
      double scale=std::sqrt(aa)*std::sqrt(bb);
      if(scale==0||std::sqrt(cc)<=1e-5*scale)continue;
      topology.counts[i]=degree;for(int k=0;k<degree;k++)topology.neighbours[i][k]=entries[k].first;
    }
  }
  void chemistry(std::string content,const std::string& file,const StartPose& pose=StartPose{}) {
    isotope_ignored=false;
    if(content.size()>1000000||content.empty())throw std::runtime_error("Molecule file is empty or too large");
    auto dot=file.find_last_of('.');
    if(dot==std::string::npos)throw std::runtime_error("Import a MOL or single-record SDF file");
    std::string suffix=file.substr(dot);
    std::transform(suffix.begin(),suffix.end(),suffix.begin(),[](unsigned char c){return std::tolower(c);});
    if(suffix!=".mol"&&suffix!=".sdf")throw std::runtime_error("Import a MOL or single-record SDF file");
    if(content.compare(0,3,"\xEF\xBB\xBF")==0)content.erase(0,3);
    content=std::regex_replace(content,std::regex("\r\n"),"\n");
    std::regex end("(^|\n)M\\s+END\\s*(\n|$)");
    if(std::distance(std::sregex_iterator(content.begin(),content.end(),end),std::sregex_iterator())!=1)
      throw std::runtime_error("MOL/SDF must contain exactly one complete molecule");
    auto delimiter=content.find("$$$$");
    if(delimiter!=std::string::npos) {
      if(content.substr(delimiter+4).find_first_not_of(" \r\n\t")!=std::string::npos)
        throw std::runtime_error("Multi-record SDF files are not supported");
      content.resize(delimiter);
    }
    molecule.reset(RDKit::MolBlockToMol(content,false,false,true));
    if(!molecule)throw std::runtime_error("Cannot parse molecule");
    RDKit::MolOps::sanitizeMol(*molecule);
    if(molecule->getNumAtoms()==0||molecule->getNumAtoms()>256)throw std::runtime_error("Molecule atom limit is 256 including hydrogen");
    std::vector<int> fragments;
    if(RDKit::MolOps::getMolFrags(*molecule,fragments)!=1)throw std::runtime_error("The car must be one connected molecule");
    for(auto* a:molecule->atoms()) {
      if(a->getNumRadicalElectrons())throw std::runtime_error("Open-shell molecules are not supported");
      if(model["effective_polarizabilities"][a->getSymbol()].isUndefined())throw std::runtime_error("Element has no game response parameters");
      isotope_ignored|=a->getIsotope()!=0;a->setIsotope(0);
    }
    for(auto* b:molecule->bonds()) {
      auto type=b->getBondType();
      if(type!=RDKit::Bond::SINGLE&&type!=RDKit::Bond::DOUBLE&&type!=RDKit::Bond::TRIPLE&&type!=RDKit::Bond::AROMATIC)
        throw std::runtime_error("Only single, double, triple and aromatic bonds are supported");
    }
    bool has3d=molecule->getNumConformers()&&molecule->getConformer().is3D();
    RDKit::MolOps::addHs(*molecule,false,has3d);
    if(molecule->getNumAtoms()>256)throw std::runtime_error("Molecule exceeds 256 atoms after hydrogen completion");
    auto types=RDKit::UFF::getAtomTypes(*molecule);
    if(!types.second)throw std::runtime_error("UFF lacks parameters for this molecule");
    for(auto* a:molecule->atoms())if(a->getSymbol()=="P"&&a->getDegree()==3)
      throw std::runtime_error("Three-coordinate phosphorus is unsupported by the current UFF validation policy");
    if(!has3d) {
      auto params=RDKit::DGeomHelpers::ETKDGv3;params.randomSeed=20261011;params.numThreads=1;
      int result=RDKit::DGeomHelpers::EmbedMolecule(*molecule,params);
      if(result!=0){params.useRandomCoords=true;result=RDKit::DGeomHelpers::EmbedMolecule(*molecule,params);}
      if(result!=0)throw std::runtime_error("Cannot create a 3D conformer; import a molecule with 3D coordinates");
    }
    uff.reset(RDKit::UFF::constructForceField(*molecule,100.0,-1,false));uff->initialize();
    if(uff->minimize(300,.001,1e-7)!=0)throw std::runtime_error("Initial UFF relaxation did not converge in 300 iterations");
    name=file;
    if(molecule->hasProp("_Name"))molecule->getProp("_Name",name);
    auto& conf=molecule->getConformer();x.resize(molecule->getNumAtoms());atoms.resize(x.size());masses.resize(x.size());
    std::vector<double> charges(x.size());RDKit::computeGasteigerCharges(*molecule,charges,12,true);
    int formal=0;double charge_sum=0;
    for(size_t i=0;i<x.size();i++) {
      auto* a=molecule->getAtomWithIdx(i);auto p=conf.getAtomPos(i);x[i]={p.x,p.y,p.z};masses[i]=a->getMass();
      if(!std::isfinite(charges[i])||!std::isfinite(masses[i])||masses[i]<=0)throw std::runtime_error("Invalid charge or mass");
      formal+=a->getFormalCharge();charge_sum+=charges[i];
      int code=a->getSymbol()=="C"?0:a->getSymbol()=="H"?1:-1;
      ForceFields::UFF::UFFVdW vd;
      if(!RDKit::UFF::getUFFVdWParams(*molecule,i,i,vd))throw std::runtime_error("Missing UFF vdW atom-type parameters");
      auto hp=model["surface"]["hetero_au"];
      double rmin=std::sqrt(vd.x_ij*number(hp,"gold_rmin_angstrom"));
      double well=std::sqrt(vd.D_ij*number(hp,"gold_well_depth_kcal_mol"));
      atoms[i]={charges[i],number(model["effective_polarizabilities"],a->getSymbol().c_str()),code,
                well*std::pow(rmin,12)/23.0605478306,-2*well*std::pow(rmin,6)/23.0605478306};
    }
    if(std::abs(charge_sum-formal)>1e-5+1e-5*std::abs(formal))throw std::runtime_error("Partial charge sum does not match formal charge");
    std::vector<int> reference;for(size_t i=0;i<x.size();i++)if(molecule->getAtomWithIdx(i)->getAtomicNum()!=1)reference.push_back(i);
    pose_placement=place_start_pose(x,masses,reference,pose);x=pose_placement.coordinates;
    for(size_t i=0;i<x.size();i++) {
      conf.setAtomPos(i,RDGeom::Point3D(x[i][0],x[i][1],x[i][2]));
    }
    for(int axis=0;axis<3;axis++) {
      double lo=x[0][axis],hi=lo;for(const auto& p:x){lo=std::min(lo,p[axis]);hi=std::max(hi,p[axis]);}
      if(hi-lo>45)throw std::runtime_error("Molecular span exceeds 45 angstrom");
    }
    uff.reset(RDKit::UFF::constructForceField(*molecule,100.0,-1,false));uff->initialize();
    bool hetero=std::any_of(atoms.begin(),atoms.end(),[](const Response& a){return a.code<0;});
    if(hetero) {
      auto flat=flatten(x);double energy=uff->calcEnergy(flat.data());std::vector<double> gradient(flat.size(),0);
      uff->calcGrad(flat.data(),gradient.data());
      if(!std::isfinite(energy)||std::abs(energy)>1e8)throw std::runtime_error("Pathological initial UFF energy");
      for(size_t i=0;i<flat.size();i++) {
        if(!std::isfinite(gradient[i]))throw std::runtime_error("Invalid UFF gradient");
        double saved=flat[i];flat[i]=saved+1e-5;double plus=uff->calcEnergy(flat.data());
        flat[i]=saved-1e-5;double minus=uff->calcEnergy(flat.data());flat[i]=saved;
        double numerical=(plus-minus)/2e-5;
        if(!std::isfinite(numerical)||std::abs(gradient[i]-numerical)>1e-3+1e-4*std::max(std::abs(gradient[i]),std::abs(numerical)))
          throw std::runtime_error("UFF energy/gradient admission check failed for this chemical environment");
      }
      for(size_t i=0;i<x.size();i++) {
        double n=0;for(int a=0;a<3;a++)n+=gradient[3*i+a]*gradient[3*i+a];
        if(std::sqrt(n)>1e5)throw std::runtime_error("Pathological initial UFF atom force");
      }
    }
    normal_topology();
  }
  double elapsed() const { return std::chrono::duration<double>(std::chrono::steady_clock::now()-prepare_started).count(); }
  void relax() {
    preparing=true;verifying=false;evaluations=verification_evaluations=0;budget_stop.clear();
    prepare_started=std::chrono::steady_clock::now();
    val stages=val::array();
    try {
      command("timer timeout off every 10");command("run 0");
      auto diagnostics=[&](){gather("f",total_force);return norm(total_force);};
      double initial=diagnostics(),residual=initial,maximum=0;
      for(const auto& f:total_force){double s=0;for(auto c:f)s+=c*c;maximum=std::max(maximum,std::sqrt(s));}
      if(maximum>1e5||std::abs(lammps_get_thermo(lmp,"pe"))>1e8)throw std::runtime_error("Pathological initial joint energy or force");
      auto protocol=model["preparation"]["stages"];
      for(unsigned index=0;index<protocol["length"].as<unsigned>();index++) {
        auto stage=protocol[index];auto style=stage[0].as<std::string>();
        if(!budget_stop.empty()||elapsed()>=wall_limit){if(budget_stop.empty())budget_stop="wall-time-budget";break;}
        if(residual<=force_tolerance)break;
        command("min_style "+style);
        if(style!="hftn")command("min_modify dmax 0.1");
        int remaining=std::max(1,evaluation_limit-evaluations);
        command("timer timeout "+format(std::max(1e-12,wall_limit-elapsed()))+" every 1");
        command("minimize 0.0 "+format(force_tolerance)+" "+std::to_string(stage[1].as<int>())+" "+std::to_string(std::min(stage[2].as<int>(),remaining)));
        command("timer timeout off every 10");verifying=true;command("run 0");verifying=false;
        residual=diagnostics();if(elapsed()>=wall_limit&&budget_stop.empty())budget_stop="wall-time-budget";
        val record=val::object();record.set("minimizer",style);record.set("global_force_norm_kcal_mol_A",residual);
        stages.call<void>("push",record);
      }
      if(residual>force_tolerance&&budget_stop.empty()) {
        if(elapsed()>=wall_limit)budget_stop="wall-time-budget";
        else if(evaluation_limit-evaluations>static_cast<int>(6*atoms.size()+1)) {
          BfgsOptions options;auto config=model["preparation"]["bfgs_rescue_options"];
          options.maximum_iterations=config["maximum_iterations"].as<int>();
          options.maximum_line_search_trials=config["maximum_line_search_trials"].as<int>();
          options.maximum_cholesky_attempts=config["maximum_cholesky_attempts"].as<int>();
          options.hessian_step=number(config,"hessian_step_angstrom");
          options.maximum_atom_displacement=number(config,"maximum_atom_displacement_angstrom");
          options.armijo=number(config,"armijo_coefficient");
          options.curvature_relative_threshold=number(config,"curvature_relative_threshold");
          options.diagonal_shift_relative_start=number(config,"diagonal_shift_relative_start");
          gather("x",x);auto checkpoint=flatten(x);int before=evaluations+verification_evaluations;
          auto stage_started=std::chrono::steady_clock::now();BfgsResult rescued;
          auto scatter=[&](const std::vector<double>& coords) {
            lammps_scatter_atoms(lmp,"x",1,3,const_cast<double*>(coords.data()));
          };
          auto finish=[&]() {
            command("timer timeout off every 10");scatter(checkpoint);verifying=true;
            try {command("run 0");}catch(...) {verifying=false;throw;}
            verifying=false;
          };
          try {
            auto objective=[&](const std::vector<double>& coords)->PreparationEvaluation {
              if(elapsed()>=wall_limit&&budget_stop.empty())budget_stop="wall-time-budget";
              if(!budget_stop.empty())throw PreparationBudgetStop(budget_stop);
              scatter(coords);command("run 0");
              if(!budget_stop.empty())throw PreparationBudgetStop(budget_stop);
              gather("f",total_force);return {lammps_get_thermo(lmp,"pe"),flatten(total_force)};
            };
            rescued=dense_bfgs(checkpoint,objective,force_tolerance,options);checkpoint=rescued.positions;
          } catch(...) {finish();throw;}
          finish();residual=diagnostics();
          if(elapsed()>=wall_limit&&budget_stop.empty())budget_stop="wall-time-budget";
          val record=val::object();record.set("minimizer",std::string("curvature-initialized-bfgs"));
          record.set("iterations",rescued.iterations);record.set("global_force_norm_kcal_mol_A",residual);
          record.set("force_callback_evaluations",evaluations+verification_evaluations-before);
          record.set("hessian_diagonal_shift",rescued.diagonal_shift);
          record.set("wall_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-stage_started).count());
          record.set("final_coordinate_force_verified",true);
          record.set("termination_reason",budget_stop.empty()?rescued.reason:budget_stop);
          stages.call<void>("push",record);
        }
      }
      if(elapsed()>=wall_limit&&budget_stop.empty())budget_stop="wall-time-budget";
      preparation=val::object();preparation.set("protocol_version",model["preparation"]["protocol_version"]);
      preparation.set("initial_force_norm_kcal_mol_A",initial);
      preparation.set("global_force_norm_kcal_mol_A",residual);preparation.set("global_force_tolerance_kcal_mol_A",force_tolerance);
      preparation.set("force_callback_evaluations",evaluations+verification_evaluations);
      preparation.set("wall_seconds",elapsed());preparation.set("stages",stages);
      preparation.set("accepted",residual<=force_tolerance&&budget_stop.empty());
      preparation.set("final_coordinate_force_verified",true);
      if(residual>force_tolerance||!budget_stop.empty())throw std::runtime_error("Joint relaxation rejected: residual force "+format(residual)+"; "+budget_stop);
    } catch(...) {
      preparing=false;verifying=false;lammps_command(lmp,"timer timeout off every 10");throw;
    }
    preparing=false;command("timer timeout off every 10");
  }
  void bath() {command("fix bath all langevin "+format(temperature)+" "+format(temperature)+" "+format(damping)+" "+std::to_string(seed)+" tally yes");}

 public:
  void install_fourier(val meta) {
    const auto schema=meta["schema"].as<std::string>();
    const bool atlas=schema=="fourier-saip-browser-v2";
    if((schema!="fourier-saip-browser-v1"&&!atlas) ||
       meta["parameter_sha256"].as<std::string>()!=model["surface"]["parameter_sha256"].as<std::string>())
      throw std::runtime_error("Fourier asset version or parameters mismatch");
    auto& p=gold.fourier;
    int nz=meta["nz"].as<int>(),na=meta["na"].as<int>(),modes=meta["modes"]["length"].as<int>();
    if(nz!=(atlas?95:111)||na!=(atlas?13:21)||modes!=12)throw std::runtime_error("Unsupported Fourier asset shape");
    if(atlas) {
      if(meta["charts"].as<int>()!=3)throw std::runtime_error("Unsupported Fourier atlas chart count");
      for(int chart=0;chart<3;chart++)for(int axis=0;axis<3;axis++)
        if(meta["chart_axes"][chart][axis].as<int>()!=(chart+axis)%3)throw std::runtime_error("Unsupported Fourier chart axes");
    }
    for(const auto& item:std::vector<std::pair<const char*,double>>{{"astep",.2},{"bound",atlas?1.2:2.},{"blend",atlas?.2:.4},{"zstep",atlas?.1:.05},{"zstart",atlas?2.1:2.65},{"zmin",atlas?2.4:2.8},{"zmax",atlas?11.2:8.},{"zblend",.3}})
      if(std::abs(number(meta,item.first)-item.second)>1e-12||!std::isfinite(number(meta,item.first)))throw std::runtime_error("Unsupported Fourier support");
    std::vector<std::array<int,2>> expected={{0,0},{0,2},{0,4},{1,-3},{1,-1},{1,1},{1,3},{2,-4},{2,-2},{2,0},{2,2},{2,4}};
    for(int k=0;k<modes;k++)if(meta["modes"][k][0].as<int>()!=expected[k][0]||meta["modes"][k][1].as<int>()!=expected[k][1])throw std::runtime_error("Fourier harmonic mismatch");
    p.enabled=false;p.loaded=0;p.nz=nz;p.na=na;p.modes=modes;p.indices=expected;p.charts=atlas?3:1;
    p.astep=.2;p.bound=atlas?1.2:2.;p.blend=atlas?.2:.4;p.zstep=atlas?.1:.05;
    p.zstart=atlas?2.1:2.65;p.zmin=atlas?2.4:2.8;p.zmax=atlas?11.2:8.;p.zblend=.3;
    p.coefficients.resize(static_cast<size_t>(p.charts)*nz*na*na*4*modes*2);
  }
  void load_fourier_species(int code,val values) {
    auto& p=gold.fourier;size_t size=static_cast<size_t>(p.charts)*p.nz*p.na*p.na*4*p.modes*2;
    if(code<0||code>1||!size||values["length"].as<size_t>()!=size)throw std::runtime_error("Invalid Fourier coefficient payload");
    if(p.coefficients.size()<(code+1)*size)p.coefficients.resize((code+1)*size);
    val(emscripten::typed_memory_view(size,p.coefficients.data()+code*size)).call<void>("set",values);
    for(size_t i=code*size;i<(code+1)*size;i++)if(!std::isfinite(p.coefficients[i]))throw std::runtime_error("Nonfinite Fourier coefficient");
    p.loaded|=1<<code;
  }
  void fourier_enabled(bool enabled) {
    if(enabled&&!(gold.fourier.loaded&1))throw std::runtime_error("Fourier carbon coefficients are not loaded");
    gold.fourier.enabled=enabled;
  }
  val proxy_counts() const {return array(std::vector<int>(gold.fourier.counts.begin(),gold.fourier.counts.end()));}
  BrowserEngine(val metadata,val tile,val saip,val hetero):model(metadata) {
    race_length=number(model["track"],"finish_x")-number(model["track"],"start_x");
    if(!std::isfinite(race_length)||race_length<=0)throw std::runtime_error("Invalid race length");
    force_tolerance=number(model["preparation"],"force_tolerance");
    wall_limit=number(model["preparation"],"wall_time_limit_seconds");
    evaluation_limit=model["preparation"]["evaluation_limit"].as<int>();
    auto p=model["field_parameters"];
    field={number(p,"ADSORPTION_DEPTH"),number(p,"ADSORPTION_HEIGHT"),number(p,"ADSORPTION_ALPHA"),
      number(p,"SURFACE_SPACING"),number(p,"CORRUGATION_DECAY"),number(p,"TIP_SOFTENING"),number(p,"FIELD_SCALE"),
      number(p,"GRADIENT_BASE_FIELD"),number(p,"GRADIENT_FIELD_SLOPE"),number(p,"GRADIENT_POLARIZATION_GAIN"),
      number(p,"BOUNDARY_WIDTH"),number(p,"BOUNDARY_STIFFNESS"),number(p,"TRACK_X_MIN"),number(p,"TRACK_X_MAX"),number(p,"LANE_HALF_WIDTH")};
    gold.tile=positions(tile);auto s=model["surface"];gold.epsilon=number(s,"normal_area_epsilon");
    for(int c=0;c<2;c++)for(int pair=0;pair<2;pair++)for(int k=0;k<8;k++)gold.parameters[c][pair][k]=s["parameters"][c][pair][k].as<double>();
    auto initialize_grid=[&](Grid& g,val v,const char* key) {
      auto shape=s["arrays"][key]["shape"];g.nz=shape[1].as<int>();g.nx=shape[2].as<int>();g.ny=shape[3].as<int>();
      g.px=s["period"][0].as<double>();g.py=s["period"][1].as<double>();g.dz=number(s,"grid_dz");g.zstart=number(s,"z_start");
      g.values=emscripten::convertJSArrayToNumberVector<double>(v);
      if(g.values.size()!=static_cast<size_t>(2*g.nz*g.nx*g.ny))throw std::runtime_error("Invalid browser potential table size");
    };
    initialize_grid(gold.saip,saip,"saip-grid");initialize_grid(gold.hetero,hetero,"hetero-grid");
  }
  ~BrowserEngine() {close();}
  void close() {if(lmp){lammps_close(lmp);lmp=nullptr;}uff.reset();molecule.reset();}
  val prepare_impl(std::string content,std::string file,val settings,const StartPose& pose) {
    close();callback_error.clear();filename=file;drive_enabled=false;sequence=0;
    thermostat=true;seed=model["seed"].as<int>();mode=1;
    strength=settings["drive_strength"].isUndefined()?number(model["drive"],"strength"):number(settings,"drive_strength");
    if(!std::isfinite(strength)||strength<0||strength>1)throw std::runtime_error("Invalid drive strength");
    temperature=number(settings,"temperature");damping=number(settings,"damping_fs");corrugation=number(settings,"corrugation");
    timestep_fs=settings["timestep_fs"].isUndefined()?.5:number(settings,"timestep_fs");
    validate_timestep(timestep_fs);
    growth=settings["field_growth"].isUndefined()?1.:number(settings,"field_growth");
    if(!std::isfinite(growth)||growth<.1||growth>2)throw std::runtime_error("Invalid field growth");
    lattice_scale=settings["lattice_scale"].isUndefined()?1.:number(settings,"lattice_scale");
    if(lattice_scale<0||lattice_scale>12)throw std::runtime_error("Lattice scale must be in [0, 12]");
    surface_mode=settings["surface"].as<std::string>();
    if(surface_mode!="periodic-au"&&surface_mode!="continuum"&&surface_mode!="fourier-au")throw std::runtime_error("Unknown browser surface");
    fourier_enabled(surface_mode=="fourier-au");
    chemistry(content,file,pose);
    char a0[]="molecular-car",a1[]="-log",a2[]="none",a3[]="-screen",a4[]="none",a5[]="-nocite";
    char* argv[]={a0,a1,a2,a3,a4,a5};lmp=lammps_open_no_mpi(6,argv,nullptr);
    if(lammps_version(lmp)!=20260930)throw std::runtime_error("Unexpected LAMMPS browser build version");
    command("units real");command("atom_style atomic");command("boundary f f f");command("atom_modify map array sort 0 0.0");
    command("region world block -110 1110 -350 350 -40 120 units box");
    std::vector<double> type_masses;std::vector<int> ids(x.size()),types(x.size());
    for(size_t i=0;i<x.size();i++) {
      auto it=std::find(type_masses.begin(),type_masses.end(),masses[i]);
      if(it==type_masses.end()){type_masses.push_back(masses[i]);types[i]=type_masses.size();}
      else types[i]=1+std::distance(type_masses.begin(),it);ids[i]=i+1;
    }
    command("create_box "+std::to_string(type_masses.size())+" world");
    for(size_t i=0;i<type_masses.size();i++)command("mass "+std::to_string(i+1)+" "+format(type_masses[i]));
    auto flat=flatten(x);
    if(lammps_create_atoms(lmp,x.size(),ids.data(),types.data(),flat.data(),nullptr,nullptr,0)!=static_cast<int>(x.size()))
      throw std::runtime_error("Failed to create all race atoms");
    command("pair_style zero 10.0");command("pair_coeff * *");command("fix molecular all external pf/callback 1 1");
    // FIRE relaxation keeps its existing integration settings for every race dt.
    command("compute internal_temperature all temp/com");command("thermo_style custom step temp pe ke etotal");command("thermo 1000000");command("timestep 0.5");
    lammps_set_fix_external_callback(lmp,"molecular",force_callback,this);relax();command("reset_timestep 0");
    command("timestep "+format(timestep_fs));
    if(temperature>0)command("velocity all create "+format(temperature)+" 17329 mom yes rot yes dist gaussian");
    else command("velocity all set 0 0 0");
    command("fix integrate all nve");if(thermostat)bath();drive_enabled=true;command("run 0");
    gather("x",x);initial_center=center_of(x);
    return info();
  }
  val prepare(std::string content,std::string file,val settings) {return prepare_impl(content,file,settings,StartPose{});}
  val prepare_with_pose(std::string content,std::string file,val settings,val pose) {return prepare_impl(content,file,settings,parse_start_pose(pose));}
  // Validation/benchmark controls, deliberately absent from game controls/UI.
  static void validate_timestep(double value) {
    if(value!=.25&&value!=.5&&value!=1.&&value!=2.)
      throw std::runtime_error("Diagnostic timestep must be 0.25, 0.5, 1 or 2 fs");
  }
  void set_timestep(double value) {
    validate_timestep(value);
    if(!lmp||lammps_get_thermo(lmp,"step")!=0)
      throw std::runtime_error("Set diagnostic timestep only at prepared step zero");
    timestep_fs=value;command("timestep "+format(value));
  }
  void surface_optimization(int level) {
    if(level<0||level>3)throw std::runtime_error("Surface optimization level must be 0..3");
    gold.optimization=level;
  }
  val force_timing(val coordinates,int repetitions) {
    if(repetitions<1||repetitions>10000)throw std::runtime_error("Invalid benchmark repeat count");
    auto p=positions(coordinates);if(p.size()!=atoms.size())throw std::runtime_error("Atom count mismatch");
    auto flat=flatten(p);std::vector<double> gradient(flat.size());std::vector<Vec3> f;
    double checksum=0;val result=val::object();
    auto started=std::chrono::steady_clock::now();
    for(int k=0;k<repetitions;k++) {
      checksum+=uff->calcEnergy(flat.data());std::fill(gradient.begin(),gradient.end(),0.);
      uff->calcGrad(flat.data(),gradient.data());checksum+=gradient[0];
    }
    auto now=std::chrono::steady_clock::now();
    result.set("uff_ms",std::chrono::duration<double,std::milli>(now-started).count());started=now;
    for(int k=0;k<repetitions;k++){checksum+=molecular_car::substrate(p,atoms,gold,topology,f,lattice_scale);checksum+=f[0][0];}
    now=std::chrono::steady_clock::now();
    result.set("surface_ms",std::chrono::duration<double,std::milli>(now-started).count());started=now;
    for(int k=0;k<repetitions;k++){checksum+=molecular_car::fields(p,atoms,field,mode,strength,Vec3{0,0,0},corrugation,false,f,growth);checksum+=f[0][0];}
    now=std::chrono::steady_clock::now();
    result.set("field_ms",std::chrono::duration<double,std::milli>(now-started).count());
    result.set("checksum",checksum);result.set("repetitions",repetitions);return result;
  }
  // Independent static evaluations never advance MD time or alter coordinates.
  val energy_forces(val coordinates,int field_mode,double field_strength,bool continuum) {
    auto p=positions(coordinates);if(p.size()!=atoms.size())throw std::runtime_error("Atom count mismatch");
    double ie,se,fe,e=evaluate(p,field_mode,field_strength,continuum,ie,se,fe);
    val result=val::object();result.set("energy",e);result.set("forces",array(flatten(total_force)));
    result.set("uff_energy",ie);result.set("uff_forces",array(flatten(internal_force)));
    result.set("surface_energy",se);result.set("surface_forces",array(flatten(surface_force)));
    result.set("field_energy",fe);result.set("field_forces",array(flatten(field_force)));
    return result;
  }
  // Test preparation exercises native RDKit without paying for minimization.
  val prepare_chemistry(std::string content,std::string file) {close();filename=file;chemistry(content,file);return info();}
  val prepare_chemistry_with_pose(std::string content,std::string file,val pose) {close();filename=file;chemistry(content,file,parse_start_pose(pose));return info();}
  val info() const {
    val result=val::object(),a=val::array(),b=val::array();
    for(size_t i=0;i<atoms.size();i++) {
      val item=val::object();item.set("id",i+1);item.set("element",molecule->getAtomWithIdx(i)->getSymbol());
      item.set("mass",masses[i]);item.set("charge",atoms[i].charge);item.set("polarizability",atoms[i].alpha);a.call<void>("push",item);
    }
    for(const auto* bond:molecule->bonds()) {
      val item=val::object();item.set("a",bond->getBeginAtomIdx()+1);item.set("b",bond->getEndAtomIdx()+1);
      item.set("order",bond->getBondTypeAsDouble());b.call<void>("push",item);
    }
    val mol=val::object();mol.set("name",name);mol.set("atoms",a);mol.set("bonds",b);
    val notes=val::array();
    if(isotope_ignored)notes.call<void>("push",std::string("比赛采用元素常规质量，忽略输入文件中的同位素标记。"));
    if(std::any_of(atoms.begin(),atoms.end(),[](const Response& a){return a.code<0;}))
      notes.call<void>("push",std::string("C/H 使用原 SAIP 游戏适配；其他元素使用尚未标定真实吸附或金属配位的范德华混合近似。"));
    mol.set("notes",notes);result.set("molecule",mol);
    val runtime=val::object();runtime.set("lammps_version",20260930);runtime.set("rdkit_version","2025.09.3");
    runtime.set("execution","browser-wasm");runtime.set("timestep_fs",timestep_fs);runtime.set("langevin_seed",seed);
    runtime.set("velocity_seed",17329);runtime.set("forcefield","RDKit C++ UFF / LAMMPS WASM integration");
    runtime.set("mass_policy","standard-element-masses-isotopes-ignored-v1");runtime.set("preparation",preparation);
    runtime.set("start_pose",start_pose_info());
    val surface=val::global("JSON").call<val>("parse",val::global("JSON").call<std::string>("stringify",model["surface"]));
    if(surface_mode=="continuum") {
      surface=val::object();surface.set("mode","continuum");surface.set("lattice_spacing_angstrom",field.spacing);
    } else if(std::none_of(atoms.begin(),atoms.end(),[](const Response& a){return a.code<0;}))
      surface.delete_("hetero_au");
    surface.set("lattice_scale",lattice_scale);surface.set("corrugation_control_version","six-phase-atomic-crystal-contrast-v3");
    surface.set("implementation",std::string("fused-phase-aligned-cache-site-stencil-v1"));
    surface.set("optimization_level",gold.optimization);
    surface.set("fourier_enabled",gold.fourier.enabled);
    if(gold.fourier.enabled) {
      surface.set("mode","fourier-au");surface.set("model","bounded Fourier/Hermite proxy with conservative exact fallback");
      surface.set("implementation",gold.fourier.charts==3?"fourier-saip-browser-v2-three-chart-atlas":"fourier-saip-browser-v1-six-phase-coefficients");
      surface.set("proxy_height_min_A",gold.fourier.zmin);surface.set("proxy_height_max_A",gold.fourier.zmax);
      surface.set("proxy_chart_bound",gold.fourier.bound);surface.set("proxy_charts",gold.fourier.charts);surface.set("fourier_modes",gold.fourier.modes);
    }
    val race=val::object();race.set("rule","longitudinal-displacement-v1");
    race.set("start_x",initial_center[0]);race.set("finish_x",initial_center[0]+race_length);
    race.set("length_angstrom",race_length);runtime.set("race",race);
    runtime.set("surface",surface);result.set("runtime",runtime);
    result.set("positions",array(flatten(x)));result.set("normal_counts",array(topology.counts));
    std::vector<int> neighbours;for(const auto& v:topology.neighbours)for(auto k:v)neighbours.push_back(k);
    result.set("normal_neighbours",array(neighbours));return result;
  }
  val controls() const {
    val d=val::object();d.set("mode","gradient");d.set("strength",strength);d.set("growth",growth);
    val c=val::object();c.set("drive",d);return c;
  }
  val frame() {
    gather("x",x);std::vector<Vec3> velocity;gather("v",velocity);
    for(const auto& p:x)for(double v:p)if(!std::isfinite(v))throw std::runtime_error("Nonfinite MD coordinate");
    for(const auto& p:velocity)for(double v:p)if(!std::isfinite(v))throw std::runtime_error("Nonfinite MD velocity");
    for(const auto& p:x)if(p[0]<=-110||p[0]>=1110||p[1]<=-350||p[1]>=350||p[2]<-5||p[2]>100)
      throw std::runtime_error("Car left the effective simulation region");
    Vec3 c=center_of(x),v=center_of(velocity);double speed=0;for(double a:v)speed+=a*a;
    int step=lammps_get_thermo(lmp,"step");val f=val::object();f.set("seq",sequence);f.set("step",step);f.set("time_fs",step*timestep_fs);f.set("timestep_fs",timestep_fs);
    f.set("positions",array(flatten(x)));f.set("center",vector3(c));f.set("speed_angstrom_ps",std::sqrt(speed)*1000);
    f.set("energy",lammps_get_thermo(lmp,"etotal"));
    auto temp=static_cast<double*>(lammps_extract_compute(lmp,"internal_temperature",0,0));f.set("internal_temperature",*temp);
    auto control=controls();f.set("drive",control["drive"]);f.set("surface_mode",surface_mode);
    f.set("fourier_counts",proxy_counts());
    val effective=val::object();effective.set("mode","gradient");effective.set("growth",growth);effective.set("strength",drive_enabled?strength:0);
    f.set("effective_drive",effective);
    double distance=c[0]-initial_center[0];
    f.set("progress",std::clamp(distance/race_length,0.0,1.0));f.set("distance_angstrom",distance);
    f.set("race_start_x",initial_center[0]);f.set("race_length_angstrom",race_length);return f;
  }
  val advance(int steps) {
    if(steps<1||steps>1000)throw std::runtime_error("Invalid MD batch length");
    command("run "+std::to_string(steps)+" post no");sequence++;return frame();
  }
  val refresh() {command("run 0");sequence++;return frame();}
  val phase_space() const {
    std::vector<Vec3> coordinates,velocity;gather("x",coordinates);gather("v",velocity);
    val result=val::object();result.set("positions",array(flatten(coordinates)));result.set("velocities",array(flatten(velocity)));
    return result;
  }
  void set_phase_space(val coordinates,val velocities) {
    auto p=positions(coordinates),v=positions(velocities);
    if(p.size()!=atoms.size()||v.size()!=atoms.size())throw std::runtime_error("Atom count mismatch in phase-space fixture");
    auto flat_x=flatten(p),flat_v=flatten(v);
    lammps_scatter_atoms(lmp,"x",1,3,flat_x.data());lammps_scatter_atoms(lmp,"v",1,3,flat_v.data());
    command("reset_timestep 0");command("run 0");sequence=0;initial_center=center_of(p);
  }
  void update(val changes) {
    auto present=[&](const char* key){return !changes[key].isUndefined();};
    if(present("drive")) {
      auto d=changes["drive"];if(!d["mode"].isUndefined()){if(d["mode"].as<std::string>()!="gradient")throw std::runtime_error("Only the fixed gradient field is supported");mode=1;}
      if(!d["strength"].isUndefined())strength=number(d,"strength");
      if(!d["growth"].isUndefined())growth=number(d,"growth");
    }
    if(present("tip")||present("autopilot"))throw std::runtime_error("Tip controls have been removed");
    if(!std::isfinite(strength)||strength<0||strength>1)throw std::runtime_error("Strength must be between zero and one");
    if(!std::isfinite(growth)||growth<.1||growth>2)throw std::runtime_error("Field growth must be in [0.1, 2]");
    bool reset_bath=present("temperature")||present("damping_fs");
    if(present("temperature"))temperature=number(changes,"temperature");
    if(present("damping_fs"))damping=number(changes,"damping_fs");
    if(present("corrugation"))corrugation=number(changes,"corrugation");
    if(present("lattice_scale"))lattice_scale=number(changes,"lattice_scale");
    if(!std::isfinite(lattice_scale)||lattice_scale<0||lattice_scale>12)throw std::runtime_error("Lattice scale must be in [0, 12]");
    if(reset_bath&&thermostat){command("unfix bath");bath();}
  }
  void conservative(bool enabled) {
    // Diagnostic only: reproduce paired NVE checks using identical x/v.
    drive_enabled=!enabled;
    if(enabled&&thermostat){command("unfix bath");thermostat=false;}
    else if(!enabled&&!thermostat){thermostat=true;bath();}
    command("run 0");
  }
  void nve_with_drive() {
    // Isolated validation only: fixed conservative field, no heat bath.
    if(thermostat){command("unfix bath");thermostat=false;}
    drive_enabled=true;command("run 0");
  }
};

EMSCRIPTEN_BINDINGS(molecular_car_browser) {
  emscripten::function("workshopRequest", &workshop::request);
  emscripten::class_<BrowserEngine>("BrowserEngine")
    .constructor<val,val,val,val>()
    .function("prepare",&BrowserEngine::prepare)
    .function("prepareWithPose",&BrowserEngine::prepare_with_pose)
    .function("prepareChemistry",&BrowserEngine::prepare_chemistry)
    .function("prepareChemistryWithPose",&BrowserEngine::prepare_chemistry_with_pose)
    .function("energyForces",&BrowserEngine::energy_forces)
    .function("surfaceOptimization",&BrowserEngine::surface_optimization)
    .function("forceTiming",&BrowserEngine::force_timing)
    .function("installFourier",&BrowserEngine::install_fourier)
    .function("loadFourierSpecies",&BrowserEngine::load_fourier_species)
    .function("fourierEnabled",&BrowserEngine::fourier_enabled)
    .function("proxyCounts",&BrowserEngine::proxy_counts)
    .function("info",&BrowserEngine::info)
    .function("frame",&BrowserEngine::frame)
    .function("advance",&BrowserEngine::advance)
    .function("refresh",&BrowserEngine::refresh)
    .function("phaseSpace",&BrowserEngine::phase_space)
    .function("setPhaseSpace",&BrowserEngine::set_phase_space)
    .function("setTimestep",&BrowserEngine::set_timestep)
    .function("controls",&BrowserEngine::controls)
    .function("update",&BrowserEngine::update)
    .function("nveWithDrive",&BrowserEngine::nve_with_drive)
    .function("conservative",&BrowserEngine::conservative)
    .function("close",&BrowserEngine::close);
}
