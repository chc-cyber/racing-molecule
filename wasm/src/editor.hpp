// Browser workshop chemistry. Separate from MD: editing never invokes the race
// engine, substrate, thermostat or its preparation minimizer.
#pragma once
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <map>
#include <set>
namespace workshop {
using emscripten::val;
inline unsigned size(val x){return x["length"].as<unsigned>();}
inline bool has(val x,const char* k){return !x[k].isUndefined()&&!x[k].isNull();}
inline std::string str(val x,const char* k,std::string fallback=""){return has(x,k)?x[k].as<std::string>():fallback;}
inline int integer(val x,const char* k,int fallback=0){return has(x,k)?x[k].as<int>():fallback;}
inline val copy(val x){return val::global("JSON").call<val>("parse",val::global("JSON").call<std::string>("stringify",x));}
inline val point(RDGeom::Point3D p,int dim){val a=val::array();a.call<void>("push",p.x,p.y);if(dim==3)a.call<void>("push",p.z);return a;}
inline std::string uid(const char* prefix){return std::string(prefix)+val::global("crypto").call<std::string>("randomUUID");}
inline RDKit::Bond::BondType bondType(double n){if(n==1)return RDKit::Bond::SINGLE;if(n==2)return RDKit::Bond::DOUBLE;if(n==3)return RDKit::Bond::TRIPLE;if(n==1.5)return RDKit::Bond::AROMATIC;throw std::runtime_error("不支持的化学键类型。");}
inline const std::vector<std::string> chiral={"CHI_UNSPECIFIED","CHI_TETRAHEDRAL_CW","CHI_TETRAHEDRAL_CCW"};
inline const std::vector<std::string> dirs={"NONE","BEGINWEDGE","BEGINDASH","ENDDOWNRIGHT","ENDUPRIGHT","EITHERDOUBLE","UNKNOWN"};
inline const std::vector<std::string> stereo={"STEREONONE","STEREOANY","STEREOZ","STEREOE","STEREOCIS","STEREOTRANS"};
inline unsigned enumeration(std::string value,const std::vector<std::string>& choices){auto it=std::find(choices.begin(),choices.end(),value);if(it==choices.end())throw std::runtime_error("不支持的立体标记。");return it-choices.begin();}
inline std::map<std::string,unsigned> indices(val p){std::map<std::string,unsigned> ids;for(unsigned i=0;i<size(p["atoms"]);++i)ids[p["atoms"][i]["id"].as<std::string>()]=i;return ids;}
inline std::unique_ptr<RDKit::RWMol> fromProject(val p,int dim,bool sanitize){
  auto m=std::make_unique<RDKit::RWMol>();auto ids=indices(p);
  if(!size(p["atoms"])||size(p["atoms"])>2000)throw std::runtime_error("请先绘制分子，编辑器最多支持 2000 个原子。");
  for(unsigned i=0;i<size(p["atoms"]);++i){auto r=p["atoms"][i];auto a=new RDKit::Atom(str(r,"element"));a->setFormalCharge(integer(r,"charge"));a->setIsotope(integer(r,"isotope"));a->setNumExplicitHs(integer(r,"explicitHs"));a->setNoImplicit(has(r,"noImplicit")&&r["noImplicit"].as<bool>());a->setNumRadicalElectrons(integer(r,"radicalElectrons"));a->setChiralTag(static_cast<RDKit::Atom::ChiralType>(enumeration(str(r,"chiralTag","CHI_UNSPECIFIED"),chiral)));m->addAtom(a,false,true);}
  for(unsigned i=0;i<size(p["bonds"]);++i){auto r=p["bonds"][i];m->addBond(ids.at(str(r,"a")),ids.at(str(r,"b")),bondType(r["order"].as<double>()));auto b=m->getBondWithIdx(i);if(b->getBondType()==RDKit::Bond::AROMATIC){b->setIsAromatic(true);b->getBeginAtom()->setIsAromatic(true);b->getEndAtom()->setIsAromatic(true);}b->setBondDir(static_cast<RDKit::Bond::BondDir>(enumeration(str(r,"direction","NONE"),dirs)));}
  m->updatePropertyCache(false);
  for(unsigned i=0;i<size(p["bonds"]);++i){auto r=p["bonds"][i];auto b=m->getBondWithIdx(i);if(has(r,"stereoAtoms")&&size(r["stereoAtoms"])==2)b->setStereoAtoms(ids.at(r["stereoAtoms"][0].as<std::string>()),ids.at(r["stereoAtoms"][1].as<std::string>()));b->setStereo(static_cast<RDKit::Bond::BondStereo>(enumeration(str(r,"stereo","STEREONONE"),stereo)));}
  auto conf=new RDKit::Conformer(m->getNumAtoms());conf->set3D(dim==3);
  for(unsigned i=0;i<m->getNumAtoms();++i){auto v=p["atoms"][i][dim==3?"pos3":"pos2"];double x=v[0].as<double>(),y=v[1].as<double>(),z=dim==3?v[2].as<double>():0;if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z))throw std::runtime_error("原子坐标无效。");conf->setAtomPos(i,{x,y,z});}m->addConformer(conf,true);
  auto name=str(p,"name","Molecule").substr(0,120);std::replace(name.begin(),name.end(),'\n',' ');std::replace(name.begin(),name.end(),'\r',' ');m->setProp("_Name",name);
  if(has(p,"metadata")&&has(p["metadata"],"ketcherGroups"))for(unsigned n=0;n<size(p["metadata"]["ketcherGroups"]);++n){auto r=p["metadata"]["ketcherGroups"][n];if(str(r,"type")!="SUP")throw std::runtime_error("只支持普通 SUP 缩写基团。");RDKit::SubstanceGroup g(m.get(),"SUP");g.setProp("LABEL",str(r,"name","Group"));if(has(r,"expanded")&&r["expanded"].as<bool>())g.setProp("ESTATE",std::string("E"));std::set<unsigned> members;for(unsigned j=0;j<size(r["atomIds"]);++j){auto index=ids.at(r["atomIds"][j].as<std::string>());g.addAtomWithIdx(index);members.insert(index);}for(auto b:m->bonds())if(members.count(b->getBeginAtomIdx())!=members.count(b->getEndAtomIdx()))g.addBondWithIdx(b->getIdx());if(has(r,"attachmentPoints"))for(unsigned j=0;j<size(r["attachmentPoints"]);++j){auto a=r["attachmentPoints"][j];g.addAttachPoint(ids.at(str(a,"atomId")),has(a,"leavingAtomId")?static_cast<int>(ids.at(str(a,"leavingAtomId"))):-1,str(a,"attachmentId","1"));}RDKit::addSubstanceGroup(*m,g);}
  if(sanitize){RDKit::MolOps::sanitizeMol(*m);RDKit::MolOps::assignChiralTypesFromBondDirs(*m,-1,false);RDKit::MolOps::assignStereochemistry(*m,false,true);}
  return m;
}
inline val atomRecord(const RDKit::Atom& a,std::string id,RDGeom::Point3D p2,RDGeom::Point3D p3){val r=val::object();r.set("id",id);r.set("element",a.getSymbol());r.set("charge",a.getFormalCharge());r.set("isotope",a.getIsotope());r.set("explicitHs",a.getNumExplicitHs());r.set("noImplicit",a.getNoImplicit());r.set("radicalElectrons",a.getNumRadicalElectrons());if(a.getChiralTag()>2)throw std::runtime_error("当前工坊不支持该立体原子类型。");r.set("chiralTag",chiral.at(a.getChiralTag()));r.set("pos2",point(p2,2));r.set("pos3",point(p3,3));return r;}
inline val bondRecord(const RDKit::Bond& b,val atoms){val r=val::object();r.set("id",uid("b-"));r.set("a",atoms[b.getBeginAtomIdx()]["id"]);r.set("b",atoms[b.getEndAtomIdx()]["id"]);r.set("order",b.getBondTypeAsDouble());r.set("direction",dirs.at(b.getBondDir()));r.set("stereo",stereo.at(b.getStereo()));val refs=val::array();for(int i:b.getStereoAtoms())refs.call<void>("push",atoms[i]["id"]);r.set("stereoAtoms",refs);return r;}
inline void layout(val p,RDKit::RWMol& mol){
  RDKit::RWMol drawing(mol);RDDepict::compute2DCoords(drawing);RDKit::WedgeMolBonds(drawing,&drawing.getConformer());
  for(unsigned i=0;i<mol.getNumAtoms();++i)p["atoms"][i].set("pos2",point(drawing.getConformer().getAtomPos(i),2));
  for(unsigned i=0;i<mol.getNumBonds();++i){auto b=drawing.getBondWithIdx(i);p["bonds"][i].set("a",p["atoms"][b->getBeginAtomIdx()]["id"]);p["bonds"][i].set("b",p["atoms"][b->getEndAtomIdx()]["id"]);p["bonds"][i].set("direction",dirs.at(b->getBondDir()));}
}
inline val importFile(val d){
  std::string text=str(d,"content"),file=str(d,"filename","molecule.mol");if(text.empty()||text.size()>2000000)throw std::runtime_error("文件为空或超过 2 MB。");
  auto end=text.find("$$$$");if(end!=std::string::npos){if(text.substr(end+4).find_first_not_of(" \r\n\t")!=std::string::npos)throw std::runtime_error("请导入单分子 SDF 文件。");text.resize(end);}
  std::regex terminator("(^|\n)M\\s+END\\s*(\n|$)");if(std::distance(std::sregex_iterator(text.begin(),text.end(),terminator),std::sregex_iterator())!=1)throw std::runtime_error("分子文件必须恰好包含一份完整结构。");
  auto m=std::unique_ptr<RDKit::RWMol>(RDKit::MolBlockToMol(text,false,false,true));if(!m||!m->getNumAtoms()||m->getNumAtoms()>2000)throw std::runtime_error("无法读取分子，或超过 2000 个原子。");
  if(!m->getStereoGroups().empty())throw std::runtime_error("当前工坊不支持增强立体组。");
  for(auto a:m->atoms())if(a->hasQuery())throw std::runtime_error("请使用具体元素，查询原子不受支持。");for(auto b:m->bonds()){if(b->hasQuery())throw std::runtime_error("查询键不受支持。");bondType(b->getBondTypeAsDouble());}
  val warnings=val::array();try{RDKit::RWMol clean(*m);RDKit::MolOps::sanitizeMol(clean);RDKit::MolOps::assignStereochemistry(clean,true,true);for(unsigned i=0;i<m->getNumAtoms();++i)m->getAtomWithIdx(i)->setChiralTag(clean.getAtomWithIdx(i)->getChiralTag());for(unsigned i=0;i<m->getNumBonds();++i){auto b=clean.getBondWithIdx(i);if(b->getStereoAtoms().size()==2)m->getBondWithIdx(i)->setStereoAtoms(b->getStereoAtoms()[0],b->getStereoAtoms()[1]);m->getBondWithIdx(i)->setStereo(b->getStereo());}}catch(const std::exception& e){warnings.call<void>("push",std::string("请检查结构价态：")+e.what());}
  bool has3=m->getNumConformers()&&m->getConformer().is3D();val p=val::object(),atoms=val::array(),bonds=val::array();std::string name=file;if(m->hasProp("_Name"))m->getProp("_Name",name);p.set("format",std::string("molecular-car-workshop-v1"));p.set("name",name);p.set("has3d",has3);p.set("fragments",val::array());
  for(auto a:m->atoms()){auto pos=m->getConformer().getAtomPos(a->getIdx());atoms.call<void>("push",atomRecord(*a,uid("a-"),pos,pos));}for(auto b:m->bonds())bonds.call<void>("push",bondRecord(*b,atoms));p.set("atoms",atoms);p.set("bonds",bonds);
  val meta=val::object(),groups=val::array();for(auto& g:RDKit::getSubstanceGroups(*m)){std::string type;g.getProp("TYPE",type);if(type!="SUP")throw std::runtime_error("当前仅支持 SUP 缩写基团。");val r=val::object(),members=val::array(),points=val::array();r.set("id",uid("g-"));r.set("type",type);std::string label="Group",estate;g.getPropIfPresent("LABEL",label);g.getPropIfPresent("ESTATE",estate);r.set("name",label);r.set("expanded",estate=="E");for(auto i:g.getAtoms())members.call<void>("push",atoms[i]["id"]);r.set("atomIds",members);for(auto& a:g.getAttachPoints()){val q=val::object();q.set("atomId",atoms[a.aIdx]["id"]);q.set("attachmentId",a.id);if(a.lvIdx>=0)q.set("leavingAtomId",atoms[a.lvIdx]["id"]);points.call<void>("push",q);}r.set("attachmentPoints",points);groups.call<void>("push",r);}meta.set("ketcherGroups",groups);meta.set("sourceFilename",file);p.set("metadata",meta);if(has3)layout(p,*m);
  val result=val::object();result.set("project",p);result.set("warnings",warnings);return result;
}
inline val exportFile(val d){auto p=d["project"];int dim=integer(d,"dimension",3);if(dim!=2&&dim!=3)throw std::runtime_error("导出维度无效。");if(dim==3&&!p["has3d"].as<bool>())throw std::runtime_error("请先生成三维构象。");auto m=fromProject(p,dim,false);auto format=str(d,"format","mol");if(format!="mol"&&format!="sdf")throw std::runtime_error("不支持的文件格式。");auto text=RDKit::MolToMolBlock(*m,true,-1,false,true);if(format=="sdf")text+="\n$$$$\n";val r=val::object();r.set("content",text);r.set("filename",str(p,"name","molecule")+"."+format);r.set("warnings",val::array());return r;}
inline void appendHydrogens(val p,RDKit::RWMol& work,RDKit::RWMol& drawing,unsigned count){
  auto atoms=p["atoms"],bonds=p["bonds"];unsigned oldBonds=size(bonds);std::map<std::string,std::string> parents;
  for(unsigned i=count;i<work.getNumAtoms();++i){auto a=work.getAtomWithIdx(i);auto r=atomRecord(*a,uid("a-"),drawing.getConformer().getAtomPos(i),work.getConformer().getAtomPos(i));r.set("generatedHydrogen",true);auto parent=*work.atomNeighbors(a).begin();parents[str(r,"id")]=str(atoms[parent->getIdx()],"id");atoms.call<void>("push",r);}
  for(unsigned i=0;i<count;++i){atoms[i].set("explicitHs",work.getAtomWithIdx(i)->getNumExplicitHs());atoms[i].set("chiralTag",chiral.at(work.getAtomWithIdx(i)->getChiralTag()));}
  for(unsigned i=oldBonds;i<work.getNumBonds();++i)bonds.call<void>("push",bondRecord(*work.getBondWithIdx(i),atoms));
  for(unsigned i=0;i<size(p["fragments"]);++i){auto members=p["fragments"][i]["atomIds"];std::set<std::string> old;for(unsigned j=0;j<size(members);++j)old.insert(members[j].as<std::string>());for(auto& [h,parent]:parents)if(old.count(parent))members.call<void>("push",h);}
}
inline val geometry(val d){
  val p=copy(d["project"]),warnings=val::array();auto op=str(d,"operation");bool has3=p["has3d"].as<bool>();auto m=fromProject(p,has3?3:2,op!="layout2d");unsigned count=m->getNumAtoms();
  if(count>600&&op!="layout2d")throw std::runtime_error("几何操作最多支持 600 个原子。");
  if(op=="layout2d")layout(p,*m);
  else {
    if(op!="prepare"&&op!="rebuild"&&op!="generate"&&op!="hydrogenate"&&op!="relax")throw std::runtime_error("未知几何操作。");
    if(op=="relax"&&!has3)throw std::runtime_error("请先生成三维构象。");
    auto ids=indices(p);std::set<unsigned> fixed;if(has(d,"fixedIds"))for(unsigned i=0;i<size(d["fixedIds"]);++i)fixed.insert(ids.at(d["fixedIds"][i].as<std::string>()));if(!has3&&!fixed.empty())throw std::runtime_error("固定原子需要已有三维坐标。");
    bool embed=op=="rebuild"||op=="generate"||(op=="prepare"&&!has3);bool pending=false;std::set<unsigned> added;
    if(op=="prepare"&&has3&&has(p,"metadata")&&has(p["metadata"],"geometryPending")){auto pendingData=p["metadata"]["geometryPending"];if(has(pendingData,"atomIds"))for(unsigned i=0;i<size(pendingData["atomIds"]);++i)added.insert(ids.at(pendingData["atomIds"][i].as<std::string>()));pending=!added.empty();}
    if(pending){embed=true;for(unsigned i=0;i<count;++i)if(!added.count(i))fixed.insert(i);}
    RDKit::RWMol drawing(*fromProject(p,2,true));RDKit::MolOps::addHs(drawing,false,true);RDKit::MolOps::addHs(*m,false,!embed);
    if(m->getNumAtoms()>2000)throw std::runtime_error("补氢后超过 2000 个原子。请减少结构大小。");
    if(embed){auto params=RDKit::DGeomHelpers::ETKDGv3;params.randomSeed=20261011;params.numThreads=1;params.timeout=20;params.maxIterations=200;std::map<int,RDGeom::Point3D> coords;for(auto i:fixed)coords[i]=m->getConformer().getAtomPos(i);if(!coords.empty()){params.coordMap=&coords;params.useRandomCoords=true;}int code=RDKit::DGeomHelpers::EmbedMolecule(*m,params);if(code!=0&&coords.empty()){params.useRandomCoords=true;code=RDKit::DGeomHelpers::EmbedMolecule(*m,params);}if(code!=0)throw std::runtime_error("无法生成三维构象，请检查结构或导入三维文件。");for(auto& [i,pos]:coords)m->getConformer().setAtomPos(i,pos);}
    if(op=="relax"||op=="generate"){if(!RDKit::UFF::getAtomTypes(*m).second)throw std::runtime_error("此结构缺少 UFF 参数。");std::unique_ptr<ForceFields::ForceField> ff(RDKit::UFF::constructForceField(*m,10.0,-1,true));for(auto i:fixed)ff->fixedPoints().push_back(i);ff->initialize();if(ff->minimize(1000,.001,1e-7))warnings.call<void>("push",std::string("UFF 尚未收敛，可继续优化。"));}
    if(op!="hydrogenate")for(unsigned i=0;i<count;++i){auto pos=m->getConformer().getAtomPos(i);if(!std::isfinite(pos.x)||!std::isfinite(pos.y)||!std::isfinite(pos.z))throw std::runtime_error("计算产生无效坐标。");p["atoms"][i].set("pos3",point(pos,3));}
    if(op!="relax"&&op!="rebuild")appendHydrogens(p,*m,drawing,count);
    if(op!="hydrogenate")p.set("has3d",true);
    if(has(p,"metadata")){p["metadata"].delete_("geometryPending");p["metadata"].delete_("geometryNotice");}
  }
  val r=val::object();r.set("project",p);r.set("revision",d["revision"]);r.set("warnings",warnings);r.set("message",std::string(op=="layout2d"?"二维排版已整理。":op=="relax"?"UFF 优化完成。":op=="hydrogenate"?"氢原子已补齐。":"三维构象已准备。"));return r;
}
inline val request(std::string op,val data){if(op=="import")return importFile(data);if(op=="export")return exportFile(data);if(op=="geometry")return geometry(data);throw std::runtime_error("未知工坊操作。");}
}
