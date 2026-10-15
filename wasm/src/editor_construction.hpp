// Included inside namespace workshop. Indexed construction mirrors the native
// editor: existing conformers and native cage seeds are rigid boundaries.
using Point = RDGeom::Point3D;
using IndexSet = std::set<unsigned>;
inline Point unit(Point p) { double n=p.length(); return n>1e-10?p/n:Point(1,0,0); }
inline Point perpendicular(Point p) { Point a=std::abs(p.x)<std::abs(p.y)?Point(1,0,0):Point(0,1,0);if(std::abs(p.z)<std::min(std::abs(p.x),std::abs(p.y)))a=Point(0,0,1);return unit(p.crossProduct(a)); }
inline Point rotate(Point p,Point source,Point target) {
  source=unit(source);target=unit(target);double c=std::clamp(source.dotProduct(target),-1.,1.);Point v=source.crossProduct(target);
  if(v.length()<1e-10){if(c>0)return p;auto axis=perpendicular(source);return axis*(2*axis.dotProduct(p))-p;}
  return p+v.crossProduct(p)+v.crossProduct(v.crossProduct(p))/(1+c);
}
inline Point center(const RDKit::RWMol& m,const IndexSet& group){Point c;for(auto i:group)c+=m.getConformer().getAtomPos(i);return c/double(group.size());}
inline double length(const RDKit::RWMol& m,unsigned a,unsigned b){auto pt=RDKit::PeriodicTable::getTable();double n=m.getBondBetweenAtoms(a,b)->getBondTypeAsDouble();return std::max(.7,(pt->getRcovalent(m.getAtomWithIdx(a)->getAtomicNum())+pt->getRcovalent(m.getAtomWithIdx(b)->getAtomicNum()))*(n==3?.8:n==2?.9:n==1.5?.92:1.));}
inline Point outward(const RDKit::RWMol& m,unsigned anchor,const IndexSet& excluded,const std::vector<IndexSet>& cages){
  auto pos=m.getConformer().getAtomPos(anchor);for(auto& g:cages)if(g.count(anchor))return unit(pos-center(m,g));
  std::vector<Point> v;for(auto a:m.atomNeighbors(m.getAtomWithIdx(anchor)))if(!excluded.count(a->getIdx())){auto d=m.getConformer().getAtomPos(a->getIdx())-pos;if(d.length()>.05)v.push_back(unit(d));}
  if(v.empty())return {1,0,0};auto a=m.getAtomWithIdx(anchor);bool planar=a->getIsAromatic()||a->getHybridization()==RDKit::Atom::SP2;
  if(a->getHybridization()==RDKit::Atom::SP)return -v[0];
  if(v.size()==1){double c=planar?-.5:-1./3;return v[0]*c+perpendicular(v[0])*std::sqrt(1-c*c);}
  if(v.size()==2&&!planar){auto bis=unit(v[0]+v[1]),normal=unit(v[0].crossProduct(v[1]));double c=std::min(.99,(1./3)/std::max(.34,bis.dotProduct(v[0])));return unit(bis*(-c)+normal*std::sqrt(1-c*c));}
  Point sum;for(auto q:v)sum+=q;return sum.length()>1e-10?unit(-sum):perpendicular(v[0]);
}
inline std::vector<std::pair<unsigned,unsigned>> boundary(const RDKit::RWMol& m,const IndexSet& g){std::vector<std::pair<unsigned,unsigned>> out;for(auto i:g)for(auto a:m.atomNeighbors(m.getAtomWithIdx(i)))if(!g.count(a->getIdx()))out.emplace_back(i,a->getIdx());return out;}
inline std::vector<IndexSet> components(const RDKit::RWMol& m,IndexSet remaining){std::vector<IndexSet> result;while(!remaining.empty()){IndexSet g;std::vector<unsigned> todo{*remaining.begin()};remaining.erase(todo[0]);for(unsigned k=0;k<todo.size();++k){auto i=todo[k];g.insert(i);for(auto a:m.atomNeighbors(m.getAtomWithIdx(i)))if(remaining.erase(a->getIdx()))todo.push_back(a->getIdx());}result.push_back(g);}return result;}
inline bool embedInitial(RDKit::RWMol& m){
  // Short ordinary attempts, then the upstream-recommended curved-aromatic
  // fallback. Only the planarity prior changes; chirality checks stay enabled.
  for(unsigned attempt=0;attempt<3;++attempt){auto params=RDKit::DGeomHelpers::ETKDGv3;params.randomSeed=20261011;params.numThreads=1;params.timeout=attempt==2?10:3;params.maxIterations=attempt==2?40:8;params.useRandomCoords=attempt>0;params.useBasicKnowledge=attempt<2;
    if(RDKit::DGeomHelpers::EmbedMolecule(m,params)==0)return attempt==2;
  }
  throw std::runtime_error("无法生成三维构象，请检查连接或导入三维文件。");
}
inline bool construct(RDKit::RWMol& m,val p,IndexSet movable,bool first){
  auto ids=indices(p);std::vector<IndexSet> cages,seeds;IndexSet seeded;bool fallback=false;
  if(has(p,"metadata")&&has(p["metadata"],"ketcherGroups"))for(unsigned k=0;k<size(p["metadata"]["ketcherGroups"]);++k){auto r=p["metadata"]["ketcherGroups"][k];if(str(r,"kind")!="c60")continue;IndexSet g;bool added=false;for(unsigned j=0;j<size(r["atomIds"]);++j){auto i=ids.at(r["atomIds"][j].as<std::string>());g.insert(i);added|=movable.count(i)>0;}cages.push_back(g);if(!first&&!added)continue;
    if(g.size()!=60)throw std::runtime_error("C₆₀ 原子映射不完整，请重新添加模板。");
    unsigned bonds=0;auto c=center(m,g);double minRadius=1e9,maxRadius=0;
    for(auto i:g){if(m.getAtomWithIdx(i)->getAtomicNum()!=6)throw std::runtime_error("C₆₀ 模板中的原子已改变，请展开检查结构。");double radius=(m.getConformer().getAtomPos(i)-c).length();minRadius=std::min(minRadius,radius);maxRadius=std::max(maxRadius,radius);}
    for(auto b:m.bonds())if(g.count(b->getBeginAtomIdx())&&g.count(b->getEndAtomIdx())){++bonds;auto d=(m.getConformer().getAtomPos(b->getBeginAtomIdx())-m.getConformer().getAtomPos(b->getEndAtomIdx())).length();if(d<.8||d>2.)throw std::runtime_error("C₆₀ 球壳坐标无效，请重新添加模板。");}
    if(bonds!=90||minRadius<2||maxRadius>5)throw std::runtime_error("C₆₀ 球壳坐标无效，请重新添加模板。");seeds.push_back(g);seeded.insert(g.begin(),g.end());
  }
  if(first&&seeds.empty()){RDKit::MolOps::addHs(m);return embedInitial(m);}
  for(auto& g:seeds){std::vector<std::pair<unsigned,unsigned>> links;for(auto link:boundary(m,g))if(!movable.count(link.second))links.push_back(link);if(links.size()>1)throw std::runtime_error("新增 C₆₀ 有多个固定连接点，请手动调整构象。");if(!links.empty()){for(auto i:g)if(!movable.count(i))throw std::runtime_error("C₆₀ 的固定原子与连接位置冲突。");auto [end,anchor]=links[0];auto dir=outward(m,anchor,movable,cages);auto origin=m.getConformer().getAtomPos(end),radial=origin-center(m,g),target=m.getConformer().getAtomPos(anchor)+dir*length(m,end,anchor);for(auto i:g)m.getConformer().setAtomPos(i,rotate(m.getConformer().getAtomPos(i)-origin,radial,-dir)+target);}}
  IndexSet rest=movable;for(auto i:seeded)rest.erase(i);
  for(auto& group:components(m,rest)){
    auto links=boundary(m,group);if(links.size()>1)throw std::runtime_error("新增片段有多个固定连接点，请手动调整或重建三维构象。");
    int anchor=links.empty()?-1:int(links[0].second);unsigned endpoint=links.empty()?*group.begin():links[0].first;
    if(group.size()==1){if(anchor>=0)m.getConformer().setAtomPos(endpoint,m.getConformer().getAtomPos(anchor)+outward(m,anchor,group,cages)*length(m,endpoint,anchor));continue;}
    IndexSet included=group;if(anchor>=0)included.insert(anchor);std::map<unsigned,unsigned> local;for(auto i:included)local[i]=local.size();
    for(auto b:m.bonds())if((group.count(b->getBeginAtomIdx())||group.count(b->getEndAtomIdx()))&&b->getStereo()>RDKit::Bond::STEREOANY)for(auto i:b->getStereoAtoms())if(!included.count(i))throw std::runtime_error("新增立体键的参考原子跨越固定片段，请手动调整或重建三维构象。");
    RDKit::RWMol fragment(m);for(int i=int(m.getNumAtoms())-1;i>=0;--i)if(!included.count(i))fragment.removeAtom(i);
    if(anchor>=0){RDKit::Atom cap(m.getAtomWithIdx(anchor)->getAtomicNum());fragment.replaceAtom(local.at(anchor),&cap);}
    fragment.clearConformers();RDKit::MolOps::sanitizeMol(fragment);RDKit::MolOps::addHs(fragment);fallback|=embedInitial(fragment);
    auto pivot=fragment.getConformer().getAtomPos(local.at(endpoint));auto originalCenter=center(m,group);Point fragmentCenter;for(auto i:group)fragmentCenter+=fragment.getConformer().getAtomPos(local.at(i));fragmentCenter/=double(group.size());
    Point dir,target,source;if(anchor>=0){dir=outward(m,anchor,group,cages);target=m.getConformer().getAtomPos(anchor)+dir*length(m,endpoint,anchor);source=pivot-fragment.getConformer().getAtomPos(local.at(anchor));}
    std::map<unsigned,Point> points;for(auto i:group){auto q=fragment.getConformer().getAtomPos(local.at(i));points[i]=anchor>=0?rotate(q-pivot,source,dir)+target:q-fragmentCenter+originalCenter;}
    // Pick a rigid roll with the best clearance, without optimizing any atom.
    if(anchor>=0){double best=-1;auto chosen=points;for(unsigned step=0;step<12;++step){double angle=step*3.141592653589793/6,c=std::cos(angle),s=std::sin(angle),score=1e9;auto candidate=points;for(auto& [i,q]:candidate){auto v=q-target;q=target+v*c+dir.crossProduct(v)*s+dir*(dir.dotProduct(v)*(1-c));if(m.getAtomWithIdx(i)->getAtomicNum()==1)continue;for(unsigned j=0;j<m.getNumAtoms();++j)if(!movable.count(j)&&m.getAtomWithIdx(j)->getAtomicNum()!=1&&!(i==endpoint&&j==unsigned(anchor)))score=std::min(score,(q-m.getConformer().getAtomPos(j)).length());}if(score>best){best=score;chosen=candidate;}}points=chosen;}
    for(auto& [i,q]:points)m.getConformer().setAtomPos(i,q);
  }
  m.getConformer().set3D(true);RDKit::MolOps::addHs(m,false,true);return fallback;
}
