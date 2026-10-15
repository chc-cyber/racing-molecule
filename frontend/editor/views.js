import * as THREE from 'three';
import {OrbitControls} from '/vendor/OrbitControls.js';
import {TransformControls} from './vendor/TransformControls.js';
import * as Model from './model.js';
import {snapBondEndpoint, pointInPolygon, segmentDistance, atomLabel} from './view-interactions.js';

const colors={C:0x536b5d,H:0xd7ddd2,N:0x467ec2,O:0xcf5e4a,F:0x8da547,Cl:0x8da547,Br:0xa17b48,I:0x906dac,S:0xcfb948,P:0xc49159,B:0xbd9b87,Si:0xc8b387};
const NS='http://www.w3.org/2000/svg';
const svg=(tag,attrs={},text)=>{const el=document.createElementNS(NS,tag);for(const [key,value] of Object.entries(attrs))el.setAttribute(key,value);if(text!==undefined)el.textContent=text;return el;};
const hitKey=hit=>`${hit?.atomId||''}:${hit?.bondId||''}`;
const hasHit=hit=>!!(hit?.atomId||hit?.bondId);

/** Transient previews never write the canonical project. Only completed gestures call the model. */
export class WorkshopViews {
  constructor(canvas2,host3,callbacks){
    this.svg=canvas2;this.host=host3;this.cb=callbacks;this.dimension=2;this.scale=48;this.center=[0,0];
    this.project={atoms:[],bonds:[],fragments:[],has3d:false};this.selected=new Set();
    this.drag=null;this.gesture=null;this.previewProject=null;this.fragmentPreview=null;this.fragmentVisual=null;this.lastPointer=null;this._hoverHit=null;this._numbers=false;
    this.renderer=new THREE.WebGLRenderer({antialias:true,alpha:true});this.renderer.setPixelRatio(Math.min(devicePixelRatio,2));host3.append(this.renderer.domElement);
    this.labels=document.createElement('div');this.labels.className='labels-3d';host3.append(this.labels);
    this.overlay=svg('svg',{class:'viewport-overlay','aria-hidden':'true'});Object.assign(this.overlay.style,{position:'absolute',inset:'0',width:'100%',height:'100%',pointerEvents:'none'});host3.append(this.overlay);
    this.scene=new THREE.Scene();this.camera=new THREE.PerspectiveCamera(38,1,.05,2000);this.camera.up.set(0,0,1);this.camera.position.set(12,-15,16);
    this.scene.add(new THREE.HemisphereLight(0xffffff,0x7c8d76,2.5));const light=new THREE.DirectionalLight(0xffffff,2);light.position.set(7,-8,12);this.scene.add(light);
    this.group=new THREE.Group();this.scene.add(this.group);this.atomMeshes=new Map();this.bondMeshes=new Map();this.materials=new Map();this.meshes=[];
    this.sphereGeometry=new THREE.SphereGeometry(1,18,12);this.bondGeometry=new THREE.CylinderGeometry(.055,.055,1,8);this.ray=new THREE.Raycaster();
    this.controls=new OrbitControls(this.camera,this.renderer.domElement);
    // Left belongs to construction; right/middle/wheel navigate in every editing mode.
    this.controls.mouseButtons.LEFT=null;this.controls.mouseButtons.MIDDLE=THREE.MOUSE.PAN;this.controls.mouseButtons.RIGHT=THREE.MOUSE.ROTATE;this.controls.enableDamping=false;
    this.controls.addEventListener('change',()=>{this.repaint3();this.cb.viewChange?.();});
    this.pivot=new THREE.Object3D();this.scene.add(this.pivot);this.transform=new TransformControls(this.camera,this.renderer.domElement);this.transform.setSize(.85);this.transform.setSpace('world');this.scene.add(this.transform.getHelper());
    this.bindTransform();this.bind2();this.bind3();this.observer=new ResizeObserver(()=>this.resize());this.observer.observe(host3.parentElement);this.resize();
  }
  get hoverHit(){return this._hoverHit;}
  display(){
    const options=this.cb.display?.()||{};
    const anchors=(options.anchors||[]).map((a,i)=>typeof a==='string'?{atomId:a,label:String.fromCharCode(65+i)}:{...a,atomId:a.atomId||a.id,label:a.label||String.fromCharCode(65+i)});
    return {...options,numbers:options.numbers??this._numbers,anchors,movingIds:new Set(options.movingIds||[]),fixedIds:new Set(options.fixedIds||[])};
  }
  setNumbers(value){this._numbers=!!value;this.paint();}
  activeProject(){return this.previewProject||this.display().previewProject||this.fragmentVisual?.project||this.project;}
  effectiveSelection(){return new Set([...this.selected,...(this.cb.moveIds?.()||[])]);}
  stateFor(id,display=this.display(),selected=this.effectiveSelection()){
    if(this.fragmentVisual?.atomIds?.includes(id))return this.fragmentVisual.error?'invalid':'ghost';
    if(display.fixedIds.has(id))return 'fixed';if(display.movingIds.has(id))return 'moving';if(selected.has(id))return 'selected';if(this._hoverHit?.atomId===id)return 'hover';return 'normal';
  }
  setDimension(d){this.cancelGesture(false);this.setHover(null);this.dimension=d;this.svg.toggleAttribute('hidden',d!==2);this.host.hidden=d!==3;this.resize();this.updateMode();}
  resize(){const r=this.host.parentElement.getBoundingClientRect();this.width=Math.max(1,r.width);this.height=Math.max(1,r.height);this.svg.setAttribute('viewBox',`0 0 ${this.width} ${this.height}`);this.overlay.setAttribute('viewBox',`0 0 ${this.width} ${this.height}`);this.renderer.setSize(this.width,this.height);this.camera.aspect=this.width/this.height;this.camera.updateProjectionMatrix();this.paint();}
  reset(){
    if(this.project.atoms.length){const a=this.project.atoms;this.center=[0,1].map(k=>(Math.min(...a.map(x=>x.pos2[k]))+Math.max(...a.map(x=>x.pos2[k])))/2);const spread=[0,1].map(k=>Math.max(...a.map(x=>x.pos2[k]))-Math.min(...a.map(x=>x.pos2[k])));this.scale=Math.max(8,Math.min(55,(this.width-100)/Math.max(3,spread[0]),(this.height-100)/Math.max(3,spread[1])));const box=new THREE.Box3();a.forEach(x=>box.expandByPoint(new THREE.Vector3(...x.pos3)));const c=box.getCenter(new THREE.Vector3()),size=box.getSize(new THREE.Vector3()).length();this.controls.target.copy(c);this.camera.position.copy(c).add(new THREE.Vector3(.7,-1,1).normalize().multiplyScalar(Math.max(10,size*1.6)));this.camera.lookAt(c);this.controls.update();}else{this.center=[0,0];this.scale=48;}this.paint();
  }
  render(project,selected){if(this.drag&&this.cb.revision?.()!==this.drag.baseRevision)this.cancelGesture(false);this.project=project;this.selected=new Set(selected);if(this._hoverHit&&((this._hoverHit.atomId&&!project.atoms.some(a=>a.id===this._hoverHit.atomId))||(this._hoverHit.bondId&&!project.bonds.some(b=>b.id===this._hoverHit.bondId)))){this._hoverHit=null;this.cb.hover?.(null);}if(this.getFragmentConfig()&&this.lastPointer)this.updateFragment(this.lastPointer);else this.paint();this.syncGizmo();}
  showPreview(project){this.previewProject=project;this.paint();}
  clearPreview(){this.previewProject=null;this.paint();}
  paint(){if(this.dimension===2)this.draw2();else this.draw3();}
  screen2(p){return [(p[0]-this.center[0])*this.scale+this.width/2,this.height/2-(p[1]-this.center[1])*this.scale];}
  point2(e){const r=this.svg.getBoundingClientRect();return [(e.clientX-r.left-this.width/2)/this.scale+this.center[0],-(e.clientY-r.top-this.height/2)/this.scale+this.center[1]];}
  localPoint(e,d=this.dimension){const r=(d===2?this.svg:this.renderer.domElement).getBoundingClientRect();return [e.clientX-r.left,e.clientY-r.top];}
  screen3(p){const q=new THREE.Vector3(...p).project(this.camera);return [(q.x+1)*this.width/2,(1-q.y)*this.height/2,q.z];}

  draw2(){
    const root=this.svg;root.replaceChildren();const p=this.activeProject(),display=this.display(),selected=this.effectiveSelection(),byId=new Map(p.atoms.map(a=>[a.id,a]));
    const bonded=new Set(p.bonds.flatMap(b=>[b.a,b.b])),labels=new Map(p.atoms.map(a=>[a.id,atomLabel(a,bonded.has(a.id))]));
    const stateColor={selected:'#277b50',moving:'#d48d31',fixed:'#4780b7',hover:'#2a9b76',ghost:'#4eaa85',invalid:'#cf5e4a',normal:'#314b43'};
    for(const b of p.bonds){
      const a=byId.get(b.a),c=byId.get(b.b);if(!a||!c)continue;
      let [x1,y1]=this.screen2(a.pos2),[x2,y2]=this.screen2(c.pos2);const raw=[x1,y1,x2,y2],len=Math.hypot(x2-x1,y2-y1)||1,ux=(x2-x1)/len,uy=(y2-y1)/len,nx=-uy,ny=ux;
      const gapA=labels.get(a.id)?Math.min(13,len/3):0,gapB=labels.get(c.id)?Math.min(13,len/3):0;x1+=ux*gapA;y1+=uy*gapA;x2-=ux*gapB;y2-=uy*gapB;
      const sa=this.stateFor(a.id,display,selected),sb=this.stateFor(c.id,display,selected),state=sa==='invalid'||sb==='invalid'?'invalid':sa==='ghost'||sb==='ghost'?'ghost':sa===sb?sa:this._hoverHit?.bondId===b.id?'hover':'normal',color=stateColor[state],ghost=['ghost','invalid'].includes(state);
      const g=svg('g',{class:`molecule-bond ${state}`,opacity:ghost?.55:1});
      if(b.direction==='BEGINWEDGE')g.append(svg('polygon',{points:`${x1},${y1} ${x2+nx*5},${y2+ny*5} ${x2-nx*5},${y2-ny*5}`,fill:color}));
      else if(b.direction==='BEGINDASH'){for(let j=1;j<=8;j++){const t=j/9,wx=nx*5*t,wy=ny*5*t;g.append(svg('line',{x1:x1+(x2-x1)*t-wx,y1:y1+(y2-y1)*t-wy,x2:x1+(x2-x1)*t+wx,y2:y1+(y2-y1)*t+wy,stroke:color,'stroke-width':1.5}));}}
      else{const offsets=b.order===3?[-4,0,4]:b.order===2||b.order===1.5?[-2.5,2.5]:[0];offsets.forEach((o,i)=>g.append(svg('line',{x1:x1+nx*o,y1:y1+ny*o,x2:x2+nx*o,y2:y2+ny*o,stroke:color,'stroke-width':state==='normal'?1.9:2.4,'stroke-linecap':'round','stroke-dasharray':b.order===1.5&&i===1?'4 3':'none'})));}
      if(b.stereo&&b.stereo!=='STEREONONE')g.append(svg('text',{x:(x1+x2)/2+nx*12,y:(y1+y2)/2+ny*12,fill:'#8c744a','font-size':9,'pointer-events':'none'},b.stereo.replace('STEREO','')));root.append(g);
      if(!ghost)root.append(svg('line',{x1:raw[0],y1:raw[1],x2:raw[2],y2:raw[3],class:'bond-hit','data-bond':b.id,stroke:'transparent','stroke-width':14}));
    }
    for(const [i,a] of p.atoms.entries()){
      const [x,y]=this.screen2(a.pos2),state=this.stateFor(a.id,display,selected),ghost=['ghost','invalid'].includes(state),g=svg('g',{class:`molecule-atom ${state}`,...(ghost?{'pointer-events':'none'}:{'data-atom':a.id})});
      if(state!=='normal')g.append(svg('circle',{cx:x,cy:y,r:labels.get(a.id)?15:9,fill:stateColor[state],'fill-opacity':.12,stroke:stateColor[state],'stroke-width':1.5}));
      if(!ghost)g.append(svg('circle',{cx:x,cy:y,r:11,fill:'transparent',class:'atom-hit','data-atom':a.id}));const label=labels.get(a.id);
      if(label)g.append(svg('text',{x,y:y+4,'text-anchor':'middle',class:'atom-label',fill:stateColor[state],'paint-order':'stroke',stroke:'#f8faf5','stroke-width':4,'stroke-linejoin':'round','pointer-events':'none',opacity:ghost?.65:1},label));
      if(display.numbers&&!ghost)g.append(svg('text',{x:x+11,y:y-10,class:'atom-id','pointer-events':'none'},`${i+1}${a.chiralTag&&a.chiralTag!=='CHI_UNSPECIFIED'?'*':''}`));root.append(g);
    }
    this.drawAnchors(root,p,false);this.drawPorts(root,p,false);this.drawGestureOverlay(root,false);
  }
  material(element,state,bond=false){
    const key=`${bond?'bond':element}:${state}`;
    if(!this.materials.has(key)){const accents={selected:0x8dbb55,moving:0xd5a04b,fixed:0x5d9aca,hover:0x55ab91,ghost:0x6db798,invalid:0xd65d52},color=accents[state]||(bond?0x899c89:colors[element]||0xb887ab),ghost=['ghost','invalid'].includes(state);this.materials.set(key,new THREE.MeshStandardMaterial({color,roughness:.5,metalness:.04,transparent:ghost,opacity:ghost?.4:1,depthWrite:!ghost}));}return this.materials.get(key);
  }
  draw3(){
    const p=this.activeProject(),display=this.display(),selected=this.effectiveSelection(),byId=new Map(p.atoms.map(a=>[a.id,a])),atomIds=new Set(byId.keys()),bondIds=new Set(p.bonds.map(b=>b.id));
    for(const [id,mesh] of this.atomMeshes)if(!atomIds.has(id)){this.group.remove(mesh);this.atomMeshes.delete(id);}
    for(const [id,meshes] of this.bondMeshes)if(!bondIds.has(id)){meshes.forEach(mesh=>this.group.remove(mesh));this.bondMeshes.delete(id);}this.meshes=[];
    for(const a of p.atoms){let mesh=this.atomMeshes.get(a.id);if(!mesh){mesh=new THREE.Mesh(this.sphereGeometry);mesh.userData.atomId=a.id;this.group.add(mesh);this.atomMeshes.set(a.id,mesh);}const state=this.stateFor(a.id,display,selected);mesh.position.set(...a.pos3);mesh.scale.setScalar(a.element==='H'?.17:.30);mesh.material=this.material(a.element,state);if(!['ghost','invalid'].includes(state))this.meshes.push(mesh);}
    for(const b of p.bonds){
      const a=byId.get(b.a),c=byId.get(b.b);if(!a||!c)continue;const start=new THREE.Vector3(...a.pos3),end=new THREE.Vector3(...c.pos3),delta=end.clone().sub(start),len=delta.length(),offsets=b.order===3?[-.13,0,.13]:b.order===2||b.order===1.5?[-.075,.075]:[0];
      let meshes=this.bondMeshes.get(b.id);if(!meshes||meshes.length!==offsets.length){meshes?.forEach(mesh=>this.group.remove(mesh));meshes=offsets.map(()=>{const mesh=new THREE.Mesh(this.bondGeometry);mesh.userData.bondId=b.id;this.group.add(mesh);return mesh;});this.bondMeshes.set(b.id,meshes);}
      const sa=this.stateFor(a.id,display,selected),sb=this.stateFor(c.id,display,selected),state=sa==='invalid'||sb==='invalid'?'invalid':sa==='ghost'||sb==='ghost'?'ghost':sa===sb?sa:this._hoverHit?.bondId===b.id?'hover':'normal',direction=len>.00001?delta.clone().normalize():new THREE.Vector3(0,1,0),side=new THREE.Vector3().crossVectors(direction,new THREE.Vector3(0,0,1));if(side.length()<.01)side.crossVectors(direction,new THREE.Vector3(0,1,0));side.normalize();
      meshes.forEach((mesh,i)=>{mesh.visible=len>.00001;mesh.material=this.material('',state,true);mesh.position.copy(start).add(end).multiplyScalar(.5).addScaledVector(side,offsets[i]);mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0),direction);mesh.scale.set(1,Math.max(len,.00001),1);if(!['ghost','invalid'].includes(state)&&mesh.visible)this.meshes.push(mesh);});
    }
    this.repaint3();
  }
  repaint3(){this.renderer.render(this.scene,this.camera);this.drawLabels();this.overlay.replaceChildren();if(this.dimension===3){const p=this.activeProject();this.drawAnchors(this.overlay,p,true);this.drawPorts(this.overlay,p,true);this.drawGestureOverlay(this.overlay,true);}}
  drawLabels(){
    this.labels.replaceChildren();if(this.dimension!==3)return;const p=this.activeProject(),display=this.display();
    for(const [i,a] of p.atoms.entries()){if(this.fragmentVisual?.atomIds?.includes(a.id))continue;if(!display.numbers&&this._hoverHit?.atomId!==a.id&&!(this.selected.size<=4&&this.selected.has(a.id)))continue;if(display.anchors.some(x=>x.atomId===a.id))continue;const [x,y,z]=this.screen3(a.pos3);if(z>1||z< -1||x<0||x>this.width||y<0||y>this.height)continue;const label=document.createElement('span');label.textContent=`${display.numbers?i+1:''}${a.element}`;label.style.left=`${x+8}px`;label.style.top=`${y-14}px`;if(this.selected.has(a.id))label.className='selected';this.labels.append(label);}
  }
  drawAnchors(root,p,three){
    const anchors=this.display().anchors,byId=new Map(p.atoms.map(a=>[a.id,a])),points=anchors.map(x=>{const a=byId.get(x.atomId);return a?(three?this.screen3(a.pos3):this.screen2(a.pos2)):null;});
    if(points.filter(Boolean).length>1)root.append(svg('polyline',{points:points.filter(Boolean).map(x=>x.slice(0,2).join(',')).join(' '),fill:'none',stroke:'#c68932','stroke-width':1.4,'stroke-dasharray':'5 4','pointer-events':'none'}));
    anchors.forEach((a,i)=>{const point=points[i];if(!point||(three&&(point[2]>1||point[2]< -1)))return;const [x,y]=point;root.append(svg('circle',{cx:x,cy:y,r:12,fill:'none',stroke:'#c68932','stroke-width':2,class:'measurement-anchor','pointer-events':'none'}));root.append(svg('text',{x:x+13,y:y-11,class:'measurement-label',fill:'#ae701b','font-size':13,'font-weight':700,'paint-order':'stroke',stroke:'#fffefa','stroke-width':4,'pointer-events':'none'},a.label));});
  }
  drawPorts(root,p,three){
    const config=this.getFragmentConfig();if(!config)return;const ports=this.fragmentVisual?.ports||[],byId=new Map(p.atoms.map(a=>[a.id,a]));
    for(const [i,port] of ports.entries()){const a=byId.get(port.atomId);if(!a)continue;const point=three?this.screen3(a.pos3):this.screen2(a.pos2);if(three&&(point[2]>1||point[2]< -1))continue;root.append(svg('circle',{cx:point[0],cy:point[1],r:i===(config.portIndex||0)?12:6,fill:'none',stroke:this.fragmentVisual?.error?'#cf5e4a':'#238f68','stroke-width':2,'stroke-dasharray':'3 3',class:'fragment-port','pointer-events':'none'}));}
    if(this.fragmentVisual?.error&&this.lastPointer){const [x,y]=this.localPoint(this.lastPointer,this.dimension);root.append(svg('text',{x:Math.max(8,Math.min(this.width-130,x+16)),y:Math.max(20,y-20),fill:'#ad4336','font-size':11,class:'fragment-error','pointer-events':'none'},'此位点无法连接'));}
  }

  setRay(e){const r=this.renderer.domElement.getBoundingClientRect();this.ray.setFromCamera(new THREE.Vector2((e.clientX-r.left)/r.width*2-1,-(e.clientY-r.top)/r.height*2+1),this.camera);}
  hit3(e){this.setRay(e);this.scene.updateMatrixWorld(true);return this.ray.intersectObjects(this.meshes)[0]?.object.userData||{};}
  hit2(e){
    const point=this.localPoint(e,2);let best=null,nearest=13;for(const a of this.project.atoms){const q=this.screen2(a.pos2),d=Math.hypot(point[0]-q[0],point[1]-q[1]);if(d<nearest){best={atomId:a.id};nearest=d;}}if(best)return best;
    const byId=new Map(this.project.atoms.map(a=>[a.id,a]));nearest=8;for(const b of this.project.bonds){const a=byId.get(b.a),c=byId.get(b.b);if(!a||!c)continue;const d=segmentDistance(point,this.screen2(a.pos2),this.screen2(c.pos2));if(d<nearest){best={bondId:b.id};nearest=d;}}return best||{};
  }
  nearestAtom(e,dimension,excluding=null){const point=this.localPoint(e,dimension);let nearest=17,id=null;for(const a of this.project.atoms){if(a.id===excluding)continue;const q=dimension===2?this.screen2(a.pos2):this.screen3(a.pos3);if(dimension===3&&(q[2]>1||q[2]< -1))continue;const d=Math.hypot(q[0]-point[0],q[1]-point[1]);if(d<nearest){id=a.id;nearest=d;}}return id;}
  planePoint(e,origin=new THREE.Vector3()){this.setRay(e);const plane=new THREE.Plane().setFromNormalAndCoplanarPoint(this.camera.getWorldDirection(new THREE.Vector3()),origin),out=new THREE.Vector3();return this.ray.ray.intersectPlane(plane,out)?out:null;}
  point3(e,hit={},origin=this.controls.target){const a=this.project.atoms.find(a=>a.id===hit.atomId);return a?[...a.pos3]:this.planePoint(e,origin)?.toArray();}
  setHover(hit){const next=hasHit(hit)?{...hit}:null;if(hitKey(next)===hitKey(this._hoverHit))return;this._hoverHit=next;this.cb.hover?.(next);this.paint();}

  beginDrag(e,id,dimension,bondId=null){
    if(id&&!this.selected.has(id))this.cb.pick({atomId:id},e);else if(!id&&bondId)this.cb.pick({bondId},e);
    const ids=new Set(this.cb.moveIds?.()||this.selected),base=structuredClone(this.project);if(!ids.size)return false;this.cb.interaction?.();
    const a=base.atoms.find(a=>a.id===id)||base.atoms.find(a=>ids.has(a.id)),mode=this.cb.mode()==='rotate'?'rotate':'move';
    this.drag={id,ids,base,baseRevision:this.cb.revision(),dimension,start:[e.clientX,e.clientY],point2:this.point2(e),origin:[...a.pos3],mode,delta:dimension===2?[0,0]:[0,0,0],degrees:0};if(dimension===3)this.drag.point3=this.planePoint(e,new THREE.Vector3(...a.pos3));return true;
  }
  previewDrag(d){
    const p=structuredClone(d.base);
    if(d.mode==='rotate'){const axis=new THREE.Vector3(...(d.axis||this.cb.axis())).normalize(),points=p.atoms.filter(a=>d.ids.has(a.id)),center=d.pivot?new THREE.Vector3(...d.pivot):new THREE.Vector3();if(!d.pivot){points.forEach(a=>center.add(new THREE.Vector3(...a.pos3)));center.multiplyScalar(1/Math.max(1,points.length));}for(const a of points)a.pos3=new THREE.Vector3(...a.pos3).sub(center).applyAxisAngle(axis,d.degrees*Math.PI/180).add(center).toArray();d.origin=center.toArray();}
    else for(const a of p.atoms)if(d.ids.has(a.id)){const key=d.dimension===2?'pos2':'pos3';a[key]=a[key].map((v,i)=>v+d.delta[i]);}this.previewProject=p;this.paint();
  }
  updateDrag(e){const d=this.drag;if(!d||d.gizmo)return;if(d.mode==='rotate')d.degrees=(e.clientX-d.start[0])*.5;else if(d.dimension===2){const q=this.point2(e);d.delta=[q[0]-d.point2[0],q[1]-d.point2[1]];}else{const q=this.planePoint(e,new THREE.Vector3(...d.origin));if(!q||!d.point3)return;d.delta=q.sub(d.point3).toArray();}this.previewDrag(d);}
  endDrag(){const d=this.drag;if(!d)return;this.drag=null;this.previewProject=null;this.paint();if(d.baseRevision===this.cb.revision())this.cb.drag(d);this.syncGizmo();}
  syncGizmo(){
    if(this.transform.dragging||this.drag?.gizmo)return;const mode=this.cb.mode(),ids=this.effectiveSelection();
    if(this.dimension!==3||!this.project.has3d||!ids.size||!['move','rotate'].includes(mode)||this.getFragmentConfig()){this.transform.detach();return;}
    const points=this.project.atoms.filter(a=>ids.has(a.id)),center=new THREE.Vector3();points.forEach(a=>center.add(new THREE.Vector3(...a.pos3)));center.multiplyScalar(1/Math.max(1,points.length));this.pivot.position.copy(center);this.pivot.quaternion.identity();this.pivot.scale.set(1,1,1);this.transform.setMode(mode==='rotate'?'rotate':'translate');this.transform.attach(this.pivot);this.repaint3();
  }
  bindTransform(){
    this.transform.addEventListener('change',()=>this.repaint3());
    this.transform.addEventListener('mouseDown',()=>{this.gizmoConsumed=true;this.controls.enabled=false;this.cb.interaction?.();this.drag={gizmo:true,ids:this.effectiveSelection(),base:structuredClone(this.project),baseRevision:this.cb.revision(),dimension:3,mode:this.transform.getMode()==='rotate'?'rotate':'move',origin:this.pivot.position.toArray(),pivot:this.pivot.position.toArray(),delta:[0,0,0],degrees:0};});
    this.transform.addEventListener('objectChange',()=>{const d=this.drag;if(!d?.gizmo)return;if(d.mode==='rotate'){d.axis=this.transform.rotationAxis.toArray();d.degrees=this.transform.rotationAngle*180/Math.PI;}else d.delta=this.pivot.position.clone().sub(new THREE.Vector3(...d.origin)).toArray();this.previewDrag(d);});
    this.transform.addEventListener('mouseUp',()=>{this.controls.enabled=true;if(this.drag?.gizmo)this.endDrag();});
  }

  getFragmentConfig(){const value=this.cb.previewFragment?.()||this.fragmentPreview;if(!value)return null;return value.template?value:{template:value,portIndex:0,order:1,dimension:this.dimension};}
  setFragmentPreview(template){this.fragmentPreview=template;this.fragmentVisual=null;if(!template)this.fragmentPlacementId=null;if(template&&this.lastPointer)this.updateFragment(this.lastPointer);else this.paint();this.syncGizmo();}
  updateFragment(e){
    const config=this.getFragmentConfig();if(!config){if(this.fragmentVisual){this.fragmentVisual=null;this.paint();}return;}
    const hit=this.dimension===2?this.hit2(e):this.hit3(e),point=this.dimension===2?this.point2(e):this.point3(e,hit);if(!point)return;
    const opts={targetId:hit.atomId||null,portIndex:config.portIndex||0,pos2:point.slice(0,2),pos3:this.dimension===3?point:[point[0],point[1],0],dimension:this.dimension,order:config.order||1,rigid:true,placementId:this.fragmentPlacementId||(this.fragmentPlacementId=`preview-${crypto.randomUUID()}`)};
    try{if(!Model.fragmentPlacementPreview)return;this.fragmentVisual={...Model.fragmentPlacementPreview(this.project,config.template,opts),error:null};}
    catch(error){try{this.fragmentVisual={...Model.fragmentPlacementPreview(this.project,config.template,{...opts,targetId:null}),error:error.message};}catch{this.fragmentVisual={project:this.project,atomIds:[],ports:[],error:error.message};}}
    this.cb.fragmentFeedback?.(this.fragmentVisual.error||null,hit);this.paint();
  }
  cancelGesture(clearFragment=true){const wasGizmo=this.drag?.gizmo;this.drag=null;this.gesture=null;this.previewProject=null;if(wasGizmo){this.transform.reset();this.transform.dragging=false;this.transform.axis=null;this.transform.detach();}this.controls.enabled=true;if(clearFragment){this.fragmentPreview=null;this.fragmentVisual=null;this.fragmentPlacementId=null;}this.paint();this.syncGizmo();}

  beginStroke(e,hit,dimension){const a=this.project.atoms.find(a=>a.id===hit.atomId),startPoint=a?[...a[dimension===2?'pos2':'pos3']]:dimension===2?this.point2(e):this.point3(e);if(!startPoint)return;this.cb.interaction?.();this.gesture={kind:'stroke',start:[e.clientX,e.clientY],startId:a?.id||null,startPoint,endId:null,endPoint:[...startPoint],hit,dimension,moved:false};}
  updateStroke(e){
    const g=this.gesture;if(g?.kind!=='stroke')return;g.moved||=Math.hypot(e.clientX-g.start[0],e.clientY-g.start[1])>4;
    const snapId=this.nearestAtom(e,g.dimension,g.startId);g.endId=snapId;
    if(snapId)g.endPoint=[...this.project.atoms.find(a=>a.id===snapId)[g.dimension===2?'pos2':'pos3']];
    else if(g.dimension===2)g.endPoint=snapBondEndpoint(g.startPoint,this.point2(e),e.altKey);
    else{const q=this.planePoint(e,new THREE.Vector3(...g.startPoint));if(!q)return;const delta=q.sub(new THREE.Vector3(...g.startPoint));if(!e.altKey){if(delta.length()<.00001)delta.set(1,0,0);delta.setLength(1.5);}g.endPoint=new THREE.Vector3(...g.startPoint).add(delta).toArray();}this.paint();
  }
  drawGestureOverlay(root,three){
    const g=this.gesture;if(!g)return;
    if(g.kind==='stroke'){const a=three?this.screen3(g.startPoint):this.screen2(g.startPoint),b=three?this.screen3(g.endPoint):this.screen2(g.endPoint);root.append(svg('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:'#24855c','stroke-width':2.4,'stroke-dasharray':'5 3',class:'drawing-preview','pointer-events':'none'}));root.append(svg('circle',{cx:b[0],cy:b[1],r:g.endId?11:4,fill:g.endId?'none':'#24855c',stroke:'#24855c','stroke-width':1.5,'pointer-events':'none'}));}
    else if(g.kind==='region'&&g.moved){const attrs={fill:'#3b8b6b','fill-opacity':.10,stroke:'#3b8b6b','stroke-width':1.2,'stroke-dasharray':'4 3',class:'selection-region','pointer-events':'none'};if(g.lasso)root.append(svg('polygon',{points:g.points.map(p=>p.join(',')).join(' '),...attrs}));else root.append(svg('rect',{x:Math.min(g.startLocal[0],g.endLocal[0]),y:Math.min(g.startLocal[1],g.endLocal[1]),width:Math.abs(g.endLocal[0]-g.startLocal[0]),height:Math.abs(g.endLocal[1]-g.startLocal[1]),...attrs}));}
  }
  regionIds(g){const min=g.startLocal.map((v,i)=>Math.min(v,g.endLocal[i])),max=g.startLocal.map((v,i)=>Math.max(v,g.endLocal[i]));return this.project.atoms.filter(a=>{const q=g.dimension===2?this.screen2(a.pos2):this.screen3(a.pos3);if(g.dimension===3&&(q[2]< -1||q[2]>1))return false;return g.lasso?pointInPolygon(q,g.points):q[0]>=min[0]&&q[0]<=max[0]&&q[1]>=min[1]&&q[1]<=max[1];}).map(a=>a.id);}
  pickClick(hit,e){
    const now=performance.now(),last=this.lastClick;
    const selectable=['select','move'].includes(this.cb.mode())&&!e.shiftKey;
    // SVG atom nodes are redrawn after selection, so native dblclick cannot be
    // trusted to retain the same DOM target. Detect it by stable molecule ID.
    if(selectable&&last&&last.key===hitKey(hit)&&now-last.time<360&&Math.hypot(e.clientX-last.x,e.clientY-last.y)<6){this.lastClick=null;this.cb.doublePick?.(hit,e);}
    else{this.lastClick=selectable?{key:hitKey(hit),time:now,x:e.clientX,y:e.clientY}:null;this.cb.pick(hit,e);}
  }

  down(e,dimension){
    if(dimension===2&&(e.button===1||e.button===2)){e.preventDefault();this.gesture={kind:'pan',start:[e.clientX,e.clientY],center:[...this.center],dimension};this.svg.setPointerCapture(e.pointerId);return;}
    if(e.button!==0)return;if(dimension===3&&this.transform.dragging){this.gizmoConsumed=true;return;}
    const el=dimension===2?this.svg:this.renderer.domElement,hit=dimension===2?this.hit2(e):this.hit3(e),mode=this.cb.mode();this.lastPointer=e;this.setHover(hit);
    if(this.getFragmentConfig())this.gesture={kind:'fragment',dimension};
    else if(mode==='bond'&&!hit.bondId)this.beginStroke(e,hit,dimension);
    else if(['select','move','rotate'].includes(mode)&&hasHit(hit)&&!e.shiftKey)this.gesture={kind:'potential-drag',start:[e.clientX,e.clientY],hit,dimension,downEvent:e};
    else if(['select','move'].includes(mode)&&!hasHit(hit)){const q=this.localPoint(e,dimension);this.gesture={kind:'region',start:[e.clientX,e.clientY],startLocal:q,endLocal:q,points:[q],lasso:!!this.display().lasso,dimension,moved:false};}
    else this.gesture={kind:'click',start:[e.clientX,e.clientY],hit,dimension};if(this.gesture)el.setPointerCapture(e.pointerId);
  }
  move(e,dimension){
    this.lastPointer=e;if(this.drag?.gizmo)return;const g=this.gesture;
    if(g?.kind==='pan'){this.center=[g.center[0]-(e.clientX-g.start[0])/this.scale,g.center[1]+(e.clientY-g.start[1])/this.scale];this.draw2();this.cb.viewChange?.();return;}
    if(g?.kind==='potential-drag'&&Math.hypot(e.clientX-g.start[0],e.clientY-g.start[1])>4){this.beginDrag(g.downEvent,g.hit.atomId,dimension,g.hit.bondId);g.kind='drag';}
    if(this.drag){this.updateDrag(e);return;}if(g?.kind==='stroke'){this.updateStroke(e);return;}
    if(g?.kind==='region'){g.endLocal=this.localPoint(e,dimension);g.moved||=Math.hypot(e.clientX-g.start[0],e.clientY-g.start[1])>4;if(Math.hypot(g.endLocal[0]-g.points.at(-1)[0],g.endLocal[1]-g.points.at(-1)[1])>2)g.points.push(g.endLocal);this.paint();return;}
    const hit=dimension===2?this.hit2(e):this.hit3(e);this.setHover(hit);if(this.getFragmentConfig())this.updateFragment(e);
  }
  up(e,dimension){
    if(e.button!==0&&this.gesture?.kind!=='pan')return;if(this.gizmoConsumed){this.gizmoConsumed=false;this.gesture=null;this.syncGizmo();return;}if(this.drag){this.endDrag();this.gesture=null;return;}
    const g=this.gesture;this.gesture=null;if(!g)return;const hit=dimension===2?this.hit2(e):this.hit3(e);
    if(g.kind==='pan'){this.paint();return;}
    if(g.kind==='fragment'){const point=dimension===2?this.point2(e):this.point3(e,hit);if(point)this.cb.fragmentPlace?.(point,hit,e,dimension);}
    else if(g.kind==='stroke'){this.gesture=g;this.updateStroke(e);this.gesture=null;if(g.moved||!g.startId)this.cb.stroke?.({...g,event:e});else this.cb.pick(g.hit,e);}
    else if(g.kind==='region'&&g.moved)this.cb.region?.(this.regionIds(g),e);
    else if(hasHit(hit))this.pickClick(hit,e);
    else{const point=dimension===2?this.point2(e):this.point3(e);if(point)this.cb.blank(point,e,dimension);}this.paint();this.syncGizmo();
  }
  bind2(){
    this.svg.addEventListener('pointerdown',e=>this.down(e,2));this.svg.addEventListener('pointermove',e=>this.move(e,2));this.svg.addEventListener('pointerup',e=>this.up(e,2));this.svg.addEventListener('pointercancel',()=>this.cancelGesture(false));this.svg.addEventListener('pointerleave',()=>{if(!this.gesture&&!this.drag)this.setHover(null);});this.svg.addEventListener('contextmenu',e=>e.preventDefault());
    this.svg.addEventListener('wheel',e=>{e.preventDefault();const before=this.point2(e);this.scale=Math.max(8,Math.min(160,this.scale*Math.exp(-e.deltaY*.001)));const after=this.point2(e);this.center=this.center.map((v,i)=>v+before[i]-after[i]);this.draw2();this.cb.viewChange?.();},{passive:false});
  }
  bind3(){const el=this.renderer.domElement;el.addEventListener('pointerdown',e=>this.down(e,3));el.addEventListener('pointermove',e=>this.move(e,3));el.addEventListener('pointerup',e=>this.up(e,3));el.addEventListener('pointercancel',()=>this.cancelGesture(false));el.addEventListener('pointerleave',()=>{if(!this.gesture&&!this.drag)this.setHover(null);});el.addEventListener('contextmenu',e=>e.preventDefault());}
  updateMode(){this.cancelGesture(false);this.controls.enabled=true;const cursor=['atom','bond'].includes(this.cb.mode())||this.getFragmentConfig()?'crosshair':this.cb.mode()==='move'?'move':'default';this.svg.style.cursor=cursor;this.renderer.domElement.style.cursor=cursor;this.syncGizmo();}
}
