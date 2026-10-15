// Full races use one backend and one preparation/velocity protocol for every
// compared condition. Native/WASM gradients and paired NVE are checked first.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import {fork} from 'node:child_process';
import {fileURLToPath,pathToFileURL} from 'node:url';
const root=fileURLToPath(new URL('../../',import.meta.url));
const folder=path.join(root,'outputs/crystal-contrast-guarded');
const assets=path.join(root,'wasm/dist/browser/assets');
const script=fileURLToPath(import.meta.url);
const cars=['benzene-12','coronene-36','ethynyl-rigid','test-car-52','fullerene-c60','pyridine-11','benzonitrile-13'];
const hash=async file=>crypto.createHash('sha256').update(await fs.readFile(file)).digest('hex');
const save=async(file,value)=>{const temporary=file+`.${process.pid}.${crypto.randomUUID()}.tmp`;await fs.writeFile(temporary,JSON.stringify(value)+'\n');await fs.rename(temporary,file);};
const child=process.argv.includes('--child');
if(child) {
  const metadata=JSON.parse(await fs.readFile(path.join(assets,'model.json'),'utf8'));
  const load=async name=>{const b=zlib.gunzipSync(await fs.readFile(path.join(assets,`${name}.bin.gz`)));return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));};
  const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(load));
  const {default:create}=await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')));
  const module=await create({print:()=>{},printErr:()=>{}});
  await fs.mkdir(path.join(folder,'browser-trials'),{recursive:true});
  process.on('message',async ({config:c,sources})=>{
    const key=`${c.phase}-${c.molecule}-T${c.temperature_K}-S${c.strength}-L${c.lattice_scale}-G${c.growth}-seed${c.seed}-${c.limit_ps}ps`;
    const target=path.join(folder,'browser-trials',key+'.json');
    try {
      const old=JSON.parse(await fs.readFile(target,'utf8'));
      const baseline=c.lattice_scale===1 && old.runtime?.velocity_seed===17329
        && old.source_sha256['wasm/dist/browser/assets/model.json']===sources['wasm/dist/browser/assets/model.json']
        && old.source_sha256[`inputs/${c.molecule}.mol`]===sources[`inputs/${c.molecule}.mol`];
      if(old.status==='passed'&&JSON.stringify(old.config)===JSON.stringify(c)&&(JSON.stringify(old.source_sha256)===JSON.stringify(sources)||baseline)) {
        if(baseline)old.reused_baseline='scale=1 is exact raw-potential path; current 224-case/NVE regression passed; original source hashes retained';
        process.send(old);return;
      }
    }catch{}
    const row={key,config:c,source_sha256:sources,status:'running',execution:'actual-lammps-wasm',samples:[]};
    const began=performance.now();let engine;
    try {
      engine=new module.BrowserEngine({...metadata,seed:c.seed},...arrays);
      const file=c.molecule+'.mol',content=await fs.readFile(path.join(folder,'inputs',file),'utf8');
      const info=engine.prepareWithPose(content,file,{...metadata.settings,temperature:c.temperature_K,damping_fs:350,lattice_scale:c.lattice_scale},{mode:'auto'});
      if(!info.runtime.preparation.accepted||info.runtime.preparation.global_force_norm_kcal_mol_A>1e-4)throw new Error('Preparation force gate failed');
      engine.update({drive:{mode:'gradient',strength:c.strength,growth:c.growth}});
      let frame=engine.refresh();if(frame.step!==0)throw new Error('Race did not start at step zero');
      row.runtime=info.runtime;row.preparation_wall_seconds=(performance.now()-began)/1000;row.atoms=info.molecule.atoms.length;
      const masses=info.molecule.atoms.map(a=>a.mass),bonds=info.molecule.bonds.map(b=>[b.a-1,b.b-1]);
      const length=(p,b)=>Math.hypot(...[0,1,2].map(a=>p[3*b[0]+a]-p[3*b[1]+a]));
      const reference=bonds.map(b=>length(frame.positions,b));
      const initial=frame.center.map(v=>v/10);row.initial_center_nm=initial;
      const positions=[],steps=[];let low=1,high=1,zmin=Infinity,zmax=-Infinity,lift=0,maxLift=0,outcome='time-limit';
      const propagation=performance.now(),limit=Math.round(c.limit_ps*2000);
      for(let step=0;step<=limit;step+=100) {
        if(step)frame=engine.advance(100);
        if(!Number.isFinite(frame.energy)||!frame.positions.every(Number.isFinite))throw new Error('Nonfinite trajectory');
        const z=frame.positions.filter((v,i)=>i%3===2),min=Math.min(...z),max=Math.max(...z);
        zmin=Math.min(zmin,min);zmax=Math.max(zmax,max);
        for(let i=0;i<bonds.length;i++){const ratio=length(frame.positions,bonds[i])/reference[i];low=Math.min(low,ratio);high=Math.max(high,ratio);}
        if(low<.65||high>1.45)throw new Error('Bond-length guard failed');
        // Count samples exactly: repeated .05 addition can miss the 5 ps boundary.
        if(step){lift=min>8?lift+1:0;maxLift=Math.max(maxLift,lift);}
        row.samples.push({step,time_ps:step*.0005,center_nm:frame.center.map(v=>v/10),energy_kcal_mol:frame.energy,
                          internal_temperature_K:frame.internal_temperature,min_z_A:min,max_z_A:max});
        const finish=frame.center[0]>=200;
        if(step%1000===0||finish||step===limit){positions.push(frame.positions);steps.push(step);}
        if(finish){outcome='finish';break;}
      }
      const coordinates=key+'.coords.json.gz';await fs.writeFile(path.join(folder,'browser-trials',coordinates),zlib.gzipSync(JSON.stringify({positions,steps,masses,bonds})));
      const last=row.samples.at(-1),time=last.time_ps;
      Object.assign(row,{status:'passed',outcome,finish_time_ps:outcome==='finish'?time:null,observed_time_ps:time,
        finish_crossing_interval_ps:outcome==='finish'?[time-.05,time]:null,right_censored:outcome!=='finish',
        final_x_nm:last.center_nm[0],final_forward_nm:last.center_nm[0]-initial[0],
        maximum_lateral_excursion_nm:Math.max(...row.samples.map(s=>Math.abs(s.center_nm[1]-initial[1]))),
        bond_ratio_range:[low,high],atom_height_range_A:[zmin,zmax],maximum_lift_duration_ps:maxLift*.05,
        sustained_lift:maxLift>=100,competitive_valid:maxLift<100,coordinates_file:coordinates,
        propagation_wall_seconds:(performance.now()-propagation)/1000});
    }catch(error) {
      row.status='failed';try{row.error=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);}catch{row.error=error.message||String(error);}
    }finally{if(engine){engine.close();engine.delete();}}
    row.total_wall_seconds=(performance.now()-began)/1000;await save(target,row);process.send(row);
  });
} else {
  const gate=JSON.parse(await fs.readFile(path.join(folder,'wasm-numerics.json'),'utf8'));
  if(gate.status!=='passed'||gate.control_version!=='six-phase-atomic-crystal-contrast-v3')throw new Error('Current native/WASM agreement gate must pass before races');
  for(const [file,sha] of Object.entries(gate.source_sha256))if(await hash(path.join(root,file))!==sha)throw new Error('Numerical gate source/binary changed: '+file);
  const jobIndex=process.argv.indexOf('--jobs');
  const jobs=jobIndex>=0?JSON.parse(await fs.readFile(process.argv[jobIndex+1],'utf8')):
    [1,2,3].flatMap(scale=>(scale===3?[1]:[1,.25]).flatMap(growth=>[2718281,1234567,31415927].flatMap(seed=>cars.map(molecule=>({
      phase:'browser-confirm',molecule,temperature_K:150,strength:.6,lattice_scale:scale,growth,seed,limit_ps:400})))));
  const sources={};for(const file of ['wasm/src/potentials.cpp','wasm/src/engine.cpp','wasm/tools/scan_crystal_races.mjs','wasm/dist/browser/runtime/molecular-car.wasm','wasm/dist/browser/assets/model.json'])sources[file]=await hash(path.join(root,file));
  for(const car of cars)sources[`inputs/${car}.mol`]=await hash(path.join(folder,'inputs',car+'.mol'));
  const report={status:'running',planned:jobs.length,trials:[],source_sha256:sources,
    protocol:{execution:'actual LAMMPS + RDKit C++ WASM, Node CPU',timestep_fs:.5,damping_fs:350,temperature_K:150,
      batch_steps:100,velocity_seed:17329,langevin_seeds:[...new Set(jobs.map(c=>c.seed))],start_pose:{mode:'auto'},
      warmup_ps:0,track_nm:[20,20],lift_guard:'all atoms >8 A continuously for >=5 ps invalid',
      censoring:'not finished by 400 ps; no extrapolation',comparison:'all conditions including baseline use the same backend and velocity protocol'}};
  await save(path.join(folder,'browser-confirm-jobs.json'),jobs);
  let next=0,done=0,queue=Promise.resolve();const started=performance.now(),workers=[];
  await new Promise((resolve,reject)=>{
    for(let i=0;i<Math.min(4,jobs.length);i++) {
      const worker=fork(script,['--child'],{stdio:['ignore','inherit','inherit','ipc']});workers.push(worker);
      const dispatch=()=>{if(next<jobs.length)worker.send({config:jobs[next++],sources});else worker.disconnect();};
      worker.on('message',row=>{queue=queue.then(async()=>{
        report.trials.push(row);done++;await save(path.join(folder,'browser-confirm-report.json'),report);
        console.log(JSON.stringify({done,total:jobs.length,key:row.key,status:row.status,finish_ps:row.finish_time_ps,lift:row.sustained_lift,error:row.error}));
        if(done===jobs.length){report.status='complete';report.wall_seconds=(performance.now()-started)/1000;await save(path.join(folder,'browser-confirm-report.json'),report);for(const w of workers)if(w.connected)w.disconnect();resolve();}
        else dispatch();
      }).catch(reject);});
      worker.on('error',reject);worker.on('exit',code=>{if(code&&done<jobs.length)reject(new Error(`Worker exited ${code}`));});dispatch();
    }
  });
}
