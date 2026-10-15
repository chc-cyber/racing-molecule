import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import assert from 'node:assert/strict';
import {root,out,read,initialize,sources,checkSources} from './atlas_support.mjs';
const {engine:e,meta}=await initialize();
const report={status:'running',sources:await sources(),static:[],finite_difference:[],barriers:[],trajectory:[]};
const fixtures=[...await read(path.join(out,'native-fixtures.json')),...await read(path.join(out,'native-extra-fixtures.json'))];
const maxdiff=(a,b)=>Math.max(...a.map((v,i)=>Math.abs(v-b[i])));
const rmsdiff=(a,b)=>Math.sqrt(a.reduce((s,v,i)=>s+(v-b[i])**2,0)/a.length);
try {
  let car;
  for(const f of fixtures) {
    if(f.car!==car){car=f.car;e.prepareChemistry(f.content,car+'.mol');}
    e.update({lattice_scale:f.scale});e.fourierEnabled(true);
    const current=e.energyForces(f.positions,1,.6,false),counts=e.proxyCounts();
    const error=maxdiff(current.surface_forces,f.forces),energy=Math.abs(current.surface_energy-f.energy);
    assert(error<1e-7&&energy<1e-7,`${car}/${f.label}/${f.scale}: ${error}, ${energy}`);
    e.fourierEnabled(false);const original=e.energyForces(f.positions,1,.6,false);e.fourierEnabled(true);
    report.static.push({car,label:f.label,scale:f.scale,port_force_error:error,port_energy_error:energy,counts,
      reference_max_force:Math.max(...original.surface_forces.map(Math.abs)),
      approximation_force_error:maxdiff(current.surface_forces,original.surface_forces),
      approximation_energy_error:Math.abs(current.surface_energy-original.surface_energy)});
    if(f.scale===2&&(car==='benzene-12'||f.label==='start')) {
      const p=f.positions.slice(),h=1e-6;let error=0;
      for(let k=0;k<p.length;k++) {
        const old=p[k];p[k]=old+h;const a=e.energyForces(p,1,.6,false).surface_energy;
        p[k]=old-h;const b=e.energyForces(p,1,.6,false).surface_energy;p[k]=old;
        error=Math.max(error,Math.abs(-(a-b)/(2*h)-current.surface_forces[k]));
      }
      assert(error<1e-4,'finite difference '+car+' '+f.label+' '+error);
      report.finite_difference.push({car,label:f.label,error});
    }
  }
  for(const car of ['benzene-12','fullerene-c60','pyridine-11']) {
    const f=fixtures.find(f=>f.car===car&&f.label==='start'&&f.scale===2);
    e.prepareChemistry(f.content,car+'.mol');e.update({lattice_scale:2});const before=[],after=[];
    for(let i=0;i<8;i++)for(let j=0;j<8;j++) {
      const p=f.positions.map((v,k)=>v+(k%3===0?i*meta.surface.period[0]/8:k%3===1?j*meta.surface.period[1]/8:0));
      e.fourierEnabled(false);before.push(e.energyForces(p,1,.6,false).surface_energy);
      e.fourierEnabled(true);after.push(e.energyForces(p,1,.6,false).surface_energy);
    }
    const span=a=>Math.max(...a)-Math.min(...a);
    report.barriers.push({car,reference_span:span(before),proxy_span:span(after),
      max_variation_error:Math.max(...after.map((v,i)=>Math.abs(v-after[0]-before[i]+before[0])))});
  }
  const folder=path.join(root,'outputs/wasm-fourier-validation/races/browser-trials');
  const f=fixtures.find(f=>f.car==='fullerene-c60');e.prepareChemistry(f.content,'fullerene-c60.mol');
  for(const file of (await fs.readdir(folder)).filter(f=>f.includes('fullerene-c60')&&f.endsWith('.coords.json.gz'))) {
    const trajectory=JSON.parse(zlib.gunzipSync(await fs.readFile(path.join(folder,file))));
    const scale=file.includes('L1.5-')?1.5:2;e.update({lattice_scale:scale});
    for(let i=0;i<trajectory.positions.length;i+=10) {
      const p=trajectory.positions[i];e.fourierEnabled(true);const a=e.energyForces(p,1,.6,false),counts=e.proxyCounts();
      e.fourierEnabled(false);const b=e.energyForces(p,1,.6,false);
      report.trajectory.push({file,step:trajectory.steps[i],counts,energy_error:Math.abs(a.surface_energy-b.surface_energy),
        reference_max_force:Math.max(...b.surface_forces.map(Math.abs)),minimum_z:Math.min(...p.filter((_,i)=>i%3===2)),
        max_force_error:maxdiff(a.surface_forces,b.surface_forces),rms_force_error:rmsdiff(a.surface_forces,b.surface_forces)});
    }
  }
  await checkSources(report.sources);report.status='passed';
}catch(error){report.status='failed';report.error=error.stack||String(error);process.exitCode=1;}
finally{e.close();e.delete();await fs.writeFile(path.join(out,'wasm-validation.json'),JSON.stringify(report,null,2)+'\n');}
console.log(JSON.stringify({status:report.status,static:report.static.length,trajectory:report.trajectory.length,fd:report.finite_difference,barriers:report.barriers,error:report.error}));
