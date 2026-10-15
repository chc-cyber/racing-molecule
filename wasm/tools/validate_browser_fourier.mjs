import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {root,out,read,initialize,sources,checkSources} from './fourier_support.mjs';
const {engine:e}=await initialize();
const report={status:'running',sources:await sources(),static:[],finite_difference:[],barriers:[]};
const fixtures=await read(path.join(out,'native-fixtures.json'));
const maxdiff=(a,b)=>Math.max(...a.map((v,i)=>Math.abs(v-b[i])));
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
      approximation_force_error:maxdiff(current.surface_forces,original.surface_forces),
      approximation_energy_error:Math.abs(current.surface_energy-original.surface_energy)});
    if(f.scale===2&&((car==='benzene-12'&&f.label!=='start')||(f.label==='start'))) {
      const p=f.positions.slice(),h=1e-6;let error=0;
      for(let k=0;k<p.length;k++) {
        const old=p[k];p[k]=old+h;const a=e.energyForces(p,1,.6,false).surface_energy;
        p[k]=old-h;const b=e.energyForces(p,1,.6,false).surface_energy;p[k]=old;
        error=Math.max(error,Math.abs(-(a-b)/(2*h)-current.surface_forces[k]));
      }
      assert(error<5e-5,'finite difference '+car+' '+f.label+' '+error);
      report.finite_difference.push({car,label:f.label,error});
    }
  }
  // Lateral energy variations, not just a potentially large adsorption offset.
  for(const car of ['benzene-12','test-car-52','fullerene-c60','pyridine-11']) {
    const f=fixtures.find(f=>f.car===car&&f.label==='start'&&f.scale===2);
    e.prepareChemistry(f.content,car+'.mol');e.update({lattice_scale:2});const before=[],after=[];
    for(let i=0;i<6;i++)for(let j=0;j<6;j++) {
      const p=f.positions.map((v,k)=>v+(k%3===0?i*2.884995667241114/6:k%3===1?j*4.996959075277683/6:0));
      e.fourierEnabled(false);before.push(e.energyForces(p,1,.6,false).surface_energy);
      e.fourierEnabled(true);after.push(e.energyForces(p,1,.6,false).surface_energy);
    }
    const span=a=>Math.max(...a)-Math.min(...a);
    report.barriers.push({car,reference_span:span(before),proxy_span:span(after),
      max_variation_error:Math.max(...after.map((v,i)=>Math.abs(v-after[0]-before[i]+before[0])))});
  }
  await checkSources(report.sources);report.status='passed';
}catch(error){report.status='failed';report.error=error.stack||String(error);process.exitCode=1;}
finally{e.close();e.delete();await fs.writeFile(path.join(out,'numerics.json'),JSON.stringify(report,null,2)+'\n');}
console.log(JSON.stringify({status:report.status,static:report.static.length,fd:report.finite_difference,error:report.error}));
