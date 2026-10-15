import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {out,read,initialize,sources,checkSources} from './atlas_support.mjs';
const {engine:e}=await initialize(), source=await sources();
const f=(await read(path.join(out,'native-fixtures.json'))).find(r=>r.car==='benzene-12'&&r.label==='start');
const info=e.prepareChemistry(f.content,'benzene-12.mol');
const center=info.normal_counts.findIndex(n=>n===3),ids=info.normal_neighbours.slice(center*3,center*3+3);
const rows=[];
try {
  e.update({lattice_scale:2});
  for(const area of [0,1e-6,2e-5,5e-5,2e-4]) {
    const p=f.positions.slice();
    // Retain preparation-time normal topology, then flatten just one local plane.
    // Nonzero distances avoid singular molecular bonds in this diagnostic.
    for(const [k,xyz] of [[0,[1,0,3.3]],[1,[2,area,3.3]],[2,[0,0,3.3]]])
      p.splice(ids[k]*3,3,...xyz);
    e.fourierEnabled(true);const a=e.energyForces(p,1,.6,false),counts=e.proxyCounts();
    e.fourierEnabled(false);const exact=e.energyForces(p,1,.6,false);e.fourierEnabled(true);
    let maxError=0,maxForce=Math.max(...a.surface_forces.map(Math.abs));const h=1e-10;
    for(let i=0;i<p.length;i++) {
      const v=p[i];p[i]=v+h;const plus=e.energyForces(p,1,.6,false).surface_energy;
      p[i]=v-h;const minus=e.energyForces(p,1,.6,false).surface_energy;p[i]=v;
      maxError=Math.max(maxError,Math.abs(-(plus-minus)/(2*h)-a.surface_forces[i]));
    }
    assert(a.surface_forces.every(Number.isFinite));
    assert(maxError/Math.max(1,maxForce)<1e-4,'degenerate FD '+area+' '+maxError);
    rows.push({area,counts,fd_error:maxError,max_force:maxForce,relative_fd_error:maxError/Math.max(1,maxForce),
      approximation_force_error:Math.max(...a.surface_forces.map((v,i)=>Math.abs(v-exact.surface_forces[i])))});
  }
  await checkSources(source);
  await fs.writeFile(path.join(out,'degenerate-validation.json'),JSON.stringify({status:'passed',sources:source,rows},null,2)+'\n');
  console.log(JSON.stringify(rows));
}finally{e.close();e.delete();}
