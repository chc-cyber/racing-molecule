// High-temperature numerical/stability guard; this is not a timing benchmark.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {root,out,read,initialize,sources,checkSources,hash} from './atlas_support.mjs';
const {engine:e,meta,module}=await initialize();
const cars={
  'fullerene-c60':'outputs/wasm-fourier-validation/races/inputs/fullerene-c60.mol',
  'folded-eicosane':'outputs/md-release-validation/inputs/folded-eicosane.mol',
  'pyridine-11':'outputs/wasm-fourier-validation/races/inputs/pyridine-11.mol',
};
const report={status:'running',sources:await sources(),protocol:{temperature_K:600,timestep_fs:.5,
  lattice_scales:[1.5,2],nve_ps:2,driven_ps:5,sample_fs:10,damping_fs:350,drive_strength:.6,
  phase_space:'same periodic-relaxed positions and seeded 600 K velocities within each car/scale pair',
  nve_comparison:'each Hamiltonian own total-energy drift; no thermostat/drive',
  purpose:'short numerical and structural guard, not full race equivalence or performance test'},rows:[]};
for(const file of [...Object.values(cars),'wasm/tools/validate_atlas_600k.mjs'])report.sources[file]=await hash(path.join(root,file));
const save=()=>fs.writeFile(path.join(out,'temperature-600-validation.json'),JSON.stringify(report,null,2)+'\n');
try {
  const gate=await read(path.join(out,'wasm-validation.json'));assert.equal(gate.status,'passed');await checkSources(gate.sources);
  for(const [car,file] of Object.entries(cars))for(const scale of [1.5,2]) {
    const content=await fs.readFile(path.join(root,file),'utf8');
    const info=e.prepareWithPose(content,car+'.mol',{...meta.settings,surface:'periodic-au',temperature:600,damping_fs:350,lattice_scale:scale,timestep_fs:.5},{mode:'auto'});
    assert(info.runtime.preparation.accepted,'preparation convergence');
    const initial=e.phaseSpace(),bonds=info.molecule.bonds.map(b=>[b.a-1,b.b-1]);
    const length=(p,b)=>Math.hypot(...[0,1,2].map(a=>p[3*b[0]+a]-p[3*b[1]+a]));
    const reference=bonds.map(b=>length(initial.positions,b));
    for(const mode of ['periodic-nve','atlas-nve','atlas-driven']) {
      const row={car,lattice_scale:scale,mode,status:'running',samples:[]};
      try {
        e.fourierEnabled(mode!=='periodic-nve');e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(.5);
        const conservative=mode.endsWith('nve');e.conservative(conservative);
        if(!conservative)e.update({temperature:600,drive:{mode:'gradient',strength:.6,growth:1}});
        const first=e.refresh();
        row.initial_total_energy=first.energy;
        if(conservative)row.initial_kinetic_energy=first.energy-e.energyForces(first.positions,1,0,false).energy;
        let minRatio=1,maxRatio=1,maxDrift=0,minZ=Infinity,maxZ=-Infinity,lift=0,maxLift=0;
        const counts=[0,0,0],samples=(conservative?2:5)*100;
        for(let k=1;k<=samples;k++) {
          const f=e.advance(20);assert.equal(f.time_fs,k*10);
          assert(f.positions.every(Number.isFinite)&&Number.isFinite(f.energy),'nonfinite propagation');
          const ratios=bonds.map((b,i)=>length(f.positions,b)/reference[i]);
          minRatio=Math.min(minRatio,...ratios);maxRatio=Math.max(maxRatio,...ratios);
          const z=f.positions.filter((_,i)=>i%3===2),low=Math.min(...z),high=Math.max(...z);
          minZ=Math.min(minZ,low);maxZ=Math.max(maxZ,high);lift=low>8?lift+1:0;maxLift=Math.max(maxLift,lift);
          maxDrift=Math.max(maxDrift,Math.abs(f.energy-first.energy));
          for(let i=0;i<3;i++)counts[i]+=f.fourier_counts[i];
          row.samples.push({time_ps:f.time_fs/1000,energy_delta:f.energy-first.energy,internal_temperature_K:f.internal_temperature,
            center_nm:f.center.map(v=>v/10),min_z_A:low,max_z_A:high,bond_min:Math.min(...ratios),bond_max:Math.max(...ratios)});
          assert(minRatio>=.65&&maxRatio<=1.45,'bond ratio stability guard');
        }
        Object.assign(row,{status:'complete',bond_ratio_range:[minRatio,maxRatio],height_range_A:[minZ,maxZ],
          maximum_lift_duration_ps:maxLift*.01,counts,stability_guard_pass:true});
        if(conservative)Object.assign(row,{max_energy_error:maxDrift,relative_to_initial_kinetic:maxDrift/row.initial_kinetic_energy,
          numerical_screen_pass:maxDrift/row.initial_kinetic_energy<.01});
        else row.numerical_screen_pass=maxLift<500;
      }catch(error){row.status='failed';row.stability_guard_pass=false;try{row.error=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);}catch{row.error=String(error.stack||error);}}
      report.rows.push(row);await save();console.log(JSON.stringify({car,scale,mode,status:row.status,pass:row.numerical_screen_pass,drift:row.relative_to_initial_kinetic,error:row.error}));
      e.conservative(false);
    }
  }
  await checkSources(report.sources);report.status=report.rows.every(r=>r.status==='complete')?'complete':'failed';
}catch(error){report.status='failed';report.error=String(error.stack||error);process.exitCode=1;}
finally{e.close();e.delete();await save();}
