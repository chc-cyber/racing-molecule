// Paired phase-space NVE screening; the stochastic bath and drive are disabled.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {root,folder,cars,initialize,sources,checkSources,read,save,stepsFor,geometry,hash} from './acceptance_support.mjs';
const {engine:e,meta,module}=await initialize();
const selected=(process.env.TIMESTEP_CARS||Object.keys(cars).join(',')).split(',');
const duration=Number(process.env.TIMESTEP_NVE_PS||2);
const report={status:'running',sources:await sources(),protocol:{duration_ps:duration,dt_fs:[.5,1],sample_fs:"every integration step",temperatures_K:[150],lattice_scales:[1.5],surfaces:['fourier-au-v2'],phase_space:'identical periodic-relaxed coordinates and seeded velocities within each temperature/scale/car pair',thermostat:false,drive:false,screen_gate:'finite/bond ratio .65..1.45 and max |E-E0| < 1% initial kinetic energy; provisional numerical screen, not race approval'},rows:[]};
report.sources[path.join(root,'wasm/tools/acceptance_nve.mjs')]=await hash(path.join(root,'wasm/tools/acceptance_nve.mjs'));
const target=path.join(folder,'nve.json');
try{
  for(const car of selected)for(const temperature of [150])for(const scale of [1.5]){
    const content=await fs.readFile(path.join(root,cars[car]),'utf8');
    const info=e.prepareWithPose(content,car+'.mol',{...meta.settings,surface:'periodic-au',temperature,lattice_scale:scale,timestep_fs:.5},{mode:'auto'});
    assert(info.runtime.preparation.accepted);const initial=e.phaseSpace(),measure=geometry(info,initial.positions);
    for(const surface of ['fourier-au'])for(const dt of [.5,1]){
      const row={car,temperature_K:temperature,lattice_scale:scale,surface,dt_fs:dt,status:'running',samples:[]};
      try{
        e.fourierEnabled(surface==='fourier-au');e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(dt);e.conservative(true);
        const first=e.refresh(),potential=e.energyForces(first.positions,1,0,false).energy;
        row.initial_total_energy=first.energy;row.initial_kinetic_energy=first.energy-potential;
        assert(Number.isFinite(row.initial_kinetic_energy)&&row.initial_kinetic_energy>0,'invalid initial kinetic energy');
        const began=performance.now();let maxDrift=0,lo=1,hi=1;
        for(let i=1;i<=Math.round(duration*1000/dt);i++){
          const f=e.advance(1),g=measure(f);assert.equal(f.time_fs,i*dt);
          maxDrift=Math.max(maxDrift,Math.abs(f.energy-first.energy));lo=Math.min(lo,g.bond_min);hi=Math.max(hi,g.bond_max);
          row.samples.push({time_ps:f.time_fs/1000,energy_delta:f.energy-first.energy,internal_temperature_K:f.internal_temperature,...g});
          assert(lo>=.65&&hi<=1.45,'bond guard');
        }
        row.wall_seconds=(performance.now()-began)/1000;row.max_energy_error=maxDrift;row.relative_to_initial_kinetic=maxDrift/row.initial_kinetic_energy;row.bond_ratio_range=[lo,hi];
        const n=row.samples.length,meanT=row.samples.reduce((s,v)=>s+v.time_ps,0)/n,meanE=row.samples.reduce((s,v)=>s+v.energy_delta,0)/n;
        row.energy_linear_slope_kcal_mol_ps=row.samples.reduce((s,v)=>s+(v.time_ps-meanT)*(v.energy_delta-meanE),0)/row.samples.reduce((s,v)=>s+(v.time_ps-meanT)**2,0);
        row.first_half_max_error=Math.max(...row.samples.slice(0,n/2).map(v=>Math.abs(v.energy_delta)));
        row.second_half_max_error=Math.max(...row.samples.slice(n/2).map(v=>Math.abs(v.energy_delta)));
        row.stability_guard_pass=true;row.numerical_screen_pass=row.relative_to_initial_kinetic<.01;row.status='complete';
      }catch(err){row.status='failed';row.stability_guard_pass=false;row.numerical_screen_pass=false;try{row.error=module.getExceptionMessage(err);module.decrementExceptionRefcount(err);}catch{row.error=String(err.stack||err);}}
      report.rows.push(row);await save(target,report);console.log(JSON.stringify({car,temperature,scale,surface,dt,status:row.status,drift:row.relative_to_initial_kinetic,pass:row.numerical_screen_pass}));
    }
  }
  await checkSources(report.sources);report.status='complete';
}catch(err){report.status='failed';report.error=String(err.stack||err);process.exitCode=1;}
finally{e.close();e.delete();await save(target,report);}
