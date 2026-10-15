// Actual WASM execution under Node, using native-model golden fixtures.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import {fileURLToPath, pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';

const root = fileURLToPath(new URL('../../', import.meta.url));
const assets = path.join(root, 'wasm/dist/browser/assets');
const metadata = JSON.parse(await fs.readFile(path.join(assets, 'model.json'), 'utf8'));
const load = async key => {
  const data = zlib.gunzipSync(await fs.readFile(path.join(assets, `${key}.bin.gz`)));
  return new Float64Array(data.buffer.slice(data.byteOffset, data.byteOffset+data.byteLength));
};
const [tile, saip, hetero] = await Promise.all(['tile', 'saip-grid', 'hetero-grid'].map(load));
const {default: createMolecularCar} = await import(pathToFileURL(path.join(root, 'wasm/dist/browser/runtime/molecular-car.mjs')));
const module = await createMolecularCar({print: () => {}, printErr: value => console.warn(value)});
const engine = new module.BrowserEngine(metadata, tile, saip, hetero);
function text(error) {
  try {const message=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);return Array.isArray(message)?message.join(': '):JSON.stringify(message);} catch {}
  return error.message || String(error);
}

const fixtures = JSON.parse(await fs.readFile(path.join(root, 'outputs/wasm-preparation/force-fixtures.json'), 'utf8'));
const report = {static_cases: [], failures: [], precision: 'float64', timestep_fs: .5, backend: false};
let filename;
try {
  for (const fixture of fixtures) {
    if (filename !== fixture.filename) {
      filename = fixture.filename;
      const mol = await fs.readFile(path.join(assets, filename), 'utf8');
      const info = engine.prepareChemistry(mol, filename);
      assert.deepEqual(info.normal_counts, fixture.normal_counts);
      assert.deepEqual(info.normal_neighbours, fixture.normal_neighbours.flat());
      assert.deepEqual(info.molecule.atoms.map(a => a.mass), fixture.masses);
      const chargeError = Math.max(...info.molecule.atoms.map((a, i) => Math.abs(a.charge-fixture.charges[i])));
      assert(chargeError < 1e-10, `Charge mismatch: ${chargeError}`);
    }
    const result = engine.energyForces(fixture.positions, fixture.mode, .75, fixture.continuum);
    const entry = {filename, mode: fixture.mode, continuum: fixture.continuum};
    for (const component of ['', 'uff_', 'surface_', 'field_']) {
      const energyError = Math.abs(result[`${component}energy`] - fixture[`${component}energy`]);
      const forceError = Math.max(...result[`${component}forces`].map((f,i) => Math.abs(f-fixture[`${component}forces`][i])));
      entry[`${component}energy_error`] = energyError;entry[`${component}force_error`] = forceError;
      assert(energyError < 1e-7, `${filename} ${component} energy ${energyError}`);
      assert(forceError < 1e-7, `${filename} ${component} force ${forceError}`);
    }
    report.static_cases.push(entry);
  }
  report.admission=[];
  const admission = JSON.parse(await fs.readFile(path.join(root, 'outputs/wasm-preparation/admission-fixtures.json'), 'utf8'));
  for (const fixture of admission) {
    let reason;
    try {engine.prepareChemistry(fixture.content,fixture.filename);}
    catch(error) {reason=text(error);}
    assert(reason && !reason.includes('[object Object]'),`Invalid molecule accepted or error unreadable: ${fixture.filename}`);
    report.admission.push({filename:fixture.filename,rejected:true,reason});
  }
  const mol = await fs.readFile(path.join(assets, 'tetraphenylbenzene.mol'), 'utf8');
  const started = performance.now();
  const info = engine.prepare(mol, 'tetraphenylbenzene.mol', metadata.settings);
  report.preparation = info.runtime.preparation;
  report.prepare_wall_seconds = (performance.now()-started)/1000;
  assert.equal(info.molecule.atoms.length, 52);
  assert(info.runtime.preparation.accepted && info.runtime.preparation.global_force_norm_kcal_mol_A <= 1e-4);
  let frame=engine.frame();
  const propagation = performance.now();
  for(let i=0;i<100;i++)frame=engine.advance(100);
  report.propagation = {steps: frame.step, wall_seconds: (performance.now()-propagation)/1000,
    time_fs: frame.time_fs, center: frame.center, energy: frame.energy};
  assert.equal(frame.step,10000);assert.equal(frame.time_fs,5000);
  assert(frame.positions.every(Number.isFinite));
  // Preserve x/v while disabling thermostat/drive; a separate conservative gate.
  engine.conservative(true);const initial=engine.frame().energy;
  let maxDeviation=0;
  for(let i=0;i<100;i++){frame=engine.advance(100);maxDeviation=Math.max(maxDeviation,Math.abs(frame.energy-initial));}
  report.nve = {duration_fs: 5000, initial_energy: initial, maximum_absolute_energy_error: maxDeviation};
  assert(maxDeviation<1, `Unexpected conservative energy drift: ${maxDeviation}`);
  const paired = JSON.parse(await fs.readFile(path.join(root, 'outputs/wasm-preparation/native-reference.json'), 'utf8'));
  engine.setPhaseSpace(paired.nve.positions, paired.nve.velocities);
  let positionError=0,energyError=0,centerError=0;
  const nveStarted=performance.now();
  for(const reference of paired.nve.frames) {
    const current=reference.step ? engine.advance(100) : engine.frame();
    assert.equal(current.step,reference.step);
    positionError=Math.max(positionError,...current.positions.map((v,i)=>Math.abs(v-reference.positions[i])));
    centerError=Math.max(centerError,...current.center.map((v,i)=>Math.abs(v-reference.center[i])));
    energyError=Math.max(energyError,Math.abs(current.energy-reference.energy));
  }
  report.paired_nve={duration_fs:5000,maximum_position_error_A:positionError,maximum_center_error_A:centerError,
    maximum_energy_error_kcal_mol:energyError,wall_seconds:(performance.now()-nveStarted)/1000,
    native_wall_seconds:paired.nve.wall_seconds,qualification:'Identical x/v, dt and batch boundaries; no bath or drive.'};
  assert(positionError<1e-5 && centerError<1e-7 && energyError<1e-5,'Paired native/WASM NVE mismatch');
  report.passed=true;
} catch(error) {report.passed=false;report.failures.push(text(error));}
finally {engine.delete();}
const output = path.join(root, 'outputs/wasm-preparation/wasm-validation.json');
await fs.writeFile(output, JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed: report.passed, static_cases: report.static_cases.length,
  preparation: report.preparation, propagation: report.propagation, nve: report.nve, failures: report.failures, output},null,2));
if(!report.passed)process.exitCode=1;
