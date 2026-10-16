import fs from 'node:fs';import path from 'node:path';import z from 'node:zlib';import assert from 'node:assert/strict';import {pathToFileURL} from 'node:url';
const root=process.cwd(),folder=path.join(root,'outputs/brotli-validation'),assets=path.join(root,'wasm/dist/browser/assets');
const meta=JSON.parse(fs.readFileSync(path.join(assets,'model.json'))),atlas=JSON.parse(fs.readFileSync(path.join(assets,'fourier-atlas.json')));
const load=d=>{const b=z.gunzipSync(fs.readFileSync(path.join(root,'wasm/dist/browser',d.path)));return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));};
const tables=['tile','saip-grid','hetero-grid'].map(k=>load(meta.surface.arrays[k])),carbon=load(atlas.arrays[0]);
const old=await(await import(pathToFileURL(path.join(folder,'baseline/molecular-car.mjs')))).default({print:()=>{},printErr:()=>{}}),fresh=await(await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')))).default({print:()=>{},printErr:()=>{}});
const results=[];
for(const [file,hetero] of [['examples/benzene.mol',false],['examples/tetraphenylbenzene.mol',false],['tests/data/workshop/ethanol.mol',true],['outputs/wasm-fourier-validation/races/inputs/fullerene-c60.mol',false]]){
 if(!fs.existsSync(file))throw Error(file);
 const content=fs.readFileSync(file,'utf8'),a=new old.BrowserEngine(meta,...tables),b=new fresh.BrowserEngine(meta,tables[0],tables[1],hetero?tables[2]:new Float64Array());
 try{for(const e of [a,b]){e.installFourier(atlas);e.loadFourierSpecies(0,carbon);e.prepareWithPose(content,path.basename(file),{...meta.settings,timestep_fs:1,drive_strength:meta.drive.strength,field_growth:meta.drive.growth},{mode:'auto'});}
 assert.deepEqual(b.frame(),a.frame());assert.deepEqual(b.phaseSpace(),a.phaseSpace());
 const coords=a.phaseSpace().positions;assert.deepEqual(b.energyForces(coords,1,.6,false),a.energyForces(coords,1,.6,false));
 for(let i=0;i<100;i++)assert.deepEqual(b.advance(50),a.advance(50));
 let finish=null;if(file==='examples/benzene.mol'){for(let i=0;i<5000;i++){const fa=a.advance(50),fb=b.advance(50);assert.deepEqual(fb,fa);if(fa.progress>=1){finish={timeFs:fa.time_fs,progress:fa.progress};break;}}assert(finish,'benzene did not finish');}
 results.push({finish,file,hetero,preparation:'bit-identical',energyForces:'bit-identical',trajectorySteps:5000,trajectory:'bit-identical'});console.log('PASS',file);
 }finally{a.close();a.delete();b.close();b.delete();}
}
fs.writeFileSync(path.join(folder,'equivalence.json'),JSON.stringify({passed:true,results},null,2));
