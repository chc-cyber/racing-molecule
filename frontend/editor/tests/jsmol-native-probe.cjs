/* Minimal real upstream popup/source bridge probe; no server/API side effects. */
const {launch}=require('./jsmol-live.cjs');
const fs=require('node:fs'),path=require('node:path');
(async()=>{
  const {browser,page,errors}=await launch();
  try{
    await page.goto('http://127.0.0.1:8765/editor/research/jsmol-host-live.html');
    await page.waitForFunction(()=>window.host?.ready,null,{timeout:60000});
    await page.evaluate(async()=>{host.nativeActionLog=[];host.nativeModeLog=[];host.callbacks.onAction=name=>host.nativeActionLog.push(name);host.callbacks.onMode=(mode,options)=>host.nativeModeLog.push({mode,options});const p={format:'molecular-car-workshop-v1',name:'Probe',atoms:[{id:'a-c1',element:'C',charge:0,pos2:[0,0],pos3:[0,0,0]},{id:'a-c2',element:'C',charge:0,pos2:[1.5,0],pos3:[1.5,0,0]}],bonds:[{id:'b-c1c2',a:'a-c1',b:'a-c2',order:1}],fragments:[],has3d:true};model.replace(p);await host.settled();});
    await page.mouse.click(20,20);await page.waitForTimeout(150);
    const out=path.resolve(__dirname,'../../../outputs/editor-jsmol-native-probe');fs.mkdirSync(out,{recursive:true});
    const evidence=await page.evaluate(()=>({text:document.body.innerText,menus:[...document.querySelectorAll('[id]')].filter(e=>/menu|popup/i.test(e.id)).map(e=>({id:e.id,text:e.textContent.slice(0,900),html:e.outerHTML.slice(0,2500)})),modelkit:host.adapter.value('getProperty("modelkitInfo")'),modelKitMode:host.viewer.g.modelKitMode,bridge:!!host.viewer._mcarModelKitBridge,mode:host.mode,source:!!J.modelkit.ModelKitPopup.prototype.appRunScript.toString().includes('_mcarModelKitBridge'),errors}));
    const applet=await page.evaluate(()=>{const p=host.viewer.modelkit.menu.popupMenu.applet;return{same:p===host.applet,keys:Object.keys(p),wrapper:!!p._mcarModelKitBridge,viewer:!!p._applet?.viewer?._mcarModelKitBridge,id:p._id};});
    fs.writeFileSync(path.join(out,'probe.json'),JSON.stringify({...evidence,applet,pageErrors:errors},null,2));
    await page.screenshot({path:path.join(out,'popup.png'),fullPage:true});
    console.log(JSON.stringify({bridge:evidence.bridge,source:evidence.source,modelKitMode:evidence.modelKitMode,applet,pageErrors:errors},null,2));
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
