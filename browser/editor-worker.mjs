import {loadRuntime} from './runtime-loader.mjs';
let modulePromise;
self.onmessage=async({data})=>{
  let module;
  try {
    module=await(modulePromise??=loadRuntime(percent=>postMessage({id:data.id,kind:'download',percent})));
    postMessage({id:data.id,kind:'download',percent:null});
    const result=module.workshopRequest(data.operation,data.payload);
    postMessage({id:data.id,result});
  } catch(error) {
    let message=error?.message||String(error);
    if(module)try {const detail=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);message=Array.isArray(detail)?detail.join(': '):String(detail);} catch {}
    postMessage({id:data.id,error:message});
  }
};
