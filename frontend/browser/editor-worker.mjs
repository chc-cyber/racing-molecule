import createModule from './runtime/molecular-car.mjs';
let modulePromise;
self.onmessage=async({data})=>{
  let module;
  try {
    module=await(modulePromise??=createModule({print:()=>{},printErr:()=>{}}));
    const result=module.workshopRequest(data.operation,data.payload);
    postMessage({id:data.id,result});
  } catch(error) {
    let message=error?.message||String(error);
    if(module)try {const detail=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);message=Array.isArray(detail)?detail.join(': '):String(detail);} catch {}
    postMessage({id:data.id,error:message});
  }
};
