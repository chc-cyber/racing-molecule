// Only same-origin, direct-parent commands are accepted. UI operations use the
// same validation and preparation path as the race's own buttons.
export function installRaceBridge({load,request,snapshot,preferences,pause}){
  if(parent===window)return;
  const send=data=>parent.postMessage({channel:'racing-molecule',...data},location.origin);
  let loadSerial=0;
  addEventListener('message',async event=>{
    if(event.origin!==location.origin||event.source!==parent||event.data?.channel!=='racing-molecule')return;
    const {type,id,data}=event.data;
    try {
      if(type==='load-race'){
        const serial=++loadSerial;
        await load(data);
        const started=performance.now();
        while(performance.now()-started<170000){
          if(serial!==loadSerial)return;
          const state=await request('/api/state');
          if(state.status==='error')throw new Error(state.error||'赛车准备失败。');
          if(state.status==='paused'&&state.frame&&state.molecule){send({type:'race-loaded',id,generation:state.generation});return;}
          await new Promise(resolve=>setTimeout(resolve,150));
        }
        throw new Error('赛车准备超时，请检查结构后重试。');
      }
      if(type==='cancel-load'){++loadSerial;await request('/api/cancel-load');}
      if(type==='pause-race'&&snapshot()?.status==='running')await (pause ? pause() : request('/api/action',{action:'pause'}));
      if(type==='preferences')preferences(data);
      if(type==='screen-visibility')document.documentElement.dataset.appActive=String(event.data.active);
    } catch(error){send({type:'race-loaded',id,error:error.message});}
  });
  send({type:'race-ready'});
  const sendState=()=>{
    const state=snapshot();if(!state)return;
    send({type:'race-status',state:{status:state.status,generation:state.generation,started:state.started,unsavedReplay:state.unsavedReplay,competition:state.competition,frame:{step:state.frame?.step||0}}});
  };
  addEventListener('race-state-change',sendState);sendState();
  document.addEventListener('click',event=>{
    const link=event.target.closest('a[data-app-route]');if(!link)return;
    event.preventDefault();send({type:'navigate',route:link.dataset.appRoute});
  });
}
