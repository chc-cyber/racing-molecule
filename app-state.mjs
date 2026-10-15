// Display preferences only. Competition rules belong to the race screen.
const key='racing-molecule-preferences-v2';
export const defaults=Object.freeze({representation:'ball-stick',quality:'balanced',motion:true});
function displayPreferences(data={}) {
  return {representation:['ball-stick','space-filling'].includes(data.representation)?data.representation:data.representation==='space-fill'?'space-filling':defaults.representation,
    quality:['balanced','high','low'].includes(data.quality)?data.quality:defaults.quality,
    motion:typeof data.motion==='boolean'?data.motion:defaults.motion};
}
export function readPreferences(storage=globalThis.localStorage) {
  try {
    const data=JSON.parse(storage.getItem(key)||'null')||JSON.parse(storage.getItem('molecular-car-ui-preferences-v1')||'{}');
    return displayPreferences(data);
  } catch { return {...defaults}; }
}
export function writePreferences(data,storage=globalThis.localStorage) {
  const result=displayPreferences({...readPreferences(storage),...data});
  storage.setItem(key,JSON.stringify(result));return result;
}
