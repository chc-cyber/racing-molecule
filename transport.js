// The renderer has one interface; native HTTP and browser-local MD supply it.
const browserMode = document.documentElement.dataset.engine === 'wasm';
let browser;
if (browserMode) {
  const {BrowserTransport} = await import('./browser/transport.mjs');
  browser = new BrowserTransport();
}

export async function request(path, data) {
  if (browser) return browser.request(path, data);
  const response = await fetch(path, data === undefined ? {cache: 'no-store'} : {
    method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(data),
  });
  let result;
  try { result = await response.json(); }
  catch { throw new Error('本地服务没有返回可读的数据。'); }
  if (!response.ok || result.error && result.status !== 'error') throw new Error(result.error || result.message || `请求失败（${response.status}）`);
  return result;
}

export async function downloadReplay() {
  const link = document.createElement('a');
  const localURL = browser ? URL.createObjectURL(new Blob([JSON.stringify(await browser.replay())], {type: 'application/json'})) : null;
  link.href = localURL || '/api/replay?download=1';
  link.download = 'molecular-car-replay.json';
  document.body.appendChild(link); link.click(); link.remove();
  if (localURL) setTimeout(() => URL.revokeObjectURL(localURL), 1000);
}

export const execution = browserMode ? 'browser-wasm' : 'native-http';
