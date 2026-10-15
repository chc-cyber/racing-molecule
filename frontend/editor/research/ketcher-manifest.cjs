/* Generates or checks the offline vendor inventory; performs no network calls. */
const fs=require('node:fs');
const path=require('node:path');
const crypto=require('node:crypto');
const assert=require('node:assert/strict');
const root=path.resolve(__dirname,'../../..');
const vendor=path.join(root,'frontend/editor/vendor/ketcher');
const manifestPath=path.join(root,'frontend/editor/vendor/ketcher-manifest.json');
const sha=file=>crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
function inventory(dir){
  return fs.readdirSync(dir,{withFileTypes:true}).flatMap(entry=>{
    const absolute=path.join(dir,entry.name);
    return entry.isDirectory()?inventory(absolute):[{path:path.relative(vendor,absolute).split(path.sep).join('/'),bytes:fs.statSync(absolute).size,sha256:sha(absolute)}];
  }).sort((a,b)=>a.path.localeCompare(b.path));
}
if(process.argv.includes('--write')){
  const archivePath=path.join(root,'outputs/ketcher-deploy/ketcher-standalone-3.18.0.zip');
  const archiveSha='484e7f10a0e74808ae5f43f0e6448db6f4d884bdd2189aa271108474167a3e7e';
  assert.equal(sha(archivePath),archiveSha,'Official release archive digest must match');
  const files=inventory(vendor);
  const manifest={
    component:'Ketcher',version:'3.18.0',distribution:'official prebuilt standalone',
    vendoredAt:new Date().toISOString(),license:'Apache-2.0; upstream dependency notices retained',
    source:{repository:'https://github.com/epam/ketcher',release:'https://github.com/epam/ketcher/releases/tag/v3.18.0',tag:'v3.18.0',tagObject:'5f0a0cdd3a1564f3a43f387d3a76e571ae7d05f1',commit:'f4f013e1855b71b0cbd37fc79d55b0cbe51c2371',publishedAt:'2026-09-02T13:38:21Z',prerelease:false},
    archive:{url:'https://github.com/epam/ketcher/releases/download/v3.18.0/ketcher-standalone-3.18.0.zip',bytes:35057618,sha256:archiveSha,officialGitHubDigest:'sha256:'+archiveSha},
    entry:'/editor/vendor/ketcher/index.html',
    deployment:{npmInstallRequired:false,externalCDNRequired:false,originalReleaseFilesModified:false,originalReleaseFileCount:30,unpackedOriginalBytes:120764113,notes:'Full prebuilt ZIP extracted unchanged. LICENSE and NOTICE added unchanged from the matching upstream tag because release ZIP omitted root copies. Bundled static/js/*.LICENSE.txt preserved unchanged. Indigo in-browser code is included in the official standalone JavaScript chunks.'},
    licenses:[{path:'LICENSE',url:'https://raw.githubusercontent.com/epam/ketcher/v3.18.0/LICENSE'},{path:'NOTICE',url:'https://raw.githubusercontent.com/epam/ketcher/v3.18.0/NOTICE'}],
    verifiedRuntime:{report:'outputs/ketcher-deploy/browser-probe.json',probe:'frontend/editor/research/probe-ketcher-browser.cjs',sameOriginIframe:true,externalRequests:0,pageErrors:0,changeSubscription:"ketcher.editor.subscribe('change', callback)",unsubscribe:'ketcher.editor.unsubscribe(eventName, tokenReturnedBySubscribe)',selectionSubscription:"editor.subscribe('selectionChange', () => editor.selection())",structure:'ketcher.editor.struct(); atoms/bonds are Pool instances extending Map',coordinateConvention:'live pp.x = KET location[0]; live pp.y = -KET location[1]',importOptions:{preserveCanvasPosition:true,needZoom:false},identityLimit:'Pool IDs are session IDs, may be nonzero and are not preserved across rebuilding. Single mol0 KET import ordinal and explicit bond endpoint references verified for benzene and C60; workshop adapter must bind ordinals to stable IDs explicitly.'},
    files,
  };
  fs.writeFileSync(manifestPath,JSON.stringify(manifest,null,2)+'\n');
  console.log(`Wrote fixed vendor manifest: ${files.length} files`);
}else{
  const manifest=JSON.parse(fs.readFileSync(manifestPath,'utf8'));
  assert.deepEqual(inventory(vendor),manifest.files,'Vendored bytes must match the pinned manifest');
  console.log(`Verified Ketcher ${manifest.version}: ${manifest.files.length} local file hashes`);
}
