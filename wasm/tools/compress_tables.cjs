const fs=require('node:fs'),path=require('node:path'),z=require('node:zlib'),crypto=require('node:crypto');
const root=path.resolve(__dirname,'../dist/browser'),hash=b=>crypto.createHash('sha256').update(b).digest('hex');
for(const name of ['model.json','fourier-atlas.json']){
  const file=path.join(root,'assets',name),meta=JSON.parse(fs.readFileSync(file));
  for(const d of name==='model.json'?Object.values(meta.surface.arrays):meta.arrays){
    if(d.download_bytes<1024*1024)continue;
    const raw=z.gunzipSync(fs.readFileSync(path.join(root,d.path)));
    const packed=z.brotliCompressSync(raw,{params:{[z.constants.BROTLI_PARAM_QUALITY]:9}});
    if(!z.brotliDecompressSync(packed).equals(raw)||hash(raw)!==d.sha256)throw Error('Table mismatch');
    const chunks=[];
    for(let offset=0;offset<packed.length;offset+=1024*1024){
      const bytes=packed.subarray(offset,offset+1024*1024),sha256=hash(bytes),relative='assets/chunks/'+sha256+'.bin';
      fs.mkdirSync(path.dirname(path.join(root,relative)),{recursive:true});fs.writeFileSync(path.join(root,relative),bytes);
      chunks.push({path:relative,bytes:bytes.length,download_bytes:bytes.length,sha256});
    }
    d.brotli={path:d.path.replace(/\.gz$/,'.br'),download_bytes:packed.length,chunks};
    fs.writeFileSync(path.join(root,d.brotli.path),packed);
  }
  fs.writeFileSync(file,JSON.stringify(meta,null,2)+'\n');
}
