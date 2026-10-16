export async function downloadTool(kind,container){
  container.querySelector(':scope > span')?.setAttribute('hidden','');
  for(const node of [...container.childNodes])if(node.nodeType===3)node.remove();
  container.querySelector('.tool-download')?.remove();
  const label=document.createElement('span'),bar=document.createElement('progress');bar.max=100;bar.value=0;
  label.textContent=kind==='2d'?'下载二维工具 0%':'下载三维工具 0%';bar.setAttribute('aria-label',label.textContent);
  const status=document.createElement('div');status.className='tool-download';status.append(label,bar);container.prepend(status);
  try{
    const response=await fetch(new URL('./tool-assets.json',import.meta.url));if(!response.ok)return;
    const files=(await response.json())[kind],counts=files.map(()=>0),total=files.reduce((n,f)=>n+f.bytes,0);
    let cursor=0;
    const next=async()=>{while(cursor<files.length){const index=cursor++,file=files[index];
      const r=await fetch(new URL(file.path,import.meta.url),{cache:'force-cache'});if(!r.ok)throw Error('工具下载失败，请重试。');
      const reader=r.body.getReader();for(;;){const {done,value}=await reader.read();if(done)break;counts[index]+=value.byteLength;
        const percent=Math.min(100,Math.floor(counts.reduce((a,b)=>a+b,0)/total*100));bar.value=percent;label.textContent=`${kind==='2d'?'下载二维工具':'下载三维工具'} ${percent}%`;
      }
    }};
    await Promise.all([next(),next()]);bar.value=100;label.textContent='下载完成，正在初始化…';
  }catch(error){status.remove();throw error;}
}
