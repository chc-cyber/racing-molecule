// Unknown/unsupported records conservatively request the full tables.
export function needsHetero(content){
  const lines=String(content).split(/\r?\n/);
  if(lines[3]?.includes('V2000')){
    const count=Number(lines[3].slice(0,3));
    if(!Number.isInteger(count)||count<1||lines.length<4+count)return true;
    return lines.slice(4,4+count).some(line=>!['C','H'].includes(line.slice(31,34).trim()));
  }
  if(lines[3]?.includes('V3000')){
    const start=lines.findIndex(l=>l.trim()==='M  V30 BEGIN ATOM'),end=lines.findIndex(l=>l.trim()==='M  V30 END ATOM');
    if(start<0||end<=start+1)return true;
    return lines.slice(start+1,end).some(line=>!/^M  V30\s+\d+\s+[CH]\s/.test(line));
  }
  return true;
}
