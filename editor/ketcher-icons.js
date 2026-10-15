/* Ketcher visual assets, not Ketcher's editor state or action pipeline.
 * Source snapshots, copyrights and SHA-256: vendor/ketcher-icons/manifest.json.
 * Ketcher originals Copyright EPAM Systems, Apache-2.0; see LICENSE / NOTICE.
 */
import {iconSources} from './vendor/ketcher-icons/icons.js';

const cache=new Map();
const provenance=name=>name.startsWith('workshop-')?'workshop':name==='template-8'?'molecular-car-c60-extension':'ketcher';

/** Update only an icon placeholder, retaining the owning control and listeners. */
export function setKetcherIcon(element,name){
  if(!element||!iconSources[name])return false;
  element.dataset.ketcherIcon=name;
  if(element.firstElementChild?.tagName.toLowerCase()==='svg'&&element.firstElementChild.dataset.ketcherIcon===name)return false;
  if(!cache.has(name)){
    const parsed=new DOMParser().parseFromString(iconSources[name],'image/svg+xml');
    const svg=parsed.documentElement;
    if(svg.localName!=='svg'||parsed.querySelector('parsererror'))throw new Error(`Invalid local icon: ${name}`);
    // All files are fixed local assets. Normalize solid monochrome colors only;
    // preserve none, geometry, viewBox, fill rules and stroke widths.
    for(const node of [svg,...svg.querySelectorAll('*')])for(const attr of ['fill','stroke']){
      const value=node.getAttribute(attr);
      if(value&&value!=='none'&&!value.startsWith('url('))node.setAttribute(attr,'currentColor');
    }
    svg.setAttribute('width','24');svg.setAttribute('height','24');svg.setAttribute('aria-hidden','true');svg.setAttribute('focusable','false');
    svg.classList.add('ketcher-svg');svg.dataset.ketcherIcon=name;svg.dataset.iconSource=provenance(name);
    cache.set(name,svg);
  }
  element.replaceChildren(element.ownerDocument.importNode(cache.get(name),true));
  return true;
}

/** Synchronous, idempotent, offline. Call after the workshop's state render. */
export function renderKetcherIcons(root=document){
  let updated=0;
  const placeholders=[...(root.matches?.('[data-ketcher-icon]')?[root]:[]),...root.querySelectorAll('[data-ketcher-icon]:not(svg)')];
  for(const node of placeholders)if(setKetcherIcon(node,node.dataset.ketcherIcon))updated++;
  const find=id=>root.id===id?root:root.querySelector(`#${id}`);
  const bond=find('bond-tool-symbol'),order=find('bond-order');
  if(bond&&order&&setKetcherIcon(bond,{'1':'bond-single','2':'bond-double','3':'bond-triple'}[order.value]||'bond-single'))updated++;
  const current=find('element');
  if(current)current.classList.toggle('common-element',['H','C','N','O','S','P','F','Cl','Br','I'].includes(current.value));
  return updated;
}
