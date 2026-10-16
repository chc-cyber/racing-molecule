import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {needsHetero} from '../frontend/browser/needs-hetero.mjs';
test('carbon/hydrogen tables selected from MOL atoms, not titles or headers',()=>{
 const benzene=fs.readFileSync('examples/benzene.mol','utf8');assert.equal(needsHetero(benzene),false);assert.equal(needsHetero(benzene.replace('Benzene','N O molecule')),false);
 assert.equal(needsHetero(fs.readFileSync('tests/data/workshop/ethanol.mol','utf8')),true);
 assert.equal(needsHetero('malformed'),true);
});
test('V3000 selection is conservative for unusual records',()=>{
 const mol=s=>`test\n\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN ATOM\nM  V30 1 ${s} 0 0 0 0\nM  V30 END ATOM`;
 assert.equal(needsHetero(mol('C')),false);assert.equal(needsHetero(mol('O')),true);assert.equal(needsHetero(mol('[C,N]')),true);
});
