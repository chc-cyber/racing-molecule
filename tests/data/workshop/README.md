# Curved-aromatic regression fixtures

`hexahelicene.mol` is [6]helicene (C26H16; PubChem CID 98863), generated as a 2D RDKit drawing from `C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=C4C6=CC=CC=C6C=C5`.

`spiral.mol` is the 11-ring example in [RDKit discussion 5841](https://github.com/rdkit/rdkit/discussions/5841), converted to a 2D MOL with RDKit. The discussion documents failure of the planar-sp2 prior and the `useBasicKnowledge=False` fallback. No source code from the discussion is bundled. Neither fixture specifies a P/M helicity; generated handedness is therefore not an assigned stereochemical constraint.
