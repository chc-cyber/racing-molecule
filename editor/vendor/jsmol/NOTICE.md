# JSmol local runtime — 16.3.31

This directory contains upstream files, selected from the official
Jmol 16.3.31 JSmol distribution released on 2025-07-22. The tested runtime reports
`_version = 1603031` and `Jmol version 16.3.31 2025-07-21 09:29`.

- Official release: https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/
- Archive: https://downloads.sourceforge.net/project/jmol/Jmol/Version%2016.3/Jmol%2016.3.31/jsmol.zip
- Archive SHA-256: `4147ba69d31a83ed6c276a1feaf3bd2b6472e82075ec5c96dd03317d0e568655`
- Archive size: 48,459,556 bytes.
- Included upstream files: `JSmol.min.js`, the complete `j2s/` tree,
  `js/JSmol.full.js` (uncompressed JavaScript bootstrap source), and `README.TXT`.
  The original upstream selection totals 1,773 files / 49,756,122 bytes, excluding these added license/notice
  and manifest files. Java binaries, PHP relay, demonstrations and unrelated
  packages are omitted. No runtime CDN is used.
- `manifest.json` records the source URL, selection and per-file SHA-256 values.
- Five library source files carry the documented `molecular-car-modelkit-1`
  extension dated 2026-10-13. Original hashes and installed hashes are recorded
  separately. Readable source, pristine originals and a unified patch are retained
  in [jsmol-source](../../jsmol-source/README.md). All other upstream files remain
  unmodified. The native menu changes keep project geometry and history with the
  application owner; atom/bond tools and mouse transforms are still JSmol's.

## Copyright and licenses

Jmol / JSmol is developed by The Jmol Development Team and contributors.
The Jmol / ModelKit sources license the library under **GNU LGPL version 2.1 or
later**. The complete LGPL 2.1 text is retained in `LICENSE-LGPL-2.1.txt`.
Existing headers and copyright notices in every copied upstream file are retained.

The bootstrap includes **jQuery 3.7.1**, copyright OpenJS Foundation and other
contributors, under the **MIT license**. Its original license is retained in
`LICENSE-jQuery-MIT.txt`, obtained from
https://raw.githubusercontent.com/jquery/jquery/3.7.1/LICENSE.txt.
The uncompressed bootstrap source also preserves upstream's jQuery code and
JSmol adaptations. The bootstrap has not been modified or minified again. The
packed Swing runtime has one conditional menu-timer expression change matching
its retained readable source; it has not been re-minified.

## Source and replacement

The corresponding full upstream source release is available at:
https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/Jmol-16.3.31-full.tar.gz/download

The library remains separate from the project adapter. Its JavaScript files are
loaded by URL and can be inspected and replaced. The independent live validation
page is `../../research/jsmol-live.html`. Its project-specific adapter is
`../../jsmol-adapter.js`; the adapter is not an upstream JSmol modification.
The native ModelKit library extension is documented separately in
`../../jsmol-source/README.md`; its changes remain LGPL-2.1-or-later, matching the
upstream library. Its complete patch and original sources can be inspected and
replaced independently of the application's adapter.

This notice documents this dependency and does not choose a license for the
Molecular Car project. Keep these notices and the library's license/source
availability and replacement terms when redistributing it.
