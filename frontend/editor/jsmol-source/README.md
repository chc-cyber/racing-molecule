# JSmol native ModelKit extension

This directory preserves the original JSmol 16.3.31 files and an auditable
unified source patch for the Molecular Car ModelKit menu bridge, dated
2026-10-13. These library changes remain under LGPL-2.1-or-later, matching
upstream. Jmol Development Team copyright and the original license notices are
preserved; the full library license is in `../vendor/jsmol/LICENSE-LGPL-2.1.txt`.
This does not select a license for the overall Molecular Car project.

The main workshop now sets the host option `nativeMenus:false`: native popup
menus and their corner/logo entry points are disabled, while ModelKit stays
open for actual coordinate/chemistry editing through the sidebars. The isolated
research page retains the menu bridge for regression checks. This switch lives
in `../jsmol-host.js`; it does not further modify the LGPL library files.

The installed readable runtime source is in `../vendor/jsmol/j2s/`. The patch
modifies five files; originals under `upstream/` are verified against the zip
manifest hashes. No core molecular simulation or renderer has been changed.

| Source | Change |
|---|---|
| `J/modelkit/ModelKitPopupResourceBundle.js` | Genuine native menu entries for project operations; retains native atom, bond, charge and geometry picking. Replaces empty custom-element `??` slots with one periodic-table owner action. Removes independent state restoration, automatic minimization and crystal commands. |
| `J/modelkit/ModelKitPopup.js` | Native menu action bridge, current shared-history enablement before popup rendering, and molecular-only menu handling. |
| `JS/JMenu.js` | Native menu actions dispatch without the original 500 ms wait when attached to a Molecular Car owner bridge. |
| `core/coreswing.js` | The same JMenu timing source change in upstream's combined runtime file. |
| `core/coreswing.z.js` | The same conditional timer expression in the actually loaded packed runtime; no re-minification. |

`JSmol.min.js`, `js/JSmol.full.js`, `core/corejmol.z.js` and the molecular
ModelKit implementation are unchanged. This is a native source extension, not a
DOM toolbar injection or an imitation of ModelKit controls.

The original distribution and corresponding upstream source are available at:

- https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/jsmol.zip/download
- https://sourceforge.net/projects/jmol/files/Jmol/Version%2016.3/Jmol%2016.3.31/Jmol-16.3.31-full.tar.gz/download

The original zip SHA-256 remains
`4147ba69d31a83ed6c276a1feaf3bd2b6472e82075ec5c96dd03317d0e568655`.
`../vendor/jsmol/manifest.json` records original and installed hashes separately.
After an intentional library source change, run
`python frontend/editor/jsmol-source/update-manifest.py` from the project root
to regenerate the patch and installed hashes. Consumers can inspect, replace
or restore these URL-loaded library sources independently of the application.

The current bridge is `molecular-car-modelkit-3`. The native
`periodicTable` menu entry dispatches `mcar:periodic-table` to the host's
existing 118-element picker; it does not evaluate a separate native element
prompt or keep another selection/history store. Original `??` definitions
remain in the preserved upstream code, but are absent from the active menu.

Real browser acceptance and the callback contract are documented in
`../research/jsmol-native.md`; no public service needs to be restarted for the
static browser tests.

The independent `../jsmol-host.js` also guards the applet instance's mouse
dispatcher against a release with no active press becoming a click. Upstream
JSmol's canvas mouseenter sends a synthetic 502 after 504; if its coordinates
equal the last press, upstream Mouse.processEvent would otherwise click again.
The host preserves native released recovery and every real press/drag/release
path. The vendor Mouse source and selection bindings are unchanged. Actual
hover, selection and shared periodic-table replacement acceptance is recorded
in `../VALIDATION.md`.

Version 3 gives independent native menu entries for coordinate rebuilding (`mcar:rebuild`), adding real H (`mcar:hydrogenate`) and UFF optimization (`mcar:relax`). The host combines ordinary bond-order picking and atom-to-atom/blank dragging in one UI tool, using declared targets and shared canonical transactions. These gestures retain native ModelKit rendering and editor IDs.
