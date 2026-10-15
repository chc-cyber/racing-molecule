/* Compatibility entry for the current workshop's actual browser acceptance.
   Use an independent test service; this clicks Send to track on that instance.
   Usage: node scripts/verify_editor_browser.cjs http://127.0.0.1:8778 */
require('./verify_editor_mature_browser.cjs');
