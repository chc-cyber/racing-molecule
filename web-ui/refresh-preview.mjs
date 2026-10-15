import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const uiRoot = path.dirname(fileURLToPath(import.meta.url));
const source = path.join(uiRoot, 'dist');
const preview = path.resolve(uiRoot, '..', 'wasm', 'dist', 'ui');
await fs.mkdir(preview, { recursive: true });
for (const file of ['index.html', 'app.js', 'style.css', 'minimal.css', 'integration.css', 'scene.js', 'brand-preview.html']) {
  await fs.copyFile(path.join(source, file), path.join(preview, file));
}
await fs.cp(path.join(source, 'assets'), path.join(preview, 'assets'), { recursive: true });
console.log(JSON.stringify({ preview, entry: '/ui/', files: 'modular interface and visual assets' }));
