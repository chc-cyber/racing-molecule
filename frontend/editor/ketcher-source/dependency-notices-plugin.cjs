// This expression is embedded into the generated Vite configuration.
// Actual included module IDs, not the complete development install, define the
// dependency notice inventory. Original upstream notices are retained too.
({
  name: 'molecular-car-dependency-notices',
  generateBundle(_options, bundle) {
    const packageRoots = new Set();
    for (const chunk of Object.values(bundle)) {
      if (chunk.type !== 'chunk') continue;
      for (const rawId of Object.keys(chunk.modules || {})) {
        const id = rawId.replaceAll('\\', '/').replace(/^\u0000/, '');
        const at = id.lastIndexOf('/node_modules/');
        if (at < 0) continue;
        const after = id.slice(at + '/node_modules/'.length).split('/');
        const name = after[0].startsWith('@') ? after.slice(0, 2).join('/') : after[0];
        const root = id.slice(0, at + '/node_modules/'.length) + name;
        if (fs.existsSync(path.join(root, 'package.json'))) packageRoots.add(root);
      }
    }
    const dependencies = [];
    for (const root of [...packageRoots].sort()) {
      const pkg = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
      const files = fs.readdirSync(root).filter(name => /^(licen[sc]e|notice|copying)([.-].*)?$/i.test(name) && fs.statSync(path.join(root, name)).isFile());
      const noticePaths = [];
      for (const name of files) {
        const filename = `licenses/${pkg.name.replaceAll('/', '_')}-${name}`;
        this.emitFile({type: 'asset', fileName: filename, source: fs.readFileSync(path.join(root, name))});
        noticePaths.push(filename);
      }
      dependencies.push({name: pkg.name, version: pkg.version, license: pkg.license || pkg.licenses || 'See upstream NOTICE', repository: pkg.repository || pkg.homepage || null, noticePaths});
    }
    this.emitFile({type: 'asset', fileName: 'BUILD-DEPENDENCIES.json', source: JSON.stringify({method: 'Package roots of modules included in emitted Vite chunks; worker notices also covered by preserved official upstream LICENSE.txt files', dependencies}, null, 2)+'\n'});
  },
})
