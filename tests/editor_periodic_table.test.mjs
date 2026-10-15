import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';

// Browser ES modules remain standalone; no project-wide Node configuration.
const source = readFileSync(new URL('../frontend/editor/periodic-table.js', import.meta.url), 'utf8');
const {ELEMENTS, findElement, searchElements, nextTableElement} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const modelSource = readFileSync(new URL('../frontend/editor/model.js', import.meta.url), 'utf8');
const {createProject, validateProject} = await import(`data:text/javascript;base64,${Buffer.from(modelSource).toString('base64')}`);

test('all 118 symbols, Chinese names and atomic numbers are complete and unique', () => {
  assert.equal(ELEMENTS.length, 118);
  assert.equal(new Set(ELEMENTS.map(element => element.symbol)).size, 118);
  assert.equal(new Set(ELEMENTS.map(element => element.atomicNumber)).size, 118);
  assert.deepEqual(ELEMENTS.map(element => element.atomicNumber), Array.from({length: 118}, (_, index) => index + 1));
  for (const element of ELEMENTS) {
    assert.match(element.symbol, /^[A-Z][a-z]?$/);
    assert.ok(element.name && /\p{Script=Han}/u.test(element.name));
  }
  assert.deepEqual(ELEMENTS.slice(-6).map(element => element.symbol), ['Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']);
});

test('element picker supports the complete canonical editor element set', () => {
  for (const element of ELEMENTS) {
    const project = createProject('element test');
    project.atoms.push({id: 'stable-atom', element: element.symbol, pos2: [0, 0], pos3: [0, 0, 0], charge: 0});
    assert.equal(validateProject(project).atoms[0].element, element.symbol);
  }
});

test('18 columns and detached 15-element lanthanoid/actinoid rows never overlap', () => {
  assert.equal(new Set(ELEMENTS.map(element => `${element.row}:${element.column}`)).size, 118);
  for (const element of ELEMENTS) assert.ok(element.column >= 1 && element.column <= 18);
  assert.deepEqual(ELEMENTS.filter(element => element.row === 10).map(element => element.atomicNumber), Array.from({length: 15}, (_, index) => index + 57));
  assert.deepEqual(ELEMENTS.filter(element => element.row === 11).map(element => element.atomicNumber), Array.from({length: 15}, (_, index) => index + 89));
  assert.equal(findElement('He').column, 18);
  assert.equal(findElement('B').column, 13);
  assert.equal(findElement('Hf').column, 4);
  assert.equal(findElement('Rf').column, 4);
  assert.equal(findElement('Og').column, 18);
});

test('search handles symbol, Chinese name and exact atomic number without simulator filtering', () => {
  assert.equal(searchElements(' c ')[0].symbol, 'C');
  assert.equal(searchElements('CL')[0].symbol, 'Cl');
  assert.deepEqual(searchElements('金').filter(element => element.name === '金').map(element => element.symbol), ['Au']);
  assert.deepEqual(searchElements('氧').map(element => element.symbol), ['O']);
  assert.deepEqual(searchElements('118').map(element => element.symbol), ['Og']);
  assert.deepEqual(searchElements('1').map(element => element.symbol), ['H']);
  assert.deepEqual(searchElements('钅尔').map(element => element.symbol), ['Nh']);
  assert.equal(searchElements('').length, 118);
  assert.equal(searchElements('119').length, 0);
  assert.equal(searchElements('unknown').length, 0);
  assert.equal(findElement('invalid'), undefined);
});

test('directional navigation skips gaps and reaches detached rows', () => {
  assert.equal(nextTableElement('H', 'ArrowRight').symbol, 'He');
  assert.equal(nextTableElement('H', 'ArrowDown').symbol, 'Li');
  assert.equal(nextTableElement('He', 'ArrowDown').symbol, 'Ne');
  assert.equal(nextTableElement('Ca', 'ArrowRight').symbol, 'Sc');
  assert.equal(nextTableElement('Ba', 'ArrowRight').symbol, 'Hf');
  assert.equal(nextTableElement('Ra', 'ArrowDown').symbol, 'La');
  assert.equal(nextTableElement('La', 'ArrowDown').symbol, 'Ac');
  assert.equal(nextTableElement('Og', 'Home').symbol, 'Fr');
  assert.equal(nextTableElement('C', 'Home', true).symbol, 'H');
  assert.equal(nextTableElement('C', 'End', true).symbol, 'Lr');
  assert.equal(nextTableElement('H', 'ArrowUp').symbol, 'H');
  assert.equal(nextTableElement('H', 'Enter'), undefined);
});

test('data is immutable and importing the component has no DOM or network side effects', () => {
  assert.ok(Object.isFrozen(ELEMENTS));
  assert.ok(ELEMENTS.every(Object.isFrozen));
  assert.throws(() => {ELEMENTS[0].symbol = 'X';}, TypeError);
  assert.doesNotMatch(source, /\bfetch\s*\(|XMLHttpRequest|https?:\/\//);
});

// Optional real-browser check using a local installed browser and Playwright.
// Supply EDITOR_PLAYWRIGHT_MODULE (the package index.mjs) and
// EDITOR_BROWSER_EXECUTABLE; no browser download or application server is used.
if (process.env.EDITOR_PLAYWRIGHT_MODULE && process.env.EDITOR_BROWSER_EXECUTABLE) {
  test('browser picker selects by mouse/keyboard, replaces explicitly, cancels and rejects a stale selection', async () => {
    const {chromium} = await import(pathToFileURL(process.env.EDITOR_PLAYWRIGHT_MODULE));
    const browser = await chromium.launch({executablePath: process.env.EDITOR_BROWSER_EXECUTABLE, headless: true});
    try {
      const page = await browser.newPage({viewport: {width: 1200, height: 900}});
      await page.setContent('<button id="opener">元素</button><dialog></dialog><style>.periodic-table{display:grid;grid-template-columns:repeat(18,42px);gap:3px}.periodic-element{height:40px}dialog{max-height:90vh} [hidden]{display:none!important}</style>');
      await page.evaluate(async code => {
        const module = await import(`data:text/javascript;base64,${btoa(unescape(encodeURIComponent(code)))}`);
        window.choices = []; window.selection = {count: 0};
        window.picker = module.mountPeriodicTable({dialog: document.querySelector('dialog'),
          getSelection: () => window.selection, onChoose: (symbol, action) => window.choices.push({symbol, ...action})});
        document.querySelector('#opener').focus(); window.picker.open();
      }, source);
      assert.equal(await page.locator('.periodic-element').count(), 118);
      assert.equal(await page.locator('.periodic-element:disabled').count(), 0);
      assert.equal(await page.locator('.periodic-replace').isVisible(), false);
      assert.equal(await page.evaluate(() => document.activeElement.className), 'periodic-search');
      await page.locator('.periodic-search').fill('氧');
      await page.locator('.periodic-search').press('Enter');
      assert.deepEqual(await page.evaluate(() => window.choices), [{symbol: 'O', replaceSelection: false}]);
      assert.equal(await page.evaluate(() => document.querySelector('dialog').open), false);
      await page.waitForFunction(() => document.activeElement.id === 'opener');

      await page.evaluate(() => {window.choices = []; window.selection = {count: 1, element: 'C'}; window.picker.open({replaceSelection: true});});
      await page.locator('.periodic-element[data-symbol="Og"]').click();
      assert.deepEqual(await page.evaluate(() => window.choices), [{symbol: 'Og', replaceSelection: true}]);

      await page.evaluate(() => {window.choices = []; window.picker.setCurrent('C'); window.picker.open();});
      await page.locator('.periodic-search').press('ArrowDown');
      assert.equal(await page.evaluate(() => document.activeElement.dataset.symbol), 'C');
      await page.keyboard.press('ArrowRight');
      assert.equal(await page.evaluate(() => document.activeElement.dataset.symbol), 'N');
      await page.keyboard.press('Enter');
      assert.deepEqual(await page.evaluate(() => window.choices), [{symbol: 'N', replaceSelection: false}]);

      await page.evaluate(() => {window.choices = []; window.picker.setCurrent('Au'); window.picker.open();});
      await page.locator('.periodic-replace').click();
      assert.deepEqual(await page.evaluate(() => window.choices), [{symbol: 'Au', replaceSelection: true}]);
      await page.evaluate(() => {window.choices = []; window.picker.open();});
      await page.keyboard.press('Escape');
      assert.deepEqual(await page.evaluate(() => window.choices), []);
      assert.equal(await page.evaluate(() => document.querySelector('dialog').open), false);

      await page.evaluate(() => {window.picker.open({replaceSelection: true}); window.selection = {count: 0};});
      await page.locator('.periodic-element[data-symbol="C"]').click();
      assert.deepEqual(await page.evaluate(() => window.choices), []);
      assert.equal(await page.evaluate(() => document.querySelector('dialog').open), true);
      assert.equal(await page.locator('.periodic-results-status').textContent(), '请先选择 1 个原子');
      await page.evaluate(() => window.picker.destroy());
      assert.equal(await page.locator('dialog > *').count(), 0);
    } finally {await browser.close();}
  });
}
