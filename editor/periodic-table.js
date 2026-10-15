// Preserve useful errors when the JSmol Java runtime exports its Error class.
const NativeError=globalThis.Clazz?._Error??globalThis.Error;
/**
 * Local element picker. Element facts are compiled from IUPAC and the Chinese
 * Chemical Society; see research/reaxys-interaction.md for source attribution.
 * This module contains no Reaxys/Marvin code, images, runtime, or network calls.
 * It chooses an element; the host owns atom IDs, topology changes and history.
 */
const SYMBOLS = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split(' ');
const NAMES = '氢 氦 锂 铍 硼 碳 氮 氧 氟 氖 钠 镁 铝 硅 磷 硫 氯 氩 钾 钙 钪 钛 钒 铬 锰 铁 钴 镍 铜 锌 镓 锗 砷 硒 溴 氪 铷 锶 钇 锆 铌 钼 锝 钌 铑 钯 银 镉 铟 锡 锑 碲 碘 氙 铯 钡 镧 铈 镨 钕 钷 钐 铕 钆 铽 镝 钬 铒 铥 镱 镥 铪 钽 钨 铼 锇 铱 铂 金 汞 铊 铅 铋 钋 砹 氡 钫 镭 锕 钍 镤 铀 镎 钚 镅 锔 锫 锎 锿 镄 钔 锘 铹 𬬻 𬭊 𬭳 𬭛 𬭶 鿏 𫟼 𬬭 鿔 鿭 𫓧 镆 𫟷 鿬 鿫'.split(' ');
// Decomposed aliases help users search and screen readers name uncommon glyphs.
const ALIASES = {
  Rf: '钅卢 金卢', Db: '钅杜 金杜', Sg: '钅喜 金喜', Bh: '钅波 金波',
  Hs: '钅黑 金黑', Mt: '钅麦 金麦', Ds: '钅达 金达', Rg: '钅仑 金仑',
  Cn: '钅哥 金哥', Nh: '钅尔 金尔', Fl: '钅夫 金夫', Lv: '钅立 金立',
  Ts: '石田', Og: '气奥',
};
const GROUPS = {
  alkali: new Set(['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']),
  alkaline: new Set(['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']),
  nonmetal: new Set(['H', 'C', 'N', 'O', 'P', 'S', 'Se']),
  metalloid: new Set(['B', 'Si', 'Ge', 'As', 'Sb', 'Te']),
  halogen: new Set(['F', 'Cl', 'Br', 'I', 'At', 'Ts']),
  noble: new Set(['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og']),
  post: new Set(['Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po', 'Nh', 'Fl', 'Mc', 'Lv']),
};

function position(atomicNumber) {
  if (atomicNumber <= 2) return {period: 1, column: atomicNumber === 1 ? 1 : 18, row: 2};
  if (atomicNumber <= 10) return {period: 2, column: atomicNumber < 5 ? atomicNumber - 2 : atomicNumber + 8, row: 3};
  if (atomicNumber <= 18) return {period: 3, column: atomicNumber < 13 ? atomicNumber - 10 : atomicNumber, row: 4};
  if (atomicNumber <= 36) return {period: 4, column: atomicNumber - 18, row: 5};
  if (atomicNumber <= 54) return {period: 5, column: atomicNumber - 36, row: 6};
  if (atomicNumber <= 56) return {period: 6, column: atomicNumber - 54, row: 7};
  if (atomicNumber <= 71) return {period: 6, column: atomicNumber - 54, row: 10};
  if (atomicNumber <= 86) return {period: 6, column: atomicNumber - 68, row: 7};
  if (atomicNumber <= 88) return {period: 7, column: atomicNumber - 86, row: 8};
  if (atomicNumber <= 103) return {period: 7, column: atomicNumber - 86, row: 11};
  return {period: 7, column: atomicNumber - 100, row: 8};
}

export const ELEMENTS = Object.freeze(SYMBOLS.map((symbol, index) => {
  const atomicNumber = index + 1;
  let category = atomicNumber >= 57 && atomicNumber <= 71 ? 'lanthanoid'
    : atomicNumber >= 89 && atomicNumber <= 103 ? 'actinoid' : 'transition';
  for (const [name, symbols] of Object.entries(GROUPS)) if (symbols.has(symbol)) category = name;
  return Object.freeze({symbol, atomicNumber, name: NAMES[index], aliases: ALIASES[symbol] ?? '', category, ...position(atomicNumber)});
}));
const BY_SYMBOL = new Map(ELEMENTS.map(element => [element.symbol.toLowerCase(), element]));
export function findElement(symbol) { return BY_SYMBOL.get(String(symbol ?? '').trim().toLowerCase()); }

/** Numeric search is exact; symbol matches rank above partial matches. */
export function searchElements(query) {
  const value = String(query ?? '').trim().toLowerCase();
  if (!value) return [...ELEMENTS];
  if (/^\d+$/.test(value)) return ELEMENTS.filter(element => element.atomicNumber === Number(value));
  return ELEMENTS.filter(element => `${element.symbol} ${element.name} ${element.aliases}`.toLowerCase().includes(value))
    .sort((a, b) => Number(b.symbol.toLowerCase() === value) - Number(a.symbol.toLowerCase() === value) || a.atomicNumber - b.atomicNumber);
}

/** Spatial keyboard navigation skips the empty cells, keeping the closest column. */
export function nextTableElement(symbol, key, control = false) {
  const current = findElement(symbol);
  if (!current) return undefined;
  const visual = [...ELEMENTS].sort((a, b) => a.row - b.row || a.column - b.column);
  if (control && key === 'Home') return visual[0];
  if (control && key === 'End') return visual.at(-1);
  const row = visual.filter(element => element.row === current.row);
  const index = row.indexOf(current);
  if (key === 'Home') return row[0];
  if (key === 'End') return row.at(-1);
  if (key === 'ArrowLeft') return row[Math.max(0, index - 1)];
  if (key === 'ArrowRight') return row[Math.min(row.length - 1, index + 1)];
  if (key !== 'ArrowUp' && key !== 'ArrowDown') return undefined;
  const direction = key === 'ArrowDown' ? 1 : -1;
  const rows = [...new Set(visual.map(element => element.row))];
  const nextRow = rows[rows.indexOf(current.row) + direction];
  if (nextRow === undefined) return current;
  return visual.filter(element => element.row === nextRow).sort((a, b) => Math.abs(a.column - current.column) - Math.abs(b.column - current.column) || a.column - b.column)[0];
}

let dialogSequence = 0;

/**
 * @param {{dialog: HTMLDialogElement, onChoose: Function, getSelection?: Function}} options
 * onChoose(symbol, {replaceSelection}) is called once per committed choice.
 * open({replaceSelection:true}) lets an atom-property button replace directly.
 * Regular opening chooses the drawing element. The explicit footer action also
 * replaces the single selection with the current drawing element.
 */
export function mountPeriodicTable({dialog, onChoose, getSelection = () => ({count: 0})}) {
  if (!dialog?.ownerDocument || dialog.tagName !== 'DIALOG') throw new TypeError('周期表需要一个 HTMLDialogElement。');
  if (typeof onChoose !== 'function') throw new TypeError('周期表需要 onChoose 回调。');
  if (dialog.dataset.periodicMounted) throw new NativeError('此周期表弹窗已经初始化。');
  const doc = dialog.ownerDocument;
  const labelId = `periodic-title-${++dialogSequence}`;
  const hintId = `periodic-hint-${dialogSequence}`;
  const events = new AbortController();
  let current = findElement('C'), replaceMode = false, opener = null, destroyed = false;
  let resultElements = [], resultButtons = [], focusSymbol = current.symbol;
  const buttons = new Map();
  const listen = (node, type, fn) => node.addEventListener(type, fn, {signal: events.signal});
  const make = (tag, className, text) => {
    const node = doc.createElement(tag); node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  };
  const header = make('div', 'periodic-header');
  const title = make('h2', 'periodic-title', '元素周期表'); title.id = labelId;
  const dismiss = make('button', 'periodic-close', '×'); dismiss.type = 'button';
  dismiss.setAttribute('aria-label', '关闭元素周期表'); dismiss.title = '关闭（Esc）';
  header.append(title, dismiss);
  const currentLine = make('div', 'periodic-current');
  const hint = make('p', 'periodic-hint'); hint.id = hintId;
  const searchRow = make('div', 'periodic-search-row');
  const input = make('input', 'periodic-search'); input.type = 'search';
  input.placeholder = '元素符号 / 中文名 / 原子序数'; input.autocomplete = 'off';
  input.setAttribute('aria-label', '搜索元素符号、中文名或原子序数');
  const resultsStatus = make('span', 'periodic-results-status'); resultsStatus.setAttribute('role', 'status');
  searchRow.append(input, resultsStatus);
  const preview = make('div', 'periodic-preview'); preview.setAttribute('aria-live', 'polite');
  const scroll = make('div', 'periodic-table-scroll');
  const table = make('div', 'periodic-table'); table.setAttribute('role', 'group'); table.setAttribute('aria-label', '118 种元素，方向键移动，Enter 选择');
  for (let group = 1; group <= 18; group++) {
    const label = make('span', 'periodic-group', String(group)); label.setAttribute('aria-hidden', 'true');
    label.style.gridColumn = String(group); label.style.gridRow = '1'; table.append(label);
  }
  for (const [text, row, column, span] of [['57–71', 7, 3, 1], ['89–103', 8, 3, 1], ['镧系', 10, 1, 2], ['锕系', 11, 1, 2]]) {
    const label = make('span', 'periodic-series', text); label.setAttribute('aria-hidden', 'true');
    label.style.gridColumn = `${column} / span ${span}`; label.style.gridRow = String(row); table.append(label);
  }
  const describe = element => `${element.symbol} · ${element.name}${element.aliases ? `（${element.aliases.split(' ')[0]}）` : ''} · 原子序数 ${element.atomicNumber}`;
  function makeElementButton(element, result = false) {
    const button = make('button', result ? 'periodic-result' : 'periodic-element'); button.type = 'button';
    button.dataset.symbol = element.symbol; button.dataset.category = element.category;
    button.setAttribute('aria-label', describe(element)); button.title = describe(element);
    const number = make('span', 'periodic-atomic-number', String(element.atomicNumber));
    const symbol = make('span', 'periodic-symbol', element.symbol);
    button.append(number, symbol);
    if (result) button.append(make('span', 'periodic-name', element.name));
    else {button.style.gridColumn = String(element.column); button.style.gridRow = String(element.row);}
    listen(button, 'click', () => choose(element));
    listen(button, 'focus', () => {preview.textContent = describe(element); if (!result) setGridFocus(element.symbol);});
    listen(button, 'mouseenter', () => {preview.textContent = describe(element);});
    listen(button, 'keydown', event => {
      if (result) {
        const index = resultButtons.indexOf(button);
        const delta = ['ArrowRight', 'ArrowDown'].includes(event.key) ? 1 : ['ArrowLeft', 'ArrowUp'].includes(event.key) ? -1 : 0;
        const target = event.key === 'Home' ? 0 : event.key === 'End' ? resultButtons.length - 1 : index + delta;
        if (delta || event.key === 'Home' || event.key === 'End') {
          event.preventDefault(); resultButtons[Math.max(0, Math.min(resultButtons.length - 1, target))]?.focus();
        }
      } else {
        const next = nextTableElement(element.symbol, event.key, event.ctrlKey || event.metaKey);
        if (next) {event.preventDefault(); buttons.get(next.symbol)?.focus();}
      }
    });
    return button;
  }
  for (const element of ELEMENTS) {const button = makeElementButton(element); buttons.set(element.symbol, button); table.append(button);}
  scroll.append(table);
  const results = make('div', 'periodic-results'); results.hidden = true; results.setAttribute('role', 'group'); results.setAttribute('aria-label', '匹配的元素');
  const footer = make('div', 'periodic-footer');
  const selection = make('span', 'periodic-selection');
  const replace = make('button', 'periodic-replace'); replace.type = 'button';
  footer.append(selection, replace);
  dialog.replaceChildren(header, currentLine, hint, searchRow, preview, scroll, results, footer);
  dialog.classList.add('periodic-dialog'); dialog.dataset.periodicMounted = 'true';
  dialog.setAttribute('aria-labelledby', labelId); dialog.setAttribute('aria-describedby', hintId);

  function setGridFocus(symbol) {
    focusSymbol = symbol;
    for (const [value, button] of buttons) button.tabIndex = value === symbol ? 0 : -1;
  }
  function selectionState() {
    const state = getSelection() ?? {};
    return {count: Number(state.count) || 0, element: findElement(state.element)};
  }
  function updateHeader() {
    currentLine.textContent = `当前元素 ${current.symbol} · ${current.name} · ${current.atomicNumber}`;
    hint.textContent = replaceMode ? '选择元素，替换当前选中原子' : '选元素后点击空白添加 / 点击原子替换';
    const state = selectionState();
    selection.textContent = state.count === 1 ? `已选 1 个原子${state.element ? ` · ${state.element.symbol}` : ''}` : state.count ? `已选 ${state.count} 个原子` : '未选中原子';
    replace.textContent = `替换当前选中原子为 ${current.symbol}`;
    replace.disabled = state.count !== 1; replace.hidden = state.count !== 1;
    for (const [symbol, button] of buttons) button.setAttribute('aria-pressed', String(symbol === current.symbol));
  }
  function choose(element, forceReplace = replaceMode) {
    if (forceReplace && selectionState().count !== 1) {
      resultsStatus.textContent = '请先选择 1 个原子'; return;
    }
    // Invoke before closing: host rejection leaves the picker available.
    onChoose(element.symbol, {replaceSelection: forceReplace});
    current = element; updateHeader(); close();
  }
  function updateSearch() {
    const active = Boolean(input.value.trim());
    resultElements = active ? searchElements(input.value) : [];
    resultButtons = resultElements.map(element => makeElementButton(element, true));
    results.replaceChildren(...resultButtons);
    if (active && !resultElements.length) results.append(make('p', 'periodic-no-results', '没有匹配的元素'));
    scroll.hidden = active; results.hidden = !active;
    resultsStatus.textContent = active ? `${resultElements.length} 个结果` : '';
  }
  function close() {
    if (dialog.open) {
      if (typeof dialog.close === 'function') dialog.close(); else dialog.removeAttribute('open');
    }
  }
  function setCurrent(symbol) {
    const element = findElement(symbol);
    if (!element) throw new RangeError(`未知元素：${symbol}`);
    current = element; setGridFocus(element.symbol); updateHeader(); preview.textContent = describe(element);
  }
  function open({replaceSelection = false} = {}) {
    if (destroyed) throw new NativeError('周期表组件已经销毁。');
    replaceMode = Boolean(replaceSelection);
    opener = doc.activeElement; input.value = ''; updateSearch(); updateHeader();
    preview.textContent = describe(current); setGridFocus(current.symbol);
    if (!dialog.open) {
      if (typeof dialog.showModal === 'function') dialog.showModal(); else dialog.setAttribute('open', '');
    }
    input.focus();
  }
  listen(dismiss, 'click', close);
  listen(replace, 'click', () => choose(current, true));
  listen(input, 'input', updateSearch);
  listen(input, 'keydown', event => {
    if (event.isComposing) return;
    if (event.key === 'ArrowDown') {event.preventDefault(); (input.value.trim() ? resultButtons[0] : buttons.get(focusSymbol))?.focus();}
    if (event.key === 'Enter' && input.value.trim() && resultElements.length) {event.preventDefault(); choose(resultElements[0]);}
  });
  listen(dialog, 'keydown', event => {
    if (event.key === 'Escape') {event.preventDefault(); close();}
  });
  listen(dialog, 'click', event => {
    if (event.target !== dialog) return;
    const rect = dialog.getBoundingClientRect();
    if (event.clientX < rect.left || event.clientX > rect.right || event.clientY < rect.top || event.clientY > rect.bottom) close();
  });
  listen(dialog, 'close', () => {opener?.isConnected && opener.focus();});
  setCurrent('C');
  return Object.freeze({open, close, setCurrent, destroy() {
    close(); events.abort(); dialog.replaceChildren(); delete dialog.dataset.periodicMounted; destroyed = true;
  }});
}
