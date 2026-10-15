/** Pure presentation/gesture helpers. They never infer or mutate topology. */
export const BOND_LAYOUT_LENGTH = 1.5;

export function snapBondEndpoint(start, pointer, free = false) {
  if (free) return [...pointer];
  const dx = pointer[0] - start[0], dy = pointer[1] - start[1];
  const angle = Math.round(Math.atan2(dy, dx) / (Math.PI / 6)) * Math.PI / 6;
  return [start[0] + BOND_LAYOUT_LENGTH * Math.cos(angle), start[1] + BOND_LAYOUT_LENGTH * Math.sin(angle)];
}

export function pointInPolygon(point, polygon) {
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const a = polygon[i], b = polygon[j];
    if ((a[1] > point[1]) !== (b[1] > point[1]) && point[0] < (b[0] - a[0]) * (point[1] - a[1]) / (b[1] - a[1]) + a[0]) inside = !inside;
  }
  return inside;
}

export function segmentDistance(point, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1], len2 = dx * dx + dy * dy;
  const t = len2 ? Math.max(0, Math.min(1, ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / len2)) : 0;
  return Math.hypot(point[0] - a[0] - t * dx, point[1] - a[1] - t * dy);
}

export function atomLabel(atom, bonded = true) {
  const charge = atom.charge ? `${Math.abs(atom.charge) === 1 ? '' : Math.abs(atom.charge)}${atom.charge > 0 ? '+' : '−'}` : '';
  const hydrogens = atom.element !== 'H' && atom.explicitHs ? `H${atom.explicitHs > 1 ? atom.explicitHs : ''}` : '';
  const ordinaryCarbon = atom.element === 'C' && !atom.isotope && !atom.charge && !atom.explicitHs && !atom.radicalElectrons;
  if (ordinaryCarbon && bonded) return '';
  return `${atom.isotope || ''}${atom.element}${hydrogens}${charge}${atom.radicalElectrons ? '•'.repeat(Math.min(atom.radicalElectrons, 2)) : ''}`;
}
