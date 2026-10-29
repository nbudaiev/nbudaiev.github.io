/* viz.js: canvas rendering (images, plots) and the live "CASA call" control panel. */
(function () {
'use strict';
const V = {};
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
V.css = css;
const color = c => (c && c.startsWith('--') ? css(c) : c);

/* ---------- colour maps ---------- */
function rgb(h) { const n = parseInt(h.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function lut(stops) {
  const out = new Uint8ClampedArray(256 * 3), c = stops.map(rgb);
  for (let i = 0; i < 256; i++) {
    const t = i / 255 * (c.length - 1), k = Math.min(c.length - 2, Math.floor(t)), f = t - k;
    for (let j = 0; j < 3; j++) out[i * 3 + j] = Math.round(c[k][j] * (1 - f) + c[k + 1][j] * f);
  }
  return out;
}
const CMAPS = {
  inferno: lut(['#000004', '#160b39', '#420a68', '#6a176e', '#932667', '#bc3754', '#dd513a', '#f37819', '#fca50a', '#f6d746', '#fcffa4']),
  viridis: lut(['#440154', '#482475', '#414487', '#355f8d', '#2a788e', '#21918c', '#22a884', '#44bf70', '#7ad151', '#bddf26', '#fde725']),
  rdbu: lut(['#173f7a', '#2f6db3', '#86afe0', '#eeeeea', '#e8907a', '#c0422e', '#7f1d14']),
  gray: lut(['#0b0d12', '#f4f5f7']),
  blues: lut(['#eef4fc', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
};
V.lut = name => CMAPS[name];
V.cmapCSS = function (name) {               // CSS gradient for legends
  const l = CMAPS[name], s = [];
  for (let i = 0; i <= 8; i++) { const k = Math.round(i / 8 * 255) * 3; s.push(`rgb(${l[k]},${l[k + 1]},${l[k + 2]}) ${i * 12.5}%`); }
  return `linear-gradient(90deg, ${s.join(',')})`;
};

/* ---------- canvas plumbing: DPR-aware, re-rendered on resize and theme change ---------- */
const live = new Set();
const ro = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(es => es.forEach(e => e.target._render && e.target._render())) : null;
function register(cv, fn) {
  cv._render = fn;
  if (!live.has(cv)) { live.add(cv); if (ro) ro.observe(cv); }
  fn();
}
V.rerenderAll = () => live.forEach(cv => cv.isConnected && cv._render && cv._render());
function setup(cv) {
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  if (!w || !h) return null;
  const W = Math.round(w * dpr), H = Math.round(h * dpr);
  if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h, fs: parseFloat(getComputedStyle(cv).fontSize) || 14 };
}
const offs = {};
function offscreen(M) {
  if (!offs[M]) { const cv = document.createElement('canvas'); cv.width = cv.height = M; offs[M] = { cv, ctx: cv.getContext('2d') }; }
  return offs[M];
}

/* ---------- ticks ---------- */
function niceStep(range, n) {
  const raw = range / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
  return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
}
function ticks(min, max, n) {
  const s = niceStep(max - min || 1, n || 5), out = [];
  for (let v = Math.ceil(min / s - 1e-9) * s; v <= max + s * 1e-9; v += s) out.push(Math.abs(v) < s * 1e-9 ? 0 : v);
  return { ticks: out, step: s };
}
function fmt(v, step) {
  const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return v.toFixed(Math.min(d, 6)).replace('-', '−');
}
V.ticks = ticks;
function logTicks(min, max) {
  const out = [], a = Math.floor(Math.log10(min)), b = Math.ceil(Math.log10(max)), few = b - a <= 2;
  for (let e = a; e <= b; e++) for (const m of few ? [1, 2, 5] : [1]) { const v = m * Math.pow(10, e); if (v >= min * 0.999 && v <= max * 1.001) out.push(v); }
  return out;
}
const fmtLog = v => (v >= 1000 || v < 0.01 ? v.toExponential(0).replace('e+', 'e') : String(+v.toPrecision(2)));

/* ---------- images ----------
   data: Float64Array N*N (x = East, y = North). o: {cell, crop, vmin, vmax, cmap, beam, mask, overlay(ctx,map),
   cbar: {label, scale}, axes, title, smooth, beamColor, symmetric} */
V.image = function (cv, data, N, o) { register(cv, () => drawImage(cv, data, N, o || {})); };
function drawImage(cv, data, N, o) {
  const s = setup(cv); if (!s) return;
  const { ctx, w, h, fs } = s, lab = fs * 0.62, M = Math.min(o.crop || N, N), cell = o.cell || 1, axes = o.axes !== false;
  const ml = axes ? lab * 3.3 : 2, mb = axes ? lab * 2.7 : 2, mt = 2, mr = o.cbar ? lab * 4.8 : 2;
  const size = Math.max(10, Math.floor(Math.min(w - ml - mr, h - mt - mb)));
  const ox = Math.round(ml + (w - ml - mr - size) / 2), oy = mt;
  const h0 = N / 2, m2 = M / 2, fx = o.flipX === false ? 1 : -1;
  const colX = c => (h0 + fx * (c - m2) + N) % N;
  let vmin = o.vmin, vmax = o.vmax;
  if (vmin == null || vmax == null) {
    let a = Infinity, b = -Infinity;
    for (let r = 0; r < M; r++) for (let c = 0; c < M; c++) { const v = data[((h0 + m2 - r + N) % N) * N + colX(c)]; if (v < a) a = v; if (v > b) b = v; }
    if (o.symmetric) { const m = Math.max(Math.abs(a), Math.abs(b)) || 1; a = -m; b = m; }
    if (vmin == null) vmin = a; if (vmax == null) vmax = b;
  }
  if (vmax <= vmin) vmax = vmin + 1e-12;
  const off = offscreen(M), id = off.ctx.createImageData(M, M), L = CMAPS[o.cmap || 'inferno'], scale = 1 / (vmax - vmin);
  const sq = o.stretch === 'sqrt', fwd = t => (sq ? Math.sqrt(t) : t);           // display stretch; colour bar follows it
  for (let r = 0; r < M; r++) {
    const y = (h0 + m2 - r + N) % N;
    for (let c = 0; c < M; c++) {
      const val = data[y * N + colX(c)], p = (r * M + c) * 4;
      if (val !== val) { id.data[p + 3] = 0; continue; }         // NaN: transparent (masked)
      let t = (val - vmin) * scale; t = t < 0 ? 0 : t > 1 ? 1 : t;
      const li = ((sq ? Math.sqrt(t) : t) * 255 | 0) * 3;
      id.data[p] = L[li]; id.data[p + 1] = L[li + 1]; id.data[p + 2] = L[li + 2]; id.data[p + 3] = 255;
    }
  }
  off.ctx.putImageData(id, 0, 0);
  ctx.imageSmoothingEnabled = !!o.smooth;
  ctx.drawImage(off.cv, 0, 0, M, M, ox, oy, size, size);
  const px = size / M;
  const map = { X: l => ox + (m2 + fx * l / cell + 0.5) * px, Y: m => oy + (m2 - m / cell + 0.5) * px, px, cell, ox, oy, size, M, N, k: px / cell, vmin, vmax };
  if (o.mask) maskOutline(ctx, o.mask, N, M, map, color(o.maskColor || '--mask'), colX);
  if (o.beam && o.beam.maj > 0) {
    const b = o.beam, rx = b.maj / 2 * map.k, ry = b.min / 2 * map.k, pad = Math.max(rx, 4) + lab * 0.7;
    ctx.beginPath(); ctx.ellipse(ox + pad, oy + size - pad, rx, ry, Math.PI / 2 - b.pa * Math.PI / 180, 0, 2 * Math.PI);
    ctx.fillStyle = o.beamColor || 'rgba(255,255,255,0.88)'; ctx.fill();
    ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(0,0,0,0.55)'; ctx.stroke();
  }
  if (o.frame !== false) { ctx.strokeStyle = css('--rule'); ctx.lineWidth = 1; ctx.strokeRect(ox + 0.5, oy + 0.5, size - 1, size - 1); }
  ctx.font = `${lab}px ${css('--f-body')}`; ctx.fillStyle = css('--muted');
  if (o.extent) {                                 // [x_left, x_right, y_bottom, y_top] in data units
    const e = o.extent;
    map.EX = v => ox + (v - e[0]) / (e[1] - e[0]) * size; map.EY = v => oy + (1 - (v - e[2]) / (e[3] - e[2])) * size;
  }
  if (axes && o.extent) {
    const e = o.extent, tx = ticks(Math.min(e[0], e[1]), Math.max(e[0], e[1]), size / (lab * 6)), ty = ticks(Math.min(e[2], e[3]), Math.max(e[2], e[3]), size / (lab * 3));
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.strokeStyle = css('--muted');
    for (const t of tx.ticks) { const xx = map.EX(t); if (xx < ox + 2 || xx > ox + size - 2) continue; ctx.fillText(fmt(t, tx.step), xx, oy + size + lab * 0.35); ctx.beginPath(); ctx.moveTo(xx, oy + size); ctx.lineTo(xx, oy + size - lab * 0.4); ctx.stroke(); }
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const t of ty.ticks) { const yy = map.EY(t); if (yy < oy + 2 || yy > oy + size - 2) continue; ctx.fillText(fmt(t, ty.step), ox - lab * 0.35, yy); ctx.beginPath(); ctx.moveTo(ox, yy); ctx.lineTo(ox + lab * 0.4, yy); ctx.stroke(); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
    ctx.fillText(o.xlabel || '', ox + size / 2, oy + size + lab * 2.45);
    ctx.save(); ctx.translate(ox - lab * 2.55, oy + size / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
  } else if (axes) {
    const ext = m2 * cell, tk = ticks(-ext, ext, size / (lab * 6));
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.strokeStyle = css('--muted');
    for (const t of tk.ticks) {
      const xx = map.X(t);
      if (xx < ox + 2 || xx > ox + size - 2) continue;
      ctx.fillText(fmt(t, tk.step), xx, oy + size + lab * 0.35);
      ctx.beginPath(); ctx.moveTo(xx, oy + size); ctx.lineTo(xx, oy + size - lab * 0.4); ctx.stroke();
    }
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const t of tk.ticks) {
      const yy = map.Y(t);
      if (yy < oy + 2 || yy > oy + size - 2) continue;
      ctx.fillText(fmt(t, tk.step), ox - lab * 0.35, yy);
      ctx.beginPath(); ctx.moveTo(ox, yy); ctx.lineTo(ox + lab * 0.4, yy); ctx.stroke();
    }
    ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
    ctx.fillText(o.xlabel || 'ΔRA (″)', ox + size / 2, oy + size + lab * 2.45);
    ctx.save(); ctx.translate(ox - lab * 2.55, oy + size / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top';
    ctx.fillText(o.ylabel || 'ΔDec (″)', 0, 0); ctx.restore();
  }
  if (o.overlay) { ctx.save(); ctx.beginPath(); ctx.rect(ox, oy, size, size); ctx.clip(); o.overlay(ctx, map); ctx.restore(); }
  if (o.cbar) {
    const bx = ox + size + lab * 0.6, bw = lab * 0.7, sc = o.cbar.scale || 1;
    for (let i = 0; i < size; i++) {
      const li = Math.round((1 - i / size) * 255) * 3;
      ctx.fillStyle = `rgb(${L[li]},${L[li + 1]},${L[li + 2]})`; ctx.fillRect(bx, oy + i, bw, 1.5);
    }
    ctx.strokeStyle = css('--rule'); ctx.strokeRect(bx + 0.5, oy + 0.5, bw - 1, size - 1);
    const tk = ticks(vmin * sc, vmax * sc, 4);
    ctx.fillStyle = css('--muted'); ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    for (const t of tk.ticks) { const f = (t / sc - vmin) / (vmax - vmin); if (f < 0 || f > 1) continue; const yy = oy + size - fwd(f) * size; ctx.fillText(fmt(t, tk.step), bx + bw + lab * 0.3, yy); }
    if (o.cbar.label) { ctx.save(); ctx.translate(w - lab * 0.2, oy + size / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.cbar.label, 0, 0); ctx.restore(); }
  }
  map.colX = colX; cv._map = map; cv._data = data; cv._N = N;
}
function maskOutline(ctx, mask, N, M, map, col, colX) {
  const h0 = N / 2, m2 = M / 2, at = (r, c) => {
    if (r < 0 || c < 0 || r >= M || c >= M) return 0;
    return mask[((h0 + m2 - r + N) % N) * N + colX(c)];
  };
  ctx.beginPath();
  for (let r = 0; r < M; r++) for (let c = 0; c < M; c++) {
    if (!at(r, c)) continue;
    const x0 = map.ox + c * map.px, y0 = map.oy + r * map.px, x1 = x0 + map.px, y1 = y0 + map.px;
    if (!at(r - 1, c)) { ctx.moveTo(x0, y0); ctx.lineTo(x1, y0); }
    if (!at(r + 1, c)) { ctx.moveTo(x0, y1); ctx.lineTo(x1, y1); }
    if (!at(r, c - 1)) { ctx.moveTo(x0, y0); ctx.lineTo(x0, y1); }
    if (!at(r, c + 1)) { ctx.moveTo(x1, y0); ctx.lineTo(x1, y1); }
  }
  ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.stroke();
}
// pixel value under the pointer, written into a caption element
V.hover = function (cv, out, f) {
  const init = out.textContent;
  cv.addEventListener('pointermove', e => {
    const m = cv._map; if (!m) return;
    const r = cv.getBoundingClientRect(), c = Math.floor((e.clientX - r.left - m.ox) / m.px), rr = Math.floor((e.clientY - r.top - m.oy) / m.px);
    if (c < 0 || rr < 0 || c >= m.M || rr >= m.M) { out.textContent = init; return; }
    const N = cv._N, x = m.colX(c), y = (N / 2 + m.M / 2 - rr + N) % N;
    out.textContent = f(cv._data[y * N + x], (x - N / 2) * m.cell, (y - N / 2) * m.cell);
  });
  cv.addEventListener('pointerleave', () => { out.textContent = init; });
};

/* ---------- line / scatter plots ----------
   o: {x: {min, max, label, log, fmt}, y: {...}, series: [{type: 'line'|'scatter'|'step'|'area', x, y, color,
   width, size, alpha, label, c (0..1 per point), cmap}], hlines: [{y, color, label}], vlines: [{x, color, label}],
   spans: [{x0, x1, color}], legend: true, equal: false} */
V.plot = function (cv, o) { register(cv, () => drawPlot(cv, o)); };
function drawPlot(cv, o) {
  const s = setup(cv); if (!s) return;
  const { ctx, w, h, fs } = s, lab = fs * 0.62;
  const ml = lab * (o.y.label ? 4.6 : 3.4), mb = lab * (o.x.label ? 3.0 : 1.9), mt = lab * 0.9, mr = lab * (o.mr || 1.1);
  let pw = w - ml - mr, ph = h - mt - mb;
  const X = Object.assign({}, o.x), Y = Object.assign({}, o.y);
  let ox = ml, oy = mt;
  if (o.equal) {                           // same units per pixel on both axes
    const kx = pw / (X.max - X.min), ky = ph / (Y.max - Y.min), k = Math.min(kx, ky);
    const nw = (X.max - X.min) * k, nh = (Y.max - Y.min) * k;
    ox += (pw - nw) / 2; oy += (ph - nh) / 2; pw = nw; ph = nh;
  }
  const sx = mk(X, ox, ox + pw), sy = mk(Y, oy + ph, oy);
  ctx.font = `${lab}px ${css('--f-body')}`;
  // grid & ticks
  const lo = a => Math.min(a.min, a.max), hi = a => Math.max(a.min, a.max);
  const gx = X.ticks ? X.ticks.slice() : X.log ? logTicks(X.min, X.max) : ticks(lo(X), hi(X), X.n || Math.max(2, pw / (lab * 7))).ticks;
  const gy = Y.ticks ? Y.ticks.slice() : Y.log ? logTicks(Y.min, Y.max) : ticks(lo(Y), hi(Y), Y.n || Math.max(2, ph / (lab * 3.4))).ticks;
  const sxStep = X.log ? 0 : niceStep(hi(X) - lo(X), X.n || Math.max(2, pw / (lab * 7))), syStep = Y.log ? 0 : niceStep(hi(Y) - lo(Y), Y.n || Math.max(2, ph / (lab * 3.4)));
  const bare = o.axes === false;
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1;
  ctx.beginPath();
  if (!bare) for (const t of gx) { const p = Math.round(sx(t)) + 0.5; ctx.moveTo(p, oy); ctx.lineTo(p, oy + ph); }
  if (!bare) for (const t of gy) { const p = Math.round(sy(t)) + 0.5; ctx.moveTo(ox, p); ctx.lineTo(ox + pw, p); }
  ctx.stroke();
  if (bare) { gx.length = 0; gy.length = 0; }
  ctx.fillStyle = css('--muted');
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (const t of gx) ctx.fillText(X.fmt ? X.fmt(t) : X.log ? fmtLog(t) : fmt(t, sxStep), sx(t), oy + ph + lab * 0.4);
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (const t of gy) ctx.fillText(Y.fmt ? Y.fmt(t) : Y.log ? fmtLog(t) : fmt(t, syStep), ox - lab * 0.4, sy(t));
  ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
  if (X.label) ctx.fillText(X.label, ox + pw / 2, oy + ph + lab * 2.6);
  if (Y.label) { ctx.save(); ctx.translate(lab * 0.9, oy + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(Y.label, 0, -lab * 0.55); ctx.restore(); }
  if (!bare) { ctx.strokeStyle = css('--axis'); ctx.beginPath(); ctx.moveTo(ox, oy + ph + 0.5); ctx.lineTo(ox + pw, oy + ph + 0.5); ctx.stroke(); }
  // data
  ctx.save(); ctx.beginPath(); ctx.rect(ox, oy, pw, ph); ctx.clip();
  for (const sp of o.spans || []) { ctx.fillStyle = color(sp.color || '--accent-soft'); ctx.fillRect(sx(sp.x0), oy, sx(sp.x1) - sx(sp.x0), ph); }
  (o.series || []).forEach((se, i) => drawSeries(ctx, se, sx, sy, i, oy + ph));
  ctx.lineWidth = 1.5;
  for (const hl of o.hlines || []) {
    ctx.strokeStyle = color(hl.color || '--ink-2'); ctx.setLineDash(hl.dash || []); ctx.beginPath();
    ctx.moveTo(ox, sy(hl.y)); ctx.lineTo(ox + pw, sy(hl.y)); ctx.stroke(); ctx.setLineDash([]);
  }
  for (const vl of o.vlines || []) {
    ctx.strokeStyle = color(vl.color || '--ink-2'); ctx.setLineDash(vl.dash || []); ctx.beginPath();
    ctx.moveTo(sx(vl.x), oy); ctx.lineTo(sx(vl.x), oy + ph); ctx.stroke(); ctx.setLineDash([]);
  }
  ctx.restore();
  ctx.fillStyle = css('--ink-2'); ctx.textBaseline = 'bottom'; ctx.textAlign = 'left';
  for (const hl of o.hlines || []) if (hl.label) ctx.fillText(hl.label, ox + lab * 0.4, sy(hl.y) - lab * 0.2);
  ctx.textBaseline = 'top';
  for (const vl of o.vlines || []) if (vl.label) ctx.fillText(vl.label, sx(vl.x) + lab * 0.3, oy + lab * 0.2);
  // legend
  const items = (o.series || []).map((se, i) => [se, i]).filter(([se]) => se.label);
  if (o.legend !== false && items.length) {
    const lw = Math.max.apply(null, items.map(([se]) => ctx.measureText(se.label).width)) + lab * 2.6, lh = lab * 1.35;
    let lx = o.legendPos === 'left' ? ox + lab * 0.5 : ox + pw - lw - lab * 0.4, ly = oy + lab * 0.4;
    ctx.fillStyle = css('--panel'); ctx.globalAlpha = 0.86; ctx.fillRect(lx - lab * 0.3, ly - lab * 0.2, lw + lab * 0.3, lh * items.length + lab * 0.3); ctx.globalAlpha = 1;
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    items.forEach(([se, i], j) => {
      const yy = ly + lh * j + lh / 2, c = color(se.color || `--s${(i % 8) + 1}`);
      ctx.strokeStyle = c; ctx.fillStyle = c; ctx.lineWidth = 2.5;
      if (se.type === 'scatter') { ctx.beginPath(); ctx.arc(lx + lab * 0.8, yy, lab * 0.3, 0, 2 * Math.PI); ctx.fill(); }
      else { ctx.setLineDash(se.dash || []); ctx.beginPath(); ctx.moveTo(lx, yy); ctx.lineTo(lx + lab * 1.6, yy); ctx.stroke(); ctx.setLineDash([]); }
      ctx.fillStyle = css('--ink-2'); ctx.fillText(se.label, lx + lab * 2.1, yy);
    });
  }
  cv._sx = sx; cv._sy = sy; cv._box = { ox, oy, pw, ph };
}
function mk(ax, a, b) {
  if (ax.log) { const l0 = Math.log10(ax.min), l1 = Math.log10(ax.max); return v => a + (Math.log10(Math.max(v, 1e-30)) - l0) / (l1 - l0) * (b - a); }
  return v => a + (v - ax.min) / (ax.max - ax.min) * (b - a);
}
function drawSeries(ctx, se, sx, sy, i, base) {
  const c = color(se.color || `--s${(i % 8) + 1}`), x = se.x, y = se.y, n = y.length;
  ctx.globalAlpha = se.alpha == null ? 1 : se.alpha;
  if (se.type === 'bars') {
    const w = Math.min(24, se.barWidth || 18), b0 = sy(0);
    ctx.fillStyle = c;
    for (let k = 0; k < n; k++) {
      const X = sx(x[k]), Y = sy(y[k]), top = Math.min(Y, b0), h = Math.abs(b0 - Y), r = Math.min(4, h), dn = Y > b0 ? -1 : 1;   // dn = -1: negative bar, rounded at the bottom
      ctx.beginPath(); ctx.moveTo(X - w / 2, b0); ctx.lineTo(X - w / 2, Y + dn * r); ctx.quadraticCurveTo(X - w / 2, Y, X - w / 2 + r, Y);
      ctx.lineTo(X + w / 2 - r, Y); ctx.quadraticCurveTo(X + w / 2, Y, X + w / 2, Y + dn * r); ctx.lineTo(X + w / 2, b0); ctx.closePath(); ctx.fill();
      if (se.labels) { ctx.globalAlpha = 1; ctx.fillStyle = css('--ink-2'); ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(se.labels[k], X, top - 3); ctx.fillStyle = c; }
    }
  } else if (se.type === 'scatter') {
    const r = se.size || 2.2, L = se.c ? CMAPS[se.cmap || 'blues'] : null;
    ctx.fillStyle = c;
    for (let k = 0; k < n; k++) {
      if (y[k] !== y[k]) continue;
      if (L) { const li = Math.round(Math.max(0, Math.min(1, se.c[k])) * 255) * 3; ctx.fillStyle = `rgb(${L[li]},${L[li + 1]},${L[li + 2]})`; }
      if (r >= 6) {                                   // highlighted markers: filled circle with a surface ring
        ctx.beginPath(); ctx.arc(sx(x[k]), sy(y[k]), r / 2, 0, 2 * Math.PI); ctx.fill();
        ctx.lineWidth = 2; ctx.strokeStyle = css('--panel'); ctx.stroke();
      } else ctx.fillRect(sx(x[k]) - r / 2, sy(y[k]) - r / 2, r, r);
    }
  } else {
    ctx.strokeStyle = c; ctx.lineWidth = se.width || 2; ctx.lineJoin = ctx.lineCap = 'round';
    ctx.setLineDash(se.dash || []);
    ctx.beginPath();
    let pen = false;
    for (let k = 0; k < n; k++) {
      if (y[k] !== y[k] || y[k] == null) { pen = false; continue; }
      const X = sx(x[k]), Y = sy(y[k]);
      if (se.type === 'step') {
        const xl = k > 0 ? sx((x[k - 1] + x[k]) / 2) : X - (sx(x[1]) - X) / 2, xr = k < n - 1 ? sx((x[k] + x[k + 1]) / 2) : X + (X - sx(x[n - 2])) / 2;
        if (!pen) { ctx.moveTo(xl, Y); pen = true; } else ctx.lineTo(xl, Y);
        ctx.lineTo(xr, Y);
      } else if (!pen) { ctx.moveTo(X, Y); pen = true; } else ctx.lineTo(X, Y);
    }
    if (se.type === 'area') {
      ctx.stroke();
      ctx.lineTo(sx(x[n - 1]), base); ctx.lineTo(sx(x[0]), base); ctx.closePath();
      ctx.globalAlpha = 0.12; ctx.fillStyle = c; ctx.fill();
    } else ctx.stroke();
    ctx.setLineDash([]);
    if (se.dots) { ctx.fillStyle = c; for (let k = 0; k < n; k++) if (y[k] === y[k]) { ctx.beginPath(); ctx.arc(sx(x[k]), sy(y[k]), se.dots, 0, 2 * Math.PI); ctx.fill(); } }
  }
  ctx.globalAlpha = 1;
}

/* ---------- the live CASA call: code you can turn into knobs ----------
   params: [{k, v} static (or vf(vals) for a computed value) | {k, type: 'range', min, max, step, value, fmt, scale: 'log', live}
            | {k, type: 'choice', options, value, labels, quote, fmt} | {k, type: 'toggle', value}], each with optional show(vals), note.
   scenario: true puts a knob below the call as a comment: it sets up the simulation and is not a CASA parameter. */
let uid = 0;
V.call = function (el, fn, params, onChange) {
  el.classList.add('call');
  el.innerHTML = '';
  const vals = {}, rows = {}, pad = ' '.repeat(fn.length + 1);
  const real = params.filter(p => !p.scenario), scen = params.filter(p => p.scenario);
  const lineText = p => {
    const v = vals[p.k];
    if (p.type && p.fmt) return p.fmt(v);
    if (p.type === 'choice') return p.quote === false ? String(v) : `'${v}'`;
    if (p.type === 'toggle') return v ? 'True' : 'False';
    if (p.type === 'range') return String(v);
    return p.vf ? p.vf(vals) : p.v;
  };
  const fire = (k, final) => {
    params.forEach(p => { if (p.show) rows[p.k].hidden = !p.show(vals); });
    // the closing parenthesis goes on the last visible parameter, so the call stays valid Python
    const shown = real.filter(p => !rows[p.k].hidden);
    shown.forEach((p, j) => { rows[p.k]._p.textContent = j === shown.length - 1 ? ')' : ','; });
    params.forEach(p => { if (rows[p.k]._val) rows[p.k]._val.textContent = lineText(p); });
    if (k) { const r = rows[k]; r.classList.remove('changed'); void r.offsetWidth; r.classList.add('changed'); }
    onChange && onChange(Object.assign({}, vals), k, final);
  };
  const span = (cls, t) => { const e = document.createElement('span'); if (cls) e.className = cls; e.textContent = t; return e; };
  const build = (p, i, isScen) => {
    const row = document.createElement('div'); row.className = 'call-row' + (p.type ? ' live' : '') + (isScen ? ' scen' : '');
    const code = document.createElement('code'), val = span('v', '');
    if (isScen) code.append(span('c', '# '), span('k', p.label || p.k), document.createTextNode(' = '), val);
    else { row._p = span('p', i === real.length - 1 ? ')' : ','); code.append(span('lead', i === 0 ? fn + '(' : pad), span('k', p.k), document.createTextNode('='), val, row._p); }
    if (p.note) code.appendChild(span('c', '  # ' + p.note));
    row.appendChild(code);
    row._val = val;
    if (p.type === 'range') {
      const inp = document.createElement('input');
      inp.type = 'range'; inp.id = 'k' + (++uid); inp.setAttribute('aria-label', p.label || p.k);
      const toPos = v => (p.scale === 'log' ? Math.log10(v) : v), fromPos = x => (p.scale === 'log' ? Math.pow(10, +x) : +x);
      inp.min = toPos(p.min); inp.max = toPos(p.max); inp.step = p.scale === 'log' ? (toPos(p.max) - toPos(p.min)) / 200 : p.step; inp.value = toPos(p.value);
      vals[p.k] = p.value;
      const snap = v => (p.snap ? p.snap(v) : p.scale === 'log' ? +v.toPrecision(2) : v);
      inp.addEventListener('input', () => { vals[p.k] = snap(fromPos(inp.value)); fire(p.k, p.live !== false); });
      inp.addEventListener('change', () => { vals[p.k] = snap(fromPos(inp.value)); fire(p.k, true); });
      row._input = inp; row.appendChild(inp);
    } else if (p.type === 'choice' || p.type === 'toggle') {
      const opts = p.type === 'toggle' ? [true, false] : p.options, labels = p.type === 'toggle' ? ['True', 'False'] : p.labels;
      vals[p.k] = p.type === 'toggle' ? !!p.value : p.value;
      const seg = document.createElement('div'); seg.className = 'seg'; seg.setAttribute('role', 'group'); seg.setAttribute('aria-label', p.label || p.k);
      opts.forEach((opt, j) => {
        const b = document.createElement('button'); b.type = 'button'; b.textContent = labels ? labels[j] : String(opt);
        b.setAttribute('aria-pressed', String(opt === vals[p.k]));
        b.addEventListener('click', () => { vals[p.k] = opt; seg.querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', String(x === b))); fire(p.k, true); });
        seg.appendChild(b);
      });
      row.appendChild(seg);
    } else if (p.vf || p.scenario) vals[p.k] = p.v;
    rows[p.k] = row; el.appendChild(row);
  };
  real.forEach((p, i) => build(p, i, false));
  if (scen.length) { const h = document.createElement('div'); h.className = 'call-row'; const c = document.createElement('code'); c.appendChild(span('c', '# scenario, not CASA parameters:')); h.appendChild(c); el.appendChild(h); }
  scen.forEach((p, i) => build(p, i, true));
  fire(null, true);
  return { vals, fire };
};

/* ---------- readout tiles ---------- */
V.readouts = function (el, items) {
  el.classList.add('readouts'); el.innerHTML = '';
  const dd = {};
  items.forEach(it => {
    const d = document.createElement('div'), t = document.createElement('dt'), v = document.createElement('dd');
    t.textContent = it.label; v.textContent = '…'; d.append(t, v); el.appendChild(d); dd[it.k] = v;
    if (it.hint) d.title = it.hint;
  });
  return { set(k, text, state) { const v = dd[k]; if (!v) return; v.textContent = text; v.dataset.state = state || ''; } };
};

/* ---------- small helpers ---------- */
V.rafThrottle = function (fn) {
  let pending = null;
  return function (...args) { const first = pending === null; pending = args; if (first) requestAnimationFrame(() => { const a = pending; pending = null; fn.apply(null, a); }); };
};
V.fmtFlux = function (jy) {                      // Jy -> readable string
  const a = Math.abs(jy);
  if (a >= 1) return jy.toFixed(a >= 10 ? 1 : 2) + ' Jy';
  if (a >= 1e-3) return (jy * 1e3).toFixed(a >= 1e-2 ? 1 : 2) + ' mJy';
  return (jy * 1e6).toFixed(a >= 1e-5 ? 1 : 2) + ' μJy';
};
V.$ = (root, sel) => root.querySelector(sel);
V.$$ = (root, sel) => Array.from(root.querySelectorAll(sel));

window.Viz = V;
})();
