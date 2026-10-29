/* sim.js: a compact radio-interferometry simulator for the workshop decks. No dependencies.
   Conventions: sky offsets in arcsec (x = +East, y = +North); u, v in wavelengths; flux in Jy;
   time in seconds from transit; images are row-major [y*N + x] with x -> East, y -> North,
   the image centre at pixel (N/2, N/2). Display code flips x so that East is on the left. */
(function (root) {
'use strict';
const S = {};
const C = 299792458, AS = Math.PI / 648000, DEG = Math.PI / 180, LAT = -23.0229 * DEG, KB = 1.380649e-23;
const SID = 2 * Math.PI / 86164.09;          // Earth rotation, rad per second
S.C = C; S.AS = AS; S.DEG = DEG; S.LAT = LAT;

/* ---------- seeded random numbers (demos must look the same every time) ---------- */
S.rng = function (seed) {
  let a = (seed >>> 0) || 1, spare = null;
  const r = function () {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  r.normal = function () {
    if (spare !== null) { const s = spare; spare = null; return s; }
    let u = 0; while (u === 0) u = r();
    const v = r(), m = Math.sqrt(-2 * Math.log(u));
    spare = m * Math.sin(2 * Math.PI * v);
    return m * Math.cos(2 * Math.PI * v);
  };
  return r;
};

/* ---------- FFT (radix 2, in place, unnormalised). sign -1: e^{-i}, +1: e^{+i} ---------- */
const plans = {};
function plan(n) {
  if (plans[n]) return plans[n];
  const bits = Math.round(Math.log2(n)), rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) { let r = 0; for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b); rev[i] = r; }
  const cos = new Float64Array(n / 2), sin = new Float64Array(n / 2);
  for (let i = 0; i < n / 2; i++) { cos[i] = Math.cos(2 * Math.PI * i / n); sin[i] = Math.sin(2 * Math.PI * i / n); }
  return (plans[n] = { rev, cos, sin });
}
function fft1(re, im, n, sign) {
  const p = plan(n), rev = p.rev, cos = p.cos, sin = p.sin;
  for (let i = 0; i < n; i++) {
    const j = rev[i];
    if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (let size = 2; size <= n; size <<= 1) {
    const half = size >> 1, step = n / size;
    for (let i = 0; i < n; i += size) {
      for (let j = 0, k = 0; j < half; j++, k += step) {
        const wr = cos[k], wi = sign * sin[k], a = i + j, b = a + half;
        const xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
        re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
      }
    }
  }
}
S.fft1 = fft1;
S.fft2 = function (re, im, nx, ny, sign) {
  const m = Math.max(nx, ny), tr = new Float64Array(m), ti = new Float64Array(m);
  for (let y = 0; y < ny; y++) {
    const o = y * nx;
    for (let x = 0; x < nx; x++) { tr[x] = re[o + x]; ti[x] = im[o + x]; }
    fft1(tr, ti, nx, sign);
    for (let x = 0; x < nx; x++) { re[o + x] = tr[x]; im[o + x] = ti[x]; }
  }
  for (let x = 0; x < nx; x++) {
    for (let y = 0; y < ny; y++) { tr[y] = re[y * nx + x]; ti[y] = im[y * nx + x]; }
    fft1(tr, ti, ny, sign);
    for (let y = 0; y < ny; y++) { re[y * nx + x] = tr[y]; im[y * nx + x] = ti[y]; }
  }
};
function checker(re, im, nx, ny) {
  for (let y = 0; y < ny; y++) for (let x = 1 - (y & 1); x < nx; x += 2) { const i = y * nx + x; re[i] = -re[i]; im[i] = -im[i]; }
}
// centred transform: the origin of both domains sits at index (n/2, n/2); n divisible by 4
S.fft2c = function (re, im, nx, ny, sign) { checker(re, im, nx, ny); S.fft2(re, im, nx, ny, sign); checker(re, im, nx, ny); };

/* ---------- Bessel functions (Numerical Recipes rational approximations) ---------- */
function j0(x) {
  const ax = Math.abs(x);
  if (ax < 8) {
    const y = x * x;
    return (57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))))) /
           (57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y)))));
  }
  const z = 8 / ax, y = z * z, xx = ax - 0.785398164;
  const p = 1 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
  const q = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934945152e-7)));
  return Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * p - z * Math.sin(xx) * q);
}
function j1(x) {
  const ax = Math.abs(x);
  if (ax < 8) {
    const y = x * x;
    return x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606)))))) /
           (144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y)))));
  }
  const z = 8 / ax, y = z * z, xx = ax - 2.356194491;
  const p = 1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
  const q = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
  const a = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * p - z * Math.sin(xx) * q);
  return x < 0 ? -a : a;
}
S.j0 = j0; S.j1 = j1;

/* ---------- ALMA-like antenna layouts (metres East, North), tuned to match the
   5th / 80th percentile baselines of the real configurations ---------- */
const CONFIGS = {"C-2":{"e":[0.2,-13.9,-7.6,15.3,-29.6,6.9,35.2,-20.4,-18.4,17.1,-49.4,-2.6,35.5,-46.9,23.0,23.8,-39.1,11.0,2.1,-72.4,62.3,51.1,-68.0,28.5,4.5,-67.6,45.5,-10.9,-96.2,83.5,38.3,-88.9,65.9,-9.9,-110.9,23.0,-57.8,-44.0,120.7,95.9,-115.8,123.3,57.7],"n":[-1.5,4.2,-16.5,9.8,9.0,-20.7,10.7,25.9,-31.6,29.5,-8.0,-36.3,40.4,8.3,-45.7,52.6,39.5,-56.6,48.8,4.5,-48.7,48.6,-13.9,-59.9,73.1,-38.2,-84.0,76.8,10.8,-38.3,67.1,-47.8,-64.2,99.3,-11.3,-115.5,128.4,-135.4,-35.1,114.1,-68.0,-83.6,159.5]},"C-4":{"e":[15.3,-3.2,20.1,23.0,1.8,36.5,2.7,38.1,64.2,-27.0,11.5,7.4,-62.5,16.8,39.9,-61.9,126.6,-17.7,27.6,110.1,-88.2,40.5,118.2,-178.2,12.5,175.1,-187.3,149.6,47.5,-228.4,49.7,-22.3,-247.6,176.7,-260.0,-247.3,327.1,-220.5,-127.2,336.6,-306.6,-34.3,345.4],"n":[10.2,17.9,-4.9,27.7,-10.9,-8.3,48.8,-37.9,12.7,48.2,-47.1,102.0,36.3,-69.3,83.9,-43.9,-54.0,123.2,-128.3,-45.3,123.6,-188.9,166.6,40.1,-193.3,14.2,-56.7,-183.4,270.4,-54.6,-239.3,283.6,-88.2,-208.3,263.2,-186.9,54.5,310.7,-274.3,253.6,283.6,-330.3,-120.6]},"C-6":{"e":[17.3,3.8,14.0,31.7,-40.6,64.5,27.0,-84.7,132.9,-50.5,71.3,164.2,-184.8,121.4,209.0,-213.1,345.1,59.4,-191.2,460.5,0.6,-163.9,502.8,-382.7,-365.1,-189.9,-515.4,276.0,-143.8,-741.5,245.8,371.8,-527.8,1065.5,-527.3,-893.2,841.6,-123.1,-901.3,1173.8,-685.3,48.2,676.9],"n":[-41.0,-29.8,-66.4,-6.7,-23.6,-73.9,34.8,-96.5,-105.7,85.3,-234.0,102.1,12.0,-220.6,80.7,-53.1,-290.9,300.3,-422.3,20.0,265.7,-509.0,38.6,255.1,-259.6,618.8,249.5,-644.5,653.5,-128.3,-735.5,636.0,-821.9,82.8,867.4,-50.0,-503.6,925.2,-226.0,-136.3,892.3,-1420.5,979.8]},"ACA":{"e":[0.3,9.6,3.1,-11.0,15.4,-6.3,4.5,15.3,-14.2,-6.0,9.3,-20.0],"n":[-2.8,-3.4,8.7,-7.1,11.2,11.7,-15.1,-11.8,6.4,-16.4,26.3,-7.6]}};
const CONFIG_INFO = {
  'C-2': { label: 'C-2-like (≈0.3 km)', D: 12 },
  'C-4': { label: 'C-4-like (≈0.8 km)', D: 12 },
  'C-6': { label: 'C-6-like (≈2.5 km)', D: 12 },
  'ACA': { label: '7-m ACA', D: 7 }
};
S.config = function (name) {
  const c = CONFIGS[name];
  return { name, label: CONFIG_INFO[name].label, D: CONFIG_INFO[name].D, e: c.e, n: c.n, nAnt: c.e.length,
           names: c.e.map((_, i) => (CONFIG_INFO[name].D === 7 ? 'CM' + String(i + 1).padStart(2, '0') : i < 25 ? 'DA' + (41 + i) : 'DV' + String(i - 24).padStart(2, '0'))) };
};
S.mergeConfigs = function (a, b) {
  return { name: a.name + '+' + b.name, label: a.label + ' + ' + b.label, e: a.e.concat(b.e), n: a.n.concat(b.n),
           nAnt: a.nAnt + b.nAnt, names: a.names.concat(b.names), dish: a.e.map(() => a.D).concat(b.e.map(() => b.D)), D: a.D };
};
S.baselineStats = function (cfg) {
  const L = [];
  for (let i = 0; i < cfg.nAnt; i++) for (let j = i + 1; j < cfg.nAnt; j++) L.push(Math.hypot(cfg.e[j] - cfg.e[i], cfg.n[j] - cfg.n[i]));
  L.sort((a, b) => a - b);
  const pc = p => L[Math.min(L.length - 1, Math.floor(p / 100 * L.length))];
  return { lmin: L[0], l5: pc(5), l80: pc(80), lmax: L[L.length - 1], lengths: L };
};
// ALMA Technical Handbook rules of thumb: resolution ~0.574 lambda/L80, MRS ~0.983 lambda/L5 (arcsec)
S.resolution = (cfg, ghz) => 0.574 * (C / ghz / 1e9) / S.baselineStats(cfg).l80 / AS;
S.mrs = (cfg, ghz) => 0.983 * (C / ghz / 1e9) / S.baselineStats(cfg).l5 / AS;
S.pbFWHM = (ghz, D) => 1.13 * (C / ghz / 1e9) / (D || 12) / AS;            // arcsec
S.pb = (r, ghz, D) => { const f = S.pbFWHM(ghz, D); return Math.exp(-4 * Math.LN2 * r * r / (f * f)); };

/* ---------- receivers & noise (round numbers, good enough to compare cases) ---------- */
S.BANDS = {
  3: { ghz: 100, tsys: 75, eta: 0.71 }, 4: { ghz: 145, tsys: 85, eta: 0.70 },
  6: { ghz: 230, tsys: 100, eta: 0.68 }, 7: { ghz: 345, tsys: 150, eta: 0.63 },
  9: { ghz: 680, tsys: 900, eta: 0.45 }
};
S.sefd = (band, D) => { const b = S.BANDS[band], A = Math.PI * Math.pow((D || 12) / 2, 2) * b.eta; return 2 * KB * b.tsys / A * 1e26; };
// noise per real/imag component of one Stokes-I visibility (dual pol), Jy
S.sigmaVis = (band, bwHz, tint, D1, D2) => {
  const s1 = S.sefd(band, D1 || 12), s2 = S.sefd(band, D2 || D1 || 12);
  return Math.sqrt(s1 * s2) / (0.88 * Math.sqrt(2 * bwHz * tint)) / Math.SQRT2;
};
// image noise for natural weighting: sigma / sqrt(n_vis)
S.sigmaImage = (band, bwHz, tTot, nAnt, D) => S.sefd(band, D) / (0.88 * Math.sqrt(2 * nAnt * (nAnt - 1) * bwHz * tTot));

/* ---------- uv tracks ---------- */
S.uvTracks = function (cfg, o) {
  const lam = C / (o.freq * 1e9);
  const sd = Math.sin(o.dec * DEG), cd = Math.cos(o.dec * DEG), sl = Math.sin(LAT), cl = Math.cos(LAT);
  const elmin = (o.elmin == null ? 15 : o.elmin) * DEG;
  let times = o.times;
  if (!times) { times = []; const dt = o.dt || 60; for (let t = o.ha0 * 3600; t <= o.ha1 * 3600 + 1e-6; t += dt) times.push(t); }
  const keep = [], scanOf = [];
  times.forEach((t, i) => {
    if (Math.asin(sl * sd + cl * cd * Math.cos(t * SID)) >= elmin) { keep.push(t); scanOf.push(o.scans ? o.scans[i] : 0); }
  });
  const na = cfg.nAnt, nb = na * (na - 1) / 2, nt = keep.length, n = nb * nt;
  const bx = new Float64Array(nb), by = new Float64Array(nb), bz = new Float64Array(nb), b1 = new Uint8Array(nb), b2 = new Uint8Array(nb);
  let b = 0;
  for (let i = 0; i < na; i++) for (let j = i + 1; j < na; j++) {
    const dE = cfg.e[j] - cfg.e[i], dN = cfg.n[j] - cfg.n[i];
    bx[b] = -sl * dN; by[b] = dE; bz[b] = cl * dN; b1[b] = i; b2[b] = j; b++;
  }
  const u = new Float64Array(n), v = new Float64Array(n), t = new Float64Array(n), a1 = new Uint8Array(n), a2 = new Uint8Array(n),
        ti = new Uint32Array(n), scan = new Uint16Array(n);
  let k = 0;
  for (let it = 0; it < nt; it++) {
    const H = keep[it] * SID, sH = Math.sin(H), cH = Math.cos(H);
    for (b = 0; b < nb; b++, k++) {
      u[k] = (sH * bx[b] + cH * by[b]) / lam;
      v[k] = (-sd * cH * bx[b] + sd * sH * by[b] + cd * bz[b]) / lam;
      t[k] = keep[it]; a1[k] = b1[b]; a2[k] = b2[b]; ti[k] = it; scan[k] = scanOf[it];
    }
  }
  return { n, u, v, t, a1, a2, ti, scan, nt, nb, nAnt: na, times: keep, lam, freq: o.freq, dec: o.dec, cfg };
};

/* ---------- analytic sky models ----------
   component: {type: 'point'|'gauss'|'disk'|'ring', flux, x, y, fwhm (gauss FWHM / disk diameter),
               r (ring radius), width (ring FWHM), ratio (minor/major = cos i), pa (deg E of N)} */
S.visModel = function (comps, u, v, re, im) {
  const n = u.length, k = Math.PI * Math.PI / (4 * Math.LN2);
  for (const c of comps) {
    const l = (c.x || 0) * AS, m = (c.y || 0) * AS, pa = (c.pa || 0) * DEG, sp = Math.sin(pa), cp = Math.cos(pa);
    const q2 = c.ratio == null ? 1 : c.ratio * c.ratio, a = (c.fwhm || 0) * AS, R = (c.r || 0) * AS, w = (c.width || 0) * AS;
    for (let i = 0; i < n; i++) {
      const uu = u[i], vv = v[i];
      let amp = c.flux;
      if (c.type !== 'point') {
        const qa = uu * sp + vv * cp, qb = uu * cp - vv * sp, q = Math.sqrt(qa * qa + q2 * qb * qb);
        if (c.type === 'gauss') amp *= Math.exp(-k * a * a * q * q);
        else if (c.type === 'disk') { const z = Math.PI * a * q; amp *= z < 1e-9 ? 1 : 2 * j1(z) / z; }
        else amp *= j0(2 * Math.PI * R * q) * Math.exp(-k * w * w * q * q);
      }
      const ph = -2 * Math.PI * (uu * l + vv * m);
      re[i] += amp * Math.cos(ph); im[i] += amp * Math.sin(ph);
    }
  }
};
// the same components rendered as an image (Jy/pixel), for "true sky" panels
S.skyImage = function (comps, N, cell) {
  const img = new Float64Array(N * N), h = N / 2, sub = 3;
  for (const c of comps) {
    const pa = (c.pa || 0) * DEG, sp = Math.sin(pa), cp = Math.cos(pa), ratio = c.ratio == null ? 1 : c.ratio;
    if (c.type === 'point') {
      const x = Math.round(h + (c.x || 0) / cell), y = Math.round(h + (c.y || 0) / cell);
      if (x >= 0 && x < N && y >= 0 && y < N) img[y * N + x] += c.flux;
      continue;
    }
    const ext = c.type === 'ring' ? c.r + 2 * c.width : (c.type === 'disk' ? c.fwhm / 2 : 1.5 * c.fwhm);
    const r = Math.ceil(ext / cell) + 1, x0 = Math.round(h + (c.x || 0) / cell), y0 = Math.round(h + (c.y || 0) / cell);
    let sum = 0;
    const tmp = [];
    for (let y = y0 - r; y <= y0 + r; y++) for (let x = x0 - r; x <= x0 + r; x++) {
      if (x < 0 || x >= N || y < 0 || y >= N) continue;
      let val = 0;
      for (let sy = 0; sy < sub; sy++) for (let sx = 0; sx < sub; sx++) {
        const dl = (x - h + (sx + 0.5) / sub - 0.5) * cell - (c.x || 0), dm = (y - h + (sy + 0.5) / sub - 0.5) * cell - (c.y || 0);
        const pa_ = dl * sp + dm * cp, pb_ = (dl * cp - dm * sp) / ratio, rr = Math.hypot(pa_, pb_);
        if (c.type === 'gauss') val += Math.exp(-4 * Math.LN2 * rr * rr / (c.fwhm * c.fwhm));
        else if (c.type === 'disk') val += rr <= c.fwhm / 2 ? 1 : 0;
        else val += Math.exp(-4 * Math.LN2 * Math.pow(rr - c.r, 2) / (c.width * c.width));
      }
      tmp.push(y * N + x, val); sum += val;
    }
    if (sum > 0) for (let i = 0; i < tmp.length; i += 2) img[tmp[i]] += c.flux * tmp[i + 1] / sum;
  }
  return img;
};

/* ---------- observe: model visibilities x antenna gains + thermal noise ----------
   o.gains: {re, im} Float64Array(nt*nAnt) indexed [ti*nAnt + ant]; o.sigma: Jy per component */
S.observe = function (tr, comps, o) {
  o = o || {};
  const n = tr.n, re = new Float64Array(n), im = new Float64Array(n), wt = new Float64Array(n);
  S.visModel(comps, tr.u, tr.v, re, im);
  const d = Object.assign({}, tr, { re, im, wt, flag: new Uint8Array(n) });
  if (o.keepModel) { d.mre = Float64Array.from(re); d.mim = Float64Array.from(im); }
  if (o.gains) S.corrupt(d, o.gains);
  const r = S.rng(o.seed || 1), sig = o.sigma || 0;
  for (let k = 0; k < n; k++) {
    const s = typeof sig === 'number' ? sig : sig[k];
    if (s) { re[k] += s * r.normal(); im[k] += s * r.normal(); wt[k] = 1 / (s * s); } else wt[k] = 1;
  }
  return d;
};
// multiply visibilities by g_a1 * conj(g_a2) (or divide, inverse=true)
S.corrupt = function (d, g, inverse) {
  const na = d.nAnt;
  for (let k = 0; k < d.n; k++) {
    const i1 = d.ti[k] * na + d.a1[k], i2 = d.ti[k] * na + d.a2[k];
    let gr = g.re[i1] * g.re[i2] + g.im[i1] * g.im[i2], gi = g.im[i1] * g.re[i2] - g.re[i1] * g.im[i2];
    if (inverse) { const m = gr * gr + gi * gi; if (m < 1e-12) { d.flag[k] = 1; continue; } gr /= m; gi = -gi / m; }
    const vr = d.re[k], vi = d.im[k];
    d.re[k] = gr * vr - gi * vi; d.im[k] = gr * vi + gi * vr;
  }
};
S.copyData = function (d) {
  return Object.assign({}, d, { re: Float64Array.from(d.re), im: Float64Array.from(d.im), wt: Float64Array.from(d.wt), flag: Uint8Array.from(d.flag) });
};

/* ---------- imaging: weighting, gridding, FFT, beam fit ----------
   o: {N, cell (arcsec), weighting: 'natural'|'uniform'|'briggs', robust, taper (arcsec FWHM), keepGrid} */
S.image = function (d, o) {
  const N = o.N, h = N / 2, nv = d.n, du = 1 / (N * o.cell * AS);
  const cell1 = new Int32Array(nv), cell2 = new Int32Array(nv), Wg = new Float64Array(N * N);
  let dropped = 0;
  for (let k = 0; k < nv; k++) {
    if (d.flag[k]) { cell1[k] = -1; continue; }
    const iu = Math.round(d.u[k] / du), iv = Math.round(d.v[k] / du);
    if (iu <= -h || iu >= h || iv <= -h || iv >= h) { cell1[k] = -1; dropped++; continue; }
    const c1 = (iv + h) * N + iu + h, c2 = (h - iv) * N + h - iu;
    cell1[k] = c1; cell2[k] = c2; Wg[c1] += d.wt[k]; Wg[c2] += d.wt[k];
  }
  const wmode = o.weighting || 'natural';
  let f2 = 0;
  if (wmode === 'briggs') {
    let sW2 = 0, sw = 0;
    for (let c = 0; c < N * N; c++) sW2 += Wg[c] * Wg[c];
    for (let k = 0; k < nv; k++) if (cell1[k] >= 0) sw += 2 * d.wt[k];
    f2 = Math.pow(5 * Math.pow(10, -(o.robust || 0)), 2) / (sW2 / sw);
  }
  const tap = o.taper ? Math.PI * Math.PI * Math.pow(o.taper * AS, 2) / (4 * Math.LN2) : 0;
  const wi = new Float64Array(nv);
  let sumw = 0, sumw2s = 0;
  const Gr = new Float64Array(N * N), Gi = new Float64Array(N * N), Pr = new Float64Array(N * N), Pi = new Float64Array(N * N);
  for (let k = 0; k < nv; k++) {
    const c1 = cell1[k]; if (c1 < 0) continue;
    let w = d.wt[k];
    if (wmode === 'uniform') w /= Wg[c1];
    else if (wmode === 'briggs') w /= 1 + Wg[c1] * f2;
    if (tap) w *= Math.exp(-tap * (d.u[k] * d.u[k] + d.v[k] * d.v[k]));
    wi[k] = w; sumw += w; sumw2s += w * w / d.wt[k];
    const c2 = cell2[k], vr = w * d.re[k], vi = w * d.im[k];
    Gr[c1] += vr; Gi[c1] += vi; Gr[c2] += vr; Gi[c2] -= vi; Pr[c1] += w; Pr[c2] += w;
  }
  const grid = o.keepGrid ? Float64Array.from(Pr) : null;
  S.fft2c(Gr, Gi, N, N, 1); S.fft2c(Pr, Pi, N, N, 1);
  // no grid correction: noise stays flat; nearest-cell gridding dims sources near the edge (<10% inside N/4)
  const norm = sumw > 0 ? 1 / (2 * sumw) : 0;
  for (let i = 0; i < N * N; i++) { Gr[i] *= norm; Pr[i] *= norm; }
  const beam = S.fitBeam(Pr, N, o.cell);
  return { dirty: Gr, psf: Pr, beam, sumw, sigma: sumw > 0 ? Math.sqrt(sumw2s) / sumw : 0, dropped, grid, wvis: wi, N, cell: o.cell, du };
};

// Gaussian fit to the PSF main lobe (like tclean's psfcutoff=0.35): weighted least squares on
// ln(psf) = -(A x^2 + B xy + C y^2) over the connected region above the cutoff
S.fitBeam = function (psf, N, cell, cutoff) {
  const h = N / 2, c0 = h * N + h, seen = new Uint8Array(N * N), stack = [c0], cut = cutoff || 0.35;
  let s40 = 0, s31 = 0, s22 = 0, s13 = 0, s04 = 0, r20 = 0, r11 = 0, r02 = 0, npix = 0;
  seen[c0] = 1;
  while (stack.length) {
    const i = stack.pop(), x = i % N, y = (i / N) | 0, dx = x - h, dy = y - h, p = psf[i];
    if (i !== c0) {
      const w = p * p, L = -Math.log(p), xx = dx * dx, yy = dy * dy, xy = dx * dy;
      s40 += w * xx * xx; s31 += w * xx * xy; s22 += w * xx * yy; s13 += w * xy * yy; s04 += w * yy * yy;
      r20 += w * xx * L; r11 += w * xy * L; r02 += w * yy * L;
    }
    npix++;
    for (const j of [i - 1, i + 1, i - N, i + N])
      if (j >= 0 && j < N * N && !seen[j] && psf[j] >= cut && Math.abs((j % N) - x) <= 1) { seen[j] = 1; stack.push(j); }
  }
  // normal equations for (A, B, C): [[s40,s31,s22],[s31,s22,s13],[s22,s13,s04]] . (A,B,C) = (r20,r11,r02)
  const M = [[s40, s31, s22], [s31, s22, s13], [s22, s13, s04]], r = [r20, r11, r02];
  const det3 = m => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  const D = det3(M);
  let sol = [0, 1, 2].map(k => det3(M.map((row, i) => row.map((v, j) => (j === k ? r[i] : v)))) / D);
  if (!sol.every(isFinite) && cut > 0.06) return S.fitBeam(psf, N, cell, cut / 2);   // undersampled PSF: include more of the main lobe
  if (!sol.every(isFinite)) {                                    // still only on-axis pixels (very coarse cell): no cross term, solve for A and C alone
    const d2 = s40 * s04 - s22 * s22;
    sol = [(r20 * s04 - s22 * r02) / d2, 0, (s40 * r02 - s22 * r20) / d2];
  }
  const A = sol[0], B = sol[1] / 2, Cc = sol[2];                 // quadratic form [[A,B],[B,C]] in pixels
  const tr = (A + Cc) / 2, dd = Math.sqrt(((A - Cc) / 2) * ((A - Cc) / 2) + B * B);
  const lmin = Math.max(tr - dd, 1e-9), lmax = Math.max(tr + dd, 1e-9);   // smallest curvature = major axis
  let ex = Math.abs(B) > 1e-12 ? B : (A <= Cc ? 1 : 0), ey = Math.abs(B) > 1e-12 ? lmin - A : (A <= Cc ? 0 : 1);
  let pa = Math.atan2(ex, ey) / DEG;
  if (pa > 90) pa -= 180; if (pa <= -90) pa += 180;
  const f = Math.sqrt(4 * Math.LN2), maj = f / Math.sqrt(lmin) * cell, min = f / Math.sqrt(lmax) * cell;
  // undersampled PSF (coarse cell): too few pixels above the cutoff for a fit; retry including more of the main lobe
  if (!(isFinite(maj) && isFinite(min) && isFinite(pa) && lmin > 1e-6) && cut > 0.06) return S.fitBeam(psf, N, cell, cut / 2);
  return { maj, min, pa, npix };
};
S.beamArea = (beam, cell) => Math.PI * beam.maj * beam.min / (4 * Math.LN2) / (cell * cell);   // pixels

/* ---------- elliptical Gaussian kernel (peak 1) and restoration ---------- */
S.gaussStamp = function (beam, cell) {
  const R = Math.ceil(1.6 * beam.maj / cell) + 1, W = 2 * R + 1, k = new Float64Array(W * W);
  const pa = beam.pa * DEG, sp = Math.sin(pa), cp = Math.cos(pa);
  for (let y = -R; y <= R; y++) for (let x = -R; x <= R; x++) {
    const dl = x * cell, dm = y * cell, a = dl * sp + dm * cp, b = dl * cp - dm * sp;
    k[(y + R) * W + x + R] = Math.exp(-4 * Math.LN2 * (a * a / (beam.maj * beam.maj) + b * b / (beam.min * beam.min)));
  }
  return { k, R, W };
};
S.smoothModel = function (model, N, cell, beam) {
  const out = new Float64Array(N * N), st = S.gaussStamp(beam, cell), R = st.R, W = st.W, k = st.k;
  for (let i = 0; i < N * N; i++) {
    const f = model[i]; if (f === 0) continue;
    const px = i % N, py = (i / N) | 0;
    for (let y = Math.max(0, py - R); y <= Math.min(N - 1, py + R); y++) {
      const ko = (y - py + R) * W - px + R, ro = y * N;
      for (let x = Math.max(0, px - R); x <= Math.min(N - 1, px + R); x++) out[ro + x] += f * k[ko + x];
    }
  }
  return out;
};
S.restore = function (model, residual, N, cell, beam) {
  const out = S.smoothModel(model, N, cell, beam);
  for (let i = 0; i < N * N; i++) out[i] += residual[i];
  return out;
};

/* ---------- Hogbom CLEAN (steppable) ---------- */
S.Hogbom = function (dirty, psf, N, o) {
  o = o || {};
  this.N = N; this.res = Float64Array.from(dirty); this.psf = psf; this.model = new Float64Array(N * N);
  this.gain = o.gain || 0.1; this.mask = o.mask || null; this.iter = 0; this.hist = []; this.flux = []; this.total = 0; this.stop = '';
};
S.Hogbom.prototype.peak = function () {
  const r = this.res, m = this.mask, n = r.length;
  let best = -1, bv = 0;
  for (let i = 0; i < n; i++) { if (m && !m[i]) continue; const a = r[i] < 0 ? -r[i] : r[i]; if (a > bv) { bv = a; best = i; } }
  return { i: best, v: bv };
};
S.Hogbom.prototype.run = function (niter, threshold) {
  const N = this.N, h = N / 2, res = this.res, psf = this.psf;
  for (let it = 0; it < niter; it++) {
    const p = this.peak();
    if (p.i < 0 || p.v <= (threshold || 0)) { this.stop = 'threshold'; return this; }
    const val = this.gain * res[p.i], px = p.i % N, py = (p.i / N) | 0;
    this.model[p.i] += val; this.total += val;
    const x0 = Math.max(0, px - h), x1 = Math.min(N, px + h), y0 = Math.max(0, py - h), y1 = Math.min(N, py + h);
    for (let y = y0; y < y1; y++) {
      const ro = y * N, po = (y - py + h) * N - px + h;
      for (let x = x0; x < x1; x++) res[ro + x] -= val * psf[po + x];
    }
    this.iter++; this.hist.push(p.v); this.flux.push(this.total);
  }
  this.stop = 'niter';
  return this;
};

/* ---------- FFT convolution helpers (kernels stored with their centre at N/2) ---------- */
function toOrigin(img, N) {     // swap quadrants so that pixel (N/2,N/2) moves to (0,0)
  const out = new Float64Array(N * N), h = N / 2;
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) out[((y + h) % N) * N + (x + h) % N] = img[y * N + x];
  return out;
}
function fwd(img, N) { const re = Float64Array.from(img), im = new Float64Array(N * N); S.fft2(re, im, N, N, -1); return { re, im }; }
function inv(F, N) {
  const re = Float64Array.from(F.re), im = Float64Array.from(F.im); S.fft2(re, im, N, N, 1);
  const s = 1 / (N * N); for (let i = 0; i < N * N; i++) re[i] *= s; return re;
}
function mul(A, B) {
  const n = A.re.length, re = new Float64Array(n), im = new Float64Array(n);
  for (let i = 0; i < n; i++) { re[i] = A.re[i] * B.re[i] - A.im[i] * B.im[i]; im[i] = A.re[i] * B.im[i] + A.im[i] * B.re[i]; }
  return { re, im };
}
// convolve img with a centred kernel (circular)
S.convolve = function (img, kernel, N) { return inv(mul(fwd(img, N), fwd(toOrigin(kernel, N), N)), N); };
S.gaussKernel = function (N, cell, beam, unitSum) {
  const k = new Float64Array(N * N), h = N / 2, pa = beam.pa * DEG, sp = Math.sin(pa), cp = Math.cos(pa);
  let s = 0;
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
    const dl = (x - h) * cell, dm = (y - h) * cell, a = dl * sp + dm * cp, b = dl * cp - dm * sp;
    s += (k[y * N + x] = Math.exp(-4 * Math.LN2 * (a * a / (beam.maj * beam.maj) + b * b / (beam.min * beam.min))));
  }
  if (unitSum) for (let i = 0; i < N * N; i++) k[i] /= s;
  return k;
};
// Fourier transform of a centred kernel, for repeated convolutions (cubes)
S.kernelFT = (kernel, N) => fwd(toOrigin(kernel, N), N);
S.convolveFT = (img, KF, N) => inv(mul(fwd(img, N), KF), N);

/* ---------- multi-scale CLEAN (Cornwell 2008 style, steppable) ---------- */
function scaleKernel(N, s) {
  const k = new Float64Array(N * N), h = N / 2;
  if (s === 0) { k[h * N + h] = 1; return k; }
  const R = 0.75 * s; let sum = 0;
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
    const r2 = ((x - h) * (x - h) + (y - h) * (y - h)) / (R * R);
    if (r2 < 1) sum += (k[y * N + x] = 1 - r2);
  }
  for (let i = 0; i < N * N; i++) k[i] /= sum;
  return k;
}
S.MSClean = function (dirty, psf, N, o) {
  this.N = N; this.scales = o.scales; this.gain = o.gain || 0.1; this.mask = o.mask || null;
  const ns = this.scales.length, bias = o.bias == null ? 0.6 : o.bias, smax = Math.max.apply(null, this.scales) || 1;
  this.K = this.scales.map(s => scaleKernel(N, s));
  const FK = this.K.map(k => fwd(toOrigin(k, N), N)), FP = fwd(toOrigin(psf, N), N), FD = fwd(dirty, N);
  this.R = FK.map(F => inv(mul(FD, F), N));
  this.P = []; this.norm = []; this.bias = [];
  for (let a = 0; a < ns; a++) {
    this.P[a] = [];
    for (let b = 0; b < ns; b++) this.P[a][b] = b < a ? this.P[b][a] : toOrigin(inv(mul(mul(FP, FK[a]), FK[b]), N), N);
  }
  for (let a = 0; a < ns; a++) { this.norm[a] = this.P[a][a][(N / 2) * N + N / 2]; this.bias[a] = 1 - bias * this.scales[a] / smax; }
  this.model = new Float64Array(N * N); this.iter = 0; this.hist = []; this.flux = []; this.total = 0; this.used = this.scales.map(() => 0);
};
S.MSClean.prototype.run = function (niter, threshold) {
  const N = this.N, h = N / 2, m = this.mask, ns = this.scales.length;
  for (let it = 0; it < niter; it++) {
    let bs = -1, bi = -1, bc = 0, peak0 = 0;
    for (let a = 0; a < ns; a++) {
      const R = this.R[a]; let vi = -1, vv = 0;
      for (let i = 0; i < N * N; i++) { if (m && !m[i]) continue; const x = R[i] < 0 ? -R[i] : R[i]; if (x > vv) { vv = x; vi = i; } }
      if (a === 0) peak0 = vv;
      const crit = this.bias[a] * vv / Math.sqrt(this.norm[a]);
      if (crit > bc) { bc = crit; bs = a; bi = vi; }
    }
    if (bs < 0 || peak0 <= (threshold || 0)) { this.stop = 'threshold'; return this; }
    const amp = this.gain * this.R[bs][bi] / this.norm[bs], px = bi % N, py = (bi / N) | 0;
    const K = this.K[bs];
    for (let y = 0; y < N; y++) {
      const ky = y - py + h; if (ky < 0 || ky >= N) continue;
      for (let x = 0; x < N; x++) { const kx = x - px + h; if (kx < 0 || kx >= N) continue; const kv = K[ky * N + kx]; if (kv) this.model[y * N + x] += amp * kv; }
    }
    for (let b = 0; b < ns; b++) {
      const R = this.R[b], P = this.P[bs][b];
      const x0 = Math.max(0, px - h), x1 = Math.min(N, px + h), y0 = Math.max(0, py - h), y1 = Math.min(N, py + h);
      for (let y = y0; y < y1; y++) { const ro = y * N, po = (y - py + h) * N - px + h; for (let x = x0; x < x1; x++) R[ro + x] -= amp * P[po + x]; }
    }
    this.total += amp; this.used[bs]++; this.iter++; this.hist.push(peak0); this.flux.push(this.total);
  }
  this.stop = 'niter';
  return this;
};
Object.defineProperty(S.MSClean.prototype, 'res', { get() { return this.R[0]; } });

/* ---------- predict model visibilities from an image (FFT + bilinear interpolation) ---------- */
S.predict = function (img, N, cell, u, v) {
  const h = N / 2, re = new Float64Array(N * N), im = new Float64Array(N * N);
  const gc = new Float64Array(N);
  for (let x = 0; x < N; x++) { const z = Math.PI * (x - h) / N; gc[x] = z === 0 ? 1 : Math.pow(Math.sin(z) / z, 2); }
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) re[y * N + x] = img[y * N + x] / (gc[x] * gc[y]);
  S.fft2c(re, im, N, N, -1);
  const du = 1 / (N * cell * AS), n = u.length, or = new Float64Array(n), oi = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const fx = u[k] / du + h, fy = v[k] / du + h, x0 = Math.floor(fx), y0 = Math.floor(fy), tx = fx - x0, ty = fy - y0;
    if (x0 < 0 || y0 < 0 || x0 + 1 >= N || y0 + 1 >= N) continue;
    const i = y0 * N + x0, w00 = (1 - tx) * (1 - ty), w10 = tx * (1 - ty), w01 = (1 - tx) * ty, w11 = tx * ty;
    or[k] = w00 * re[i] + w10 * re[i + 1] + w01 * re[i + N] + w11 * re[i + N + 1];
    oi[k] = w00 * im[i] + w10 * im[i + 1] + w01 * im[i + N] + w11 * im[i + N + 1];
  }
  return { re: or, im: oi };
};

/* ---------- statistics ---------- */
S.madRms = function (a, mask) {
  const v = [];
  for (let i = 0; i < a.length; i++) if (!mask || mask[i]) v.push(a[i]);
  if (!v.length) return 0;
  v.sort((x, y) => x - y);
  const med = v[v.length >> 1];
  for (let i = 0; i < v.length; i++) v[i] = Math.abs(v[i] - med);
  v.sort((x, y) => x - y);
  return 1.4826 * v[v.length >> 1];
};
S.maxAbs = function (a) { let m = 0, i0 = 0; for (let i = 0; i < a.length; i++) { const x = Math.abs(a[i]); if (x > m) { m = x; i0 = i; } } return { v: m, i: i0 }; };
S.max = function (a) { let m = -Infinity, i0 = 0; for (let i = 0; i < a.length; i++) if (a[i] > m) { m = a[i]; i0 = i; } return { v: m, i: i0 }; };
S.min = function (a) { let m = Infinity, i0 = 0; for (let i = 0; i < a.length; i++) if (a[i] < m) { m = a[i]; i0 = i; } return { v: m, i: i0 }; };
S.circleMask = function (N, cell, radius, x0, y0) {
  const m = new Uint8Array(N * N), h = N / 2;
  for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) m[y * N + x] = Math.hypot((x - h) * cell - (x0 || 0), (y - h) * cell - (y0 || 0)) <= radius ? 1 : 0;
  return m;
};
S.sumIn = function (img, mask) { let s = 0; for (let i = 0; i < img.length; i++) if (mask[i]) s += img[i]; return s; };

/* ---------- atmosphere: frozen Kolmogorov (von Karman) phase screen ----------
   returns excess path in micrometres on an nx*ny grid of spacing dx (m), normalised so that
   the rms path difference between two points 1 km apart equals rms1km */
S.phaseScreen = function (o) {
  const nx = o.nx || 2048, ny = o.ny || 256, dx = o.dx || 16, L0 = o.L0 || 3000, r = S.rng(o.seed || 7);
  const re = new Float64Array(nx * ny), im = new Float64Array(nx * ny);
  for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
    const kx = (x < nx / 2 ? x : x - nx) / (nx * dx), ky = (y < ny / 2 ? y : y - ny) / (ny * dx);
    const amp = Math.pow(kx * kx + ky * ky + 1 / (L0 * L0), -11 / 12);
    re[y * nx + x] = amp * r.normal(); im[y * nx + x] = amp * r.normal();
  }
  re[0] = im[0] = 0;
  S.fft2(re, im, nx, ny, 1);
  const lag = Math.round(1000 / dx);
  let d = 0, c = 0;
  for (let y = 0; y < ny; y += 4) for (let x = 0; x < nx; x += 2) { const a = re[y * nx + x] - re[y * nx + (x + lag) % nx]; d += a * a; c++; }
  const s = (o.rms1km || 100) / Math.sqrt(d / c);
  for (let i = 0; i < nx * ny; i++) re[i] *= s;
  return { path: re, nx, ny, dx };
};
S.screenAt = function (sc, xm, ym) {          // bilinear, periodic
  let fx = xm / sc.dx, fy = ym / sc.dx;
  fx = ((fx % sc.nx) + sc.nx) % sc.nx; fy = ((fy % sc.ny) + sc.ny) % sc.ny;
  const x0 = Math.floor(fx), y0 = Math.floor(fy), x1 = (x0 + 1) % sc.nx, y1 = (y0 + 1) % sc.ny, tx = fx - x0, ty = fy - y0, p = sc.path, n = sc.nx;
  return (1 - ty) * ((1 - tx) * p[y0 * n + x0] + tx * p[y0 * n + x1]) + ty * ((1 - tx) * p[y1 * n + x0] + tx * p[y1 * n + x1]);
};
// antenna phases (rad) on a set of times: phi[ti*nAnt + a]; offset (m) shifts the line of sight
S.antennaPhases = function (cfg, sc, times, o) {
  const lam = C / (o.freq * 1e9), wind = o.wind || 10, na = cfg.nAnt, out = new Float64Array(times.length * na);
  const ox = o.offsetX || 0, oy = o.offsetY || 0, cx = o.x0 || 0, cy = o.y0 || (sc.ny * sc.dx / 2);
  const scale = (o.factor == null ? 1 : o.factor) * 2 * Math.PI * 1e-6 / lam;
  for (let it = 0; it < times.length; it++) for (let a = 0; a < na; a++)
    out[it * na + a] = scale * S.screenAt(sc, cx + cfg.e[a] + ox - wind * times[it], cy + cfg.n[a] + oy);
  return out;
};

/* ---------- gain solver (StefCal-style alternating least squares) ----------
   intervals: 'int' (each integration), 'inf' (each scan), 'all' (whole data set) or seconds.
   Returns solutions per interval and the map vis -> interval. */
S.solveGains = function (d, mre, mim, o) {
  const na = d.nAnt, n = d.n, mode = o.mode || 'p', refant = o.refant || 0, minsnr = o.minsnr == null ? 3 : o.minsnr;
  const iv = new Int32Array(n), keys = new Map(), ivTime = [];
  const t0 = d.times[0];
  for (let k = 0; k < n; k++) {
    let key;
    if (o.solint === 'int') key = d.ti[k];
    else if (o.solint === 'inf') key = d.scan[k];
    else if (o.solint === 'all') key = 0;
    else key = d.scan[k] * 100000 + Math.floor((d.t[k] - t0) / o.solint);
    let id = keys.get(key);
    if (id === undefined) { id = keys.size; keys.set(key, id); ivTime.push([0, 0]); }
    iv[k] = id; ivTime[id][0] += d.t[k]; ivTime[id][1]++;
  }
  const ns = keys.size, gre = new Float64Array(ns * na).fill(1), gim = new Float64Array(ns * na), snr = new Float64Array(ns * na), flag = new Uint8Array(ns * na);
  const nre = new Float64Array(ns * na), nim = new Float64Array(ns * na), den = new Float64Array(ns * na);
  for (let pass = 0; pass < (o.iter || 40); pass++) {
    nre.fill(0); nim.fill(0); den.fill(0);
    for (let k = 0; k < n; k++) {
      if (d.flag[k]) continue;
      const s = iv[k] * na, i = s + d.a1[k], j = s + d.a2[k], w = d.wt[k];
      const vr = d.re[k], vi = d.im[k], mr = mre[k], mi = mim[k], m2 = mr * mr + mi * mi;
      // antenna a1: V g_j conj(M)
      let ar = gre[j] * mr + gim[j] * mi, ai = gim[j] * mr - gre[j] * mi;        // g_j * conj(M)
      nre[i] += w * (vr * ar - vi * ai); nim[i] += w * (vr * ai + vi * ar); den[i] += w * (gre[j] * gre[j] + gim[j] * gim[j]) * m2;
      // antenna a2: conj(V) g_i M
      ar = gre[i] * mr - gim[i] * mi; ai = gre[i] * mi + gim[i] * mr;            // g_i * M
      nre[j] += w * (vr * ar + vi * ai); nim[j] += w * (vr * ai - vi * ar); den[j] += w * (gre[i] * gre[i] + gim[i] * gim[i]) * m2;
    }
    for (let q = 0; q < ns * na; q++) {
      if (den[q] <= 0) continue;
      let r = nre[q] / den[q], im = nim[q] / den[q];
      if (mode === 'p') { const a = Math.hypot(r, im) || 1; r /= a; im /= a; }
      if (pass % 2 === 1) { r = 0.5 * (r + gre[q]); im = 0.5 * (im + gim[q]); if (mode === 'p') { const a = Math.hypot(r, im) || 1; r /= a; im /= a; } }
      gre[q] = r; gim[q] = im;
    }
  }
  for (let s = 0; s < ns; s++) {                              // reference phases to refant, SNR, flags
    const rr = gre[s * na + refant], ri = gim[s * na + refant], ra = Math.hypot(rr, ri) || 1;
    const cr = rr / ra, ci = -ri / ra;
    for (let a = 0; a < na; a++) {
      const q = s * na + a, r = gre[q] * cr - gim[q] * ci, im = gre[q] * ci + gim[q] * cr;
      gre[q] = r; gim[q] = im;
      snr[q] = Math.hypot(r, im) * Math.sqrt(den[q]);
      if (snr[q] < minsnr || den[q] <= 0) flag[q] = 1;
    }
  }
  if (mode === 'ap' && o.solnorm) {
    let s = 0, c = 0;
    for (let q = 0; q < ns * na; q++) if (!flag[q]) { s += Math.hypot(gre[q], gim[q]); c++; }
    const f = c ? c / s : 1;
    for (let q = 0; q < ns * na; q++) { gre[q] *= f; gim[q] *= f; }
  }
  return { ns, na, gre, gim, snr, flag, iv, time: ivTime.map(x => x[0] / x[1]) };
};
// apply solutions: V / (g1 conj g2); 'calflag' flags data whose solutions failed, 'calonly' leaves them
S.applyGains = function (d, sol, mode) {
  const out = S.copyData(d), na = sol.na;
  for (let k = 0; k < d.n; k++) {
    const s = sol.iv[k] * na, i = s + d.a1[k], j = s + d.a2[k];
    if (sol.flag[i] || sol.flag[j]) { if (mode !== 'calonly') out.flag[k] = 1; continue; }
    let gr = sol.gre[i] * sol.gre[j] + sol.gim[i] * sol.gim[j], gi = sol.gim[i] * sol.gre[j] - sol.gre[i] * sol.gim[j];
    const m = gr * gr + gi * gi; gr /= m; gi = -gi / m;
    const vr = d.re[k], vi = d.im[k];
    out.re[k] = gr * vr - gi * vi; out.im[k] = gr * vi + gi * vr;
    out.wt[k] = d.wt[k] * m;                                  // calwt: weights scale with |g|^2
  }
  return out;
};

root.Sim = S;
})(typeof window !== 'undefined' ? window : globalThis);
