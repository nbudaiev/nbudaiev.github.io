/* kit.js: shared scenes and helpers for the decks (observations, test sources, clean + restore). */
(function () {
'use strict';
const S = window.Sim, K = {};
const cache = {};
K.memo = (key, f) => (key in cache ? cache[key] : (cache[key] = f()));

// a protoplanetary-disk-like test source: core + three rings, inclined 40 deg, plus a faint background point
K.disk = function (scale, opts) {
  const f = scale || 1, o = opts || {}, ratio = Math.cos((o.inc == null ? 35 : o.inc) * Math.PI / 180), pa = o.pa == null ? 120 : o.pa;
  return [
    { type: 'gauss', flux: 0.5e-3 * f, fwhm: 0.15, ratio, pa },
    { type: 'ring', flux: 2.0e-3 * f, r: 0.45, width: 0.10, ratio, pa },
    { type: 'ring', flux: 4.5e-3 * f, r: 1.10, width: 0.13, ratio, pa },
    { type: 'ring', flux: 2.7e-3 * f, r: 1.75, width: 0.20, ratio, pa },
    { type: 'point', flux: 0.4e-3 * f, x: 2.5, y: -1.7 }
  ];
};
K.fluxOf = comps => comps.reduce((s, c) => s + c.flux, 0);

// observation: configuration, frequency, track; noise from the radiometer equation (band, bandwidth, dt)
K.obs = function (o) {
  const key = JSON.stringify(o);
  return K.memo('obs' + key, () => {
    const cfg = o.aca ? S.mergeConfigs(S.config(o.cfg), S.config('ACA')) : S.config(o.cfg);
    const tr = S.uvTracks(cfg, { freq: o.freq || 230, dec: o.dec == null ? -30 : o.dec, ha0: o.ha0 == null ? -0.5 : o.ha0, ha1: o.ha1 == null ? 0.5 : o.ha1, dt: o.dt || 60 });
    const band = o.band || 6, bw = o.bw || 7.5e9, dt = o.dt || 60;
    let sigma = S.sigmaVis(band, bw, dt) * (o.noiseScale || 1);
    if (o.aca) {                                     // per-baseline noise for a mixed 12 m + 7 m array
      const s = new Float64Array(tr.n), D = cfg.dish;
      for (let k = 0; k < tr.n; k++) s[k] = S.sigmaVis(band, bw, dt, D[tr.a1[k]], D[tr.a2[k]]) * (o.noiseScale || 1);
      sigma = s;
    }
    return { cfg, tr, sigma };
  });
};

// clean an image result from Sim.image; returns restored image and bookkeeping
K.clean = function (img, o) {
  o = o || {};
  const N = img.N, cl = o.deconvolver === 'multiscale'
    ? new S.MSClean(img.dirty, img.psf, N, { scales: o.scales || [0, 4, 12], gain: o.gain || 0.1, mask: o.mask })
    : new S.Hogbom(img.dirty, img.psf, N, { gain: o.gain || 0.1, mask: o.mask });
  cl.run(o.niter == null ? 1000 : o.niter, o.threshold || 0);
  const restored = S.restore(cl.model, cl.res, N, img.cell, img.beam);
  return { cl, restored, model: cl.model, residual: cl.res, iter: cl.iter, stop: cl.stop };
};

// highest |PSF| outside the main lobe (main lobe = positive region around the centre, within 2 beams)
K.sidelobe = function (img) {
  const N = img.N, cell = img.cell, b = img.beam, pa = b.pa * Math.PI / 180, sp = Math.sin(pa), cp = Math.cos(pa), h = N / 2, psf = img.psf;
  const inside = (x, y) => { const dl = (x - h) * cell, dm = (y - h) * cell, a = (dl * sp + dm * cp) / b.maj, c = (dl * cp - dm * sp) / b.min; return a * a + c * c < 4; };
  const lobe = new Uint8Array(N * N), st = [h * N + h]; lobe[h * N + h] = 1;
  while (st.length) {
    const i = st.pop(), x = i % N, y = (i / N) | 0;
    for (const j of [i - 1, i + 1, i - N, i + N]) {
      const xx = j % N, yy = (j / N) | 0;
      if (j < 0 || j >= N * N || lobe[j] || Math.abs(xx - x) > 1 || psf[j] <= 0 || !inside(xx, yy)) continue;
      lobe[j] = 1; st.push(j);
    }
  }
  let m = 0;
  for (let i = 0; i < N * N; i++) if (!lobe[i]) { const a = Math.abs(psf[i]); if (a > m) m = a; }
  return m;
};

// connected regions of a boolean image (4-neighbour); returns label image and sizes
K.regions = function (on, N) {
  const lab = new Int32Array(N * N).fill(-1), sizes = [];
  for (let s = 0; s < N * N; s++) {
    if (!on[s] || lab[s] >= 0) continue;
    const id = sizes.length, st = [s]; lab[s] = id; let n = 0;
    while (st.length) {
      const i = st.pop(), x = i % N; n++;
      for (const j of [i - 1, i + 1, i - N, i + N]) if (j >= 0 && j < N * N && on[j] && lab[j] < 0 && Math.abs((j % N) - x) <= 1) { lab[j] = id; st.push(j); }
    }
    sizes.push(n);
  }
  return { lab, sizes };
};

// a teaching version of tclean's auto-multithresh: each major cycle, threshold the residual at
// max(noisethreshold*rms, sidelobethreshold*sidelobe*peak), prune small islands, smooth by the beam,
// grow the existing mask into emission above lownoisethreshold*rms, then clean inside the mask
K.autoClean = function (img, o) {
  o = Object.assign({ sidelobethreshold: 3, noisethreshold: 5, lownoisethreshold: 1.5, minbeamfrac: 0.3, gain: 0.1, niter: 3000, threshold: 0 }, o);
  const S = window.Sim, N = img.N, cell = img.cell, sl = K.sidelobe(img), beamPix = S.beamArea(img.beam, cell);
  const cl = new S.Hogbom(img.dirty, img.psf, N, { gain: o.gain }), mask = new Uint8Array(N * N);
  cl.mask = mask;
  const kern = S.gaussKernel(N, cell, img.beam), KF = S.kernelFT(kern, N);
  const cycles = [];
  for (let cyc = 0; cyc < 40 && cl.iter < o.niter; cyc++) {
    const rms = S.madRms(cl.res), pk = S.maxAbs(cl.res).v;
    const T = Math.max(o.noisethreshold * rms, o.sidelobethreshold * sl * pk);
    const on = new Uint8Array(N * N);
    for (let i = 0; i < N * N; i++) on[i] = cl.res[i] > T ? 1 : 0;
    const reg = K.regions(on, N);
    const keep = new Float64Array(N * N);
    for (let i = 0; i < N * N; i++) if (on[i] && reg.sizes[reg.lab[i]] >= o.minbeamfrac * beamPix) keep[i] = 1;
    const sm = S.convolveFT(keep, KF, N);
    let smax = 0; for (let i = 0; i < N * N; i++) if (sm[i] > smax) smax = sm[i];
    let added = 0;
    if (smax > 0) for (let i = 0; i < N * N; i++) if (!mask[i] && sm[i] > 0.01 * smax) { mask[i] = 1; added++; }
    for (let g = 0; g < 8; g++) {                    // grow into low-level emission next to the mask
      const grow = [];
      for (let i = 0; i < N * N; i++) {
        if (mask[i] || cl.res[i] < o.lownoisethreshold * rms) continue;
        const x = i % N;
        if ((x > 0 && mask[i - 1]) || (x < N - 1 && mask[i + 1]) || (i >= N && mask[i - N]) || (i < N * N - N && mask[i + N])) grow.push(i);
      }
      if (!grow.length) break;
      grow.forEach(i => { mask[i] = 1; }); added += grow.length;
    }
    let pin = 0; for (let i = 0; i < N * N; i++) if (mask[i]) { const a = Math.abs(cl.res[i]); if (a > pin) pin = a; }
    const stopAt = Math.max(o.threshold, sl * pin);
    cycles.push({ T, rms, added });
    if (pin <= o.threshold && !added) break;
    if (!pin) break;
    cl.run(Math.min(o.niter - cl.iter, 400), stopAt);
  }
  const restored = S.restore(cl.model, cl.res, N, cell, img.beam);
  return { cl, mask, restored, model: cl.model, residual: cl.res, iter: cl.iter, cycles, sidelobe: sl };
};

// cut through an image along position angle pa (deg E of N), +-half arcsec
K.cut = function (img, N, cell, pa, half, n) {
  n = n || 161;
  const x = new Float64Array(n), y = new Float64Array(n), s = Math.sin(pa * Math.PI / 180), c = Math.cos(pa * Math.PI / 180), h = N / 2;
  for (let i = 0; i < n; i++) {
    const r = -half + 2 * half * i / (n - 1), fx = h + r * s / cell, fy = h + r * c / cell, x0 = Math.floor(fx), y0 = Math.floor(fy), tx = fx - x0, ty = fy - y0;
    x[i] = r;
    if (x0 < 0 || y0 < 0 || x0 + 1 >= N || y0 + 1 >= N) { y[i] = NaN; continue; }
    const p = y0 * N + x0;
    y[i] = (1 - tx) * (1 - ty) * img[p] + tx * (1 - ty) * img[p + 1] + (1 - tx) * ty * img[p + N] + tx * ty * img[p + N + 1];
  }
  return { x, y };
};

// gridded imaging weights for display: sum of per-visibility weights on an M x M grid spanning +-ext (lambda)
K.uvDensity = function (d, w, ext, M) {
  const g = new Float64Array(M * M), h = M / 2, du = 2 * ext / M;
  for (let k = 0; k < d.n; k++) {
    if (!w[k]) continue;
    const iu = Math.round(d.u[k] / du), iv = Math.round(d.v[k] / du);
    if (Math.abs(iu) >= h || Math.abs(iv) >= h) continue;
    g[(iv + h) * M + iu + h] += w[k]; g[(h - iv) * M + h - iu] += w[k];
  }
  return g;
};
// log stretch for display (0..1), ignoring empty cells
K.logStretch = function (g, decades) {
  let mx = 0; for (let i = 0; i < g.length; i++) if (g[i] > mx) mx = g[i];
  const out = new Float64Array(g.length), lo = Math.log10(mx) - (decades || 2.5);
  for (let i = 0; i < g.length; i++) out[i] = g[i] > 0 ? Math.max(0.04, (Math.log10(g[i]) - lo) / (decades || 2.5)) : 0;
  return out;
};

// atmosphere + phase referencing: one shared scene for the calibration and self-calibration decks.
// o: {cfg, band, rms1km (um of path at 1 km), wind (m/s), seed, dur (min), dt (s), cycle (min), calScan (s),
//     sep (deg), h (m, screen height), wvr (bool), calFlux (Jy), dec}
K.atmScene = function (o) {
  o = Object.assign({ cfg: 'C-6', band: 6, rms1km: 120, wind: 10, seed: 7, dur: 60, dt: 10, cycle: 8, calScan: 30, sep: 3, h: 1000, wvr: false, calFlux: 0.5, dec: -30, posErr: 0.25 }, o);
  const S = window.Sim, ghz = S.BANDS[o.band].ghz, cfg = S.config(o.cfg), na = cfg.nAnt, T = o.dur * 60;
  const sc = K.memo('screen:' + o.rms1km + ':' + o.seed, () => S.phaseScreen({ nx: 2048, ny: 256, dx: 16, rms1km: o.rms1km, seed: o.seed }));
  const times = [], isCal = [], scans = [];
  let t = 0, id = 0;
  const gap = Math.min(10, 0.08 * o.cycle * 60);                  // slew + settle
  while (t < T) {
    scans.push({ id: id++, kind: 'cal', t0: t, t1: t + o.calScan }); for (let u = t; u < t + o.calScan; u += o.dt) { times.push(u); isCal.push(1); }
    t += o.calScan + gap;
    const tEnd = Math.min(T, t + o.cycle * 60 - o.calScan - 2 * gap);
    if (tEnd - t >= o.dt) { scans.push({ id: id++, kind: 'tgt', t0: t, t1: tEnd }); for (let u = t; u < tEnd; u += o.dt) { times.push(u); isCal.push(0); } }
    t = tEnd + gap;
  }
  scans.push({ id: id++, kind: 'cal', t0: t, t1: t + o.calScan }); for (let u = t; u < t + o.calScan; u += o.dt) { times.push(u); isCal.push(1); }
  const f = o.wvr ? 0.45 : 1, off = o.h * Math.tan(o.sep * S.DEG);
  const pT = S.antennaPhases(cfg, sc, times, { freq: ghz, wind: o.wind, factor: f });
  const pC = S.antennaPhases(cfg, sc, times, { freq: ghz, wind: o.wind, factor: f, offsetX: 0.6 * off, offsetY: 0.8 * off });
  const r = S.rng(o.seed + 100), inst = new Float64Array(times.length * na);
  // direction-dependent term seen only by the target: antenna position errors (mm) x calibrator offset (rad)
  const lam = 299.792458 / ghz, sx = 0.6 * o.sep * S.DEG, sy = 0.8 * o.sep * S.DEG;
  for (let a = 0; a < na; a++) {
    const dd = 2 * Math.PI / lam * o.posErr * (r.normal() * sx + r.normal() * sy);
    for (let k = 0; k < times.length; k++) pT[k * na + a] += dd;
  }
  for (let a = 0; a < na; a++) {
    const p0 = 2 * Math.PI * r(), A = 0.6 * r(), P = 1800 + 2400 * r(), q = 2 * Math.PI * r(), sl = (r() - 0.5) * 0.8 / 3600;
    for (let k = 0; k < times.length; k++) inst[k * na + a] = p0 + A * Math.sin(2 * Math.PI * times[k] / P + q) + sl * times[k];
  }
  // calibrator solution per cal scan: scan-averaged phase + thermal noise of an antenna-based solution
  const sigVis = S.sigmaVis(o.band, 7.5e9, o.dt), nint = Math.max(1, Math.round(o.calScan / o.dt));
  const snrCal = o.calFlux / sigVis * Math.sqrt((na - 1) * nint), sigSol = 1 / snrCal;
  const sol = [];
  scans.filter(s => s.kind === 'cal').forEach(s => {
    const idx = []; times.forEach((u, k) => { if (u >= s.t0 && u < s.t1) idx.push(k); });
    const ph = new Float64Array(na);
    for (let a = 0; a < na; a++) {
      let cr = 0, ci = 0; idx.forEach(k => { const p = pC[k * na + a] + inst[k * na + a]; cr += Math.cos(p); ci += Math.sin(p); });
      ph[a] = Math.atan2(ci, cr) + sigSol * r.normal();
    }
    for (let a = na - 1; a >= 0; a--) ph[a] -= ph[0];              // reference antenna 0
    sol.push({ t: (s.t0 + s.t1) / 2, ph });
  });
  const wrap = x => Math.atan2(Math.sin(x), Math.cos(x));
  const interp = (u, a) => {
    let i = 0; while (i < sol.length - 2 && sol[i + 1].t < u) i++;
    const s0 = sol[i], s1 = sol[i + 1], w = Math.max(0, Math.min(1, (u - s0.t) / (s1.t - s0.t)));
    return s0.ph[a] + w * wrap(s1.ph[a] - s0.ph[a]);
  };
  const truth = new Float64Array(times.length * na), resid = new Float64Array(times.length * na);
  for (let k = 0; k < times.length; k++) for (let a = 0; a < na; a++) {
    truth[k * na + a] = pT[k * na + a] + inst[k * na + a] - (pT[k * na] + inst[k * na]);
    resid[k * na + a] = wrap(truth[k * na + a] - interp(times[k], a));
  }
  return { o, cfg, ghz, times, isCal, scans, sol, truth, resid, interp, pT, pC, inst, na, snrCal, sigSol, screen: sc };
};
// target-only visibilities of a scene with its residual phases applied (gains on the target samples)
K.targetData = function (A, comps, opt) {
  opt = opt || {};
  const S = window.Sim, idx = [], tt = [];
  A.times.forEach((u, k) => { if (!A.isCal[k]) { idx.push(k); tt.push(u - A.o.dur * 30); } });
  const scanOf = idx.map(k => A.scans.findIndex(s => A.times[k] >= s.t0 && A.times[k] < s.t1));
  const tr = S.uvTracks(A.cfg, { freq: A.ghz, dec: A.o.dec, times: tt, scans: scanOf, elmin: 0 });
  const na = A.na, g = { re: new Float64Array(tr.nt * na), im: new Float64Array(tr.nt * na) };
  for (let j = 0; j < tr.nt; j++) for (let a = 0; a < na; a++) {
    const p = opt.perfect ? 0 : A.resid[idx[j] * na + a], am = opt.amp ? opt.amp(a, j) : 1;
    g.re[j * na + a] = am * Math.cos(p); g.im[j * na + a] = am * Math.sin(p);
  }
  const sigma = S.sigmaVis(A.o.band, 7.5e9, A.o.dt) * (opt.noiseScale || 1);
  return S.observe(tr, comps, { sigma, seed: opt.seed || 5, gains: g, keepModel: true });
};

// formatting
K.as = v => (v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v.toFixed(3).replace(/0$/, '')) + '″';
K.beamStr = b => `${b.maj.toFixed(2)}″ × ${b.min.toFixed(2)}″`;
K.uJy = v => (v * 1e6 >= 100 ? (v * 1e6).toFixed(0) : (v * 1e6).toFixed(1)) + ' μJy/bm';
K.mJy = v => (v * 1e3).toFixed(v * 1e3 >= 10 ? 1 : 2) + ' mJy';
K.pct = v => (v * 100).toFixed(0) + '%';
// brightness temperature of S (Jy/beam) for beam (arcsec) at nu (GHz), Rayleigh-Jeans
K.tb = (S, beam, ghz) => 1.222e6 * S / (ghz * ghz * beam.maj * beam.min);

window.Kit = K;
})();
