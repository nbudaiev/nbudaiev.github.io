/* deck.js: slide navigation. Arrow keys / PageUp / PageDown / space, Home, End;
   F = fullscreen, T = light/dark. Demos (window.DEMOS[name]) start the first time their slide is shown. */
(function () {
'use strict';
const slides = Array.from(document.querySelectorAll('.slide'));
const mq = matchMedia('(min-width: 700px) and (min-height: 380px) and (min-aspect-ratio: 1/1)');
const deckTitle = (document.querySelector('meta[name="deck-title"]') || {}).content || document.title;
const started = new WeakSet();
let cur = 0;

slides.forEach((s, i) => {
  if (!s.id) s.id = 's' + (i + 1);
  if (s.classList.contains('title-slide')) return;
  const f = document.createElement('div'); f.className = 'foot';
  const a = document.createElement('span'); a.textContent = `${deckTitle} · ${i + 1} / ${slides.length}`;
  f.append(a); s.appendChild(f);
});
const bar = document.createElement('div'); bar.className = 'progress'; document.body.appendChild(bar);
const nav = document.createElement('nav'); nav.className = 'nav'; nav.setAttribute('aria-label', 'Slides');
const mkBtn = (t, label, fn) => { const b = document.createElement('button'); b.type = 'button'; b.textContent = t; b.setAttribute('aria-label', label); b.addEventListener('click', fn); nav.appendChild(b); return b; };
mkBtn('‹', 'Previous slide', () => go(cur - 1));
const count = document.createElement('button'); count.type = 'button'; count.setAttribute('aria-label', 'Back to the workshop index');
count.addEventListener('click', () => { location.href = 'index.html'; }); nav.appendChild(count);
mkBtn('›', 'Next slide', () => go(cur + 1));
mkBtn('◐', 'Switch light or dark', toggleTheme);
document.body.appendChild(nav);
let navTimer = 0;
const wake = () => { nav.classList.add('show'); clearTimeout(navTimer); navTimer = setTimeout(() => nav.classList.remove('show'), 2200); };
document.addEventListener('pointermove', e => { if (e.clientY > innerHeight - 90 && e.clientX > innerWidth - 260) wake(); });
document.addEventListener('touchstart', wake, { passive: true });

function start(s) {
  const k = s.dataset.demo;
  if (!k || started.has(s) || !window.DEMOS || !window.DEMOS[k]) return;
  started.add(s);
  try { window.DEMOS[k](s); } catch (e) { console.error('demo ' + k, e); }
}
function go(i) {
  i = Math.max(0, Math.min(slides.length - 1, i));
  if (mq.matches) {
    if (document.activeElement && slides[cur].contains(document.activeElement)) document.activeElement.blur();
    slides.forEach((s, j) => s.classList.toggle('active', j === i));
  } else slides[i].scrollIntoView({ block: 'start' });
  cur = i;
  start(slides[i]);
  bar.style.width = (100 * (i + 1) / slides.length) + '%';
  count.textContent = `${i + 1} / ${slides.length}`;
  try { history.replaceState(null, '', '#' + slides[i].id); } catch (e) { /* sandboxed frames may refuse */ }
}
function toggleTheme() {
  const r = document.documentElement, dark = getComputedStyle(r).colorScheme === 'dark';
  r.setAttribute('data-theme', dark ? 'light' : 'dark');
  try { localStorage.setItem('deck-theme', r.getAttribute('data-theme')); } catch (e) { /* storage optional */ }
  requestAnimationFrame(() => window.Viz && window.Viz.rerenderAll());
}
try { const t = localStorage.getItem('deck-theme'); if (t) document.documentElement.setAttribute('data-theme', t); } catch (e) { /* storage optional */ }
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => window.Viz && window.Viz.rerenderAll());

document.addEventListener('keydown', e => {
  if (e.metaKey || e.ctrlKey || e.altKey) return;
  const t = e.target, tag = t && t.tagName, typing = tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA';
  const k = e.key;
  if (k === 'PageDown' || (!typing && (k === 'ArrowRight' || k === 'ArrowDown' || (k === ' ' && tag !== 'BUTTON')))) { e.preventDefault(); go(cur + 1); }
  else if (k === 'PageUp' || (!typing && (k === 'ArrowLeft' || k === 'ArrowUp'))) { e.preventDefault(); go(cur - 1); }
  else if (!typing && k === 'Home') { e.preventDefault(); go(0); }
  else if (!typing && k === 'End') { e.preventDefault(); go(slides.length - 1); }
  else if (!typing && (k === 'f' || k === 'F')) {
    const d = document;
    try { if (d.fullscreenElement) d.exitFullscreen(); else d.documentElement.requestFullscreen().catch(() => {}); } catch (err) { /* optional */ }
  } else if (!typing && (k === 't' || k === 'T')) toggleTheme();
});

// in scrolling (phone) mode demos start as they come into view
const io = new IntersectionObserver(es => es.forEach(en => { if (en.isIntersecting && !mq.matches) start(en.target); }), { rootMargin: '300px 0px' });
slides.forEach(s => io.observe(s));
mq.addEventListener('change', () => go(cur));

const fromHash = () => {
  const h = decodeURIComponent((location.hash || '').slice(1));
  let i = slides.findIndex(s => s.id === h);
  if (i < 0 && /^\d+$/.test(h)) i = +h - 1;
  return i;
};
window.addEventListener('hashchange', () => { const i = fromHash(); if (i >= 0 && i !== cur) go(i); });
const first = fromHash();
go(first < 0 ? 0 : first);
})();
