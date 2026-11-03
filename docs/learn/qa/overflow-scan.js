/* Paste into the browser console on docs/learn/index.html or simulator.html, at each width to check.
   It opens every folded section, then reports text that spills out of the bordered or shaded box it sits in,
   text past the page edge, SVG labels past their chart, and SVG labels that overlap. Text inside a box that
   scrolls sideways is skipped, since it can be scrolled into view. It measures the rendered text, not the
   element's box, because text that overflows a fixed column keeps its element's box inside the card.
   Checked on 2026-11-03: it caught the guide's example column at its old 150 px width (6 px past the card). */
(() => {
  document.querySelectorAll("details").forEach((d) => { d.open = true; });
  const W = document.documentElement.clientWidth, out = [];
  const boxed = (e) => { const s = getComputedStyle(e); return parseFloat(s.borderLeftWidth) > 0 || parseFloat(s.borderRightWidth) > 0 || (s.backgroundColor !== "rgba(0, 0, 0, 0)" && s.backgroundColor !== "transparent"); };
  const scrolls = (e) => /(auto|scroll|hidden|clip)/.test(getComputedStyle(e).overflowX) && e.tagName !== "BODY";
  const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  let n;
  while ((n = walker.nextNode())) {
    if (!n.textContent.trim()) continue;
    const el = n.parentElement;
    if (!el || el.closest("svg, script, style")) continue;
    const range = document.createRange(); range.selectNodeContents(n);
    const r = range.getBoundingClientRect();
    if (!r.width) continue;
    let p = el, scroller = false;
    while (p && p !== document.body) {
      if (p !== el && scrolls(p)) { scroller = true; break; }
      if (p !== el && boxed(p)) break;
      if (p === el && boxed(p) && getComputedStyle(p).display !== "inline") break;
      p = p.parentElement;
    }
    if (scroller || /(auto|scroll)/.test(getComputedStyle(el).overflowX)) continue;
    const label = n.textContent.trim().slice(0, 40);
    if (r.right > W + 1) out.push(`past the page edge: "${label}"`);
    if (p && p !== document.body) {
      const b = p.getBoundingClientRect(), cs = getComputedStyle(p);
      const innerR = b.right - parseFloat(cs.borderRightWidth), innerL = b.left + parseFloat(cs.borderLeftWidth);
      if (r.right > innerR + 1 || r.left < innerL - 1) out.push(`spills out of ${p.tagName.toLowerCase()}.${(p.className || "").toString().split(" ")[0]} by ${Math.round(Math.max(r.right - innerR, innerL - r.left))} px: "${label}"`);
    }
  }
  for (const svg of document.querySelectorAll("svg")) {
    const sb = svg.getBoundingClientRect();
    if (!sb.width) continue;
    const texts = [...svg.querySelectorAll("text")].map((t) => ({ t, r: t.getBoundingClientRect() })).filter((x) => x.r.width);
    for (const { t, r } of texts) if (r.right > sb.right + 1 || r.left < sb.left - 1) out.push(`svg text past its chart: "${t.textContent.slice(0, 40)}"`);
    for (let i = 0; i < texts.length; i++) for (let j = i + 1; j < texts.length; j++) {
      const a = texts[i].r, b = texts[j].r;
      if (a.left < b.right - 1 && b.left < a.right - 1 && a.top < b.bottom - 1 && b.top < a.bottom - 1) out.push(`svg labels overlap: "${texts[i].t.textContent.slice(0, 24)}" / "${texts[j].t.textContent.slice(0, 24)}"`);
    }
  }
  return JSON.stringify({ width: W, docWidth: document.documentElement.scrollWidth, problems: [...new Set(out)] });
})()
