/* Runs the guide's method picker against a stub page and writes every R snippet it can produce, with the
   answers that produce it. Usage: node docs/learn/qa/picker-snippets.js docs/learn/index.html SNIPPETS.json
   tools/learn-picker-check.R then runs each snippet in R. */
const fs = require("fs");
const html = fs.readFileSync(process.argv[2], "utf8");
const data = html.match(/<script type="application\/json" id="picker-text">([\s\S]*?)<\/script>/)[1];
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map((m) => m[1]);
const main = scripts.find((s) => s.includes("picker-text"));
const state = {};
const els = {};
const make = () => ({ textContent: "", className: "", dataset: {}, children: [], appendChild(c) { this.children.push(c); }, addEventListener(type, fn) { this.handler = fn; } });
const document = {
  getElementById: (id) => { if (id === "picker-text") return { textContent: data }; return (els[id] = els[id] || make()); },
  querySelector: (sel) => ({ value: state[sel.match(/name=(\w+)/)[1]] }),
  createElement: () => make(),
};
const names = { role: ["predictor", "outcome"], chosen: ["random", "targeted", "truth"], sites: ["no", "yes"], errors: ["small", "gross", "unsure"], checked: ["exact", "known", "unknown"] };
Object.keys(names).forEach((k) => { state[k] = names[k][0]; });
new Function("document", "navigator", "getSelection", "setTimeout", main)(document, {}, () => ({}), setTimeout);
const out = {};
const keys = Object.keys(names);
(function rec(i) {
  if (i === keys.length) { els["pick-form"].handler(); out[els["pick-code"].textContent] = { method: els["pick-result"].dataset.method, ...state }; return; }
  for (const v of names[keys[i]]) { state[keys[i]] = v; rec(i + 1); }
})(0);
const list = Object.entries(out).map(([code, s]) => ({ code, ...s }));
const digest = require("crypto").createHash("sha256").update(list.map((s) => s.code).join("\n---\n")).digest("hex");
if (process.argv[3]) fs.writeFileSync(process.argv[3], JSON.stringify(list, null, 1));
console.log(`${list.length} distinct snippets, sha256 ${digest}`);
