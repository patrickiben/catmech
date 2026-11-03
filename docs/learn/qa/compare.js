/* Checks the audit simulator's engine against the R package.
   Usage: node docs/learn/qa/compare.js [OUTPUT.txt]   (exits 1 on any failure)
   Exact checks: least-squares intervals on every study; one imputation draw with fixed deviates through the
   package's own fitting and drawing steps; the analysis and pooling of the package's own completed copies; and
   every study where the audit found no disagreement, so the recorded values were kept.
   Distribution checks: the browser's imputations against the package's on the same studies, which differ only
   by random draws, so they are compared on average. Mutated engines and pooling rules must be caught. */
"use strict";
const fs = require("fs");
const path = require("path");
const vm = require("vm");
const MI = require(path.join(__dirname, "..", "linear-mi.js"));
const S = require(path.join(__dirname, "..", "audit-sim.js"));
const { studies } = require(path.join(__dirname, "make-datasets.js"));

const lines = [];
const say = (s = "") => lines.push(s);
let failures = 0;
const check = (ok, text) => { say(`${ok ? "PASS" : "FAIL"}  ${text}`); if (!ok) failures += 1; };
const rel = (a, b) => Math.abs(a - b) / Math.max(1e-12, Math.abs(b));
const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const zOf = (a) => mean(a) / (Math.sqrt(a.reduce((s, v) => s + (v - mean(a)) ** 2, 0) / (a.length - 1)) / Math.sqrt(a.length));

const csv = fs.readFileSync(path.join(__dirname, "r-results.csv"), "utf8").trim().split("\n");
const head = csv[0].split(",").map((h) => h.replace(/"/g, ""));
const R = csv.slice(1).map((line) => {
  const v = line.split(",").map((x) => x.replace(/"/g, ""));
  return Object.fromEntries(head.map((h, i) => [h, ["error", "selection", "method"].includes(h) ? v[i] : Number(v[i])]));
});
const rOf = (id, method) => R.find((r) => r.id === id && r.method === method);
const all = studies();
say("catmech audit simulator: browser engine against the R package");
say(`Studies: ${all.length}, from docs/learn/qa/make-datasets.js; R results from tools/learn-simulator-check.R`);
say();

// 1. Exact: least squares for trusting the recorded value, audited rows only, and every record checked.
let worst = 0;
all.forEach(({ id, study, audited }) => {
  const rows = Array.from({ length: study.n }, (_, i) => i), aud = rows.filter((i) => audited[i]);
  const ols = (use, value) => {
    const a = MI.analyze(use.map((i) => [1, value[i], study.age[i]]), use.map((i) => study.sbp[i]), 1);
    return MI.interval(a.estimate, a.variance, a.df);
  };
  [["trust", ols(rows, study.recorded)], ["audited", ols(aud, study.truth)], ["oracle", ols(rows, study.truth)]].forEach(([m, js]) => {
    const r = rOf(id, m);
    ["estimate", "se", "lower", "upper"].forEach((k) => { worst = Math.max(worst, rel(js[k], r[k])); });
    worst = Math.max(worst, Math.abs(js.df - r.df));
  });
});
check(worst < 1e-9, `least squares on ${all.length} studies, three methods: largest relative difference ${worst.toExponential(2)}`);

// 2. Exact: one imputation draw with fixed deviates, against the package's catmech_model_design,
// catmech_model_fit and parameter draw. Householder signs may differ by row, so the coefficient deviates
// are flipped to match, which leaves the draw unchanged.
const draw = JSON.parse(fs.readFileSync(path.join(__dirname, "r-draw.json"), "utf8"));
function drawWith(engine) {
  const { study, audited } = all[draw.id];
  const design = (i) => [study.age[i], study.sbp[i], study.recorded[i]];
  const rows = audited.map((a, i) => (a ? i : -1)).filter((i) => i >= 0);
  const fit = engine.qrFit(rows.map((i) => [1].concat(design(i))), rows.map((i) => study.truth[i]));
  const flip = fit.R.map((row, k) => Math.sign(row[k]) * Math.sign(draw.R[k][k]));
  const deviates = draw.z.map((z, k) => z * flip[k]).concat(draw.e);
  let next = 0;
  // The stub's chi-square scales with the df it is asked for, so a draw on the wrong df gives a different value.
  const stub = { chisq: (df) => draw.chisq * df / draw.df, norm: () => deviates[next++], unif: () => 0.5 };
  const checked = study.truth.map((v, i) => (audited[i] ? v : NaN));
  const copy = engine.fillIn(checked, audited, design, 1, stub)[0];
  const filled = copy.filter((v, i) => !audited[i]);
  return { fit, filled };
}
const ours = drawWith(MI);
const fitDiff = Math.max(...draw.beta_hat.map((b, k) => rel(ours.fit.beta[k], b)), rel(ours.fit.sse, draw.sse), Math.abs(ours.fit.df - draw.df));
const drawDiff = Math.max(...draw.expected.map((v, k) => rel(ours.filled[k], v)));
check(draw.pivot.every((p, k) => p === k + 1), "the package's fit used no column pivoting, so coefficient order matches");
check(fitDiff < 1e-9, `the fit on audited rows (coefficients, residual sum of squares, df): largest relative difference ${fitDiff.toExponential(2)}`);
check(drawDiff < 1e-9, `one draw of ${draw.expected.length} filled values with fixed deviates: largest relative difference ${drawDiff.toExponential(2)}`);

// Negative controls: mutated engines must fail the exact draw.
const source = fs.readFileSync(path.join(__dirname, "..", "linear-mi.js"), "utf8");
const mutants = [
  ["sigma without the chi-square draw", "Math.sqrt(fit.sse / random.chisq(fit.df))", "Math.sqrt(fit.sse / fit.df)"],
  ["chi-square on the wrong df", "random.chisq(fit.df)", "random.chisq(fit.df - 1)"],
  ["no coefficient draw", "bk + sigma * noise[k]", "bk + 0 * noise[k]"],
  ["no noise on filled rows", "mean + sigma * random.norm()", "mean + 0 * random.norm()"]
];
mutants.forEach(([name, from, to]) => {
  if (!source.includes(from)) { check(false, `negative control "${name}": pattern not found in linear-mi.js`); return; }
  const sandbox = { module: { exports: {} }, globalThis: {} };
  vm.runInNewContext(source.replace(from, to), sandbox);
  let caught;
  try { const m = drawWith(sandbox.module.exports); caught = Math.max(...draw.expected.map((v, k) => rel(m.filled[k], v))) > 1e-6; }
  catch (e) { caught = true; }
  check(caught, `negative control: a mutated engine (${name}) fails the exact draw`);
});

// 3. Exact: analysis and pooling of the package's completed copies (m = 5), Barnard-Rubin df included.
const exact = JSON.parse(fs.readFileSync(path.join(__dirname, "r-exact.json"), "utf8"));
const poolCase = (c, poolFn) => {
  const study = all[c.id].study;
  return poolFn(c.copies, study.sbp, (i, v) => [1, v, study.age[i]]);
};
let worstPool = 0, cases = 0;
exact.forEach((pair) => ["default", "linear"].forEach((method) => {
  const c = pair[method];
  const js = poolCase(c, (copies, y, row) => MI.poolCopies(copies, y, row, 1));
  ["estimate", "se", "df", "lower", "upper"].forEach((k) => { worstPool = Math.max(worstPool, rel(js[k], c.pooled[k])); });
  cases += 1;
}));
check(worstPool < 1e-8, `pooling the package's own completed copies, ${cases} cases: largest relative difference ${worstPool.toExponential(2)}`);
const wrongPool = (copies, y, row) => {
  const Q = [], U = [];
  copies.forEach((copy) => { const a = MI.analyze(copy.map((v, i) => row(i, v)), y, 1); Q.push(a.estimate); U.push(a.variance); });
  const m = Q.length, qbar = mean(Q), ubar = mean(U), b = Q.reduce((s, q) => s + (q - qbar) ** 2, 0) / (m - 1);
  return MI.interval(qbar, ubar + b, 50);
};
check(exact.some((pair) => rel(poolCase(pair.linear, wrongPool).se, pair.linear.pooled.se) > 1e-6),
  "negative control: a pooling rule without the (1 + 1/m) factor is detected");

// 4. Exact: studies where every audited record agreed, so the recorded values were kept and m identical
// copies were pooled (b = 0).
let keptCount = 0, worstKept = 0;
all.forEach(({ id, study, audited }, k) => {
  const js = S.estimate(study, audited, 20, MI.rng(5000 + k)).linear;
  if (!js.kept) return;
  keptCount += 1;
  const r = rOf(id, "linear");
  ["estimate", "se", "df", "lower", "upper"].forEach((key) => { worstKept = Math.max(worstKept, rel(js[key], r[key])); });
});
check(keptCount >= 10, `studies where the audit found no disagreement: ${keptCount}`);
check(worstKept < 1e-8, `kept recorded values, against the package on those ${keptCount} studies: largest relative difference ${worstKept.toExponential(2)}`);
const floored = (() => {
  const c = all.find(({ study, audited }, k) => S.estimate(study, audited, 20, MI.rng(5000 + k)).linear.kept);
  const a = MI.analyze(c.study.truth.map((v, i) => [1, c.audited[i] ? v : c.study.recorded[i], c.study.age[i]]), c.study.sbp, 1);
  const lambda = 1e-4, m = 20, dfcom = a.df;
  const dfold = (m - 1) / (lambda * lambda), dfobs = (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda);
  return rel(dfold * dfobs / (dfold + dfobs), rOf(c.id, "linear").df) > 1e-6;
})();
check(floored, "negative control: the older pooling rule with a floor on lambda is detected on a kept study");

// 5. Distribution: the browser's imputation against the package's, 20 copies each.
say();
say(`Imputation, browser against package, same studies, m = 20 each. Differences come from random draws, so these use z-tests (|z| < 3.5).`);
["default", "linear"].forEach((method) => {
  const d = [], logRatio = [];
  let coverJ = 0, coverR = 0, n = 0;
  all.forEach(({ id, study, audited }, k) => {
    const js = S.estimate(study, audited, 20, MI.rng(5000 + k))[method], r = rOf(id, method);
    coverJ += js.lower <= S.TRUTH && S.TRUTH <= js.upper ? 1 : 0;
    coverR += r.lower <= S.TRUTH && S.TRUTH <= r.upper ? 1 : 0;
    if (js.kept) return;
    n += 1;
    d.push(js.estimate - r.estimate);
    logRatio.push(Math.log(js.se / r.se));
  });
  const ratio = Math.exp(mean(logRatio));
  say(`  ${method}, ${n} studies with imputation: mean difference in estimate ${mean(d).toExponential(2)} (z = ${zOf(d).toFixed(2)}); geometric mean SE ratio ${ratio.toFixed(4)} (z = ${zOf(logRatio).toFixed(2)}); intervals covering the true slope ${coverJ} and ${coverR} of ${all.length}`);
  check(Math.abs(zOf(d)) < 3.5, `${method}: no systematic difference in estimates`);
  check(Math.abs(zOf(logRatio)) < 3.5, `${method}: no systematic difference in standard errors`);
});
const noDraws = { chisq: (df) => df, norm: () => 0, unif: Math.random };
const control = all.filter(({ study, audited }, k) => !S.estimate(study, audited, 20, MI.rng(5000 + k)).linear.kept)
  .map(({ id, study, audited }) => Math.log(S.estimate(study, audited, 20, noDraws).linear.se / rOf(id, "linear").se));
check(Math.abs(zOf(control)) >= 3.5, `negative control: imputation without random draws is detected (z = ${zOf(control).toFixed(2)})`);

say();
say("Notes: a first version required average standard errors within 3%. It flagged the default method at 1.0353,");
say("but two browser runs with different seeds differed by as much, so the standard-error check is a z-test, which");
say("misses smaller systematic differences; the exact draw above is what pins down the imputation step. The first");
say("version also floored lambda at 1e-4, as older mice did; mice 3.19.0, which the package ran with, does not.");
say();
say(failures ? `RESULT: ${failures} FAILED` : "RESULT: PASS");
const text = lines.join("\n") + "\n";
if (process.argv[2]) fs.writeFileSync(process.argv[2], text);
process.stdout.write(text);
process.exit(failures ? 1 : 0);
