#!/usr/bin/env node
/* Runs the exact engine text embedded in ../index.html under Node and compares it with:
   (1) the three R test vectors in validation/audit-plan-study/planner-data.json, exact values;
   (2) the pool values, recomputed here from the raw frozen pools with the same covariates as tools/verify-artifacts.py;
   (3) R's bootstrap limits on the vectors, as a reasonableness check only, because the random streams differ.
   Writes qa/k-port-check.json. Exit 1 if (1) or (2) fails. */
const fs = require('fs'), path = require('path'), vm = require('vm');
const here = __dirname, root = path.resolve(here, '../../..');
const html = fs.readFileSync(path.join(here, '../index.html'), 'utf8');
const engine = html.match(/<script id="k-engine">([\s\S]*?)<\/script>/)[1];
const box = { setTimeout, console }; box.globalThis = box; vm.createContext(box); vm.runInContext(engine, box);
const K = box.AuditK;
const planner = JSON.parse(fs.readFileSync(path.join(root, 'validation/audit-plan-study/planner-data.json'), 'utf8'));
function readCsv(file) {
  const lines = fs.readFileSync(file, 'utf8').trim().split(/\r?\n/), split = l => { const out = []; let cur = '', q = false; for (const c of l) { if (c === '"') q = !q; else if (c === ',' && !q) { out.push(cur); cur = ''; } else cur += c; } out.push(cur); return out; };
  const head = split(lines[0]); return lines.slice(1).map(l => { const v = split(l), r = {}; head.forEach((h, i) => r[h] = v[i]); return r; });
}
function poolInputs(name) {
  const pool = readCsv(path.join(root, `validation/real-data/${name}-pool.csv`)).filter(r => r.split === 'evaluation');
  if (name === 'ctgov') return { covariates: pool.map(r => [+(r.phase === 'early'), +(r.phase === 'phase2'), +(r.phase === 'phase3_4'), +(r.sponsor === 'industry'), Math.log(+r.sites), +r.start_year - 2010]), outcome: pool.map(r => Math.log(+r.duration_months)), reference: pool.map(r => Math.log(+r.reference_n)), recorded: pool.map(r => Math.log(+r.recorded_n)) };
  const levels = ['black', 'mexican_american', 'other_hispanic', 'other'];
  return { covariates: pool.map(r => [+r.age, +(r.sex === 'male')].concat(levels.map(v => +(r.race === v)), [+r.height])), outcome: pool.map(r => +r.sbp), reference: pool.map(r => +r.reference_kg), recorded: pool.map(r => r.recorded_kg === 'NA' ? null : +r.recorded_kg) };
}
(async () => {
  const report = { generated_from: 'docs/audit-planner/index.html, script#k-engine', node: process.version, vectors: [], pools: [], bootstrap: [], normal_functions: {}, records_rule: [] };
  let ok = true;
  for (const v of planner.js_test_vectors) { const r = K.checkVector(v); report.vectors.push({ name: r.name, pass: r.pass, worst_relative_difference: r.worstRelativeDifference, k: r.k, discard_k: r.discardK, classical: r.classical, fails: r.fails }); ok = ok && r.pass; }
  for (const name of ['ctgov', 'nhanes']) {
    const b = planner.datasets[name], ctx = K.prepare(poolInputs(name)), pt = K.point(ctx);
    const diffs = { k: pt.k / b.pool_k - 1, discard_k: pt.discardK / b.pool_discard_k - 1, classical: pt.classical / b.pool_classical - 1, discard_classical: pt.discardClassical / b.pool_discard_classical - 1 };
    const pass = Math.abs(diffs.k) < 1e-7 && Math.abs(diffs.discard_k) < 1e-7 && Math.abs(diffs.classical) < 1e-6 && Math.abs(diffs.discard_classical) < 1e-6;
    ok = ok && pass; report.pools.push({ name, rows: ctx.n, js_k: pt.k, json_k: b.pool_k, js_discard_k: pt.discardK, json_discard_k: b.pool_discard_k, js_classical: pt.classical, json_classical: b.pool_classical, relative_differences: diffs, pass });
    report.records_rule.push({ name, k: b.pool_k, discard_k: b.pool_discard_k, fraction_keep: K.required(b.pool_k, 1.25), fraction_discard: K.required(b.pool_discard_k, 1.25), exact_records_keep: K.required(b.pool_k, 1.25) * 600, exact_records_discard: K.required(b.pool_discard_k, 1.25) * 600, records_keep_ceiling: K.recordsFor(b.pool_k, 1.25, 600), records_discard_ceiling: K.recordsFor(b.pool_discard_k, 1.25, 600), ratio_at_rounded_down: K.relative(b.pool_k, Math.floor(K.required(b.pool_k, 1.25) * 600) / 600) });
  }
  for (const v of planner.js_test_vectors) {
    const ctx = K.prepare(v), pt = K.point(ctx), t0 = Date.now(), iv = await K.interval(ctx, pt, { bootstrap: 500, seed: 1 }), rb = v.expected.r_bootstrap_only;
    const pred = K.predictive(ctx, iv, [0.05, 0.1, 0.25, 0.5, 0.75], ctx.n);
    report.bootstrap.push({ name: v.name, ms: Date.now() - t0, usable: iv.usable, js: [iv.kLower, iv.kUpper, iv.discardLower, iv.discardUpper], r: [rb.k_lower, rb.k_upper, rb.discard_k_lower, rb.discard_k_upper], predictive: pred });
  }
  report.normal_functions = { qnorm_975: K.qnorm(0.975), expected_qnorm_975: 1.959963984540054, pnorm_196: K.pnorm(1.959963984540054), qnorm_001: K.qnorm(0.001), expected_qnorm_001: -3.090232306167813, round_trip_error: Math.max(...[-4, -2, -1, 0.3, 1.5, 3].map(z => Math.abs(K.qnorm(K.pnorm(z)) - z))) };
  report.all_exact_checks_pass = ok;
  fs.writeFileSync(path.join(here, 'k-port-check.json'), JSON.stringify(report, null, 1));
  for (const r of report.vectors) console.log(r.name, r.pass ? 'PASS' : 'FAIL', 'worst rel diff', r.worst_relative_difference.toExponential(2), r.fails.join('; '));
  for (const r of report.pools) console.log(r.name, 'pool', r.pass ? 'PASS' : 'FAIL', JSON.stringify(r.relative_differences));
  for (const r of report.records_rule) console.log(r.name, 'exact keep', r.exact_records_keep.toFixed(3), 'discard', r.exact_records_discard.toFixed(3), 'ceil', r.records_keep_ceiling, r.records_discard_ceiling, 'ratio at floor', r.ratio_at_rounded_down.toFixed(4));
  for (const r of report.bootstrap) console.log(r.name, r.ms + 'ms', 'JS', r.js.map(x => x.toFixed(4)).join(' '), '| R', r.r.map(x => x.toFixed(4)).join(' '));
  console.log('normal', JSON.stringify(report.normal_functions));
  process.exit(ok ? 0 : 1);
})();
