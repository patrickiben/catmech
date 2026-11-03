/* The teaching model behind the catmech audit simulator: simulated people, an error in their recorded
   weight, an audit, and five ways to estimate how blood pressure changes with weight. It is a classroom
   model, loosely shaped like the README's NHANES example. It is not the package's validation studies. */
(function (root) {
  "use strict";
  const MI = typeof module === "object" && module.exports ? require("./linear-mi.js") : root.LinearMI;

  const TRUTH = 0.13;        // 0.130 mmHg of systolic blood pressure per kg of true weight, the slope to estimate
  const PEOPLE = 600;        // records per study, as in the package's cohorts of 600
  const SITES = 30;          // 20 records per site

  const ERRORS = {
    none: "No errors",
    noise: "Random errors",
    slips: "Decimal slips",
    rare: "Rare slips only",
    site: "Bad scales by site"
  };
  const SELECTION = {
    random: "Random records",
    sites: "Whole sites at random",
    recorded: "Heavier recorded weights more often"
  };
  const METHODS = [
    { key: "trust", label: "Trust the recorded value", short: "Trust recorded" },
    { key: "audited", label: "Audited records only", short: "Audited only" },
    { key: "default", label: "Fill in without the recorded value", short: "Fill in, no recorded" },
    { key: "linear", label: "Fill in from the recorded value", short: "Fill in from recorded" },
    { key: "oracle", label: "Every record checked", short: "Every record checked" }
  ];

  // One simulated study. size: error spread in kg (noise, site); share: records with an extra digit. With
  // slips the other records carry small errors of 2 kg; with rare slips the other records are exactly right.
  function makeStudy(opts, random) {
    const n = PEOPLE, site = [], age = [], truth = [], recorded = [], sbp = [];
    const offset = Array.from({ length: SITES }, () => random.norm());
    for (let i = 0; i < n; i++) {
      const s = Math.floor(i / (n / SITES));
      const a = 20 + 60 * random.unif();
      const w = 58 + 0.3 * a + 15 * random.norm();
      let r = w;
      if (opts.error === "noise") r = w - 1 + opts.size * random.norm();
      else if (opts.error === "slips") r = random.unif() < opts.share ? 10 * w : w - 1 + 2 * random.norm();
      else if (opts.error === "rare") r = random.unif() < opts.share ? 10 * w : w;
      else if (opts.error === "site") r = w + opts.size * offset[s] + 0.5 * random.norm();
      site.push(s); age.push(a); truth.push(w); recorded.push(r);
      sbp.push(96 + TRUTH * w + 0.45 * a + 14 * random.norm());
    }
    return { site, age, truth, recorded, sbp, n };
  }

  // Choose the audited records. fraction of records; whole sites when selection is "sites".
  function chooseAudit(study, opts, random) {
    const n = study.n, audited = new Array(n).fill(false);
    if (opts.selection === "sites") {
      const k = Math.max(2, Math.round(opts.fraction * SITES));
      const order = shuffle(Array.from({ length: SITES }, (_, s) => s), random).slice(0, k);
      const chosen = new Set(order);
      for (let i = 0; i < n; i++) audited[i] = chosen.has(study.site[i]);
      return audited;
    }
    const k = Math.max(8, Math.round(opts.fraction * n));
    let keys;
    if (opts.selection !== "recorded") keys = study.truth.map(() => random.unif());
    else {
      // Weighted sampling without replacement on the recorded weight's rank, so slips do not dominate:
      // the heaviest recorded weight is nine times as likely to be audited as the lightest.
      const rank = new Array(n);
      study.recorded.map((v, i) => [v, i]).sort((x, y) => x[0] - y[0]).forEach(([, i], r) => { rank[i] = r / (n - 1); });
      keys = rank.map((q) => Math.pow(random.unif(), 1 / Math.pow(9, q)));
    }
    const order = keys.map((key, i) => [key, i]).sort((x, y) => y[0] - x[0]).slice(0, k);
    order.forEach(([, i]) => { audited[i] = true; });
    return audited;
  }

  function shuffle(a, random) {
    for (let i = a.length - 1; i > 0; i--) {
      const j = Math.floor(random.unif() * (i + 1));
      [a[i], a[j]] = [a[j], a[i]];
    }
    return a;
  }

  // The five estimates of the weight slope, each with a 95% interval, or an error message.
  function estimate(study, audited, m, random) {
    const n = study.n, out = {};
    const all = Array.from({ length: n }, (_, i) => i);
    const aud = all.filter((i) => audited[i]);
    const ols = (rows, value) => {
      const a = MI.analyze(rows.map((i) => [1, value[i], study.age[i]]), rows.map((i) => study.sbp[i]), 1);
      return MI.interval(a.estimate, a.variance, a.df);
    };
    const attempt = (key, fn) => {
      try { out[key] = fn(); } catch (e) { out[key] = { error: e.message }; }
    };
    const checked = study.truth.map((v, i) => (audited[i] ? v : NaN));
    const row = (i, v) => [1, v, study.age[i]];
    attempt("trust", () => ols(all, study.recorded));
    attempt("audited", () => ols(aud, study.truth));
    attempt("default", () => MI.poolCopies(
      MI.fillIn(checked, audited, (i) => [study.age[i], study.sbp[i]], m, random), study.sbp, row, 1));
    attempt("linear", () => {
      // As in the package: if every audited record agrees (catmech_contamination_agreement with tolerance 0),
      // unaudited records keep their recorded value with no imputation uncertainty, and the result says so.
      // catmech_pool_lm still pools the m identical copies, so the df is Barnard-Rubin's with b = 0.
      if (aud.every((i) => Math.abs(study.recorded[i] - study.truth[i]) <= 1e-8 * Math.max(1, Math.abs(study.truth[i])))) {
        const kept = study.truth.map((v, i) => (audited[i] ? v : study.recorded[i]));
        const a = MI.analyze(all.map((i) => [1, kept[i], study.age[i]]), all.map((i) => study.sbp[i]), 1);
        return Object.assign(MI.pool(new Array(m).fill(a.estimate), new Array(m).fill(a.variance), a.df), { kept: true });
      }
      return MI.poolCopies(
        MI.fillIn(checked, audited, (i) => [study.age[i], study.sbp[i], study.recorded[i]], m, random), study.sbp, row, 1);
    });
    attempt("oracle", () => ols(all, study.truth));
    return out;
  }

  // Repeat whole studies and summarize each method against the true slope. A runner works in steps,
  // so a page can stay responsive; repeat() runs it to the end.
  function runner(opts, m, seed) {
    const random = MI.rng(seed);
    const tally = Object.fromEntries(METHODS.map((x) => [x.key, { cover: 0, sq: 0, width: 0, ok: 0, bias: 0, kept: 0 }]));
    let done = 0;
    function step(count) {
      for (let r = 0; r < count; r++) {
        const study = makeStudy(opts, random);
        const audited = chooseAudit(study, opts, random);
        const res = estimate(study, audited, m, random);
        METHODS.forEach(({ key }) => {
          const e = res[key];
          if (!e || e.error) return;
          const t = tally[key];
          t.ok += 1;
          t.cover += e.lower <= TRUTH && TRUTH <= e.upper ? 1 : 0;
          t.sq += (e.estimate - TRUTH) * (e.estimate - TRUTH);
          t.bias += e.estimate - TRUTH;
          t.width += e.upper - e.lower;
          t.kept += e.kept ? 1 : 0;
        });
        done += 1;
      }
      return done;
    }
    function summary() {
      return METHODS.map(({ key, label, short }) => {
        const t = tally[key];
        return { key, label, short, reps: done, ok: t.ok, kept: t.kept,
          coverage: t.ok ? 100 * t.cover / t.ok : NaN,
          rmse: t.ok ? Math.sqrt(t.sq / t.ok) : NaN,
          bias: t.ok ? t.bias / t.ok : NaN,
          width: t.ok ? t.width / t.ok : NaN };
      });
    }
    return { step, summary, done: () => done };
  }

  function repeat(opts, reps, m, seed) {
    const run = runner(opts, m, seed);
    run.step(reps);
    return run.summary();
  }

  const api = { TRUTH, PEOPLE, SITES, ERRORS, SELECTION, METHODS, makeStudy, chooseAudit, estimate, runner, repeat };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.AuditSim = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
