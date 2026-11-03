/* Teaching engine for the catmech audit simulator.
   It follows catmech_reference(use_recorded = TRUE) and catmech_pool_lm(): a least-squares fit on the
   audited rows through a QR decomposition, sigma drawn as sqrt(SSE / chi-square(df)), coefficients drawn as
   beta + sigma * R^-1 z, normal noise on each filled row, and Rubin's rules with Barnard-Rubin degrees of
   freedom (mice::pool.scalar). docs/learn/qa/ checks it against the R package. Works in a browser and in node. */
(function (root) {
  "use strict";

  // Random numbers: a seeded generator, so a run can be repeated.
  function rng(seed) {
    let a = seed >>> 0;
    let spare = null;
    const unif = () => {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return (((t ^ (t >>> 14)) >>> 0) + 0.5) / 4294967296;
    };
    const norm = () => {
      if (spare !== null) { const s = spare; spare = null; return s; }
      const u = unif(), v = unif();
      const r = Math.sqrt(-2 * Math.log(u));
      spare = r * Math.sin(2 * Math.PI * v);
      return r * Math.cos(2 * Math.PI * v);
    };
    const gamma = (shape) => {
      if (shape < 1) return gamma(shape + 1) * Math.pow(unif(), 1 / shape);
      const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
      for (;;) {
        let x, v;
        do { x = norm(); v = 1 + c * x; } while (v <= 0);
        v = v * v * v;
        const u = unif();
        if (u < 1 - 0.0331 * x * x * x * x) return d * v;
        if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
      }
    };
    return { unif, norm, chisq: (df) => 2 * gamma(df / 2) };
  }

  // Least squares through a Householder QR, as R's qr() does. X is an array of rows.
  function qrFit(X, y) {
    const n = X.length, p = X[0].length;
    const A = X.map((r) => r.slice());
    const b = y.slice();
    const R = Array.from({ length: p }, () => new Array(p).fill(0));
    for (let k = 0; k < p; k++) {
      let norm = 0;
      for (let i = k; i < n; i++) norm += A[i][k] * A[i][k];
      norm = Math.sqrt(norm);
      if (!(norm > 1e-10)) throw new Error("The audited rows cannot support this model; audit more records.");
      const alpha = A[k][k] > 0 ? -norm : norm;
      const v = new Array(n).fill(0);
      v[k] = A[k][k] - alpha;
      for (let i = k + 1; i < n; i++) v[i] = A[i][k];
      let vv = 0;
      for (let i = k; i < n; i++) vv += v[i] * v[i];
      for (let j = k; j < p; j++) {
        let s = 0;
        for (let i = k; i < n; i++) s += v[i] * A[i][j];
        const f = 2 * s / vv;
        for (let i = k; i < n; i++) A[i][j] -= f * v[i];
      }
      let s = 0;
      for (let i = k; i < n; i++) s += v[i] * b[i];
      const f = 2 * s / vv;
      for (let i = k; i < n; i++) b[i] -= f * v[i];
      for (let j = k; j < p; j++) R[k][j] = A[k][j];
    }
    const beta = backsolve(R, b.slice(0, p));
    let sse = 0;
    for (let i = 0; i < n; i++) {
      let fit = 0;
      for (let j = 0; j < p; j++) fit += X[i][j] * beta[j];
      sse += (y[i] - fit) * (y[i] - fit);
    }
    const df = n - p;
    if (df < 1) throw new Error("More audited records than model terms are needed.");
    if (!(sse > 1e-24 * n)) throw new Error("The audited rows leave no residual variance.");
    return { beta, sse, df, R, p, n };
  }

  function backsolve(R, z) {
    const p = z.length, x = new Array(p).fill(0);
    for (let i = p - 1; i >= 0; i--) {
      let s = z[i];
      for (let j = i + 1; j < p; j++) s -= R[i][j] * x[j];
      x[i] = s / R[i][i];
    }
    return x;
  }

  // Variance of coefficient j: sigma^2 times the j-th diagonal of (R'R)^-1.
  // (R'R)^-1 = R^-1 R^-T, so its j-th diagonal is the sum of squares of row j of R^-1.
  // Column k of R^-1 solves R x = e_k.
  function coefVariance(fit, j) {
    const p = fit.p;
    let total = 0;
    for (let k = 0; k < p; k++) {
      const e = new Array(p).fill(0);
      e[k] = 1;
      const col = backsolve(fit.R, e);
      total += col[j] * col[j];
    }
    return (fit.sse / fit.df) * total;
  }

  // Log gamma (Lanczos), regularized incomplete beta, Student t distribution and quantile.
  const LG = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
    x -= 1;
    let a = LG[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LG[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }
  function betacf(a, b, x) {
    const tiny = 1e-300;
    let c = 1, d = 1 - (a + b) * x / (a + 1);
    if (Math.abs(d) < tiny) d = tiny;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 10000; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((a + m2 - 1) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny;
      c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1));
      d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny;
      c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-16) break;
    }
    return h;
  }
  function ibeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  function pt(t, df) {
    const tail = 0.5 * ibeta(df / (df + t * t), df / 2, 0.5);
    return t >= 0 ? 1 - tail : tail;
  }
  function qt(p, df) {
    if (!(df > 0)) return NaN;
    if (p === 0.5) return 0;
    if (p < 0.5) return -qt(1 - p, df);
    let lo = 0, hi = 2;
    while (pt(hi, df) < p) { lo = hi; hi *= 2; if (hi > 1e300) return Infinity; }
    for (let i = 0; i < 300 && hi - lo > 1e-15 * Math.max(1, hi); i++) {
      const mid = (lo + hi) / 2;
      if (pt(mid, df) < p) lo = mid; else hi = mid;
    }
    return (lo + hi) / 2;
  }

  // One analysis: ordinary least squares, the term's estimate, variance and residual df.
  function analyze(X, y, term) {
    const fit = qrFit(X, y);
    return { estimate: fit.beta[term], variance: coefVariance(fit, term), df: fit.df, n: fit.n };
  }

  function interval(estimate, variance, df, level) {
    const se = Math.sqrt(variance), crit = qt((1 + (level || 0.95)) / 2, df);
    return { estimate, se, df, lower: estimate - crit * se, upper: estimate + crit * se };
  }

  // Rubin's rules with Barnard-Rubin df, as mice::pool.scalar(Q, U, n, k) with dfcom = n - k, written as
  // mice 3.19.0's barnard.rubin, which has no floor on lambda and so also handles identical copies (b = 0).
  function pool(Q, U, dfcom, level) {
    const m = Q.length;
    const qbar = Q.reduce((s, q) => s + q, 0) / m;
    const ubar = U.reduce((s, u) => s + u, 0) / m;
    const b = Q.reduce((s, q) => s + (q - qbar) * (q - qbar), 0) / (m - 1);
    const t = ubar + (1 + 1 / m) * b;
    const lambda = (1 + 1 / m) * b / t;
    const tmp = (1 - lambda) * (1 + dfcom) * dfcom;
    const df = (m - 1) * tmp / ((dfcom + 3) * (m - 1) + lambda * lambda * tmp);
    return Object.assign(interval(qbar, t, df, level), { within: ubar, between: b });
  }

  // Fill in the checked value on unaudited rows. design(i) gives row i's imputation predictors,
  // without the intercept. Returns m arrays of length n, with audited rows keeping their checked value.
  function fillIn(checked, audited, design, m, random) {
    const n = checked.length;
    const fitRows = [], fillRows = [];
    for (let i = 0; i < n; i++) (audited[i] ? fitRows : fillRows).push(i);
    const X = fitRows.map((i) => [1].concat(design(i)));
    const fit = qrFit(X, fitRows.map((i) => checked[i]));
    const Xfill = fillRows.map((i) => [1].concat(design(i)));
    const copies = [];
    for (let j = 0; j < m; j++) {
      const sigma = Math.sqrt(fit.sse / random.chisq(fit.df));
      const z = Array.from({ length: fit.p }, () => random.norm());
      const noise = backsolve(fit.R, z);
      const beta = fit.beta.map((bk, k) => bk + sigma * noise[k]);
      const copy = checked.slice();
      Xfill.forEach((x, r) => {
        let mean = 0;
        for (let k = 0; k < fit.p; k++) mean += x[k] * beta[k];
        copy[fillRows[r]] = mean + sigma * random.norm();
      });
      copies.push(copy);
    }
    return copies;
  }

  // Analyze each completed copy and pool. analysisRow(i, value) gives the analysis design row.
  function poolCopies(copies, outcome, analysisRow, term, level) {
    const Q = [], U = [];
    let dfcom = null;
    copies.forEach((copy) => {
      const X = copy.map((v, i) => analysisRow(i, v));
      const a = analyze(X, outcome, term);
      Q.push(a.estimate); U.push(a.variance); dfcom = a.df;
    });
    return pool(Q, U, dfcom, level);
  }

  const api = { rng, qrFit, backsolve, coefVariance, lgamma, ibeta, pt, qt, analyze, interval, pool, fillIn, poolCopies };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.LinearMI = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
