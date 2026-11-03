/* Writes the simulated studies that tools/learn-simulator-check.R runs through the R package.
   Usage: node docs/learn/qa/make-datasets.js OUTPUT.json
   The same seeds always give the same studies, so compare.js can rebuild them without this file. */
"use strict";
const fs = require("fs");
const path = require("path");
const S = require(path.join(__dirname, "..", "audit-sim.js"));
const MI = require(path.join(__dirname, "..", "linear-mi.js"));

const DESIGN = [
  { count: 40, opts: { error: "noise", size: 3, selection: "random", fraction: 0.25 } },
  { count: 40, opts: { error: "slips", share: 0.02, selection: "random", fraction: 0.25 } },
  { count: 20, opts: { error: "site", size: 4, selection: "sites", fraction: 0.25 } },
  { count: 20, opts: { error: "noise", size: 3, selection: "recorded", fraction: 0.25 } },
  { count: 10, opts: { error: "none", selection: "random", fraction: 0.25 } },
  { count: 10, opts: { error: "rare", share: 0.01, selection: "random", fraction: 0.05 } }
];
const SEED = 20261102;

function studies() {
  const out = [];
  DESIGN.forEach(({ count, opts }) => {
    for (let k = 0; k < count; k++) {
      const random = MI.rng(SEED + out.length);
      const study = S.makeStudy(opts, random);
      const audited = S.chooseAudit(study, opts, random);
      out.push({ id: out.length, opts, study, audited });
    }
  });
  return out;
}

if (require.main === module) {
  const target = process.argv[2];
  if (!target) { console.error("Usage: node make-datasets.js OUTPUT.json"); process.exit(2); }
  const data = studies().map(({ id, opts, study, audited }) => ({
    id, error: opts.error, selection: opts.selection,
    age: study.age, sbp: study.sbp, recorded: study.recorded, truth: study.truth, audited
  }));
  fs.writeFileSync(target, JSON.stringify(data));
  console.log(`Wrote ${data.length} studies to ${target}`);
}
module.exports = { studies, SEED, DESIGN };
