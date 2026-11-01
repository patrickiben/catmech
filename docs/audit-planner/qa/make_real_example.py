#!/usr/bin/env python3
"""Builds qa/example-pilot-nhanes.csv: 300 evaluation rows of the frozen NHANES pool, with the measured weight kept on a random 60 only."""
import csv, pathlib, random
ROOT = pathlib.Path(__file__).resolve().parents[3]
pool = [r for r in csv.DictReader(open(ROOT / "validation/real-data/nhanes-pool.csv")) if r["split"] == "evaluation"]
rng = random.Random(20261101)
rows = rng.sample(pool, 300)
pilot = set(rng.sample(range(300), 60))
out = ROOT / "docs/audit-planner/qa/example-pilot-nhanes.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["record_id", "age", "sex", "race", "height", "sbp", "recorded_kg", "reference_kg"])
    for i, r in enumerate(rows):
        w.writerow([f"R{i+1:03d}", r["age"], r["sex"], r["race"], r["height"], r["sbp"], r["recorded_kg"], r["reference_kg"] if i in pilot else "."])
print(out, "recorded missing:", sum(r["recorded_kg"] == "NA" for r in rows))
