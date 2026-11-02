#!/usr/bin/env python3
"""Applies the comparison's scorer with the post hoc rerun in place of the original rows. Plain Python, from saved results.
Usage: python3 -I tools/outside-methods-sensitivity-score.py [OUTPUT_FILE]"""
import csv
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT/"validation/outside-methods-study/replicates.csv.gz"
RERUN = ROOT/"validation/outside-methods-sensitivity/replicates.csv.gz"
NOTE = "POST HOC. The comparison's rule was KEEP as written; these numbers do not change it. See validation/OUTSIDE_METHODS_SENSITIVITY_PLAN.md. mice with eps = 0; raking with method = approx.\n\n"


def key(r):
    return (r["analysis"], r["fraction"], r["replicate"], r["method"])


with gzip.open(RERUN, "rt", newline="", encoding="utf-8") as stream:
    rerun = {key(r): r for r in csv.DictReader(stream)}
replaced = 0
with tempfile.TemporaryDirectory() as folder:
    with gzip.open(STUDY, "rt", newline="", encoding="utf-8") as source, gzip.open(Path(folder)/"replicates.csv.gz", "wt", newline="", encoding="utf-8") as target:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(target, reader.fieldnames)
        writer.writeheader()
        for r in reader:
            new = rerun.get(key(r))
            if new is not None:
                replaced += 1
                if r["fitted_estimate"]:
                    new = dict(new, fitted_estimate=r["fitted_estimate"])
            writer.writerow(new if new is not None else r)
    if replaced != len(rerun):
        raise SystemExit(f"Rerun rows {len(rerun)} do not all match the comparison: {replaced}")
    # Seconds for the rerun methods reflect different work, so the scorer's timing lines are read from the comparison, not from here.
    done = subprocess.run([sys.executable, "-I", str(ROOT/"tools/outside-methods-score.py"), folder], capture_output=True, text=True, check=True)
text = NOTE+f"Rows replaced by the rerun: {replaced}\n\n"+done.stdout
if len(sys.argv) > 1:
    Path(sys.argv[1]).write_text(text, encoding="utf-8")
print(text)
