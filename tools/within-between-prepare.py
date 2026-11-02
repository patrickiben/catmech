#!/usr/bin/env python3
"""Builds validation/within-between/snap-pool.csv from the public SNAP QC file, as validation/WITHIN_BETWEEN_PLAN.md requires. Writes counts only, with no analysis.
Usage: python3 -I tools/within-between-prepare.py RAW_CSV"""
import collections
import csv
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(sys.argv[1])
assert hashlib.md5(RAW.read_bytes()).hexdigest() == "d35bf400376f1774109c9f3b713458eb", "Not the published file"
OUT = ROOT/"validation/within-between"
OUT.mkdir(parents=True, exist_ok=True)
kept = []
dropped = collections.Counter()
with RAW.open(newline="", encoding="latin-1") as stream:
    for r in csv.DictReader(stream):
        status, error, benefit = r["STATUS"], float(r["AMTERR"]), float(r["RAWBEN"])
        if r["CASE"] != "1":
            dropped["case_not_1"] += 1
        elif (status == "1" and error != 0) or (status in ("2", "3") and error == 0):
            dropped["status_and_error_disagree"] += 1
        else:
            truth = benefit-error if status == "2" else benefit+error if status == "3" else benefit
            kept.append((r["STATE"], benefit/100, truth/100, int(float(r["CERTHHSZ"])), float(r["RAWGROSS"])/100))
with (OUT/"snap-pool.csv").open("w", newline="", encoding="utf-8") as stream:
    w = csv.writer(stream)
    w.writerow(["state", "recorded", "reference", "size", "gross"])
    w.writerows(kept)
sizes = collections.Counter(k[0] for k in kept)
print("rows kept", len(kept), "dropped", dict(dropped), "states", len(sizes), "smallest", min(sizes.values()), "largest", max(sizes.values()))
