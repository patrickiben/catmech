#!/usr/bin/env python3
"""Check saved evidence against code and replicate-level data; no validity claim."""
import argparse
import csv
import hashlib
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def verify_summary(summary, records):
    for row in summary:
        group = [r for r in records if r["scenario"] == row["scenario"]
                 and r["method"] == row["method"]]
        good = [r for r in group if r["status"] == "ok"]
        assert len(group) == int(row["requested"])
        assert len(good) == int(row["successful"])
        assert len(group)-len(good) == int(row["failures"])
        error = [float(r["estimate"])-float(row["estimand_truth"]) for r in good]
        cover = [float(r["lower"]) <= float(row["estimand_truth"]) <= float(r["upper"])
                 for r in good]
        calculated = {
            "bias": statistics.mean(error),
            "bias_mcse": statistics.stdev(error)/math.sqrt(len(error)),
            "rmse": math.sqrt(statistics.mean(e*e for e in error)),
            "coverage": statistics.mean(cover),
            "mean_se": statistics.mean(float(r["se"]) for r in good),
            "mean_interval_width": statistics.mean(float(r["upper"])-float(r["lower"]) for r in good),
        }
        calculated["coverage_mcse"] = math.sqrt(calculated["coverage"]*(1-calculated["coverage"])/len(good))
        for name, value in calculated.items():
            assert math.isclose(float(row[name]), value, rel_tol=1e-11, abs_tol=1e-12), (row["scenario"], row["method"], name)


def check_hashes(manifest, base):
    for row in rows(manifest):
        actual = hashlib.sha256((base/row["file"]).read_bytes()).hexdigest()
        assert actual == row["sha256"], row["file"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path)
    args = parser.parse_args()
    benchmark = ROOT/"validation/benchmark"
    check_hashes(benchmark/"source-hashes.csv", ROOT)
    records = rows(benchmark/"replicates.csv")
    summary = rows(benchmark/"summary.csv")
    assert len({(r["scenario"],r["replicate"],r["method"]) for r in records}) == len(records)
    verify_summary(summary,records)
    wrong = [dict(r) for r in summary]
    wrong[0]["bias"] = str(float(wrong[0]["bias"])+0.2)
    try:
        verify_summary(wrong,records)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control: corrupted summary went undetected")
    if args.original:
        check_hashes(ROOT/"validation/original-sha256.csv",args.original)
        print("Original draft file hashes unchanged: PASS")
    print(f"Source hashes and independently recomputed {len(summary)} summary rows: PASS")
    print("Known-wrong bias mutation detected: PASS")
    print("This verifies correspondence and repeatability, not inferential validity.")


if __name__ == "__main__":
    main()
