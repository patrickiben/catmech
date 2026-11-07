#!/usr/bin/env python3
"""G2 / G17 / G19 / G22 deterministic verifier for paper/manuscript.Rmd, paper/manuscript.tex, and paper/paper.md."""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
VAL = ROOT / "validation"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def verify_synthetic_claims(manuscript_rmd: str, manuscript_tex: str, joss_md: str) -> int:
    summary = read_csv(VAL / "contamination-study" / "summary.csv")
    lookup = {(r["scenario"], r["method"]): r for r in summary}

    # Check total request counts from summary.csv
    total_req = sum(int(r["requested"]) for r in summary)
    total_ret = sum(int(r["returned_intervals"]) for r in summary)
    total_ref = sum(int(r["refused"]) for r in summary)
    total_conv = sum(int(r["convergence_failures"]) for r in summary)
    total_err = sum(int(r["other_errors"]) for r in summary)
    assert (total_req, total_ret, total_ref, total_conv, total_err) == (52950, 52900, 16, 33, 1)

    checked = 0
    for doc_name, text in [("manuscript.Rmd", manuscript_rmd), ("manuscript.tex", manuscript_tex), ("paper.md", joss_md)]:
        for token in ["52,950", "52,900", "0.800"]:
            assert token in text, f"{doc_name} missing {token}"
            checked += 1

    # Check every row in missed_flags and error_gross_and_small
    truth = 0.8
    for scenario in ["missed_flags", "error_gross_and_small"]:
        for method in ["oracle", "v02_mask_normal", "gold_only", "reference_normal", "cole_mi", "raking", "contamination"]:
            row = lookup[(scenario, method)]
            b = float(row["bias"])
            est = f"{truth + b:.4f}"
            bias_str = f"{b:+.4f}" if b > 0 else f"{b:.4f}"
            mcse_str = f"({float(row['bias_mcse']):.4f})"
            emp_sd = f"{float(row['empirical_sd']):.4f}"
            mean_se = f"{float(row['mean_se']):.4f}"
            cov = f"{100 * float(row['coverage_conditional']):.1f}%"
            width = f"{float(row['mean_width']):.4f}"
            rmse = f"{float(row['rmse']):.4f}"
            for doc_name, text in [("manuscript.Rmd", manuscript_rmd), ("paper.md", joss_md)]:
                for val in [est, bias_str, mcse_str, emp_sd, mean_se, cov, width, rmse]:
                    assert val in text, f"{doc_name} missing {scenario}/{method} value {val}"
                    checked += 1

    # Check G17 primary screen honesty disclosure on site_shift_declared
    ssd = lookup[("site_shift_declared", "contamination")]
    assert f"{100 * float(ssd['coverage_conditional']):.1f}" == "93.7"
    assert f"{float(ssd['coverage_wilson_lower']):.4f}" == "0.9202"
    assert f"{float(ssd['coverage_wilson_upper']):.4f}" == "0.9505"
    assert ssd["approximate_nominal_criterion"] == "FALSE"
    for doc_name, text in [("manuscript.Rmd", manuscript_rmd), ("paper.md", joss_md)]:
        for token in ["93.7%", "0.9202", "0.9505", "-0.0085"]:
            assert token in text, f"{doc_name} missing G17 honesty disclosure token {token}"
            checked += 1

    return checked


def verify_real_data_claims(manuscript_rmd: str, manuscript_tex: str, joss_md: str) -> int:
    summary = read_csv(VAL / "real-data-study" / "summary.csv")
    lookup = {(r["family"], r["scenario"], r["method"]): r for r in summary}
    targets = {r["dataset"]: r for r in read_csv(VAL / "real-data-study" / "pool-targets.csv")}
    overall = {r["dataset"]: r for r in read_csv(VAL / "real-data-study" / "pool-audit-overall.csv")}

    total_req = sum(int(r["requested"]) for r in summary)
    assert total_req == 33800
    assert int(targets["ctgov"]["pool_rows"]) == 27208
    assert int(targets["nhanes"]["pool_rows"]) == 3908
    assert f"{float(targets['ctgov']['target_coefficient']):.4f}" == "0.0357"
    assert f"{float(targets['nhanes']['target_coefficient']):.3f}" == "0.127"
    assert f"{100 * float(overall['ctgov']['disagreement_rate']):.1f}%" == "14.2%"
    assert f"{100 * float(overall['nhanes']['disagreement_rate']):.1f}%" == "67.3%"

    checked = 0
    for family in ["ctgov", "nhanes"]:
        for method in ["reference_normal", "contamination", "cole_mi", "oracle"]:
            row = lookup[(family, "audit_25pct", method)]
            b = float(row["bias"])
            bias_str = f"{b:+.4f}" if b > 0 else f"{b:.4f}"
            cov_str = f"{100 * float(row['coverage_conditional']):.1f}%"
            width_str = f"{float(row['mean_width']):.3f}"
            rmse_str = f"{float(row['rmse']):.4f}"
            for doc_name, text in [("manuscript.Rmd", manuscript_rmd), ("paper.md", joss_md)]:
                for val in [bias_str, cov_str, width_str, rmse_str]:
                    assert val in text, f"{doc_name} missing {family}/{method} value {val}"
                    checked += 1
    return checked


def verify_clustered_claims(manuscript_rmd: str, joss_md: str) -> int:
    decisions = (VAL / "variance-mi-decisions.txt").read_text(encoding="utf-8")
    assert "Requests: 148500; error 1, ok 148486, refused 13" in decisions
    assert "R3: mean distance 1.47 against 3.18; holds True" in decisions
    assert "OUTCOME: cv3" in decisions

    recount_cv3 = [
        float(m.group(1))
        for m in re.finditer(r"oracle cr1s [0-9.]+, cr2 [0-9.]+, cv3 ([0-9.]+)", decisions)
    ]
    assert len(recount_cv3) == 11
    assert f"{min(recount_cv3):.1f}%" == "92.0%"
    assert f"{max(recount_cv3):.1f}%" == "97.5%"

    checked = 0
    for doc_name, text in [("manuscript.Rmd", manuscript_rmd), ("paper.md", joss_md)]:
        for token in ["148,500", "148,486", "92.0%", "97.5%"]:
            assert token in text, f"{doc_name} missing clustered token {token}"
            checked += 1
    return checked


def verify_code_surface(manuscript_rmd: str, joss_md: str) -> int:
    namespace_text = (ROOT / "NAMESPACE").read_text(encoding="utf-8")
    exported = set(re.findall(r"export\(([A-Za-z0-9_.]+)\)", namespace_text))

    checked = 0
    for doc_name, text in [("manuscript.Rmd", manuscript_rmd), ("paper.md", joss_md)]:
        fn_mentions = set(re.findall(r"\b(catmech_[a-z0-9_]+|simulate_catmech_[a-z0-9_]+)\b", text))
        for fn in sorted(fn_mentions):
            assert fn in exported, f"{doc_name} mentions non-exported function {fn}"
            checked += 1
        script_mentions = set(re.findall(r"\b([A-Za-z0-9_.-]+\.(?:py|R))\b", text))
        for name in sorted(script_mentions):
            assert (ROOT / "tools" / name).is_file() or (PAPER / name).is_file(), (
                f"{doc_name} mentions missing script {name}"
            )
            checked += 1
    return checked


def main() -> int:
    manuscript_rmd = (PAPER / "manuscript.Rmd").read_text(encoding="utf-8")
    manuscript_tex = (PAPER / "manuscript.tex").read_text(encoding="utf-8")
    joss_md = (PAPER / "paper.md").read_text(encoding="utf-8")

    c1 = verify_synthetic_claims(manuscript_rmd, manuscript_tex, joss_md)
    c2 = verify_real_data_claims(manuscript_rmd, manuscript_tex, joss_md)
    c3 = verify_clustered_claims(manuscript_rmd, joss_md)
    c4 = verify_code_surface(manuscript_rmd, joss_md)

    # Negative control: corrupting a table value must trigger AssertionError
    try:
        verify_real_data_claims(manuscript_rmd.replace("0.0286", "0.9999"), manuscript_tex, joss_md)
    except AssertionError:
        pass
    else:
        raise AssertionError("Negative control failed: corrupted RMSE went undetected")

    print(
        f"G2/G17/G19/G22 paper verification: {c1} synthetic checks, {c2} real-data checks, "
        f"{c3} clustered checks, {c4} code/script surface checks, negative control: PASS"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
