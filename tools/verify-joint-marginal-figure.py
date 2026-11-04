#!/usr/bin/env python3
"""Independently check every marginal-evidence figure point against raw records.

This verifies plotted arithmetic and file provenance, not model validity or the
visual layout. The companion study verifier checks the complete scored study.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys

MODULE = Path(__file__).with_name("verify-joint-marginal.py")
spec = importlib.util.spec_from_file_location("marginal_verifier", MODULE)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def verify(study, figures):
    settings = json.loads((study/"settings.json").read_text())
    complete = json.loads((study/"COMPLETE.json").read_text())
    provenance = json.loads((figures/"figure-provenance.json").read_text())
    v.require(settings["mode"] == complete["mode"] == "confirmation" and complete["status"] == "complete",
              "Only completed confirmation evidence is supported")
    v.require(provenance["stage"] == "confirmation" and provenance["source_fingerprint"] == settings["source_fingerprint"] == complete["source_fingerprint"], "Figure source identifier differs")
    expected_files = {"COMPLETE.json", "settings.json", "records.csv", "fit-records.csv", "summary.csv", "prior-triples.csv",
                      "prior-summary.csv", "availability.csv", "adoption.json", "source-hashes-before.csv", "source-hashes-after.csv"}
    v.require(set(provenance["input_files"]) == expected_files, "Figure input provenance is incomplete")
    for filename, sha in provenance["input_files"].items():
        v.require(hashlib.sha256((study/filename).read_bytes()).hexdigest() == sha, f"Figure used different input: {filename}")
    script = Path(__file__).with_name("plot-joint-marginal.R")
    v.require(hashlib.sha256(script.read_bytes()).hexdigest() == provenance["script_sha256"], "Plot script changed after generation")
    records, fits = v.read_csv(study/"records.csv"), v.read_csv(study/"fit-records.csv")
    rm, fm = {v.key(r):r for r in records}, {v.key(r,"fit_id"):r for r in fits}
    er, ef = v.expected_grid("confirmation")
    v.require(len(rm)==len(records) and set(rm)==er and len(fm)==len(fits) and set(fm)==ef, "Figure input request grid differs from confirmation protocol")
    actual = v.index_rows(v.read_csv(figures/"plotted-data.csv"),("panel","cell","method"),"figure point")
    expected, noise = {}, {}
    methods = ("baseline_direct","marginal_direct")
    for cell in v.CELLS:
        for method in methods:
            rows = [rm[cell,r,method] for r in range(1,401)]
            good = [r for r in rows if r["status"]=="ok"]
            n = len(good); covers = sum(v.number(r,"covered") for r in rows)
            lo,hi = v.wilson(covers,n)
            expected["A",cell,method] = dict(estimate=covers/n if n else math.nan,lower=lo,upper=hi,n=n,events=covers)
            errors = [v.number(r,"estimate")-v.number(r,"target") for r in good]
            bias,se = v.mean(errors),v.mcse(errors)
            expected["B",cell,method] = dict(estimate=bias,lower=bias-v.Z975*se,upper=bias+v.Z975*se,n=n,events=math.nan)
            family = method.replace("_direct","")
            refused = sum(fm[cell,r,family]["fit_status"]!="ok" for r in range(1,401))
            lo,hi = v.wilson(refused,400)
            expected["D",cell,method] = dict(estimate=refused/400,lower=lo,upper=hi,n=400,events=refused)
    for cell in v.PRIOR_CELLS:
        values = {m:[] for m in methods}; noises = {m:[] for m in methods}
        for r in range(1,201):
            selected = {m:[rm[cell,r,m+s] for s in ("","_prior_quarter","_prior_four")+
                          (("_regression_half","_regression_double") if m=="marginal_direct" else ())] for m in methods}
            if not all(row["status"]=="ok" for group in selected.values() for row in group): continue
            for m,group in selected.items():
                base = v.number(group[0],"estimate")
                hw = v.mean([abs(v.number(z,"estimate")-base) for z in group[1:3]])
                value = max(hw,v.mean([abs(v.number(z,"estimate")-base) for z in group[3:5]])) if m=="marginal_direct" else hw
                values[m].append(value)
                noises[m].append(v.mean([v.number(z,"mcse_mean") for z in group]))
        for method in methods:
            avg,se = v.mean(values[method]),v.mcse(values[method]); n=len(values[method])
            expected["C",cell,method] = dict(estimate=avg,lower=avg-v.Z975*se,upper=avg+v.Z975*se,n=n,events=math.nan)
            noise[cell,method] = dict(common_sets=n,mean_absolute_coefficient_shift=avg,mean_single_summary_mcse=v.mean(noises[method]))
    v.require(set(actual)==set(expected),"Figure point grid differs")
    for k, row in expected.items():
        for field,value in row.items(): v.close(v.number(actual[k],field),value,f"Figure {k}/{field}")
        v.require("Monte Carlo" in actual[k]["uncertainty"], f"Figure uncertainty mislabeled: {k}")
    actual_noise = v.index_rows(v.read_csv(figures/"prior-shift-mcse.csv"),("cell","method"),"MCSE comparison")
    v.require(set(actual_noise)==set(noise),"MCSE comparison grid differs")
    for k,row in noise.items():
        for field,value in row.items(): v.close(v.number(actual_noise[k],field),value,f"MCSE comparison {k}/{field}")
    v.require(provenance["points"]==len(expected)==62,"Figure provenance point count differs")
    png=(figures/"joint-marginal-evidence.png").read_bytes()
    v.require(png[:8]==b'\x89PNG\r\n\x1a\n' and struct.unpack('>II',png[16:24])==(3200,2400),"PNG dimensions differ")
    v.require((figures/"joint-marginal-evidence.pdf").read_bytes().startswith(b'%PDF-'),"Vector PDF missing")
    return dict(status="PASS",points=62,interval_bounds=124,denominators=62,prior_mcse_comparisons=8,
                input_files=len(expected_files),png_dimensions=[3200,2400],
                scope="Independent plotted arithmetic and provenance; visual inspection and statistical validity are separate.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("study",type=Path);p.add_argument("figures",type=Path);p.add_argument("--report",type=Path)
    a=p.parse_args()
    try:
        result=verify(a.study,a.figures);text=json.dumps(result,indent=2)+"\n"
        if a.report:
            v.require(not a.report.exists(),"Refusing to overwrite a verification receipt")
            a.report.write_text(text)
        print(text,end="")
    except (AssertionError,KeyError,ValueError,FileNotFoundError)as e:
        print(f"FAIL: {e}",file=sys.stderr);return 1
    return 0

if __name__=="__main__":sys.exit(main())
