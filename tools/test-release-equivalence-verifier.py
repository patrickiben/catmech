#!/usr/bin/env python3
"""Exercise release-verifier acceptance/refusal policy using synthetic files only.

This is not a study rerun or an assertion that real preparation equivalence holds.
The real R checker must still run on all completed confirmation cohorts.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tarfile
import tempfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(rows)


def run():
    repository = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("release_verifier", repository / "tools/verify-release-archive.py")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    checks = []
    with tempfile.TemporaryDirectory(prefix="catmech-equivalence-controls-") as temporary:
        work = Path(temporary)
        root, frozen, study = [work / x for x in ("root", "frozen", "study")]
        for directory in (root, frozen, study):
            directory.mkdir()
        sampler = "R/marginal-sampler.R"
        old_call = b"ave(x[, j], ids, FUN = mean)"
        new_call = b"stats::ave(x[, j], ids, FUN = mean)"
        scientific = {
            "DESCRIPTION": b"Package: catmech\nVersion: 0.10.0\nNeedsCompilation: yes\n",
            "NAMESPACE": b"# synthetic namespace\n",
            "validation/JOINT_MARGINAL_PLAN.md": b"synthetic protocol\n",
            "validation/run-joint-marginal.R": b"# synthetic runner\n",
            sampler: b"prepare <- function(x, ids, j) " + old_call + b"\n",
            "src/marginal.cpp": b"// synthetic C++ source\n",
            "src/init.cpp": b"// synthetic registration\n",
            "src/Makevars": b"CXX_STD = CXX17\n",
            "src/Makevars.win": b"CXX_STD = CXX17\n",
        }
        scientific.update({f"R/fixture{i}.R": f"value{i} <- {i}\n".encode() for i in range(14)})
        scientific.update({f"validation/v0.9.0-source/fixture{i}": b"baseline\n" for i in range(12)})
        assert len(scientific) == 35
        for relative, data in scientific.items():
            write(root / relative, data)
            write(frozen / relative, data)
        for relative in ("LICENSE", "README.md", "NEWS.md"):
            write(root / relative, b"synthetic archive fixture\n")
        write(root / "tools/check-marginal-release-equivalence.R",
              (repository / "tools/check-marginal-release-equivalence.R").read_bytes())
        before = [dict(file=name, bytes=len(data), sha256=sha(data)) for name, data in scientific.items()]
        for filename in ("source-hashes-before.csv", "source-hashes-after.csv"):
            write_csv(study / filename, before, ["file", "bytes", "sha256"])
        cells = [f"cell{i}" for i in range(9)]
        settings = dict(mode="confirmation", cells=cells, primary_reps=400,
                        prior_reps=200, prior_cells=cells[:4], data_seed_base=1710000000,
                        source_fingerprint="f" * 64)
        complete = dict(status="complete", mode="confirmation", cohorts=3600,
                        analysis_requests=27600, mcmc_fit_requests=12000,
                        fitted_workflow_requests=15600, source_fingerprint="f" * 64)
        write_json(study / "settings.json", settings)
        write_json(study / "COMPLETE.json", complete)
        records, cohort_checks = [], []
        for index, cell in enumerate(cells, 1):
            for replicate in range(1, 401):
                value = f"{cell}:{replicate}".encode()
                input_hash, preparation_hash = sha(value), sha(b"prepared:" + value)
                count = 17 if index <= 4 and replicate <= 200 else 5
                records.extend(dict(cell=cell, replicate=replicate, input_sha256=input_hash) for _ in range(count))
                cohort_checks.append(dict(cell=cell, replicate=replicate,
                    data_seed=settings["data_seed_base"] + index * 100000 + replicate,
                    input_sha256=input_hash, old_preparation_sha256=preparation_hash,
                    new_preparation_sha256=preparation_hash))
        write_csv(study / "records.csv", records, ["cell", "replicate", "input_sha256"])
        archive = work / "catmech.tar.gz"

        def bundle():
            with tarfile.open(archive, "w:gz") as out:
                for path in root.rglob("*"):
                    if path.is_file() and path.relative_to(root).parts[0] not in {"validation", "tools"}:
                        out.add(path, arcname="catmech/" + path.relative_to(root).as_posix())

        bundle()
        verifier.verify(archive, root, study)
        checks.append(dict(case="unchanged source passes strict default", result="pass"))
        old = scientific[sampler]
        new = old.replace(old_call, new_call, 1)
        write(root / sampler, new)
        bundle()
        release = copy.deepcopy(before)
        for row in release:
            if row["file"] == sampler:
                row.update(bytes=len(new), sha256=sha(new))
        receipt = dict(schema_version=1, status="pass", correction_id="qualify_stats_ave",
            correction_kind="post_confirmation_namespace_portability", stage="confirmation",
            frozen_source_fingerprint="f" * 64,
            study_complete_sha256=sha((study / "COMPLETE.json").read_bytes()),
            study_manifest_sha256=sha((study / "source-hashes-before.csv").read_bytes()),
            study_records_sha256=sha((study / "records.csv").read_bytes()),
            check_script_sha256=sha((root / "tools/check-marginal-release-equivalence.R").read_bytes()),
            permitted_source_delta=dict(file=sampler, replacement_count=1,
                old_text=old_call.decode(), new_text=new_call.decode(), old_bytes=len(old),
                new_bytes=len(new), old_sha256=sha(old), new_sha256=sha(new)),
            cohorts=3600, fitting_or_accuracy_evaluation_performed=False,
            frozen_source_checks=before, release_source_checks=release, cohort_checks=cohort_checks)
        for key in ("all_frozen_sources_verified", "all_other_frozen_sources_byte_identical",
                    "release_source_file_set_unchanged", "all_cohort_input_hashes_match",
                    "all_preparation_objects_byte_identical", "legacy_ave_resolved_to_stats_ave",
                    "stats_detached_old_prepare_failed", "stats_detached_current_prepare_passed",
                    "stats_detached_current_preparation_unchanged"):
            receipt[key] = True
        receipt_path = work / "equivalence.json"
        write_json(receipt_path, receipt)
        report = verifier.verify(archive, root, study, receipt_path, frozen)
        assert report["release_source_equivalence"]["identical_to_confirmation_source"] is False
        checks.append(dict(case="exact named correction passes with explicit receipt", result="pass"))

        def reject(label, callback):
            try:
                callback()
            except AssertionError:
                checks.append(dict(case=label, result="rejected"))
            else:
                raise AssertionError("Tamper accepted: " + label)

        call = lambda: verifier.verify(archive, root, study, receipt_path, frozen)
        reject("changed source rejected without waiver", lambda: verifier.verify(archive, root, study))
        reject("receipt without frozen directory rejected", lambda: verifier.verify(archive, root, study, receipt_path))

        def altered_file(label, path, replacement):
            original = path.read_bytes()
            try:
                path.write_bytes(replacement)
                reject(label, call)
            finally:
                path.write_bytes(original)

        altered_file("unrelated R change", root / "R/fixture0.R", b"value0 <- 99\n")
        altered_file("native source change", root / "src/marginal.cpp", b"// changed\n")
        altered_file("extra sampler change", root / sampler, new + b"# another edit\n")
        altered_file("archived source alteration", frozen / "R/fixture0.R", b"value0 <- 99\n")
        altered_file("study source manifest change", study / "source-hashes-after.csv", b"changed\n")
        altered_file("study input records change", study / "records.csv", b"changed\n")
        altered_file("equivalence checker change", root / "tools/check-marginal-release-equivalence.R", b"changed\n")

        def altered_receipt(label, mutate):
            changed = copy.deepcopy(receipt)
            mutate(changed)
            try:
                write_json(receipt_path, changed)
                reject(label, call)
            finally:
                write_json(receipt_path, receipt)

        altered_receipt("cohort input hash change", lambda x: x["cohort_checks"][0].update(input_sha256="0" * 64))
        altered_receipt("preparation hash mismatch", lambda x: x["cohort_checks"][0].update(new_preparation_sha256="0" * 64))
        altered_receipt("missing cohort", lambda x: x["cohort_checks"].pop())
        altered_receipt("duplicate cohort", lambda x: x["cohort_checks"].__setitem__(1, x["cohort_checks"][0]))
        altered_receipt("new source hash change", lambda x: x["permitted_source_delta"].update(new_sha256="0" * 64))
        altered_receipt("namespace test failure", lambda x: x.update(stats_detached_current_prepare_passed=False))
        altered_receipt("incorrect cohort seed", lambda x: x["cohort_checks"][0].update(data_seed=1))
        altered_receipt("unapproved correction name", lambda x: x.update(correction_id="arbitrary_change"))
    return dict(status="PASS", evidence="synthetic verifier-policy fixtures; no real study evaluated",
                positive_controls=2, rejected_controls=len(checks) - 2, checks=checks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = run()
    if args.report:
        write_json(args.report, report)
    print(json.dumps(report, indent=2))
