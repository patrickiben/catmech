#!/usr/bin/env python3
"""Check a built R archive against its checkout and a completed study source freeze."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dcf(data):
    fields = {}
    key = None
    for line in data.decode().splitlines():
        if not line.strip():
            continue
        if line[0].isspace():
            if key is None:
                raise AssertionError("Malformed DESCRIPTION continuation")
            fields[key] += " " + line.strip()
        else:
            key, value = line.split(":", 1)
            fields[key] = value.strip()
    return {k: " ".join(v.split()) for k, v in fields.items()}


def assert_source(row, root):
    data = (root / row["file"]).read_bytes()
    assert len(data) == int(row["bytes"]) and sha(data) == row["sha256"], row["file"]


def verify_equivalence(receipt_path, frozen, root, study, before):
    """Permit only one named portability correction, never an arbitrary waiver."""
    receipt = json.loads(receipt_path.read_text())
    settings = json.loads((study / "settings.json").read_text())
    complete = json.loads((study / "COMPLETE.json").read_text())
    assert receipt["schema_version"] == 1 and receipt["status"] == "pass"
    assert receipt["correction_id"] == "qualify_stats_ave"
    assert receipt["correction_kind"] == "post_confirmation_namespace_portability"
    assert receipt["stage"] == settings["mode"] == complete["mode"] == "confirmation"
    assert complete["status"] == "complete"
    assert receipt["frozen_source_fingerprint"] == settings["source_fingerprint"] == complete["source_fingerprint"]
    for key, file in (("study_complete_sha256", "COMPLETE.json"),
                      ("study_manifest_sha256", "source-hashes-before.csv"),
                      ("study_records_sha256", "records.csv")):
        assert receipt[key] == sha((study / file).read_bytes()), key
    checker = root / "tools/check-marginal-release-equivalence.R"
    assert receipt["check_script_sha256"] == sha(checker.read_bytes()), "Equivalence checker hash differs"
    for key in ("all_frozen_sources_verified", "all_other_frozen_sources_byte_identical",
                "release_source_file_set_unchanged", "all_cohort_input_hashes_match",
                "all_preparation_objects_byte_identical", "legacy_ave_resolved_to_stats_ave",
                "stats_detached_old_prepare_failed", "stats_detached_current_prepare_passed",
                "stats_detached_current_preparation_unchanged"):
        assert receipt[key] is True, key
    assert receipt["fitting_or_accuracy_evaluation_performed"] is False

    sampler = "R/marginal-sampler.R"
    old_text = b"ave(x[, j], ids, FUN = mean)"
    new_text = b"stats::ave(x[, j], ids, FUN = mean)"
    old = (frozen / sampler).read_bytes()
    new = (root / sampler).read_bytes()
    assert old.count(old_text) == 1 and new_text not in old
    assert new == old.replace(old_text, new_text, 1), "Not the exact permitted stats::ave correction"
    wanted_delta = dict(file=sampler, replacement_count=1, old_text=old_text.decode(),
                        new_text=new_text.decode(), old_bytes=len(old), new_bytes=len(new),
                        old_sha256=sha(old), new_sha256=sha(new))
    assert receipt["permitted_source_delta"] == wanted_delta, "Source-delta receipt mismatch"
    frozen_rows = {row["file"]: row for row in receipt["frozen_source_checks"]}
    release_rows = {row["file"]: row for row in receipt["release_source_checks"]}
    assert len(frozen_rows) == len(receipt["frozen_source_checks"]) == len(before)
    assert len(release_rows) == len(receipt["release_source_checks"]) == len(before)
    assert set(frozen_rows) == set(release_rows) == {row["file"] for row in before}
    for row in before:
        assert_source(row, frozen)
        expected = dict(file=row["file"], bytes=int(row["bytes"]), sha256=row["sha256"])
        assert frozen_rows[row["file"]] == expected, row["file"]
        current = (root / row["file"]).read_bytes()
        assert release_rows[row["file"]] == dict(file=row["file"], bytes=len(current), sha256=sha(current))
        if row["file"] != sampler:
            assert current == (frozen / row["file"]).read_bytes(), row["file"]

    cells = settings["cells"]
    assert len(cells) == len(set(cells)) == 9
    assert settings["primary_reps"] == 400 and settings["prior_reps"] == 200
    assert len(settings["prior_cells"]) == len(set(settings["prior_cells"])) == 4
    assert set(settings["prior_cells"]) <= set(cells)
    grid = {(cell, replicate) for cell in cells for replicate in range(1, 401)}
    assert complete["cohorts"] == receipt["cohorts"] == len(grid) == 3600
    assert complete["analysis_requests"] == 27600
    assert complete["mcmc_fit_requests"] == 12000 and complete["fitted_workflow_requests"] == 15600
    records = list(csv.DictReader((study / "records.csv").open()))
    assert len(records) == complete["analysis_requests"]
    inputs = {}
    for row in records:
        key = (row["cell"], int(row["replicate"]))
        assert key in grid and re.fullmatch(r"[0-9a-f]{64}", row["input_sha256"])
        assert inputs.setdefault(key, row["input_sha256"]) == row["input_sha256"]
    assert set(inputs) == grid
    checks = receipt["cohort_checks"]
    assert len(checks) == len(grid)
    seen = set()
    for row in checks:
        key = (row["cell"], row["replicate"])
        assert key in grid and key not in seen, "Duplicate or unexpected equivalence cohort"
        seen.add(key)
        assert row["data_seed"] == settings["data_seed_base"] + (cells.index(key[0]) + 1) * 100000 + key[1]
        assert row["input_sha256"] == inputs[key], "Equivalence input hash mismatch"
        assert row["old_preparation_sha256"] == row["new_preparation_sha256"]
        assert re.fullmatch(r"[0-9a-f]{64}", row["old_preparation_sha256"])
    assert seen == grid
    return dict(mode="explicit_namespace_correction_equivalence", identical_to_confirmation_source=False,
                correction_id="qualify_stats_ave", delta=wanted_delta, cohorts=3600,
                receipt=str(receipt_path.resolve()), receipt_sha256=sha(receipt_path.read_bytes()),
                frozen_source=str(frozen.resolve()), all_other_frozen_sources_byte_identical=True,
                input_and_preparation_equivalence="validated receipt; no new fits or accuracy scoring")


def verify(archive, root, study, equivalence=None, frozen_source=None):
    assert (study / "COMPLETE.json").is_file(), "Study has no completion receipt"
    before = list(csv.DictReader((study / "source-hashes-before.csv").open()))
    after = list(csv.DictReader((study / "source-hashes-after.csv").open()))
    assert before == after, "Study sources changed during execution"
    assert len({r["file"] for r in before}) == len(before), "Duplicate source manifest row"
    scientific = [r for r in before if r["file"].startswith(("R/", "src/"))
                  or r["file"] in {"DESCRIPTION", "NAMESPACE"}]
    actual_scientific = {"DESCRIPTION", "NAMESPACE"}
    actual_scientific.update(p.relative_to(root).as_posix() for p in (root / "R").glob("*.R"))
    actual_scientific.update(p.relative_to(root).as_posix() for p in (root / "src").iterdir()
                             if p.is_file() and (p.suffix in {".c", ".cpp", ".h"} or p.name.startswith("Makevars")))
    assert actual_scientific == {r["file"] for r in scientific}, "Scientific source-file set changed"
    assert (equivalence is None) == (frozen_source is None), "Both --equivalence and --frozen-source are required"
    if equivalence is None:
        for row in before:
            assert_source(row, root)
        source_equivalence = dict(mode="exact_source_freeze", identical_to_confirmation_source=True)
    else:
        source_equivalence = verify_equivalence(equivalence, frozen_source, root, study, before)
    expected = {"DESCRIPTION", "NAMESPACE", "LICENSE", "README.md", "NEWS.md"}
    for directory in ("R", "src", "man", "inst", "tests"):
        for path in (root / directory).rglob("*"):
            if path.is_file() and not path.name.startswith(".") and path.suffix not in {
                ".o", ".so", ".dll", ".dylib"}:
                expected.add(path.relative_to(root).as_posix())
    archived = {}
    with tarfile.open(archive, "r:gz") as bundle:
        for item in bundle.getmembers():
            assert not item.issym() and not item.islnk(), "Archive contains a link"
            if not item.isfile():
                continue
            parts = Path(item.name).parts
            assert parts[0] == "catmech" and ".." not in parts, item.name
            relative = "/".join(parts[1:])
            assert relative not in archived, "Duplicate archive member"
            assert not relative.startswith(("validation/", "tools/", "docs/", ".git/", ".claude/")), relative
            assert not relative.endswith((".tar.gz", ".o", ".so", ".dll", ".dylib")), relative
            archived[relative] = bundle.extractfile(item).read()
    missing = expected - archived.keys()
    assert not missing, f"Archive omitted source files: {sorted(missing)}"
    receipts = []
    for relative, data in sorted(archived.items()):
        if relative == "build/partial.rdb":
            # R CMD build generates this empty partial-Rd cache. Recreate it
            # exactly instead of waiving arbitrary generated binary contents.
            with tempfile.TemporaryDirectory(prefix="catmech-partial-rdb-") as temporary:
                expected_cache = Path(temporary) / "partial.rdb"
                subprocess.run(["Rscript", "--vanilla", "-e",
                    "saveRDS(structure(list(), names=character()), commandArgs(TRUE)[1L], version=2)",
                    str(expected_cache)], check=True, capture_output=True, text=True)
                assert data == expected_cache.read_bytes(), "Unexpected generated partial-Rd cache"
            mode = "R-generated empty named list; byte-matched independent version-2 serialization"
        elif relative == "DESCRIPTION":
            original = (root / relative).read_bytes()
            wanted, actual = dcf(original), dcf(data)
            assert all(actual.get(k) == value for k, value in wanted.items()), "DESCRIPTION semantic mismatch"
            assert set(actual) - set(wanted) <= {"NeedsCompilation", "Packaged", "Author", "Maintainer"}
            assert actual["NeedsCompilation"] == "yes", "Native source not recognized"
            mode = "original fields match after DCF unfolding; generated compilation metadata checked"
        else:
            original = (root / relative).read_bytes()
            assert original == data, f"Archive differs from checkout: {relative}"
            mode = "byte-identical"
        receipts.append(dict(file=relative, bytes=len(data), sha256=sha(data), check=mode))
    return dict(status="PASS", archive=str(archive.resolve()), sha256=sha(archive.read_bytes()),
                bytes=archive.stat().st_size, scientific_source_files=len(scientific),
                archive_files=len(archived), source_freeze=str(study.resolve()),
                release_source_equivalence=source_equivalence, files=receipts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("study", type=Path)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--report", type=Path)
    parser.add_argument("--equivalence", type=Path, help="Explicit receipt for the single stats::ave qualification")
    parser.add_argument("--frozen-source", type=Path, help="Original completed confirmation source archive")
    args = parser.parse_args()
    report = verify(args.archive, args.root, args.study, args.equivalence, args.frozen_source)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "files"}, indent=2))
