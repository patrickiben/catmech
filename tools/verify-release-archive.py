#!/usr/bin/env python3
"""Check a built R archive against its checkout and a completed study source freeze."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tarfile


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


def verify(archive, root, study):
    assert (study / "COMPLETE.json").is_file(), "Study has no completion receipt"
    before = list(csv.DictReader((study / "source-hashes-before.csv").open()))
    after = list(csv.DictReader((study / "source-hashes-after.csv").open()))
    assert before == after, "Study sources changed during execution"
    scientific = [r for r in before if r["file"].startswith(("R/", "src/"))
                  or r["file"] in {"DESCRIPTION", "NAMESPACE"}]
    for row in scientific:
        data = (root / row["file"]).read_bytes()
        assert len(data) == int(row["bytes"]) and sha(data) == row["sha256"], row["file"]
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
        original = (root / relative).read_bytes()
        if relative == "DESCRIPTION":
            wanted, actual = dcf(original), dcf(data)
            assert all(actual.get(k) == value for k, value in wanted.items()), "DESCRIPTION semantic mismatch"
            assert set(actual) - set(wanted) <= {"NeedsCompilation", "Packaged", "Author", "Maintainer"}
            assert actual["NeedsCompilation"] == "yes", "Native source not recognized"
            mode = "original fields match after DCF unfolding; generated compilation metadata checked"
        else:
            assert original == data, f"Archive differs from checkout: {relative}"
            mode = "byte-identical"
        receipts.append(dict(file=relative, bytes=len(data), sha256=sha(data), check=mode))
    return dict(status="PASS", archive=str(archive.resolve()), sha256=sha(archive.read_bytes()),
                bytes=archive.stat().st_size, scientific_source_files=len(scientific),
                archive_files=len(archived), source_freeze=str(study.resolve()), files=receipts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("study", type=Path)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = verify(args.archive, args.root, args.study)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "files"}, indent=2))
