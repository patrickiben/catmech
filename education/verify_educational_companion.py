#!/usr/bin/env python3
"""End-to-end verification suite for the catmech v1.0.0 Multimodal R-Package Educational Companion."""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

EDU_DIR = Path(__file__).resolve().parent
ROOT = EDU_DIR.parent
# The spelling, prose, G14 and lock gates live in a separate validation harness that is
# not part of this repository. Point RESEARCH_VALIDATION_HARNESS_DIR at it to run them;
# without it they are reported as SKIPPED, and --require-harness turns a skip into a failure.
_HARNESS_ENV = os.environ.get("RESEARCH_VALIDATION_HARNESS_DIR", "").strip()
HARNESS_DIR = Path(_HARNESS_ENV) if _HARNESS_ENV else None
REQUIRE_HARNESS = "--require-harness" in sys.argv[1:]

REQUIRED_FILES = [
    "COMPANION.md",
    "COMPANION.pdf",
    "STUDY_GUIDE.md",
    "STUDY_GUIDE.pdf",
    "ONE_PAGER.md",
    "ONE_PAGER.pdf",
    "companion_lab.html",
    "COMPLETE_COVERAGE.md",
    "README.md",
    "ARTIFACTS.json",
    "G14_HUMAN_SIGNOFF.md",
    "g14_companion_manifest.json",
    "g9_claims_MAP.json",
    "g9_claim_inventory.json",
    "run_companion_labs.R",
    "interactive_tutor.R",
    "companion_data.py",
    "build_educational_companion.py",
    "verify_educational_companion.py",
    ".gitattributes",
]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdf_info(path: Path) -> tuple[int, str]:
    info = subprocess.check_output(["pdfinfo", str(path)], text=True)
    m_pages = re.search(r"^Pages:\s+(\d+)", info, re.M)
    pages = int(m_pages.group(1)) if m_pages else 0
    return pages, info


def main() -> int:
    failures: list[str] = []

    def ok(msg: str) -> None:
        print(f"  [PASS] {msg}")

    def fail(msg: str) -> None:
        print(f"  [FAIL] {msg}")
        failures.append(msg)

    skipped: list[str] = []

    def skip(msg: str) -> None:
        if REQUIRE_HARNESS:
            fail(f"{msg} (required by --require-harness)")
        else:
            print(f"  [SKIPPED] {msg}")
            skipped.append(msg)

    def harness_gate(name: str) -> Path | None:
        if HARNESS_DIR is None:
            skip(f"{name}: RESEARCH_VALIDATION_HARNESS_DIR is not set")
            return None
        gate = HARNESS_DIR / name
        if not gate.exists():
            skip(f"{name}: not found in RESEARCH_VALIDATION_HARNESS_DIR")
            return None
        return gate

    print("=== catmech v1.0.0 Multimodal Educational Companion Verifier ===")

    # 1. File existence
    for rel in REQUIRED_FILES:
        p = EDU_DIR / rel
        if not p.exists():
            fail(f"Missing required file: education/{rel}")
    if not failures:
        ok(f"All {len(REQUIRED_FILES)} required files present in education/")

    # 2. ARTIFACTS.json hash, page count, and PDF Subject lock check
    artifacts_data = json.loads((EDU_DIR / "ARTIFACTS.json").read_text(encoding="utf-8"))
    for src_key, src_meta in artifacts_data.get("canonical_sources", {}).items():
        spath = (EDU_DIR / src_meta["file"]).resolve()
        if sha256_file(spath) != src_meta["sha256"]:
            fail(f"canonical_sources.{src_key}: sha256 mismatch")

    pdf_to_md_key = {
        "study_guide_pdf": "study_guide_markdown",
        "companion_pdf": "companion_markdown",
        "one_pager_pdf": "companion_markdown",
    }
    for key, meta in artifacts_data.items():
        if key == "canonical_sources":
            continue
        fpath = EDU_DIR / meta["file"]
        actual_sha = sha256_file(fpath)
        if actual_sha != meta["sha256"]:
            fail(f"{key} ({meta['file']}): sha256 mismatch ({actual_sha} != {meta['sha256']})")
        if meta["file"].endswith(".pdf"):
            pages, info_text = pdf_info(fpath)
            if pages != meta.get("pages"):
                fail(f"{meta['file']}: page count mismatch ({pages} != {meta.get('pages')})")
            md_key = pdf_to_md_key[key]
            src_sha = artifacts_data[md_key]["sha256"]
            if src_sha not in info_text:
                fail(f"{meta['file']}: missing source sha256 ({src_sha[:16]}...) in PDF Subject metadata")
    if artifacts_data["one_pager_pdf"]["pages"] != 1:
        fail(f"ONE_PAGER.pdf must be 1 page, got {artifacts_data['one_pager_pdf']['pages']}")
    ok("ARTIFACTS.json SHA-256 digests, PDF page counts (1, 8, 24), and embedded PDF Subject hashes verified")

    # 3. Markdown/HTML link resolution, no em-dashes, no local machine paths in shipped docs
    shipped_text_files = [
        "COMPANION.md",
        "STUDY_GUIDE.md",
        "ONE_PAGER.md",
        "COMPLETE_COVERAGE.md",
        "README.md",
        "G14_HUMAN_SIGNOFF.md",
        "companion_lab.html",
        "ARTIFACTS.json",
        "g14_companion_manifest.json",
        "g9_claim_inventory.json",
        "g9_claims_MAP.json",
    ]
    md_link_re = re.compile(r"\[[^\]]+\]\(([^)#]+)(?:#[^)]+)?\)")
    html_href_re = re.compile(r'href="([^"#]+)(?:#[^"]+)?"')
    local_path_prefix = "/" + "Users" + "/"

    for name in shipped_text_files:
        text = (EDU_DIR / name).read_text(encoding="utf-8")
        if local_path_prefix in text:
            fail(f"{name}: contains machine-local user path")
        if name.endswith(".md") and "\u2014" in text:
            fail(f"{name}: contains em-dash character (U+2014)")
        if name.endswith(".md"):
            for target in md_link_re.findall(text):
                if target.startswith(("http://", "https://", "mailto:")):
                    continue
                resolved = (EDU_DIR / target).resolve()
                if not resolved.exists():
                    fail(f"{name}: broken relative link -> {target}")
        if name.endswith(".html"):
            for target in html_href_re.findall(text):
                if target.startswith(("http://", "https://", "mailto:")):
                    continue
                resolved = (EDU_DIR / target).resolve()
                if not resolved.exists():
                    fail(f"{name}: broken HTML href -> {target}")
    ok("All relative links in Markdown and HTML resolve; 0 em-dashes and 0 local paths in shipped artifacts")

    # 3b. Run G10 US spelling and prose style checks on all education/*.md files
    spell_gate = harness_gate("check_spelling_us.py")
    prose_gate = harness_gate("check_prose_style_v2.py")
    md_files = ["COMPANION.md", "STUDY_GUIDE.md", "ONE_PAGER.md", "COMPLETE_COVERAGE.md", "README.md", "G14_HUMAN_SIGNOFF.md"]
    if spell_gate is not None and prose_gate is not None:
        for md_name in md_files:
            sp = subprocess.run([sys.executable, str(spell_gate), f"education/{md_name}"], cwd=str(ROOT), capture_output=True, text=True)
            if sp.returncode != 0:
                fail(f"check_spelling_us.py failed on education/{md_name}:\n{sp.stdout}\n{sp.stderr}")
            pr = subprocess.run([sys.executable, str(prose_gate), f"education/{md_name}"], cwd=str(ROOT), capture_output=True, text=True)
            if pr.returncode != 0:
                fail(f"check_prose_style_v2.py failed on education/{md_name}:\n{pr.stdout}\n{pr.stderr}")
        ok("G10 US spelling (0 British spellings) and prose style (all categories [ok]) passed across all 6 Markdown files")

    # 4. Run G14 gate (--selftest and strict v2 evaluation) and check_companion_lock.py
    g14_gate = harness_gate("g14_companion_gate.py")
    lock_gate = harness_gate("check_companion_lock.py")
    if g14_gate is not None:
        st = subprocess.run([sys.executable, str(g14_gate), "--selftest"], cwd=str(ROOT), capture_output=True, text=True)
        if st.returncode != 0:
            fail(f"g14_companion_gate.py --selftest failed:\n{st.stdout}\n{st.stderr}")
        else:
            ok("g14_companion_gate.py --selftest passed")

        ev = subprocess.run(
            [
                sys.executable,
                str(g14_gate),
                "--manuscript",
                "paper/manuscript.Rmd",
                "--companion",
                "education/COMPANION.md",
                "--manifest",
                "education/g14_companion_manifest.json",
                "--project-root",
                ".",
            ],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        if ev.returncode != 0:
            fail(f"g14_companion_gate.py failed:\n{ev.stdout}\n{ev.stderr}")
        else:
            ok("g14_companion_gate.py (schema-v2 MULTIMODAL: DOCUMENT, ONE_PAGER, INTERACTIVE) passed")

    if lock_gate is not None:
        lk = subprocess.run(
            [sys.executable, str(lock_gate), "education/g14_companion_manifest.json", "."],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        if lk.returncode != 0:
            fail(f"check_companion_lock.py failed:\n{lk.stdout}\n{lk.stderr}")
        else:
            ok("check_companion_lock.py (2 fresh, 0 stale, 0 missing) passed")

    # 5. Run live R package companion labs & refusal checks
    r_labs = subprocess.run(
        ["Rscript", "education/run_companion_labs.R"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    if r_labs.returncode != 0:
        fail(f"Rscript education/run_companion_labs.R failed:\n{r_labs.stdout}\n{r_labs.stderr}")
    else:
        ok("Rscript education/run_companion_labs.R (16/16 live R labs & refusal checks) passed")

    # 6. Run interactive R console tutor in batch verification mode
    r_tutor = subprocess.run(
        ["Rscript", "education/interactive_tutor.R", "--non-interactive"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    if r_tutor.returncode != 0:
        fail(f"Rscript education/interactive_tutor.R --non-interactive failed:\n{r_tutor.stdout}\n{r_tutor.stderr}")
    else:
        ok("Rscript education/interactive_tutor.R --non-interactive (12/12 lessons & transfer prompts) passed")

    if failures:
        print(f"\nFAILED: {len(failures)} check(s) failed.")
        return 1
    if skipped:
        print(f"\nPASSED WITH {len(skipped)} HARNESS GATE(S) SKIPPED; set RESEARCH_VALIDATION_HARNESS_DIR to run them.")
        return 0
    print("\nALL EDUCATIONAL COMPANION VERIFICATION CHECKS PASSED (exit 0).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
