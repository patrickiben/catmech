#!/usr/bin/env python3
"""Download the public real-data sources named in validation/REAL_DATA_PLAN.md.

Usage: python3 -I tools/real-data-fetch.py RAW_DIRECTORY OUTPUT_DIRECTORY

RAW_DIRECTORY receives the four NHANES transport files and must be new or empty.
OUTPUT_DIRECTORY receives ctgov-pool.csv and fetch-manifest.json. Each identifier
is assigned to a development or evaluation split by its hash.
"""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys
import time
import urllib.parse
import urllib.request

API = "https://clinicaltrials.gov/api/v2/studies"
QUERY = ("AREA[HasResults]true AND AREA[StudyType]INTERVENTIONAL AND "
         "AREA[CompletionDate]RANGE[2015-01-01,2022-12-31]")
FIELDS = ("NCTId|EnrollmentCount|EnrollmentType|Phase|LeadSponsorClass|LocationFacility|"
          "StartDate|PrimaryCompletionDate|ParticipantFlowModule")
NHANES = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/{}.xpt"


def split_of(identifier):
    digest = hashlib.sha256(str(identifier).encode()).hexdigest()
    return "development" if int(digest, 16) % 5 == 0 else "evaluation"


def fetch(url):
    request = urllib.request.Request(url, headers={"User-Agent": "catmech-real-data-study"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))


def months(date):
    parts = [int(p) for p in date.split("-")]
    day = parts[2] if len(parts) > 2 else 15
    return parts[0] * 12 + parts[1] + (day - 1) / 30.4


def phase_group(phases):
    if not phases or phases == ["NA"]:
        return "not_applicable"
    if "PHASE3" in phases or "PHASE4" in phases:
        return "phase3_4"
    if "PHASE2" in phases:
        return "phase2"
    return "early"


def row_for(study):
    p = study["protocolSection"]
    enrollment = p.get("designModule", {}).get("enrollmentInfo", {})
    if enrollment.get("type") != "ACTUAL" or enrollment.get("count") is None:
        return None, "enrollment_not_actual"
    periods = study.get("resultsSection", {}).get("participantFlowModule", {}).get("periods", [])
    started = None
    if periods:
        for milestone in periods[0].get("milestones", []):
            if milestone.get("type", "").upper() == "STARTED":
                try:
                    started = sum(int(a["numSubjects"]) for a in milestone.get("achievements", []))
                except (KeyError, TypeError, ValueError):
                    started = None
    if started is None:
        return None, "no_started_total"
    if enrollment["count"] <= 0 or started <= 0:
        return None, "nonpositive_count"
    status = p.get("statusModule", {})
    start = status.get("startDateStruct", {}).get("date")
    finish = status.get("primaryCompletionDateStruct", {}).get("date")
    try:
        duration = months(finish) - months(start)
    except (AttributeError, ValueError, IndexError):
        return None, "unparseable_date"
    if duration <= 0:
        return None, "nonpositive_duration"
    sites = len(p.get("contactsLocationsModule", {}).get("locations", []))
    if sites == 0:
        return None, "no_location"
    sponsor = p.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("class")
    identifier = p["identificationModule"]["nctId"]
    return dict(nct_id=identifier, split=split_of(identifier), recorded_n=enrollment["count"],
                reference_n=started, duration_months=round(duration, 3), sites=sites,
                start_year=int(start[:4]), phase=phase_group(p.get("designModule", {}).get("phases")),
                sponsor="industry" if sponsor == "INDUSTRY" else "other"), None


def main():
    raw, output = Path(sys.argv[1]), Path(sys.argv[2])
    if raw.exists() and any(raw.iterdir()):
        raise SystemExit("RAW_DIRECTORY must be new or empty.")
    raw.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    if (output/"ctgov-pool.csv").exists():
        raise SystemExit("An existing pool is never overwritten.")
    manifest = {"retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "nhanes": {}}
    for name in ["DEMO_J", "BMX_J", "WHQ_J", "BPX_J"]:
        body = fetch(NHANES.format(name))
        (raw/f"{name}.xpt").write_bytes(body)
        manifest["nhanes"][name] = {"url": NHANES.format(name), "bytes": len(body),
                                    "sha256": hashlib.sha256(body).hexdigest()}
        print(name, len(body), "bytes", flush=True)
    params = {"filter.advanced": QUERY, "fields": FIELDS, "pageSize": "1000", "format": "json", "countTotal": "true"}
    rows, excluded, read, total, token = [], {}, 0, None, None
    streamed = 0
    while True:
        query = dict(params)
        if token:
            query["pageToken"] = token
        body = fetch(API + "?" + urllib.parse.urlencode(query, safe="|[]"))
        streamed += len(body)
        page = json.loads(body)
        total = page.get("totalCount", total)
        for study in page["studies"]:
            read += 1
            row, reason = row_for(study)
            if row:
                rows.append(row)
            else:
                excluded[reason] = excluded.get(reason, 0) + 1
        token = page.get("nextPageToken")
        print(f"ctgov {read}/{total}", flush=True)
        if not token:
            break
        time.sleep(1.3)
    rows.sort(key=lambda r: r["nct_id"])
    assert len({r["nct_id"] for r in rows}) == len(rows)
    with (output/"ctgov-pool.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest["ctgov"] = {"api": API, "filter": QUERY, "fields": FIELDS, "total_matching": total, "records_read": read,
                         "rows_retained": len(rows), "excluded": excluded, "bytes_streamed": streamed,
                         "pool_sha256": hashlib.sha256((output/"ctgov-pool.csv").read_bytes()).hexdigest()}
    (output/"fetch-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["ctgov"], indent=2))


if __name__ == "__main__":
    main()
