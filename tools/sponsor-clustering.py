#!/usr/bin/env python3
"""Sponsor clustering screen on the frozen ClinicalTrials.gov pool.

Usage:
  python3 -I tools/sponsor-clustering.py fetch RAW_DIRECTORY OUTPUT_DIRECTORY
  python3 -I tools/sponsor-clustering.py run development NEW_OUTPUT_DIRECTORY
  python3 -I tools/sponsor-clustering.py run evaluation validation/sponsor-clustering/evaluation

The design is fixed in validation/SPONSOR_CLUSTERING_PLAN.md. Run this file with
python3 -I because it reads downloaded data. Development output is engineering only.
"""
from collections import Counter
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import platform
import random
import sys
import time
import urllib.parse
import urllib.request

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
API = "https://clinicaltrials.gov/api/v2/studies"
QUERY = ("AREA[HasResults]true AND AREA[StudyType]INTERVENTIONAL AND "
         "AREA[CompletionDate]RANGE[2015-01-01,2022-12-31]")
FIELDS = "NCTId|LeadSponsorName"
POOL = ROOT/"validation/real-data/ctgov-pool.csv"
POOL_SHA256 = "c2e02ea2542b8291dff54799a9f0dfecfdbcadb71ebf28e90aa798764e0a91ca"
PLAN = ROOT/"validation/SPONSOR_CLUSTERING_PLAN.md"
PLAN_SHA256 = "d8b330b40fea6e929661e0e510b2dbc0f49f859bab392878e3e7e380bb80d222"
SPONSOR_MAP = ROOT/"validation/sponsor-clustering/sponsor-map.csv"
MAP_SHA256 = "31a223e87de063d0def6d6befc5f9a4626ddfba9a6fe55291dca56144f46ec81"
EVALUATION_OUTPUT = ROOT/"validation/sponsor-clustering/evaluation"
PHASES = ["early", "phase2", "phase3_4", "not_applicable"]
CLASSES = ["industry", "other"]
OUTCOMES = ["D", "D2"]
SEED_PERMUTE = 20261101
SEED_REVEAL = 20261102
SEED_BOOTSTRAP = 20261103
SEED_AUDIT = 20261104
SEED_BOOTSTRAP_AUDIT = 20261105
PERMUTATIONS = 2000
BOOTSTRAPS = 2000
AUDIT_SEEDS = 200


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def normalize(name):
    return name.strip().casefold()


def download(url):
    request = urllib.request.Request(url, headers={"User-Agent": "catmech-real-data-study"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))


def fetch_command(raw, output):
    raw, output = Path(raw), Path(output)
    if raw.exists() and any(raw.iterdir()):
        raise SystemExit("RAW_DIRECTORY must be new or empty.")
    raw.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    target = output/"sponsor-map.csv"
    if target.exists() or (output/"fetch-manifest.json").exists():
        raise SystemExit("An existing sponsor map is never overwritten.")
    if sha256_file(POOL) != POOL_SHA256:
        raise SystemExit("The frozen pool hash does not match.")
    retrieved = datetime.now(timezone.utc).isoformat(timespec="seconds")
    params = {"filter.advanced": QUERY, "fields": FIELDS, "pageSize": "1000", "format": "json",
              "countTotal": "true"}
    names, pages, read, total, token, streamed = {}, [], 0, None, None, 0
    missing_name, duplicate = 0, 0
    while True:
        query = dict(params)
        if token:
            query["pageToken"] = token
        body = download(API + "?" + urllib.parse.urlencode(query, safe="|[]"))
        streamed += len(body)
        if streamed > 120_000_000:
            raise SystemExit("Approved download size of about 120 MB exceeded.")
        name = f"page-{len(pages)+1:04d}.json"
        (raw/name).write_bytes(body)
        pages.append(dict(file=name, bytes=len(body), sha256=hashlib.sha256(body).hexdigest()))
        page = json.loads(body)
        total = page.get("totalCount", total)
        for study in page["studies"]:
            read += 1
            section = study.get("protocolSection", {})
            identifier = section.get("identificationModule", {}).get("nctId")
            sponsor = section.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name")
            if identifier in names:
                duplicate += 1
            if not isinstance(sponsor, str) or not sponsor.strip():
                missing_name += 1
                continue
            names[identifier] = sponsor
        token = page.get("nextPageToken")
        print(f"ctgov {read}/{total}", flush=True)
        if not token:
            break
        time.sleep(1.3)
    pool = read_csv(POOL)
    ids = [r["nct_id"] for r in pool]
    assert len(set(ids)) == len(ids)
    matched = [r for r in pool if r["nct_id"] in names]
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["nct_id", "lead_sponsor_name", "sponsor_key"])
        for r in matched:
            writer.writerow([r["nct_id"], names[r["nct_id"]], normalize(names[r["nct_id"]])])
    by_split = {}
    for split in ["development", "evaluation"]:
        in_split = [r for r in pool if r["split"] == split]
        hit = [r for r in in_split if r["nct_id"] in names]
        by_split[split] = dict(pool_rows=len(in_split), matched=len(hit), rate=len(hit)/len(in_split))
    raw_names = {names[r["nct_id"]] for r in matched}
    keys = {normalize(n) for n in raw_names}
    manifest = dict(
        retrieved_utc=retrieved, api=API, filter=QUERY, fields=FIELDS, page_size=1000,
        total_matching=total, records_read=read, records_without_name=missing_name,
        duplicate_records=duplicate, bytes_streamed=streamed, pages=len(pages),
        raw_directory_note="Raw pages were kept outside the repository; hashes are listed here.",
        pool_sha256=sha256_file(POOL), pool_rows=len(pool), pool_rows_matched=len(matched),
        match_rate=len(matched)/len(pool), match_by_split=by_split,
        api_records_not_in_pool=len(set(names)-set(ids)),
        distinct_raw_names_in_pool=len(raw_names), distinct_normalized_names_in_pool=len(keys),
        sponsor_map_sha256=sha256_file(target), raw_pages=pages)
    (output/"fetch-manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "raw_pages"}, indent=1))


def band_of(year):
    return 0 if year <= 2013 else 1 if year <= 2015 else 2 if year <= 2017 else 3


def load_trials(split):
    names = {r["nct_id"]: r["sponsor_key"] for r in read_csv(SPONSOR_MAP)}
    pool = read_csv(POOL)
    ids = [r["nct_id"] for r in pool]
    assert ids == sorted(set(ids))
    chosen = [r for r in pool if r["split"] == split]
    recorded = np.array([int(r["recorded_n"]) for r in chosen], dtype=float)
    reference = np.array([int(r["reference_n"]) for r in chosen], dtype=float)
    assert (recorded > 0).all() and (reference > 0).all()
    data = dict(ids=[r["nct_id"] for r in chosen], key=[names[r["nct_id"]] for r in chosen],
                phase=np.array([PHASES.index(r["phase"]) for r in chosen]),
                cls=np.array([CLASSES.index(r["sponsor"]) for r in chosen]),
                band=np.array([band_of(int(r["start_year"])) for r in chosen]),
                recorded=recorded, reference=reference,
                D=(recorded != reference).astype(float),
                D2=(np.abs(np.log(recorded/reference)) > 0.25).astype(float))
    classes = {}
    for key, cls in zip(data["key"], data["cls"]):
        classes.setdefault(key, set()).add(int(cls))
    assert all(len(v) == 1 for v in classes.values()), "A sponsor has two classes"
    data["class_of"] = {k: next(iter(v)) for k, v in classes.items()}
    data["tie"] = np.argsort(np.argsort([hashlib.sha256(i.encode()).hexdigest() for i in data["ids"]]))
    return data


def percentile_limits(values):
    good = np.asarray(values, dtype=float)
    good = good[np.isfinite(good)]
    low, high = np.quantile(good, [.025, .975])
    return float(low), float(high), int(len(values)-good.size)


def test_one(data, add):
    counts = Counter(data["key"])
    index = np.array([i for i, k in enumerate(data["key"]) if counts[k] >= 3])
    sponsors = sorted({data["key"][i] for i in index})
    position = {k: j for j, k in enumerate(sponsors)}
    owner = np.array([position[data["key"][i]] for i in index])
    size = data["recorded"][index]
    cuts = np.quantile(size, [.25, .5, .75])
    quartile = np.searchsorted(cuts, size, side="left")
    code = (data["phase"][index]*2 + data["cls"][index])*16 + data["band"][index]*4 + quartile
    _, strata = np.unique(code, return_inverse=True)
    m, trials = len(sponsors), index.size
    sponsor_class = np.array([data["class_of"][k] for k in sponsors])
    grouped = np.argsort(strata, kind="stable")
    outcome = {}
    for name in OUTCOMES:
        y = data[name][index]
        rate = np.bincount(strata, weights=y)/np.bincount(strata)
        expected = np.bincount(owner, weights=rate[strata], minlength=m)
        outcome[name] = dict(y=y, expected=expected, positive=expected > 0, null=np.zeros((PERMUTATIONS, 3)))

    def terms(name, y):
        state = outcome[name]
        observed = np.bincount(owner, weights=y, minlength=m)
        value = np.zeros(m)
        good = state["positive"]
        value[good] = (observed[good]-state["expected"][good])**2/state["expected"][good]
        return value, observed

    rng = np.random.Generator(np.random.PCG64(SEED_PERMUTE))
    shuffled = np.empty(trials)
    for b in range(PERMUTATIONS):
        order = np.lexsort((rng.random(trials), strata))
        for name in OUTCOMES:
            shuffled[grouped] = outcome[name]["y"][order]
            value, _ = terms(name, shuffled)
            outcome[name]["null"][b] = value.sum(), value[sponsor_class == 0].sum(), value[sponsor_class == 1].sum()
    results = {}
    for name in OUTCOMES:
        state = outcome[name]
        value, observed = terms(name, state["y"])
        state.update(term=value, observed=observed)
        parts = [("all", value.sum(), 0, np.ones(m, bool)), ("industry", value[sponsor_class == 0].sum(), 1, sponsor_class == 0),
                 ("other", value[sponsor_class == 1].sum(), 2, sponsor_class == 1)]
        for subgroup, statistic, column, chosen in parts:
            null = state["null"][:, column]
            p = (1 + int((null >= statistic).sum()))/(PERMUTATIONS + 1)
            ratio = statistic/null.mean()
            in_part = chosen[owner]
            add("test1", name, subgroup, "sponsors", int(chosen.sum()))
            add("test1", name, subgroup, "trials", int(in_part.sum()))
            add("test1", name, subgroup, "discrepancies", int(state["y"][in_part].sum()))
            add("test1", name, subgroup, "statistic", float(statistic))
            add("test1", name, subgroup, "null_mean", float(null.mean()), float(np.quantile(null, .025)), float(np.quantile(null, .975)))
            add("test1", name, subgroup, "p_value", p)
            add("test1", name, subgroup, "observed_to_null_ratio", float(ratio))
            if subgroup == "all":
                results[name] = dict(p=p, ratio=float(ratio))
    add("test1", "both", "all", "nonempty_strata", int(strata.max())+1)
    for label, value in zip(["q1", "q2", "q3"], cuts):
        add("test1", "both", "all", f"size_cut_{label}", float(value))
    table = {k: dict(trials=int(counts[k]), sponsor_class=CLASSES[data["class_of"][k]]) for k in sponsors}
    for name in OUTCOMES:
        for j, k in enumerate(sponsors):
            table[k].update({f"observed_{name}": int(outcome[name]["observed"][j]),
                             f"expected_{name}": float(outcome[name]["expected"][j]),
                             f"term_{name}": float(outcome[name]["term"][j])})
    null_rows = [dict(permutation=b+1, **{f"{part}_{name}": float(outcome[name]["null"][b][c])
                                          for name in OUTCOMES
                                          for c, part in enumerate(["statistic", "industry", "other"])})
                 for b in range(PERMUTATIONS)]
    return results, table, null_rows


def split_sponsors(data):
    counts = Counter(data["key"])
    eligible = sorted((k for k in counts if counts[k] >= 4), key=lambda k: hashlib.sha256(k.encode()).hexdigest())
    first = (len(eligible)+1)//2
    return eligible[:first], eligible[first:]


def reveal_trials(data, eligible):
    wanted = set(eligible)
    rnd = random.Random(SEED_REVEAL)
    draws = {}
    for i, key in enumerate(data["key"]):
        if key in wanted:
            draws.setdefault(key, []).append((rnd.random(), i))
    revealed = np.zeros(len(data["key"]), bool)
    for key, values in draws.items():
        for _, i in sorted(values)[:len(values)//2]:
            revealed[i] = True
    return revealed


def design(data, index, history=None):
    columns = [np.ones(index.size)]
    names = ["intercept"]
    for level in range(1, 4):
        columns.append((data["phase"][index] == level).astype(float))
        names.append(f"phase_{PHASES[level]}")
    columns.append((data["cls"][index] == 0).astype(float))
    names.append("class_industry")
    for level, label in zip(range(1, 4), ["2014_2015", "2016_2017", "2018_or_later"]):
        columns.append((data["band"][index] == level).astype(float))
        names.append(f"start_{label}")
    columns.append(np.log(data["recorded"][index]))
    names.append("log_recorded_n")
    if history is not None:
        columns.append(history[index])
        names.append("revealed_logit")
    return np.column_stack(columns), names


def fit_logistic(x, y):
    beta = np.zeros(x.shape[1])
    for step in range(1, 101):
        mu = 1/(1 + np.exp(-(x @ beta)))
        delta = np.linalg.solve((x.T*(mu*(1-mu))) @ x, x.T @ (y-mu))
        beta = beta + delta
        if np.max(np.abs(delta)) < 1e-10:
            return beta, step
    raise RuntimeError("Newton's method did not converge within 100 steps")


def tied_groups(score):
    values, inverse = np.unique(score, return_inverse=True)
    return inverse, values.size


def weighted_auc(groups, y, w):
    inverse, size = groups
    pos = np.bincount(inverse, weights=w*y, minlength=size)
    neg = np.bincount(inverse, weights=w*(1-y), minlength=size)
    p_total, n_total = pos.sum(), neg.sum()
    if p_total == 0 or n_total == 0:
        return float("nan")
    below = np.concatenate(([0.0], np.cumsum(neg)[:-1]))
    return float((pos*(below + 0.5*neg)).sum()/(p_total*n_total))


def weighted_brier(score, y, w):
    total = w.sum()
    return float((w*(score-y)**2).sum()/total) if total else float("nan")


def audit_once(chosen, y, label, first_stage, second_stage, flagged, audited):
    found_random = y[chosen].sum()
    first = chosen[:first_stage]
    hit = y[first]
    flagged[:] = False
    flagged[label[first[hit == 1]]] = True
    audited[:] = False
    audited[first] = True
    candidates = flagged[label] & ~audited
    picks = np.flatnonzero(candidates)[:second_stage]
    found_second = y[picks].sum()
    shortfall = second_stage - picks.size
    if shortfall:
        found_second += y[np.flatnonzero(~audited & ~candidates)[:shortfall]].sum()
    return found_random, hit.sum() + found_second, found_second, int(flagged.sum()), shortfall


def point_audit(data, y_all, population, prediction, name, subgroup, add, seed_rows):
    """population is in nct_id order; prediction is the baseline score of every trial."""
    population = np.asarray(population)
    n = population.size
    budget = (n + 5)//10
    first_stage = budget//2
    second_stage = budget-first_stage
    rank = np.lexsort((data["tie"][population], -prediction[population]))
    position = np.empty(n, dtype=int)
    position[rank] = np.arange(n)
    sponsors = sorted({data["key"][i] for i in population})
    numbering = {k: j for j, k in enumerate(sponsors)}
    y = y_all[population][rank]
    label = np.array([numbering[data["key"][i]] for i in population])[rank]
    found_risk = float(y[:budget].sum())
    flagged, audited = np.zeros(len(sponsors), bool), np.zeros(n, bool)
    random_found, adaptive_found, flagged_count, second_found, shortfalls = [], [], [], [], []
    for s in range(AUDIT_SEEDS):
        rnd = random.Random(SEED_AUDIT + s)
        draws = np.array([rnd.random() for _ in range(n)])
        chosen = position[np.argsort(draws)[:budget]]
        r, a, sec, f, short = audit_once(chosen, y, label, first_stage, second_stage, flagged, audited)
        random_found.append(float(r))
        adaptive_found.append(float(a))
        second_found.append(float(sec))
        flagged_count.append(f)
        shortfalls.append(short)
        seed_rows.append(dict(outcome=name, subgroup=subgroup, seed=s, budget=budget, found_risk=int(found_risk),
                              found_random=int(r), found_adaptive=int(a), found_second_stage=int(sec),
                              flagged_sponsors=f, second_stage_shortfall=short))
    result = dict(population_trials=n, sponsors=len(sponsors), budget=budget, discrepancies=float(y.sum()),
                  yield_random=100*float(np.mean(random_found))/budget, yield_risk=100*found_risk/budget,
                  yield_adaptive=100*float(np.mean(adaptive_found))/budget,
                  yield_second_stage=100*float(np.mean(second_found))/second_stage,
                  mean_flagged_sponsors=float(np.mean(flagged_count)),
                  mean_second_stage_shortfall=float(np.mean(shortfalls)),
                  yield_adaptive_seed_se=100*float(np.std(adaptive_found, ddof=1))/budget/AUDIT_SEEDS**.5)
    result["adaptive_to_risk"] = result["yield_adaptive"]/result["yield_risk"]
    result["adaptive_to_random"] = result["yield_adaptive"]/result["yield_random"]
    result["risk_to_random"] = result["yield_risk"]/result["yield_random"]
    return result, rank


def bootstrap_audit(y, sponsor_of, sizes_k, draw, replicate, flagged_pool):
    """y and sponsor_of are in rank order over test trials; draw lists resampled sponsors."""
    counts = np.bincount(draw, minlength=sizes_k)
    weight = counts[sponsor_of]
    position = np.repeat(np.arange(y.size), weight)
    starts = np.concatenate(([0], np.cumsum(weight)[:-1]))
    copy = np.arange(position.size) - np.repeat(starts, weight)
    offset = np.concatenate(([0], np.cumsum(counts)[:-1]))
    label = offset[sponsor_of[position]] + copy
    values = y[position]
    n = values.size
    budget = (n + 5)//10
    first_stage = budget//2
    second_stage = budget-first_stage
    flagged, audited = flagged_pool, np.zeros(n, bool)
    total_random = total_adaptive = 0.0
    for s in range(AUDIT_SEEDS):
        rng = np.random.Generator(np.random.PCG64([SEED_BOOTSTRAP_AUDIT, replicate, s]))
        chosen = rng.choice(n, budget, replace=False)
        r, a, _, _, _ = audit_once(chosen, values, label, first_stage, second_stage, flagged, audited)
        total_random += r
        total_adaptive += a
    risk = 100*float(values[:budget].sum())/budget
    randomized = 100*total_random/AUDIT_SEEDS/budget
    adaptive = 100*total_adaptive/AUDIT_SEEDS/budget
    return dict(budget=budget, yield_random=randomized, yield_risk=risk, yield_adaptive=adaptive,
                adaptive_to_risk=adaptive/risk if risk else float("nan"),
                adaptive_to_random=adaptive/randomized if randomized else float("nan"),
                risk_to_random=risk/randomized if randomized else float("nan"))


def verdict_of(p, auc_difference, auc_lower, ratio, ratio_lower):
    if p >= 0.01:
        return "NO-CLUSTERING"
    if auc_difference >= 0.02 and auc_lower > 0 and ratio >= 1.2 and ratio_lower > 1:
        return "CLUSTERING-SUPPORTED"
    return "WEAK"


def write_csv(path, rows_out):
    fields = list(rows_out[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)


def tests_two_three(data, add, sponsor_table, seed_rows):
    train, test = split_sponsors(data)
    roles = {k: "train" for k in train}
    roles.update({k: "test" for k in test})
    revealed = reveal_trials(data, train+test)
    n_all = len(data["key"])
    role = np.array([roles.get(k, "") for k in data["key"]])
    unrevealed_train = np.flatnonzero((role == "train") & ~revealed)
    unrevealed_test = np.flatnonzero((role == "test") & ~revealed)
    all_test = np.flatnonzero(role == "test")
    every = np.arange(n_all)
    for k in sponsor_table:
        sponsor_table[k]["role"] = roles.get(k, "not_eligible")
    for i in np.flatnonzero(revealed):
        sponsor_table[data["key"][i]]["revealed"] = sponsor_table[data["key"][i]].get("revealed", 0) + 1
    numbering = {k: j for j, k in enumerate(test)}
    owner_unrevealed = np.array([numbering[data["key"][i]] for i in unrevealed_test])
    owner_all = np.array([numbering[data["key"][i]] for i in all_test])
    class_unrevealed = data["cls"][unrevealed_test]
    rnd = random.Random(SEED_BOOTSTRAP)
    k = len(test)
    draws = [np.array([int(rnd.random()*k) for _ in range(k)]) for _ in range(BOOTSTRAPS)]
    coefficient_rows, prediction_columns, boot_rows, results = [], {}, [], {}
    for name in OUTCOMES:
        y_all = data[name]
        revealed_hits, revealed_trials = {}, {}
        for i in np.flatnonzero(revealed):
            key = data["key"][i]
            revealed_hits[key] = revealed_hits.get(key, 0.0) + y_all[i]
            revealed_trials[key] = revealed_trials.get(key, 0) + 1
        history = np.zeros(n_all)
        for i, key in enumerate(data["key"]):
            if key in revealed_trials:
                rate = (revealed_hits[key] + 0.5)/(revealed_trials[key] + 1)
                history[i] = np.log(rate/(1-rate))
        x_train_base, base_names = design(data, unrevealed_train)
        x_train_hist, hist_names = design(data, unrevealed_train, history)
        y_train = y_all[unrevealed_train]
        beta_base, steps_base = fit_logistic(x_train_base, y_train)
        beta_hist, steps_hist = fit_logistic(x_train_hist, y_train)
        for model, labels, beta, steps in [("baseline", base_names, beta_base, steps_base),
                                           ("history", hist_names, beta_hist, steps_hist)]:
            for term, value in zip(labels, beta):
                coefficient_rows.append(dict(outcome=name, model=model, term=term, estimate=float(value), newton_steps=steps))
        x_all_base, _ = design(data, every)
        score_base_all = 1/(1 + np.exp(-(x_all_base @ beta_base)))
        x_test_hist, _ = design(data, unrevealed_test, history)
        y_test = y_all[unrevealed_test]
        score_base = score_base_all[unrevealed_test]
        score_hist = 1/(1 + np.exp(-(x_test_hist @ beta_hist)))
        prediction_columns[name] = (score_base, score_hist)
        groups_base, groups_hist = tied_groups(score_base), tied_groups(score_hist)
        one = np.ones(y_test.size)

        def measures(weights):
            auc_b = weighted_auc(groups_base, y_test, weights)
            auc_h = weighted_auc(groups_hist, y_test, weights)
            brier_b = weighted_brier(score_base, y_test, weights)
            brier_h = weighted_brier(score_hist, y_test, weights)
            return dict(auc_baseline=auc_b, auc_history=auc_h, auc_difference=auc_h-auc_b,
                        brier_baseline=brier_b, brier_history=brier_h, brier_difference=brier_h-brier_b)

        point = {"all": measures(one)}
        masks = {"industry": (class_unrevealed == 0).astype(float), "other": (class_unrevealed == 1).astype(float)}
        for subgroup, mask in masks.items():
            point[subgroup] = measures(mask)
        audit_all, rank = point_audit(data, y_all, all_test, score_base_all, name, "all", add, seed_rows)
        y_ranked = y_all[all_test][rank]
        owner_ranked = owner_all[rank]
        flagged_pool = np.zeros(k, bool)
        started = time.time()
        for b, draw in enumerate(draws):
            counts = np.bincount(draw, minlength=k)
            weight = counts[owner_unrevealed].astype(float)
            row = dict(outcome=name, replicate=b+1, **measures(weight))
            for subgroup, mask in masks.items():
                part = measures(weight*mask)
                row[f"auc_difference_{subgroup}"] = part["auc_difference"]
                row[f"brier_difference_{subgroup}"] = part["brier_difference"]
            row.update(bootstrap_audit(y_ranked, owner_ranked, k, draw, b, flagged_pool))
            boot_rows.append(row)
            if (b+1) % 500 == 0:
                print(f"{name} bootstrap {b+1}/{BOOTSTRAPS}, {time.time()-started:.0f} s", flush=True)
        here = [r for r in boot_rows if r["outcome"] == name]
        for subgroup in ["all", "industry", "other"]:
            add("test2", name, subgroup, "test_trials_scored", int(one.sum() if subgroup == "all" else masks[subgroup].sum()))
            chosen = one if subgroup == "all" else masks[subgroup]
            add("test2", name, subgroup, "test_discrepancies", int((y_test*chosen).sum()))
        add("test2", name, "all", "train_sponsors", len(train))
        add("test2", name, "all", "test_sponsors", len(test))
        add("test2", name, "all", "train_trials_fit", int(unrevealed_train.size))
        add("test2", name, "all", "train_discrepancies_fit", int(y_train.sum()))
        add("test2", name, "all", "revealed_trials", int(revealed.sum()))
        for metric in ["auc_baseline", "auc_history", "auc_difference", "brier_baseline", "brier_history", "brier_difference"]:
            low, high, missing = percentile_limits([r[metric] for r in here])
            add("test2", name, "all", metric, point["all"][metric], low, high, missing)
        for subgroup in ["industry", "other"]:
            add("test2", name, subgroup, "auc_baseline", point[subgroup]["auc_baseline"])
            add("test2", name, subgroup, "auc_history", point[subgroup]["auc_history"])
            add("test2", name, subgroup, "brier_baseline", point[subgroup]["brier_baseline"])
            add("test2", name, subgroup, "brier_history", point[subgroup]["brier_history"])
            for metric in ["auc_difference", "brier_difference"]:
                low, high, missing = percentile_limits([r[f"{metric}_{subgroup}"] for r in here])
                add("test2", name, subgroup, metric, point[subgroup][metric], low, high, missing)
        for metric, value in audit_all.items():
            if metric in ["yield_random", "yield_risk", "yield_adaptive", "adaptive_to_risk", "adaptive_to_random", "risk_to_random"]:
                low, high, missing = percentile_limits([r[metric] for r in here])
                add("test3", name, "all", metric, value, low, high, missing)
            else:
                add("test3", name, "all", metric, value)
        for c, subgroup in enumerate(CLASSES):
            population = all_test[data["cls"][all_test] == c]
            audit_part, _ = point_audit(data, y_all, population, score_base_all, name, subgroup, add, seed_rows)
            for metric, value in audit_part.items():
                add("test3", name, subgroup, metric, value)
        results[name] = dict(auc_difference=point["all"]["auc_difference"],
                             auc_lower=percentile_limits([r["auc_difference"] for r in here])[0],
                             ratio=audit_all["adaptive_to_risk"],
                             ratio_lower=percentile_limits([r["adaptive_to_risk"] for r in here])[0])
    prediction_rows = []
    for j, i in enumerate(unrevealed_test):
        prediction_rows.append(dict(nct_id=data["ids"][i], D=int(data["D"][i]), D2=int(data["D2"][i]),
                                    baseline_D=float(prediction_columns["D"][0][j]), history_D=float(prediction_columns["D"][1][j]),
                                    baseline_D2=float(prediction_columns["D2"][0][j]), history_D2=float(prediction_columns["D2"][1][j])))
    return results, coefficient_rows, prediction_rows, boot_rows


def run_command(split, output):
    if not sys.flags.isolated:
        raise SystemExit("Run with python3 -I because this reads downloaded data.")
    if split not in ["development", "evaluation"]:
        raise SystemExit("Split must be development or evaluation.")
    output = Path(output).resolve()
    if split == "evaluation" and output != EVALUATION_OUTPUT:
        raise SystemExit(f"The evaluation run writes only to {EVALUATION_OUTPUT}.")
    if split == "evaluation" and EVALUATION_OUTPUT.exists():
        raise SystemExit("The evaluation run happens once; its output already exists.")
    inputs = {"plan": (PLAN, PLAN_SHA256), "pool": (POOL, POOL_SHA256), "sponsor_map": (SPONSOR_MAP, MAP_SHA256),
              "script": (Path(__file__).resolve(), None)}
    before = {}
    for name, (path, expected) in inputs.items():
        before[name] = sha256_file(path)
        if expected and before[name] != expected:
            raise SystemExit(f"Hash mismatch for {name}; stopping before any analysis.")
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    print(f"Split {split}. Hashes match.", flush=True)
    data = load_trials(split)
    summary = []

    def add(section, outcome, subgroup, metric, value, lower="", upper="", undefined=""):
        summary.append(dict(section=section, outcome=outcome, subgroup=subgroup, metric=metric, value=value,
                            lower=lower, upper=upper, undefined_replicates=undefined))

    counts = Counter(data["key"])
    add("counts", "both", "all", "trials_in_split", len(data["key"]))
    add("counts", "both", "all", "sponsors_in_split", len(counts))
    add("counts", "both", "all", "sponsors_at_least_3", sum(v >= 3 for v in counts.values()))
    add("counts", "both", "all", "sponsors_at_least_4", sum(v >= 4 for v in counts.values()))
    for name in OUTCOMES:
        add("counts", name, "all", "discrepancies_in_split", int(data[name].sum()))
    first, sponsor_table, null_rows = test_one(data, add)
    print(f"Test 1 done, {time.time()-started:.0f} s", flush=True)
    seed_rows = []
    later, coefficient_rows, prediction_rows, boot_rows = tests_two_three(data, add, sponsor_table, seed_rows)
    print(f"Tests 2 and 3 done, {time.time()-started:.0f} s", flush=True)
    verdict = verdict_of(first["D"]["p"], later["D"]["auc_difference"], later["D"]["auc_lower"],
                         later["D"]["ratio"], later["D"]["ratio_lower"])
    add("verdict", "D", "all", "test1_p_below_0.01", first["D"]["p"] < 0.01)
    add("verdict", "D", "all", "auc_difference_at_least_0.02", later["D"]["auc_difference"] >= 0.02)
    add("verdict", "D", "all", "auc_lower_above_0", later["D"]["auc_lower"] > 0)
    add("verdict", "D", "all", "ratio_at_least_1.2", later["D"]["ratio"] >= 1.2)
    add("verdict", "D", "all", "ratio_lower_above_1", later["D"]["ratio_lower"] > 1)
    add("verdict", "D", "all", "verdict", verdict)
    table_rows = [dict(sponsor_key=k, sponsor_class=v["sponsor_class"], trials=v["trials"], role=v["role"], revealed=v.get("revealed", 0),
                       observed_D=v["observed_D"], expected_D=v["expected_D"], term_D=v["term_D"],
                       observed_D2=v["observed_D2"], expected_D2=v["expected_D2"], term_D2=v["term_D2"])
                  for k, v in sponsor_table.items()]
    write_csv(output/"summary.csv", summary)
    write_csv(output/"sponsors.csv", table_rows)
    write_csv(output/"test1-null.csv", null_rows)
    write_csv(output/"coefficients.csv", coefficient_rows)
    write_csv(output/"predictions.csv", prediction_rows)
    write_csv(output/"test3-seeds.csv", seed_rows)
    write_csv(output/"bootstrap.csv", boot_rows)
    elapsed = time.time()-started
    info = dict(split=split, use="confirmatory" if split == "evaluation" else "engineering only",
                finished_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"), elapsed_seconds=round(elapsed, 1),
                python=platform.python_version(), numpy=np.__version__,
                seeds=dict(permute=SEED_PERMUTE, reveal=SEED_REVEAL, bootstrap=SEED_BOOTSTRAP, audit=SEED_AUDIT, bootstrap_audit=SEED_BOOTSTRAP_AUDIT),
                permutations=PERMUTATIONS, bootstraps=BOOTSTRAPS, audit_seeds=AUDIT_SEEDS, verdict=verdict)
    (output/"run-info.json").write_text(json.dumps(info, indent=1) + "\n")
    hash_rows, unchanged = [], True
    for name, (path, expected) in inputs.items():
        after = sha256_file(path)
        good = after == before[name] and (expected is None or after == expected)
        unchanged = unchanged and good
        hash_rows.append(dict(input=name, path=str(path.relative_to(ROOT)), sha256_before=before[name], sha256_after=after, unchanged=good))
    write_csv(output/"input-hashes.csv", hash_rows)
    print(f"Verdict on D ({info['use']}): {verdict}. Elapsed {elapsed:.0f} s.", flush=True)
    if not unchanged:
        raise SystemExit("An input changed during the run; treat this output as invalid.")


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "fetch":
        fetch_command(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "run":
        run_command(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
