"""Stage 2 — Collect full study records from ClinicalTrials.gov API v2.

Reads the CRC scope NCT list, batch-fetches full study records (all modules),
caches each raw study to data/raw_ctgov/{nct}.json (idempotent), then flattens
the fields we need into data/cohorts_crc.csv (one row per trial).
"""
import json
import time
from pathlib import Path

import requests

from cancer_config import current

CFG = current()
ROOT = Path(__file__).resolve().parent.parent
SCOPE = ROOT / "data" / f"scope_{CFG['slug']}.json"
RAW = ROOT / "data" / "raw_ctgov"                     # shared cache across cancers
OUT_CSV = ROOT / "data" / f"cohorts_{CFG['slug']}.csv"

API = "https://clinicaltrials.gov/api/v2/studies"
BATCH = 50


def fetch_missing(ncts: list[str]) -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    missing = [n for n in ncts if not (RAW / f"{n}.json").exists()]
    print(f"cached: {len(ncts) - len(missing)} | to fetch: {len(missing)}")
    for i in range(0, len(missing), BATCH):
        chunk = missing[i : i + BATCH]
        r = requests.get(
            API,
            params={"filter.ids": ",".join(chunk), "pageSize": BATCH},
            headers={"accept": "application/json"},
            timeout=60,
        )
        r.raise_for_status()
        got = r.json().get("studies", [])
        for s in got:
            nct = s["protocolSection"]["identificationModule"]["nctId"]
            (RAW / f"{nct}.json").write_text(json.dumps(s))
        print(f"  batch {i//BATCH + 1}: requested {len(chunk)}, got {len(got)}")
        time.sleep(0.4)  # polite


def flatten(nct: str) -> dict | None:
    p = RAW / f"{nct}.json"
    if not p.exists():
        return None
    ps = json.loads(p.read_text())["protocolSection"]
    idm = ps.get("identificationModule", {})
    st = ps.get("statusModule", {})
    des = ps.get("designModule", {})
    spon = ps.get("sponsorCollaboratorsModule", {})
    arms = ps.get("armsInterventionsModule", {})
    elig = ps.get("eligibilityModule", {})
    cond = ps.get("conditionsModule", {})
    desc = ps.get("descriptionModule", {})
    out = ps.get("outcomesModule", {})

    arm_groups = arms.get("armGroups", []) or []
    interventions = arms.get("interventions", []) or []
    lead = spon.get("leadSponsor", {})

    return {
        "nct_id": idm.get("nctId"),
        "brief_title": idm.get("briefTitle"),
        "official_title": idm.get("officialTitle"),
        "phase": ";".join(des.get("phases", []) or []),
        "overall_status": st.get("overallStatus"),
        "start_date": (st.get("startDateStruct") or {}).get("date"),
        "primary_completion_date": (st.get("primaryCompletionDateStruct") or {}).get("date"),
        "enrollment": (des.get("enrollmentInfo") or {}).get("count"),
        "enrollment_type": (des.get("enrollmentInfo") or {}).get("type"),
        "allocation": (des.get("designInfo") or {}).get("allocation"),
        "lead_sponsor": lead.get("name"),
        "lead_sponsor_class": lead.get("class"),
        "collaborators": ";".join(c.get("name", "") for c in spon.get("collaborators", []) or []),
        "conditions": ";".join(cond.get("conditions", []) or []),
        "n_arms": len(arm_groups),
        "arm_types": ";".join(a.get("type", "") or "" for a in arm_groups),
        # arm label ::: type ::: interventions | next arm
        "arms_detail": " || ".join(
            f"{a.get('label','')} ::: {a.get('type','')} ::: "
            + ",".join(a.get("interventionNames", []) or [])
            for a in arm_groups
        ),
        "interventions": ";".join(
            f"{iv.get('type','')}:{iv.get('name','')}" for iv in interventions
        ),
        "eligibility": (elig.get("eligibilityCriteria") or "").replace("\n", " ").replace("\r", " "),
        "primary_outcomes": ";".join(o.get("measure", "") for o in out.get("primaryOutcomes", []) or []),
        "brief_summary": (desc.get("briefSummary") or "").replace("\n", " ").replace("\r", " "),
    }


def main() -> None:
    import csv

    ncts = json.loads(SCOPE.read_text())["ncts"]
    fetch_missing(ncts)

    rows = [r for n in ncts if (r := flatten(n))]
    if not rows:
        print("no rows flattened")
        return
    cols = list(rows[0].keys())
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"\nflattened {len(rows)} trials -> {OUT_CSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
