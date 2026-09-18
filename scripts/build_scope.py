"""Stage 1 — Screen pipeline.json to the active cancer scope (env CANCER).

Filters: <CANCER> + ongoing + Phase 1/2/3, then emits the distinct NCT list
that Stage 2 (enrich_ctgov.py) will fully fetch from ClinicalTrials.gov.
"""
import json
from pathlib import Path

from cancer_config import current

CFG = current()
CANCER = CFG["name"]
ROOT = Path(__file__).resolve().parent.parent
PIPELINE = ROOT / "data" / "pipeline.json"
OUT = ROOT / "data" / f"scope_{CFG['slug']}.json"

ONGOING = {
    "RECRUITING",
    "ACTIVE_NOT_RECRUITING",
    "NOT_YET_RECRUITING",
    "ENROLLING_BY_INVITATION",
}
PHASES = {
    "EARLY_PHASE1",
    "PHASE1",
    "PHASE1/PHASE2",
    "PHASE2",
    "PHASE2/PHASE3",
    "PHASE3",
}


def main() -> None:
    data = json.loads(PIPELINE.read_text())
    drugs = [d for d in data["drugs"] if d.get("is_oncology")]

    # A drug record can span several cancer categories / NCTs. Keep it if it
    # touches Colorectal and is ongoing + early-to-pivotal phase.
    scope = [
        d
        for d in drugs
        if CANCER in (d.get("cancer_categories") or [d.get("cancer_category")])
        and d.get("overall_status") in ONGOING
        and d.get("phase") in PHASES
    ]

    ncts = sorted({n for d in scope for n in (d.get("nct_ids") or [])})

    OUT.write_text(
        json.dumps(
            {"cancer": CANCER, "n_drug_records": len(scope), "ncts": ncts},
            indent=2,
        )
    )
    print(f"drug records in {CANCER} scope : {len(scope)}")
    print(f"distinct NCT trials       : {len(ncts)}")
    print(f"wrote -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
