"""Stage 4 — Extract standard-of-care comparator B + engineer cohort features.

Input : data/cohorts_crc.csv (from enrich_ctgov.py) + data/pipeline.json (join)
Output: data/cohorts_crc_features.csv (one row per trial, analysis-ready)
"""
import json
import re
from pathlib import Path

import pandas as pd

from drug_normalize import (
    clean_name, canon, is_supportive, is_big_pharma, TME_MODALITIES,
)
from cancer_config import current

CFG = current()
ROOT = Path(__file__).resolve().parent.parent
IN_CSV = ROOT / "data" / f"cohorts_{CFG['slug']}.csv"
PIPELINE = ROOT / "data" / "pipeline.json"
OUT = ROOT / "data" / f"cohorts_{CFG['slug']}_features.csv"


def parse_arms(arms_detail: str):
    """'label ::: type ::: iv1,iv2 || ...' -> [(label, type, [ivs]), ...]."""
    arms = []
    if not isinstance(arms_detail, str) or not arms_detail.strip():
        return arms
    for chunk in arms_detail.split(" || "):
        parts = chunk.split(" ::: ")
        if len(parts) != 3:
            continue
        label, atype, ivs = parts
        drugs = [d for d in ivs.split(",") if d.strip()]
        arms.append((label.strip(), atype.strip(), drugs))
    return arms


def extract_B(arms):
    """Return (comparator_B_set, arm_structure).

    B = the standard-of-care regimen a sponsor wants TME data on.
    """
    active = [a for a in arms if a[1] == "ACTIVE_COMPARATOR"]
    placebo = [a for a in arms if a[1] == "PLACEBO_COMPARATOR"]
    exp = [a for a in arms if a[1] == "EXPERIMENTAL"]

    def drugset(armlist):
        s = set()
        for _, _, ivs in armlist:
            for iv in ivs:
                s |= {d for d in [canon(iv)] if d and not is_supportive(d)}
        return s

    if active:
        return drugset(active), "active_comparator"
    if placebo:
        # backbone shared with the placebo arm = non-placebo drugs given to all
        shared = drugset(placebo)  # placebo arm often carries the backbone SoC
        return shared, "placebo_comparator"
    if len(arms) <= 1:
        # single arm: backbone = combo partners in the (only) arm
        return drugset(arms), "single_arm"
    # multi experimental arms, no formal comparator: shared backbone = intersection
    per_arm = [drugset([a]) for a in exp or arms]
    inter = set.intersection(*per_arm) if per_arm else set()
    return inter, "multi_experimental"


# eligibility biomarker / line patterns come from the active cancer config
BIOMARKERS = CFG["biomarkers"]
LINE = CFG["line"]


def flags(text: str, patterns: dict) -> dict:
    t = (text or "").lower()
    return {k: int(bool(re.search(p, t))) for k, p in patterns.items()}


# ── biomarker positivity labeling (HER2 -> HER2+, EGFR_mut -> EGFRm, …) ───────
# (base_label, mode): "posneg" = infer +/- from eligibility context; "fixed:X" = use X
BM_LABEL = {
    "HER2": ("HER2", "posneg"), "PD-L1": ("PD-L1", "posneg"), "PD-L1_CPS": ("PD-L1", "posneg"),
    "HR/ER": ("HR", "posneg"), "CLDN18.2": ("CLDN18.2", "posneg"), "AR": ("AR", "posneg"),
    "HPV_p16": ("HPV/p16", "posneg"), "PSMA": ("PSMA", "posneg"), "EGFR": ("EGFR", "posneg"),
    "FGFR": ("FGFR", "fixed:FGFR-alt"), "SSTR": ("SSTR", "fixed:SSTR+"),
    "ALK": ("ALK", "fixed:ALK+"), "ROS1": ("ROS1", "fixed:ROS1+"),
    "RET": ("RET", "fixed:RET+"), "MET": ("MET", "fixed:MET+"),
    "EGFR_mut": ("EGFR", "fixed:EGFRm"), "BRAF_mut": ("BRAF", "fixed:BRAFm"),
    "BRAF": ("BRAF", "fixed:BRAFm"), "KRAS_G12C": ("KRAS", "fixed:KRAS G12C"),
    "KRAS": ("KRAS", "fixed:KRASm"), "NRAS": ("NRAS", "fixed:NRASm"),
    "PIK3CA": ("PIK3CA", "fixed:PIK3CAm"), "IDH": ("IDH", "fixed:IDHm"),
    "RAS/RAF_wt": ("RAS/RAF", "fixed:RAS/RAF wt"), "MSS/pMMR": ("MSS", "fixed:MSS/pMMR"),
    "MSI-H/dMMR": ("MSI", "fixed:MSI-H/dMMR"),
    "BRCA/HRD": ("BRCA/HRD", "fixed:BRCA/HRD+"), "BRCA/HRR": ("BRCA/HRR", "fixed:BRCA/HRR+"),
    "MGMT": ("MGMT", "fixed:MGMT"), "1p19q": ("1p19q", "fixed:1p19q-codel"),
    "EBV": ("EBV", "fixed:EBV+"), "AFP": ("AFP", "fixed:AFP"), "POLE": ("POLE", "fixed:POLEm"),
    "Ki-67": ("Ki-67", "fixed:Ki-67"), "translocation": ("fusion", "fixed:fusion+"),
    "platinum_status": ("platinum", "fixed:platinum-status"),
    "sarcomatoid": ("sarcomatoid", "fixed:sarcomatoid"), "BAP1": ("BAP1", "fixed:BAP1"),
}
_NEG = re.compile(r"negative|\bneg\b|(?<![a-z])low(?![a-z])|absent|loss of|non[- ]?amplif|proficient|wild[\s-]?type", re.I)
_POS = re.compile(r"positive|\bpos\b|high|over[- ]?express|amplif|expressing|\bcps\b|\btps\b|selected|\+", re.I)


def polarize(key: str, text: str) -> str:
    base, mode = BM_LABEL.get(key, (key, "posneg"))
    if mode.startswith("fixed:"):
        return mode[6:]
    pat = BIOMARKERS.get(key)
    neg = pos = False
    if pat:
        for m in re.finditer(pat, text, re.I):
            win = text[max(0, m.start() - 25): m.end() + 35]
            if _NEG.search(win):
                neg = True
            if _POS.search(win):
                pos = True
    return base + ("−" if (neg and not pos) else "+")


# ── rule-based eligibility one-line summary ──────────────────────────────────
_SETTING = [(r"metasta", "metastatic"), (r"unresectable", "unresectable"),
            (r"locally advanced", "locally advanced"), (r"\brecurrent\b", "recurrent"),
            (r"stage\s*iv|stage\s*4", "stage IV"), (r"\badvanced\b", "advanced"),
            (r"neoadjuvant", "neoadjuvant"), (r"adjuvant", "adjuvant"), (r"resectable", "resectable")]
_HISTO = [(r"adenocarcinoma", "adenoca."), (r"squamous", "squamous"),
          (r"small[\s-]?cell", "small-cell"), (r"ductal", "ductal")]


def summarize_elig(text: str) -> str:
    t = (text or "").lower()
    if not t.strip():
        return ""
    f = []
    for pat, lab in _SETTING:
        if re.search(pat, t):
            f.append(lab); break
    if re.search(r"first[\s-]?line|1st[\s-]?line|treatment[\s-]?na[iï]ve|previously untreated", t):
        f.append("1st-line")
    elif re.search(r"second[\s-]?line|2nd[\s-]?line|third[\s-]?line|refractory|relapsed|previously treated", t):
        f.append("2nd-line+")
    for pat, lab in _HISTO:
        if re.search(pat, t):
            f.append(lab); break
    if re.search(r"measurable disease|recist", t):
        f.append("measurable (RECIST)")
    m = re.search(r"ecog[^0-9]{0,12}(0\s*(?:to|-|–|,)?\s*1|0|1|2)", t)
    if m:
        f.append("ECOG " + re.sub(r"\s*(to|,)\s*", "-", m.group(1)))
    return " · ".join(dict.fromkeys(f))


def build_nct_pipeline_map():
    """From pipeline.json build:
    - nct_map: nct_id -> {targets, modalities, biomarkers, company}
    - scale:   company_normalized -> distinct oncology-trial count (R&D-scale proxy)
    """
    data = json.loads(PIPELINE.read_text())
    nct_map = {}
    company_ncts: dict[str, set] = {}
    for d in data["drugs"]:
        if not d.get("is_oncology"):
            continue
        comp = d.get("company_normalized") or d.get("company") or ""
        for nct in d.get("nct_ids") or []:
            e = nct_map.setdefault(
                nct, {"targets": set(), "modalities": set(), "biomarkers": set(), "company": comp}
            )
            if comp and not e["company"]:
                e["company"] = comp
            if d.get("target") and d["target"] != "Unknown":
                e["targets"].add(d["target"])
            if d.get("modality") and d["modality"] != "Unknown":
                e["modalities"].add(d["modality"])
            for b in d.get("biomarker_list") or []:
                e["biomarkers"].add(b)
            if comp:
                company_ncts.setdefault(comp, set()).add(nct)
    scale = {c: len(v) for c, v in company_ncts.items()}
    return nct_map, scale


def main() -> None:
    df = pd.read_csv(IN_CSV)
    pipe, scale = build_nct_pipeline_map()

    recs = []
    for _, r in df.iterrows():
        arms = parse_arms(r.get("arms_detail"))
        B, structure = extract_B(arms)
        # experimental (drug A) set
        A = set()
        for _, atype, ivs in arms:
            if atype == "EXPERIMENTAL":
                for iv in ivs:
                    A |= {d for d in [canon(iv)] if d and not is_supportive(d)}
        pm = pipe.get(r["nct_id"], {})
        modalities = pm.get("modalities", set())
        company = pm.get("company") or r.get("lead_sponsor") or ""
        rec = {
            "nct_id": r["nct_id"],
            "brief_title": r["brief_title"],
            "phase": r["phase"],
            "phase_simple": simplify_phase(r["phase"]),
            "overall_status": r["overall_status"],
            "enrollment": r["enrollment"],
            "lead_sponsor": r["lead_sponsor"],
            "company": company,
            "sponsor_scale": scale.get(company, 1),  # # oncology trials by this sponsor (R&D-scale proxy)
            "is_industry": int(r.get("lead_sponsor_class") == "INDUSTRY"),
            "is_big_pharma": int(is_big_pharma(str(r.get("lead_sponsor", "")))
                                 or any(is_big_pharma(c) for c in str(r.get("collaborators", "")).split(";"))),
            "arm_structure": structure,
            "n_arms": r["n_arms"],
            "comparator_B": ";".join(sorted(B)),
            "experimental_A": ";".join(sorted(A)),
            "targets": ";".join(sorted(pm.get("targets", set()))),
            "modalities": ";".join(sorted(modalities)),
            "biomarkers_pipeline": ";".join(sorted(pm.get("biomarkers", set()))),
            "tme_relevant": int(bool(modalities & TME_MODALITIES)),
            "is_combination": int(len(A | B) > 1 or len(B) >= 1 and len(A) >= 1),
        }
        rec.update(flags(r.get("eligibility"), BIOMARKERS))
        rec.update(flags(r.get("eligibility"), LINE))
        rec["biomarker_selected"] = int(any(rec[k] for k in BIOMARKERS))
        elig = str(r.get("eligibility") or "")
        rec["biomarker_detail"] = ", ".join(polarize(k, elig) for k in BIOMARKERS if rec[k] == 1)
        rec["elig_summary"] = summarize_elig(elig)
        recs.append(rec)

    out = pd.DataFrame(recs)
    out.to_csv(OUT, index=False)
    print(f"wrote {len(out)} trials -> {OUT.relative_to(ROOT)}")
    print("\narm_structure:\n", out.arm_structure.value_counts().to_string())
    print("\ntop comparator_B tokens:")
    from collections import Counter
    cb = Counter(x for s in out.comparator_B.dropna() for x in s.split(";") if x)
    for k, v in cb.most_common(15):
        print(f"  {v:4d}  {k}")


def simplify_phase(p: str) -> str:
    p = str(p)
    if "PHASE3" in p:
        return "Phase 3"
    if "PHASE2" in p:
        return "Phase 2"
    if "PHASE1" in p or "EARLY_PHASE1" in p:
        return "Phase 1"
    return "NA"


if __name__ == "__main__":
    main()
