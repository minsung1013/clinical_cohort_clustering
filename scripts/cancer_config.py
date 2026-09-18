"""Per-cancer configuration (scope, slug, eligibility biomarker patterns).

The active cancer is chosen by the CANCER environment variable (default Colorectal),
so every stage of the pipeline (build_scope -> enrich -> extract -> notebook ->
dashboard) shares one source of truth.
"""
import os

# line-of-therapy patterns are cancer-agnostic
LINE = {
    "first_line": r"first[\s-]?line|1st[\s-]?line|treatment[\s-]?na[iï]ve|previously untreated",
    "later_line": r"second[\s-]?line|2nd[\s-]?line|third[\s-]?line|refractory|relapsed|previously treated",
}

CONFIGS = {
    "Colorectal": {
        "slug": "crc",
        "title": "대장암(Colorectal)",
        "biomarkers": {
            "RAS/RAF_wt": r"\b(k?ras|nras|braf)\b.{0,30}\bwild[\s-]?type|\bwild[\s-]?type\b.{0,30}\b(k?ras|nras|braf)\b",
            "MSI-H/dMMR": r"\bmsi[\s-]?h\b|microsatellite instab|\bdmmr\b|mismatch repair defic",
            "MSS/pMMR": r"\bmss\b|microsatellite stable|\bpmmr\b|proficient mismatch",
            "HER2": r"\bher2\b|erbb2",
            "BRAF_mut": r"braf\s*v?600|braf mutat",
            "KRAS_G12C": r"kras\s*g12c",
        },
    },
    "Gastric": {
        "slug": "gastric",
        "title": "위암(Gastric)",
        "biomarkers": {
            "HER2": r"\bher2\b|erbb2",
            "CLDN18.2": r"claudin\s*18\.?2|cldn\s*18\.?2|zolbetuximab",
            "MSI-H/dMMR": r"\bmsi[\s-]?h\b|microsatellite instab|\bdmmr\b|mismatch repair defic",
            "PD-L1_CPS": r"pd[\s-]?l1|\bcps\b|combined positive score",
            "EBV": r"epstein[\s-]?barr|\bebv\b",
        },
    },
}


def current() -> dict:
    name = os.environ.get("CANCER", "Colorectal")
    if name not in CONFIGS:
        raise SystemExit(f"Unknown CANCER={name!r}. Options: {list(CONFIGS)}")
    cfg = dict(CONFIGS[name])
    cfg["name"] = name
    cfg["line"] = LINE
    cfg["bm_cols"] = list(cfg["biomarkers"].keys())
    return cfg
