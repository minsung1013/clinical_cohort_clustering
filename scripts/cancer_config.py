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
    "Head & Neck": {
        "slug": "hnc",
        "title": "두경부암(Head & Neck)",
        "biomarkers": {
            "PD-L1_CPS": r"pd[\s-]?l1|\bcps\b|combined positive score",
            "HPV_p16": r"\bhpv\b|human papilloma|\bp16\b",
            "EGFR": r"\begfr\b|epidermal growth factor",
            "MSI-H/dMMR": r"\bmsi[\s-]?h\b|microsatellite instab|\bdmmr\b|mismatch repair defic",
        },
    },
    "Lung": {
        "slug": "lung",
        "title": "폐암(Lung)",
        "biomarkers": {
            "EGFR_mut": r"egfr[\s-]*(mut|exon\s*19|exon\s*20|l858r|del\b|del19|t790m|activating|sensiti)|exon\s*19\s*del|l858r|t790m",
            "ALK": r"\balk\b[\s-]*(positive|rearrang|fusion|translocation|\+)|anaplastic lymphoma kinase",
            "ROS1": r"\bros1\b",
            "KRAS_G12C": r"kras\s*g12c",
            "MET": r"met\s*exon\s*14|metex14|\bmet\b[\s-]*(amplif|exon|skipping)",
            "PD-L1": r"pd[\s-]?l1|\btps\b|tumou?r proportion score",
        },
    },
    "Breast": {"slug": "breast", "title": "유방암(Breast)", "biomarkers": {
        "HR/ER": r"estrogen receptor|hormone receptor|\ber[\s\+/]|\bhr[\s\+]|progesterone receptor",
        "HER2": r"\bher2\b|erbb2",
        "PD-L1": r"pd[\s-]?l1",
        "BRCA/HRD": r"\bbrca\b|\bhrd\b|homologous recombination",
        "PIK3CA": r"pik3ca|\bakt\b|pten"}},
    "Prostate": {"slug": "prostate", "title": "전립선암(Prostate)", "biomarkers": {
        "AR": r"androgen receptor|\bar[\s-]?v7\b",
        "BRCA/HRR": r"\bbrca\b|\bhrr\b|homologous recombination",
        "PSMA": r"\bpsma\b",
        "MSI-H/dMMR": r"\bmsi[\s-]?h\b|\bdmmr\b|mismatch repair defic"}},
    "Renal": {"slug": "renal", "title": "신장암(Renal)", "biomarkers": {
        "PD-L1": r"pd[\s-]?l1",
        "sarcomatoid": r"sarcomatoid"}},
    "Liver": {"slug": "liver", "title": "간암(Liver)", "biomarkers": {
        "AFP": r"\bafp\b|alpha[\s-]?fetoprotein",
        "PD-L1": r"pd[\s-]?l1"}},
    "Pancreatic": {"slug": "panc", "title": "췌장암(Pancreatic)", "biomarkers": {
        "BRCA/HRD": r"\bbrca\b|\bhrd\b|homologous recombination",
        "KRAS": r"\bkras\b",
        "MSI-H/dMMR": r"\bmsi[\s-]?h\b|\bdmmr\b|mismatch repair defic"}},
    "Ovarian": {"slug": "ovarian", "title": "난소암(Ovarian)", "biomarkers": {
        "BRCA/HRD": r"\bbrca\b|\bhrd\b|homologous recombination",
        "platinum_status": r"platinum[\s-]?(sensitive|resistant|refractory)",
        "PD-L1": r"pd[\s-]?l1"}},
    "Melanoma": {"slug": "melanoma", "title": "흑색종(Melanoma)", "biomarkers": {
        "BRAF": r"braf\s*v?600|braf mutat",
        "NRAS": r"\bnras\b",
        "PD-L1": r"pd[\s-]?l1"}},
    "Bladder": {"slug": "bladder", "title": "방광암(Bladder)", "biomarkers": {
        "PD-L1": r"pd[\s-]?l1",
        "FGFR": r"\bfgfr\b",
        "HER2": r"\bher2\b|erbb2",
        "MSI-H/dMMR": r"\bmsi[\s-]?h\b|\bdmmr\b|mismatch repair defic"}},
    "Endometrial": {"slug": "endometrial", "title": "자궁내막암(Endometrial)", "biomarkers": {
        "MSI-H/dMMR": r"\bmsi[\s-]?h\b|microsatellite instab|\bdmmr\b|mismatch repair defic",
        "HER2": r"\bher2\b|erbb2",
        "PD-L1": r"pd[\s-]?l1",
        "POLE": r"\bpole\b"}},
    "Esophageal": {"slug": "esoph", "title": "식도암(Esophageal)", "biomarkers": {
        "HER2": r"\bher2\b|erbb2",
        "PD-L1_CPS": r"pd[\s-]?l1|\bcps\b|combined positive score",
        "CLDN18.2": r"claudin\s*18\.?2|cldn\s*18\.?2|zolbetuximab",
        "MSI-H/dMMR": r"\bmsi[\s-]?h\b|\bdmmr\b|mismatch repair defic"}},
    "Glioma": {"slug": "glioma", "title": "신경교종(Glioma)", "biomarkers": {
        "IDH": r"\bidh1?\b|isocitrate dehydrogenase",
        "MGMT": r"\bmgmt\b",
        "EGFR": r"\begfr\b|egfrviii",
        "1p19q": r"1p/?19q"}},
    "Neuroendocrine": {"slug": "net", "title": "신경내분비종양(NET)", "biomarkers": {
        "SSTR": r"\bsstr\b|somatostatin receptor",
        "Ki-67": r"ki[\s-]?67"}},
    "Sarcoma": {"slug": "sarcoma", "title": "육종(Sarcoma)", "biomarkers": {
        "PD-L1": r"pd[\s-]?l1",
        "translocation": r"translocation|fusion|rearrang"}},
    "Mesothelioma": {"slug": "meso", "title": "중피종(Mesothelioma)", "biomarkers": {
        "PD-L1": r"pd[\s-]?l1",
        "BAP1": r"\bbap1\b"}},
    "Skin": {"slug": "skin", "title": "피부암(Skin)", "biomarkers": {
        "PD-L1": r"pd[\s-]?l1",
        "BRAF": r"braf\s*v?600"}},
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
