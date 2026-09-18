"""Drug-name normalization, regimen expansion, big-pharma set, TME modality.

Reusable helpers imported by extract_comparator.py and the analysis notebook.
"""
import re

# ── synonym map (lower-case surface -> canonical) ────────────────────────────
SYNONYMS = {
    "5-fu": "fluorouracil",
    "5 fu": "fluorouracil",
    "5-fluorouracil": "fluorouracil",
    "5fu": "fluorouracil",
    "fu": "fluorouracil",
    "leucovorin calcium": "leucovorin",
    "folinic acid": "leucovorin",
    "levoleucovorin": "leucovorin",
    "oxaliplatin injection": "oxaliplatin",
    "xeloda": "capecitabine",
    "avastin": "bevacizumab",
    "erbitux": "cetuximab",
    "vectibix": "panitumumab",
    "keytruda": "pembrolizumab",
    "opdivo": "nivolumab",
    "camptosar": "irinotecan",
    "eloxatin": "oxaliplatin",
    "calcium folinate": "leucovorin",
    "folinate": "leucovorin",
    "trifluridine and tipiracil hydrochloride": "trifluridine/tipiracil",
    "trifluridine/tipiracil hydrochloride": "trifluridine/tipiracil",
    "tas-102": "trifluridine/tipiracil",
    "lonsurf": "trifluridine/tipiracil",
    # gastric
    "ts-1": "s-1", "tegafur": "s-1", "tegafur/gimeracil/oteracil": "s-1",
    "tegafur gimeracil oteracil": "s-1", "tegafur-gimeracil-oteracil": "s-1",
    "tegafur/gimeracil/oteracil potassium": "s-1", "teysuno": "s-1",
    "herceptin": "trastuzumab", "cyramza": "ramucirumab", "enhertu": "trastuzumab deruxtecan",
    "vyloy": "zolbetuximab",
}

# ── established standard-of-care agents (recognized comparators/backbones) ────
# Used to keep clustering "standard-therapy-centered": novel first-in-class
# experimental drugs are excluded from the B axis.
STANDARD_OF_CARE = {
    # chemo backbones
    "fluorouracil", "capecitabine", "oxaliplatin", "irinotecan", "leucovorin",
    "trifluridine/tipiracil", "gemcitabine", "cisplatin", "carboplatin",
    "paclitaxel", "docetaxel", "pemetrexed", "mitomycin",
    "s-1", "epirubicin",                                    # gastric
    # anti-angiogenic
    "bevacizumab", "aflibercept", "ramucirumab", "fruquintinib", "regorafenib",
    # anti-EGFR
    "cetuximab", "panitumumab",
    # immuno-oncology (approved / standard comparators)
    "pembrolizumab", "nivolumab", "ipilimumab", "dostarlimab",
    "atezolizumab", "durvalumab", "tislelizumab", "sintilimab",  # gastric IO
    # HER2
    "trastuzumab", "pertuzumab", "tucatinib", "trastuzumab deruxtecan", "lapatinib",
    # CLDN18.2 (gastric)
    "zolbetuximab",
    # BRAF / MEK
    "encorafenib", "binimetinib",
}

# standard-of-care therapeutic classes (for interpretable, treatment-centered clusters)
SOC_CLASS = {
    "chemo": {
        "fluorouracil", "capecitabine", "oxaliplatin", "irinotecan", "leucovorin",
        "trifluridine/tipiracil", "gemcitabine", "cisplatin", "carboplatin",
        "paclitaxel", "docetaxel", "pemetrexed", "mitomycin", "s-1", "epirubicin",
        # regimen-level tokens (kept atomic, not decomposed)
        "FOLFOX", "FOLFIRI", "FOLFOXIRI", "FOLFIRINOX", "CAPOX", "CAPIRI",
        "FLOT", "SOX", "FP", "XP", "ECF", "ECX",
    },
    "anti_VEGF": {"bevacizumab", "aflibercept", "ramucirumab", "fruquintinib", "regorafenib"},
    "anti_EGFR": {"cetuximab", "panitumumab"},
    "IO": {"pembrolizumab", "nivolumab", "ipilimumab", "dostarlimab", "atezolizumab",
           "durvalumab", "tislelizumab", "sintilimab"},
    "HER2": {"trastuzumab", "pertuzumab", "tucatinib", "trastuzumab deruxtecan", "lapatinib"},
    "CLDN18.2": {"zolbetuximab"},
    "BRAF_MEK": {"encorafenib", "binimetinib"},
    "SoC_generic": {"standard_of_care"},  # comparator named as a category, not a specific drug
}


def soc_classes(drugs) -> set[str]:
    """Map a set/list of canonical drug tokens to their SoC classes."""
    ds = set(drugs)
    return {cls for cls, members in SOC_CLASS.items() if ds & members}


# regimen label -> component drugs (inverse of reconstruction).
# Used to build the CLUSTERING matrix at component level so partial regimens
# (e.g. 5FU+oxaliplatin) sit close to FOLFOX (graded similarity); regimen
# labels themselves are kept only for display/interpretation.
REGIMEN_COMPONENTS = {
    "FOLFOX": {"fluorouracil", "leucovorin", "oxaliplatin"},
    "FOLFIRI": {"fluorouracil", "leucovorin", "irinotecan"},
    "FOLFOXIRI": {"fluorouracil", "leucovorin", "oxaliplatin", "irinotecan"},
    "FOLFIRINOX": {"fluorouracil", "leucovorin", "oxaliplatin", "irinotecan"},
    "CAPOX": {"capecitabine", "oxaliplatin"},
    "CAPIRI": {"capecitabine", "irinotecan"},
    # gastric
    "FLOT": {"fluorouracil", "leucovorin", "oxaliplatin", "docetaxel"},
    "SOX": {"s-1", "oxaliplatin"},
    "FP": {"cisplatin", "fluorouracil"},
    "XP": {"capecitabine", "cisplatin"},
    "ECF": {"epirubicin", "cisplatin", "fluorouracil"},
    "ECX": {"epirubicin", "cisplatin", "capecitabine"},
}


def to_components(tokens) -> set[str]:
    """Expand regimen labels to component drugs; pass other tokens through."""
    out: set[str] = set()
    for t in tokens:
        if t:
            out |= REGIMEN_COMPONENTS.get(t, {t})
    return out


# generic standard-of-care phrases (comparator named as a category, not a drug)
GENERIC_SOC = [
    "standard of care", "standard chemotherapy", "soc therapy", "soc chemotherapy",
    "investigator's choice", "investigators choice", "physician's choice", "physicians choice",
    "treatment of physician", "doublets chemotherapy", "doublet chemotherapy",
]
SOC_GENERIC_TOKEN = "standard_of_care"

# non-drug comparators (not a treatment cohort for TME analysis)
NONDRUG_HINTS = [
    "surgery", "resection", "device", "imaging", "diagnostic", "assay",
    "procedure", "radiation therapy", "observation", "placebo", "best supportive",
]


def soc_from_text(raw: str) -> set[str]:
    """Substring-scan a raw comparator string for established SoC.

    Robust to regimen suffixes ('FOLFIRI regimen'), compound arm strings
    ('bevacizumab + FOLFOXIRI') and generic labels ('standard of care').
    Returns canonical SoC drug tokens (regimens expanded); generic → {SOC_GENERIC_TOKEN}.
    """
    t = str(raw).lower()
    found: set[str] = set()
    work = t
    # regimens as ATOMIC tokens (longest variant first so FOLFOXIRI≠FOLFOX/FOLFIRI)
    for var in sorted(REGIMEN_LABEL, key=len, reverse=True):
        if re.search(r"(?<![a-z])" + re.escape(var) + r"(?![a-z0-9])", work):
            found.add(REGIMEN_LABEL[var])
            work = re.sub(re.escape(var), " ", work)
    # remaining individual SoC drugs (regimen labels excluded)
    for d in STANDARD_OF_CARE:
        if d in work:
            found.add(d)
    for syn, canon_ in SYNONYMS.items():
        if syn in work and canon_ in STANDARD_OF_CARE:
            found.add(canon_)
    if any(g in t for g in GENERIC_SOC):
        found.add(SOC_GENERIC_TOKEN)
    # fold separately-listed chemo components into a named regimen (most specific first)
    for comps, label in _RECON:
        if comps <= found:
            found -= comps
            found.add(label)
    return found


# component sets -> regimen label (larger/more-specific supersets first)
_RECON = [
    ({"fluorouracil", "leucovorin", "oxaliplatin", "docetaxel"}, "FLOT"),
    ({"fluorouracil", "leucovorin", "oxaliplatin", "irinotecan"}, "FOLFOXIRI"),
    ({"epirubicin", "cisplatin", "fluorouracil"}, "ECF"),
    ({"epirubicin", "cisplatin", "capecitabine"}, "ECX"),
    ({"fluorouracil", "leucovorin", "oxaliplatin"}, "FOLFOX"),
    ({"fluorouracil", "leucovorin", "irinotecan"}, "FOLFIRI"),
    ({"s-1", "oxaliplatin"}, "SOX"),
    ({"cisplatin", "fluorouracil"}, "FP"),
    ({"capecitabine", "cisplatin"}, "XP"),
    ({"capecitabine", "oxaliplatin"}, "CAPOX"),
    ({"capecitabine", "irinotecan"}, "CAPIRI"),
]

# ── multi-drug regimens -> component set (canonical) ─────────────────────────
REGIMENS = {
    "folfox": {"fluorouracil", "leucovorin", "oxaliplatin"},
    "mfolfox6": {"fluorouracil", "leucovorin", "oxaliplatin"},
    "folfox6": {"fluorouracil", "leucovorin", "oxaliplatin"},
    "folfox4": {"fluorouracil", "leucovorin", "oxaliplatin"},
    "folfiri": {"fluorouracil", "leucovorin", "irinotecan"},
    "folfoxiri": {"fluorouracil", "leucovorin", "oxaliplatin", "irinotecan"},
    "folfirinox": {"fluorouracil", "leucovorin", "oxaliplatin", "irinotecan"},
    "capox": {"capecitabine", "oxaliplatin"},
    "xelox": {"capecitabine", "oxaliplatin"},
    "capiri": {"capecitabine", "irinotecan"},
    "xeliri": {"capecitabine", "irinotecan"},
    "capeox": {"capecitabine", "oxaliplatin"},
    # gastric
    "flot": {"fluorouracil", "leucovorin", "oxaliplatin", "docetaxel"},
    "sox": {"s-1", "oxaliplatin"},
    "fp": {"cisplatin", "fluorouracil"},
    "cf": {"cisplatin", "fluorouracil"},
    "xp": {"capecitabine", "cisplatin"},
    "ecf": {"epirubicin", "cisplatin", "fluorouracil"},
    "ecx": {"epirubicin", "cisplatin", "capecitabine"},
}

# regimen name variants -> canonical atomic label (kept as a single feature)
REGIMEN_LABEL = {
    "folfoxiri": "FOLFOXIRI", "folfirinox": "FOLFIRINOX", "mfolfirinox": "FOLFIRINOX",
    "mfolfox6": "FOLFOX", "mfolfox": "FOLFOX", "folfox6": "FOLFOX", "folfox4": "FOLFOX", "folfox": "FOLFOX",
    "mfolfiri": "FOLFIRI", "folfiri": "FOLFIRI",
    "capox": "CAPOX", "capeox": "CAPOX", "xelox": "CAPOX",
    "capiri": "CAPIRI", "xeliri": "CAPIRI",
    # gastric
    "flot": "FLOT", "sox": "SOX", "fp": "FP", "cf": "FP", "xp": "XP", "ecf": "ECF", "ecx": "ECX",
}

# supportive / non-treatment agents — flagged, excluded from "standard therapy B"
SUPPORTIVE = {
    "placebo", "dexamethasone", "prednisone", "prednisolone", "methylprednisolone",
    "saline", "normal saline", "vitamin", "folic acid", "vitamin b12",
    "ondansetron", "loperamide", "atropine", "best supportive care", "bsc",
}

# ── big pharma (canonical company_normalized surfaces) ───────────────────────
BIG_PHARMA = {
    "roche", "genentech", "merck", "merck sharp & dohme", "msd", "pfizer",
    "astrazeneca", "bristol myers squibb", "bristol-myers squibb", "novartis",
    "johnson & johnson", "janssen", "sanofi", "gsk", "glaxosmithkline", "abbvie",
    "amgen", "gilead", "gilead sciences", "takeda", "daiichi sankyo",
    "eli lilly", "lilly", "bayer", "astellas", "emd serono", "merck kgaa",
    "beigene", "regeneron", "boehringer ingelheim", "novo nordisk",
    "genmab", "servier",
}

# TME-dependent modalities (spatial/TME analysis adds most value)
TME_MODALITIES = {
    "ADC", "Bispecific Antibody", "CAR-T", "Cell Therapy",
    "Monoclonal Antibody", "Oncolytic Virus", "Vaccine",
}

_STRIP = re.compile(r"^(drug|biological|combination product|other|genetic|radiation|procedure)\s*:\s*", re.I)
_DOSE = re.compile(r"\b\d+\s*(mg|mcg|g|ml|%|iu)\b.*$", re.I)


def clean_name(raw: str) -> str:
    """Lower-cased canonical drug token from a raw intervention string."""
    s = _STRIP.sub("", str(raw)).strip()
    s = _DOSE.sub("", s).strip()
    s = re.sub(r"[®™]", "", s)
    s = re.sub(r"\s+", " ", s).strip().lower().strip(" -,")
    return SYNONYMS.get(s, s)


def expand(name: str) -> set[str]:
    """Return the component drug set for a name (regimen -> components)."""
    c = clean_name(name)
    if c in REGIMENS:
        return set(REGIMENS[c])
    return {c} if c else set()


def canon(name: str) -> str:
    """Canonical ATOMIC token: regimen kept whole (mFOLFOX6 -> FOLFOX), else cleaned drug."""
    s = clean_name(name)
    s = re.sub(r"\s*(regimen|chemotherapy)\s*$", "", s).strip()
    if s in REGIMEN_LABEL:
        return REGIMEN_LABEL[s]
    return SYNONYMS.get(s, s)


def is_supportive(name: str) -> bool:
    return clean_name(name) in SUPPORTIVE


def is_big_pharma(company: str) -> bool:
    if not company:
        return False
    c = company.lower().strip()
    return any(bp in c for bp in BIG_PHARMA)
