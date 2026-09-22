"""Build data/nct_links.json — map each dashboard NCT to linked papers/abstracts.

Source: the oncology_pipeline_dashboard R2 corpus (publications + AACR/ASCO/ESMO
abstracts). Scans each record's nct_ids, keeps only NCTs shown in our dashboards,
and links: publications -> PubMed (pmid), conference abstracts -> DOI (source.doi).

Usage:  python3 scripts/build_nct_links.py
Re-run whenever the oncology corpus updates or new cohorts are added.
"""
import ast
import json
import re
import urllib.request
from pathlib import Path

DATA_BASE = "https://pub-38ddc55a3aa34cf0988d355d9a0abe74.r2.dev"
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "nct_links.json"
NCT_RE = re.compile(r"NCT\d{8}")


def get_json(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.loads(r.read())


def as_list(v):
    if not v:
        return []
    if isinstance(v, list):
        return v
    s = str(v)
    if s.startswith("["):
        try:
            return ast.literal_eval(s)
        except Exception:
            pass
    return NCT_RE.findall(s)


def source_doi(rec):
    src = rec.get("source")
    if isinstance(src, str) and src.startswith("{"):
        try:
            src = ast.literal_eval(src)
        except Exception:
            src = {}
    return src.get("doi") if isinstance(src, dict) else None


def dashboard_ncts():
    ncts = set()
    for f in ROOT.glob("outputs/cohort_map_*.html"):
        h = f.read_text(encoding="utf-8")
        line = next(l for l in h.splitlines() if l.startswith("const TRIALS="))
        trials = json.loads(line[len("const TRIALS="):-1].split(", CLUSTERS=", 1)[0])
        ncts.update(t["nct"] for t in trials)
    return ncts


def main():
    keep = dashboard_ncts()
    print(f"dashboard NCTs: {len(keep)}")
    links = {}

    def add(nct, d):
        links.setdefault(nct, []).append(d)

    # publications -> PubMed
    npub = 0
    for p in get_json(f"{DATA_BASE}/pub_index.json")["publications"]:
        recs = get_json(f"{DATA_BASE}/{p['file']}").get("abstracts", [])
        for r in recs:
            hits = [n for n in as_list(r.get("nct_ids")) if n in keep]
            if not hits:
                continue
            pmid = r.get("pmid")
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None
            for n in hits:
                add(n, {"type": "pub", "title": (r.get("title") or "").strip(),
                        "venue": r.get("journal") or "", "year": str(r.get("year") or ""), "url": url})
                npub += 1

    # conference abstracts -> DOI
    nabs = 0
    for a in get_json(f"{DATA_BASE}/index.json")["abstracts"]:
        recs = get_json(f"{DATA_BASE}/{a['file']}").get("abstracts", [])
        for r in recs:
            hits = [n for n in as_list(r.get("nct_ids")) if n in keep]
            if not hits:
                continue
            doi = source_doi(r)
            url = f"https://doi.org/{doi}" if doi else None
            venue = f"{r.get('conference','')} {r.get('year','')}".strip()
            for n in hits:
                add(n, {"type": "abs", "title": (r.get("title") or "").strip(),
                        "venue": venue, "year": str(r.get("year") or ""), "url": url})
                nabs += 1

    # dedupe + sort (papers first, newest first)
    for n, ls in links.items():
        seen, uniq = set(), []
        for d in ls:
            k = (d["type"], d["title"][:60], d["url"])
            if k not in seen:
                seen.add(k)
                uniq.append(d)
        uniq.sort(key=lambda d: (d["type"] != "pub", -(int(d["year"]) if d["year"].isdigit() else 0)))
        links[n] = uniq

    OUT.write_text(json.dumps(links, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"publications: {npub} | abstracts: {nabs} | NCTs with links: {len(links)}/{len(keep)} "
          f"({100*len(links)/max(1,len(keep)):.1f}%)")
    print(f"wrote {OUT} ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
