"""Assemble a static site/ folder: a navigation landing page + per-cancer dashboards.

Scans outputs/ for cohort_map_<slug>.html (dashboards) and copies them, their
report and thumbnail into site/, then generates site/index.html (nav page).
Deploy: GitHub Pages (main branch /docs) or drag docs/ to https://app.netlify.com/drop.
"""
import shutil
from datetime import date
from pathlib import Path

from cancer_config import CONFIGS

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
SITE = ROOT / "docs"  # GitHub Pages serves main branch /docs

# slug -> display title
SLUG_TITLE = {c["slug"]: c["title"] for c in CONFIGS.values()}


def cards():
    SITE.mkdir(exist_ok=True)
    items = []
    for dash in sorted(OUT.glob("cohort_map_*.html")):
        slug = dash.stem.replace("cohort_map_", "")
        title = SLUG_TITLE.get(slug, slug)
        shutil.copy(dash, SITE / f"{slug}.html")
        report = OUT / f"report_{slug}.html"
        thumb = OUT / f"cohort_map_{slug}.png"
        rep_link = ""
        if report.exists():
            shutil.copy(report, SITE / f"report_{slug}.html")
            rep_link = f'<a class="rep" href="report_{slug}.html">📄 분석 리포트</a>'
        img = ""
        if thumb.exists():
            shutil.copy(thumb, SITE / f"{slug}.png")
            img = f'<img src="{slug}.png" alt="{title} map">'
        items.append(f"""
      <a class="card" href="{slug}.html">
        <div class="thumb">{img}</div>
        <div class="body">
          <h2>{title}</h2>
          <p>진행 중 임상 코호트 · 표준치료(B) 데이터 자산 · 클러스터 랜드스케이프</p>
          <div class="links"><span class="open">대시보드 열기 →</span>{rep_link}</div>
        </div>
      </a>""")
    return "\n".join(items)


def main():
    body = cards()
    html = _PAGE.replace("__CARDS__", body).replace("__DATE__", date.today().isoformat())
    (SITE / "index.html").write_text(html, encoding="utf-8")
    print(f"wrote {SITE/'index.html'} with dashboards:",
          ", ".join(sorted(p.stem for p in SITE.glob('*.html') if p.stem not in ('index',) and not p.stem.startswith('report'))))


_PAGE = r"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>임상 코호트 인텔리전스</title>
<style>
  *{box-sizing:border-box}
  body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Noto Sans KR",sans-serif;
       color:#0f172a;background:#f1f5f9}
  header{background:#0f172a;color:#fff;padding:26px 22px;text-align:center}
  header h1{margin:0;font-size:22px}
  header p{margin:8px auto 0;color:#94a3b8;font-size:13px;max-width:680px;line-height:1.6}
  .container{max-width:900px;margin:0 auto;padding:26px 22px 40px}
  .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px}
  .card{display:flex;flex-direction:column;background:#fff;border:1px solid #e5e7eb;border-radius:14px;
        overflow:hidden;text-decoration:none;color:inherit;transition:.15s;box-shadow:0 1px 2px rgba(0,0,0,.04)}
  .card:hover{transform:translateY(-3px);box-shadow:0 8px 24px rgba(2,6,23,.12);border-color:#cbd5e1}
  .thumb{aspect-ratio:16/10;background:#f8fafc;overflow:hidden;border-bottom:1px solid #eef2f7}
  .thumb img{width:100%;height:100%;object-fit:cover;object-position:center;display:block}
  .body{padding:15px 17px;display:flex;flex-direction:column;flex:1}
  .body h2{margin:0 0 6px;font-size:17px}
  .body p{margin:0 0 13px;color:#64748b;font-size:12.5px;line-height:1.5;flex:1}
  .links{display:flex;align-items:center;justify-content:space-between;gap:12px}
  .open{color:#2563eb;font-weight:600;font-size:13px}
  .rep{color:#475569;font-size:12px;text-decoration:none}.rep:hover{text-decoration:underline}
  footer{text-align:center;color:#94a3b8;font-size:12px;padding:0 22px 34px}
</style></head>
<body>
<header>
  <h1>임상 코호트 인텔리전스</h1>
  <p>진행 중인 항암 임상 코호트를 분석해, 공간전사체 + AI 종양미세환경(TME) 분석 관점에서
     표준치료(B) 데이터 자산의 가치를 우선순위화합니다. 암종을 선택하세요.</p>
</header>
<div class="container">
  <div class="grid">__CARDS__</div>
</div>
<footer>업데이트 __DATE__ · 데이터 출처: ClinicalTrials.gov</footer>
</body></html>"""


if __name__ == "__main__":
    main()
