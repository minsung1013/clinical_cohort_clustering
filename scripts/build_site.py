"""Assemble a static docs/ folder: a bilingual navigation page + per-cancer dashboards.

Scans outputs/ for cohort_map_<slug>.html (dashboards), copies them + report +
thumbnail into docs/, then generates docs/index.html (KO/EN nav with a language
toggle). Deploy: GitHub Pages (main branch /docs) or drag docs/ to Netlify Drop.
"""
import shutil
from datetime import date
from pathlib import Path

from cancer_config import CONFIGS

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
SITE = ROOT / "docs"  # GitHub Pages serves main branch /docs

# slug -> (korean, english) display names
SLUG_NAMES = {c["slug"]: (c["title"].split("(")[0].strip(), name) for name, c in CONFIGS.items()}


def cards():
    SITE.mkdir(exist_ok=True)
    items = []
    for dash in sorted(OUT.glob("cohort_map_*.html")):
        slug = dash.stem.replace("cohort_map_", "")
        ko, en = SLUG_NAMES.get(slug, (slug, slug))
        shutil.copy(dash, SITE / f"{slug}.html")
        report = OUT / f"report_{slug}.html"
        thumb = OUT / f"cohort_map_{slug}.png"
        rep_link = ""
        if report.exists():
            shutil.copy(report, SITE / f"report_{slug}.html")
            rep_link = f'<a class="rep" href="report_{slug}.html"></a>'
        img = ""
        if thumb.exists():
            shutil.copy(thumb, SITE / f"{slug}.png")
            img = f'<img src="{slug}.png" alt="{en} map">'
        items.append(f"""
      <a class="card" href="{slug}.html">
        <div class="thumb">{img}</div>
        <div class="body">
          <h2 data-ko="{ko}" data-en="{en}">{ko}</h2>
          <p class="cdesc"></p>
          <div class="links"><span class="open"></span>{rep_link}</div>
        </div>
      </a>""")
    return "\n".join(items)


def main():
    html = _PAGE.replace("__CARDS__", cards()).replace("__DATE__", date.today().isoformat())
    (SITE / "index.html").write_text(html, encoding="utf-8")
    n = len([p for p in SITE.glob('*.html') if p.stem != 'index' and not p.stem.startswith('report')])
    print(f"wrote {SITE/'index.html'} with {n} dashboards")


_PAGE = r"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title id="ttl">Clinical Cohort Intelligence</title>
<style>
  *{box-sizing:border-box}
  body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Noto Sans KR",sans-serif;color:#0f172a;background:#f1f5f9}
  header{background:#0f172a;color:#fff;padding:22px;text-align:center;position:relative}
  header h1{margin:0;font-size:22px}
  header p{margin:8px auto 0;color:#94a3b8;font-size:13px;max-width:680px;line-height:1.6}
  .lang{position:absolute;top:16px;right:18px;display:flex;border:1px solid #334155;border-radius:8px;overflow:hidden}
  .lang button{background:#1e293b;color:#94a3b8;border:0;padding:6px 11px;font-size:12px;cursor:pointer}
  .lang button.on{background:#2563eb;color:#fff}
  .container{max-width:900px;margin:0 auto;padding:26px 22px 40px}
  .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px}
  .card{display:flex;flex-direction:column;background:#fff;border:1px solid #e5e7eb;border-radius:14px;
        overflow:hidden;text-decoration:none;color:inherit;transition:.15s;box-shadow:0 1px 2px rgba(0,0,0,.04)}
  .card:hover{transform:translateY(-3px);box-shadow:0 8px 24px rgba(2,6,23,.12);border-color:#cbd5e1}
  .thumb{aspect-ratio:16/10;background:#f8fafc;overflow:hidden;border-bottom:1px solid #eef2f7}
  .thumb img{width:100%;height:100%;object-fit:cover;object-position:center;display:block}
  .body{padding:15px 17px;display:flex;flex-direction:column;flex:1}
  .body h2{margin:0 0 6px;font-size:17px}
  .cdesc{margin:0 0 13px;color:#64748b;font-size:12.5px;line-height:1.5;flex:1}
  .links{display:flex;align-items:center;justify-content:space-between;gap:12px}
  .open{color:#2563eb;font-weight:600;font-size:13px}
  .rep{color:#475569;font-size:12px;text-decoration:none}.rep:hover{text-decoration:underline}
  footer{text-align:center;color:#94a3b8;font-size:12px;padding:0 22px 34px}
</style></head>
<body>
<header>
  <div class="lang"><button id="ko" class="on">KO</button><button id="en">EN</button></div>
  <h1 id="h1"></h1><p id="desc"></p>
</header>
<div class="container"><div class="grid">__CARDS__</div></div>
<footer id="foot"></footer>
<script>
const DATE="__DATE__";
const I18N={
 ko:{title:"임상 코호트 인텔리전스",
     desc:"진행 중인 항암 임상 코호트를 분석해, 공간전사체 + AI 종양미세환경(TME) 분석 관점에서 표준치료(B) 데이터 자산의 가치를 우선순위화합니다. 암종을 선택하세요.",
     cdesc:"진행 중 임상 코호트 · 표준치료(B) 데이터 자산 · 클러스터 랜드스케이프",
     open:"대시보드 열기 →", report:"📄 분석 리포트", foot:d=>`업데이트 ${d} · 데이터 출처: ClinicalTrials.gov`},
 en:{title:"Clinical Cohort Intelligence",
     desc:"Analyze ongoing oncology clinical cohorts and prioritize standard-of-care (B) data assets for spatial-transcriptomics + AI tumor-microenvironment analysis. Choose a cancer type.",
     cdesc:"Ongoing cohorts · standard-of-care (B) data assets · cluster landscape",
     open:"Open dashboard →", report:"📄 Analysis report", foot:d=>`Updated ${d} · Source: ClinicalTrials.gov`}
};
let lang=localStorage.getItem('lang')||'ko';
function applyLang(l){lang=l;const T=I18N[l];localStorage.setItem('lang',l);document.documentElement.lang=l;
  document.getElementById('ko').classList.toggle('on',l==='ko');
  document.getElementById('en').classList.toggle('on',l==='en');
  document.getElementById('ttl').textContent=T.title;
  document.getElementById('h1').textContent=T.title;
  document.getElementById('desc').textContent=T.desc;
  document.getElementById('foot').textContent=T.foot(DATE);
  document.querySelectorAll('.card h2').forEach(h=>h.textContent=l==='ko'?h.dataset.ko:h.dataset.en);
  document.querySelectorAll('.cdesc').forEach(e=>e.textContent=T.cdesc);
  document.querySelectorAll('.open').forEach(e=>e.textContent=T.open);
  document.querySelectorAll('.rep').forEach(e=>e.textContent=T.report);}
document.getElementById('ko').onclick=()=>applyLang('ko');
document.getElementById('en').onclick=()=>applyLang('en');
applyLang(lang);
</script>
</body></html>"""


if __name__ == "__main__":
    main()
