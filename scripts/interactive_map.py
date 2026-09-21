"""Build a standalone bilingual (KO/EN) interactive HTML dashboard per cancer.

Header has a language toggle + sponsor filter dropdown; left = plotly map;
right = cluster overview + detail cards. Selecting a sponsor highlights its
cohorts on the map (red rings) and lists its trials (CT.gov link + full info).
Clicking a point opens that trial. All UI text switches KO<->EN client-side.
"""
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.colors as mcolors

from cancer_config import current

CFG = current()
PHASE_ALPHA = {"Phase 1": 0.35, "Phase 2": 0.65, "Phase 3": 1.0}


def _load_extra():
    """brief_summary + eligibility from the raw collected table, keyed by NCT."""
    p = Path(__file__).resolve().parent.parent / "data" / f"cohorts_{CFG['slug']}.csv"
    if not p.exists():
        return {}
    raw = pd.read_csv(p)
    ex = {}
    for _, r in raw.iterrows():
        ex[r["nct_id"]] = {
            "summary": str(r.get("brief_summary") or "")[:600],
            "eligibility": str(r.get("eligibility") or "")[:700],
            "start": str(r.get("start_date") or ""),
            "pcd": str(r.get("primary_completion_date") or ""),
        }
    return ex


def _build_fig(dfp, hexcol):
    scale = dfp.sponsor_scale.fillna(1).clip(lower=1)
    smin, smax = np.sqrt(scale.min()), np.sqrt(scale.max())
    dfp = dfp.assign(msize=np.interp(np.sqrt(scale), [smin, smax], [7, 40]))

    def _reg_disp(r):
        add = str(r.get("regimen_add", "") or "")
        return str(r.get("regimen", "-")) + (" + " + add if add and add not in ("chemo only", "-", "—", "backbone only") else "")

    fig = go.Figure()
    for ph, a in PHASE_ALPHA.items():
        s = dfp[dfp.phase_simple == ph]
        if not len(s):
            continue
        cd = np.column_stack([
            s.nct_id.values, s.lead_sponsor.fillna("-").values,
            s.apply(_reg_disp, axis=1).values, s.phase_simple.values,
            s.overall_status.fillna("-").values,
        ])
        fig.add_trace(go.Scatter(
            x=s.x, y=s.y, mode="markers", name=ph,
            marker=dict(size=s.msize, opacity=a,
                        color=[hexcol[int(c)] for c in s.cluster],
                        line=dict(width=[2 if b else 0.4 for b in s.is_big_pharma],
                                  color=["black" if b else "lightgray" for b in s.is_big_pharma])),
            customdata=cd,
            hovertemplate=("<b>%{customdata[1]}</b><br>%{customdata[2]}<br>"
                           "%{customdata[3]} · %{customdata[4]}<br>"
                           "<span style='color:#94a3b8'>%{customdata[0]}</span><extra></extra>")))
    fig.add_trace(go.Scatter(x=[], y=[], mode="markers", name="sel",
                             marker=dict(size=26, color="rgba(0,0,0,0)",
                                         line=dict(width=3, color="#e11d48")),
                             hoverinfo="skip", showlegend=False))
    hl_index = len(fig.data) - 1
    fig.update_layout(template="plotly_white", autosize=True,
                      legend_title="Phase", margin=dict(l=8, r=8, t=10, b=8),
                      legend=dict(orientation="h", y=1.02, x=0))
    div = fig.to_html(full_html=False, include_plotlyjs="cdn", div_id="map",
                      default_height="100%", default_width="100%",
                      config={"responsive": True, "displaylogo": False})
    return div, hl_index


def build_interactive(df, rep_tbl, palette, bm_cols, out_path):
    dfp = df.copy()
    hexcol = {int(c): mcolors.to_hex(palette[c]) for c in sorted(dfp.cluster.unique())}
    fig_div, hl_index = _build_fig(dfp, hexcol)
    extra = _load_extra()

    trials = []
    for _, r in dfp.iterrows():
        e = extra.get(r.nct_id, {})
        trials.append({
            "nct": r.nct_id, "title": str(r.brief_title), "phase": r.phase_simple,
            "status": str(r.get("overall_status") or ""),
            "cluster": int(r.cluster), "x": float(r.x), "y": float(r.y),
            "regimen": r.get("regimen", "-"), "add": r.get("regimen_add", ""),
            "B": r.B_soc or "-", "A": r.get("experimental_A") or "-",
            "targets": r.get("targets") or "-", "modalities": r.get("modalities") or "-",
            "tme": int(r.get("tme_relevant", 0)),
            "bio": (r.get("biomarker_detail") if pd.notna(r.get("biomarker_detail")) and r.get("biomarker_detail") else "-"),
            "elig_sum": (r.get("elig_summary") if pd.notna(r.get("elig_summary")) and r.get("elig_summary") else ""),
            "sponsor": r.lead_sponsor or "-",
            "scale": int(r.sponsor_scale) if pd.notna(r.sponsor_scale) else 0,
            "enroll": int(r.enrollment) if pd.notna(r.enrollment) else None,
            "big": int(r.is_big_pharma), "start": e.get("start", ""), "pcd": e.get("pcd", ""),
            "summary": e.get("summary", ""), "eligibility": e.get("eligibility", ""),
        })

    NEAR = " (빅파마 없음·최근접)"
    clusters = []
    for c in sorted(dfp.cluster.unique()):
        sub = dfp[dfp.cluster == c]
        cbt = Counter(x for s in sub.B_soc for x in s.split(";") if x)
        rep = rep_tbl[rep_tbl.cluster == c]
        rep_name = (rep.iloc[0]["대표 제약사"] if len(rep) else "-")
        near = rep_name.endswith(NEAR)
        clusters.append({
            "id": int(c), "n": int(len(sub)), "color": hexcol[int(c)],
            "name": ", ".join(k for k, _ in cbt.most_common(3)) or "-",
            "rep": rep_name.replace(NEAR, ""), "near": int(near),
            "B": (rep.iloc[0]["대조군 치료 B"] if len(rep) else "-"),
        })

    comp = (dfp.groupby("lead_sponsor")
            .agg(scale=("sponsor_scale", "max"), n=("nct_id", "count"), big=("is_big_pharma", "max"))
            .reset_index().sort_values("lead_sponsor", key=lambda s: s.str.lower()))
    companies = [{"name": r.lead_sponsor, "scale": int(r.scale), "n": int(r.n), "big": int(r.big)}
                 for _, r in comp.iterrows() if pd.notna(r.lead_sponsor)]

    ko = CFG["title"].split("(")[0].strip()
    en = CFG["name"]
    html = (_TEMPLATE
            .replace("__CANCER_KO__", ko).replace("__CANCER_EN__", en)
            .replace("__FIG__", fig_div).replace("__HL__", str(hl_index))
            .replace("__TRIALS__", json.dumps(trials, ensure_ascii=False))
            .replace("__CLUSTERS__", json.dumps(clusters, ensure_ascii=False))
            .replace("__COMPANIES__", json.dumps(companies, ensure_ascii=False)))
    Path(out_path).write_text(html, encoding="utf-8")
    return out_path


_TEMPLATE = r"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title id="ttl"></title>
<script src="https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js"></script>
<style>
  *{box-sizing:border-box}
  html,body{height:100%;margin:0}
  body{display:flex;flex-direction:column;overflow:hidden;color:#1f2937;background:#f1f5f9;
       font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Noto Sans KR",sans-serif}
  header{flex:0 0 auto;background:#0f172a;color:#fff;padding:9px 16px;display:flex;align-items:center;gap:12px;flex-wrap:wrap}
  header h1{font-size:15px;margin:0;font-weight:700}
  header .sub{color:#94a3b8;font-size:11px}
  .ctrl{margin-left:auto;display:flex;align-items:center;gap:8px}
  .ctrl label{font-size:12px;color:#cbd5e1}
  .lang{display:flex;border:1px solid #334155;border-radius:8px;overflow:hidden}
  .lang button{background:#1e293b;color:#94a3b8;border:0;padding:6px 10px;font-size:12px;cursor:pointer}
  .lang button.on{background:#2563eb;color:#fff}
  select{padding:7px 10px;border:1px solid #334155;border-radius:8px;font-size:13px;background:#fff;min-width:300px}
  #reset{padding:7px 12px;border:1px solid #475569;background:#1e293b;color:#fff;border-radius:8px;font-size:12px;cursor:pointer}
  .wrap{flex:1 1 auto;min-height:0;display:flex;gap:12px;padding:12px;overflow:hidden}
  .map{flex:1 1 60%;min-width:0;height:100%;background:#fff;border:1px solid #e5e7eb;border-radius:12px;padding:6px;overflow:hidden}
  .map .js-plotly-plot,.map .plotly-graph-div{height:100%!important;width:100%!important}
  .side{flex:1 1 40%;min-width:340px;max-width:560px;height:100%;display:flex;flex-direction:column;gap:12px;overflow:hidden}
  .panel{background:#fff;border:1px solid #e5e7eb;border-radius:12px;padding:10px 13px}
  .panel.clusters{flex:0 0 auto;max-height:26vh;overflow:auto}
  .panel.cards{flex:1 1 auto;min-height:0;display:flex;flex-direction:column;overflow:hidden}
  .muted{color:#64748b;font-size:12px}
  details summary{cursor:pointer;font-size:13px;font-weight:600;margin-bottom:6px}
  .clu{display:flex;gap:8px;padding:5px;border-radius:7px;font-size:12px}
  .dot{display:inline-block;width:11px;height:11px;border-radius:50%;margin-top:3px;flex:0 0 auto}
  .clu .sub{color:#64748b;font-size:11px}
  .cnt{flex:0 0 auto;color:#334155;font-size:13px;font-weight:600;margin:0 0 8px}
  #cards{flex:1 1 auto;min-height:0;overflow:auto}
  .card{border:1px solid #e5e7eb;border-left-width:6px;border-radius:9px;padding:10px 12px;margin-bottom:10px;background:#fff}
  .card a{color:#2563eb;text-decoration:none;font-size:13px}.card a:hover{text-decoration:underline}
  .title{font-size:12.5px;margin:4px 0 6px;color:#111827;line-height:1.35}
  .row{font-size:12px;color:#374151;margin:3px 0}.row b{color:#0f172a}
  .chip{display:inline-block;border-radius:5px;font-size:11px;padding:1px 7px;margin:0 4px 3px 0;background:#eef2ff;color:#3730a3}
  .chip.ph{background:#ecfeff;color:#155e75}.chip.st{background:#f1f5f9;color:#334155}.chip.tme{background:#fef2f2;color:#b91c1c}
  .big{background:#111827;color:#fff;border-radius:4px;font-size:10px;padding:1px 6px;margin-left:6px}
  details.sum{margin-top:6px}details.sum summary{font-size:11.5px;color:#2563eb;font-weight:500}
  details.sum p{font-size:11.5px;color:#475569;margin:5px 0 0;line-height:1.45}
  .exbar{flex:0 0 auto;display:flex;align-items:center;gap:8px;padding:1px 2px 9px;border-bottom:1px solid #eef2f7;margin-bottom:8px;flex-wrap:wrap}
  .exbar .wlc{font-size:12px;color:#92400e;font-weight:600}
  .exbar button{border:1px solid #cbd5e1;background:#f8fafc;color:#0f172a;border-radius:7px;font-size:12px;padding:5px 10px;cursor:pointer}
  .exbar button:hover{background:#eef2ff;border-color:#a5b4fc}
  .exbtns{margin-left:auto;display:flex;gap:6px}
  .card{position:relative;padding-right:32px}
  .card.wl{background:#fffbeb;border-color:#fcd34d}
  .card .star{position:absolute;top:8px;right:9px;font-size:18px;line-height:1;cursor:pointer;color:#cbd5e1;background:none;border:0;padding:2px}
  .card .star.on{color:#f59e0b}
  .card .star:hover{color:#f59e0b}
</style></head>
<body>
<header>
  <h1 id="h1"></h1><span class="sub" id="sub"></span>
  <div class="ctrl">
    <div class="lang"><button id="ko" class="on">KO</button><button id="en">EN</button></div>
    <label id="lbl_sponsor"></label>
    <select id="company"></select>
    <button id="reset"></button>
  </div>
</header>
<div class="wrap">
  <div class="map">__FIG__</div>
  <div class="side">
    <div class="panel clusters"><details open><summary id="cl_head"></summary><div id="clusters"></div></details></div>
    <div class="panel cards">
      <div class="exbar">
        <span class="wlc" id="wlcount"></span>
        <button id="wlView"></button>
        <button id="wlClear"></button>
        <span class="exbtns"><button id="exCsv"></button><button id="exXlsx"></button></span>
      </div>
      <div id="cards"></div>
    </div>
  </div>
</div>
<script>
const TRIALS=__TRIALS__, CLUSTERS=__CLUSTERS__, COMPANIES=__COMPANIES__;
const HL=__HL__, CK="__CANCER_KO__", CE="__CANCER_EN__";
const CC={}; CLUSTERS.forEach(c=>CC[c.id]=c.color);
const gd=document.getElementById('map');

const I18N={
 ko:{cancer:CK, title:c=>`${c} 임상 코호트 대시보드`,
     sub:c=>`진행 중 ${c} · 제약사 주도 · 표준치료 코호트 (색=클러스터 · 투명도=Phase · 크기=R&D 규모 · 테두리=빅파마)`,
     sponsor:"제약사", all:"— 전체 보기 —", reset:"초기화",
     clusters:"클러스터 개요 (표준치료 아키타입)", rep:"대표", near:"(빅파마 없음·최근접)",
     hint:"상단에서 제약사를 선택하거나 맵의 점을 클릭하세요. 카드를 우클릭(또는 ☆ 클릭)하면 관심 목록에 담깁니다.",
     reg:"표준치료 레지멘(B)", A:"실험약(A)", target:"타겟", modality:"모달리티",
     bio:"바이오마커(선택)", elig_sum:"선정기준(요약)", spon:"스폰서", start:"시작", pcd:"1차완료(예정)",
     summary:"연구 요약", elig:"선정기준(발췌)", trials:n=>`${n}건`, scale:"R&D 규모",
     wln:n=>`⭐ 관심 목록 ${n}개`, wlview:"관심목록 보기", wlclear:"비우기",
     csv:"CSV 내려받기", xlsx:"Excel 내려받기",
     wlhead:n=>`⭐ 관심 목록 — ${n}건`,
     nosel:"관심 목록이 비어 있습니다. 카드를 우클릭하거나 ☆를 눌러 담으세요.",
     noxlsx:"Excel 라이브러리를 불러오지 못했습니다(인터넷 연결 확인). CSV로 받아주세요."},
 en:{cancer:CE, title:c=>`${c} clinical-cohort dashboard`,
     sub:c=>`Ongoing ${c} · industry-sponsored · standard-of-care cohorts (color=cluster · opacity=phase · size=R&D scale · outline=big pharma)`,
     sponsor:"Sponsor", all:"— Show all —", reset:"Reset",
     clusters:"Cluster overview (standard-of-care archetypes)", rep:"Rep", near:"(no big pharma · nearest)",
     hint:"Select a sponsor above, or click a point on the map. Right-click a card (or click ☆) to add it to your watchlist.",
     reg:"Standard-of-care regimen (B)", A:"Experimental (A)", target:"Target", modality:"Modality",
     bio:"Biomarker (selection)", elig_sum:"Eligibility (summary)", spon:"Sponsor", start:"Start", pcd:"Primary completion (est.)",
     summary:"Study summary", elig:"Eligibility (excerpt)", trials:n=>`${n} trial${n>1?'s':''}`, scale:"R&D",
     wln:n=>`⭐ Watchlist: ${n}`, wlview:"View watchlist", wlclear:"Clear",
     csv:"Download CSV", xlsx:"Download Excel",
     wlhead:n=>`⭐ Watchlist — ${n} trial${n>1?'s':''}`,
     nosel:"Watchlist is empty. Right-click a card or click ☆ to add trials.",
     noxlsx:"Excel library failed to load (check your connection). Please use CSV."}
};
let lang=localStorage.getItem('lang')||'ko', T=I18N[lang];
let view={type:'none'};
const WKEY='watch_'+CE;   // watchlist (관심 목록) persisted per cancer
const WATCH=new Set(JSON.parse(localStorage.getItem(WKEY)||'[]'));
function saveWatch(){localStorage.setItem(WKEY,JSON.stringify([...WATCH]));}

function fmt(s){return s&&s!=='nan'?s:'—';}
function card(t){const url='https://clinicaltrials.gov/study/'+t.nct;
  const reg=t.regimen+(t.add&&!['chemo only','—','-','backbone only'].includes(t.add)?' + '+t.add:'');
  const on=WATCH.has(t.nct);
  return `<div class="card${on?' wl':''}" data-nct="${t.nct}" style="border-left-color:${CC[t.cluster]||'#ccc'}">
    <button class="star${on?' on':''}" data-nct="${t.nct}" title="관심 목록 / watchlist">${on?'★':'☆'}</button>
    <a href="${url}" target="_blank"><b>${t.nct}</b> ↗</a>${t.big?'<span class="big">★ big pharma</span>':''}
    <div class="title">${t.title}</div>
    <div class="row"><span class="chip ph">${t.phase}</span><span class="chip st">${t.status}</span>
      <span class="chip">C${t.cluster}</span>${t.tme?'<span class="chip tme">TME</span>':''}<span class="chip">enroll ${t.enroll??'—'}</span></div>
    <div class="row"><b>${T.reg}:</b> ${reg}</div>
    <div class="row"><b>${T.A}:</b> ${t.A}</div>
    <div class="row"><b>${T.target}:</b> ${t.targets}　<b>${T.modality}:</b> ${t.modalities}</div>
    <div class="row"><b>${T.bio}:</b> ${t.bio}</div>
    ${t.elig_sum?`<div class="row"><b>${T.elig_sum}:</b> ${t.elig_sum}</div>`:''}
    <div class="row"><b>${T.spon}:</b> ${t.sponsor} <span class="muted">(${T.scale} ${t.scale})</span></div>
    <div class="row muted">${T.start} ${fmt(t.start)} · ${T.pcd} ${fmt(t.pcd)}</div>
    ${t.summary?`<details class="sum"><summary>${T.summary}</summary><p>${t.summary}</p></details>`:''}
    ${t.eligibility?`<details class="sum"><summary>${T.elig}</summary><p>${t.eligibility}</p></details>`:''}
  </div>`;}
function renderCards(ts,head){document.getElementById('cards').innerHTML=`<div class="cnt">${head}</div>`+(ts.length?ts.map(card).join(''):`<div class="muted">—</div>`);}
function syncWatchUI(){document.getElementById('wlcount').textContent=T.wln(WATCH.size);}

// --- watchlist (관심 목록 / cart) ---
let WLI=null;   // plotly trace index for watchlist star markers
function initWatchTrace(){
  if(!window.Plotly||!gd||!gd.data){setTimeout(initWatchTrace,200);return;}
  WLI=gd.data.length;
  Plotly.addTraces('map',{x:[],y:[],mode:'markers',name:'watch',
    marker:{symbol:'star',size:15,color:'#f59e0b',line:{width:1.3,color:'#7c2d12'}},
    hoverinfo:'skip',showlegend:false}).then(updateWatchMap);
}
function updateWatchMap(){if(WLI==null)return;const p=TRIALS.filter(t=>WATCH.has(t.nct));
  Plotly.restyle('map',{x:[p.map(t=>t.x)],y:[p.map(t=>t.y)]},[WLI]);}
function showWatch(){view={type:'watch'};const ts=TRIALS.filter(t=>WATCH.has(t.nct)).sort((a,b)=>a.cluster-b.cluster);
  highlight(ts.map(t=>[t.x,t.y]));renderCards(ts,T.wlhead(ts.length));}
function toggleWatch(nct){if(WATCH.has(nct))WATCH.delete(nct);else WATCH.add(nct);
  saveWatch();updateWatchMap();syncWatchUI();
  const c=document.querySelector(`#cards .card[data-nct="${nct}"]`);
  if(c){const on=WATCH.has(nct);c.classList.toggle('wl',on);
    const b=c.querySelector('.star');if(b){b.classList.toggle('on',on);b.textContent=on?'★':'☆';}}
  if(view.type==='watch')showWatch();}
function clearWatch(){if(!WATCH.size)return;WATCH.clear();saveWatch();updateWatchMap();syncWatchUI();
  if(view.type==='watch')showWatch();
  else document.querySelectorAll('#cards .card').forEach(c=>{c.classList.remove('wl');
    const b=c.querySelector('.star');if(b){b.classList.remove('on');b.textContent='☆';}});}

const EXCOLS=[["nct","NCT"],["title","Title"],["phase","Phase"],["status","Status"],["cluster","Cluster"],
  ["B","StandardOfCare_B"],["A","Experimental_A"],["regimen","Regimen"],["add","Regimen_add"],
  ["targets","Targets"],["modalities","Modalities"],["tme","TME_relevant"],["bio","Biomarker"],
  ["elig_sum","Eligibility_summary"],["sponsor","Sponsor"],["big","BigPharma"],["scale","Sponsor_RnD_scale"],
  ["enroll","Enrollment"],["start","Start"],["pcd","PrimaryCompletion"],["url","URL"]];
function watchRows(){const rows=[];TRIALS.forEach(t=>{if(!WATCH.has(t.nct))return;const o={};
  EXCOLS.forEach(([k,h])=>{o[h]=k==='url'?('https://clinicaltrials.gov/study/'+t.nct):(t[k]??'');});rows.push(o);});return rows;}
function toCSV(rows){const hs=Object.keys(rows[0]);
  const esc=v=>{v=String(v==null?'':v);return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;};
  return hs.join(',')+'\n'+rows.map(r=>hs.map(h=>esc(r[h])).join(',')).join('\n');}
function fname(ext){return `${CE.replace(/[^A-Za-z0-9]+/g,'_')}_cohorts_${new Date().toISOString().slice(0,10)}.${ext}`;}
function dl(blob,name){const u=URL.createObjectURL(blob),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),1000);}
function exportCSV(){const rows=watchRows();if(!rows.length){alert(T.nosel);return;}
  dl(new Blob(['﻿'+toCSV(rows)],{type:'text/csv;charset=utf-8;'}),fname('csv'));}
function exportXLSX(){const rows=watchRows();if(!rows.length){alert(T.nosel);return;}
  if(typeof XLSX==='undefined'){alert(T.noxlsx);return;}
  const wb=XLSX.utils.book_new();XLSX.utils.book_append_sheet(wb,XLSX.utils.json_to_sheet(rows),'cohorts');
  XLSX.writeFile(wb,fname('xlsx'));}
function renderClusters(){document.getElementById('clusters').innerHTML='';
  CLUSTERS.forEach(c=>{const d=document.createElement('div');d.className='clu';
    d.innerHTML=`<span class="dot" style="background:${c.color}"></span>
      <div><b>C${c.id}</b> · ${c.name} <span class="sub">(n=${c.n})</span>
      <div class="sub">${T.rep}: ${c.rep}${c.near?' '+T.near:''} · B: ${c.B}</div></div>`;
    d.onclick=()=>{view={type:'cluster',id:c.id};
      const ts=TRIALS.filter(x=>x.cluster===c.id).sort((a,b)=>b.scale-a.scale);
      highlight(ts.map(x=>[x.x,x.y]));renderCards(ts,`C${c.id} — ${T.trials(ts.length)}`);};
    document.getElementById('clusters').appendChild(d);});}
function highlight(coords){Plotly.restyle('map',{x:[coords.map(c=>c[0])],y:[coords.map(c=>c[1])]},[HL]);}
function rebuildSelect(){const sel=document.getElementById('company');const cur=sel.value;
  sel.innerHTML=`<option value="__ALL__">${T.all}</option>`+COMPANIES.map(c=>
    `<option value="${c.name.replace(/"/g,'&quot;')}">${c.name} — ${T.scale} ${c.scale} · ${T.trials(c.n)}${c.big?'  ★':''}</option>`).join('');
  if(cur)sel.value=cur;}
function showCompany(name){view={type:'company',name};
  if(name==='__ALL__'){highlight([]);renderCards([],T.hint);return;}
  const ts=TRIALS.filter(t=>t.sponsor===name).sort((a,b)=>a.cluster-b.cluster);
  highlight(ts.map(t=>[t.x,t.y]));renderCards(ts,`${name} — ${T.trials(ts.length)}`);}
function rerender(){if(view.type==='company')showCompany(view.name);
  else if(view.type==='cluster'){const ts=TRIALS.filter(x=>x.cluster===view.id).sort((a,b)=>b.scale-a.scale);renderCards(ts,`C${view.id} — ${T.trials(ts.length)}`);}
  else if(view.type==='trial'){const t=TRIALS.find(x=>x.nct===view.nct);if(t)renderCards([t],`${t.nct} · ${t.sponsor}`);}
  else if(view.type==='watch')showWatch();
  else renderCards([],T.hint);}

function applyLang(l){lang=l;T=I18N[l];localStorage.setItem('lang',l);
  document.documentElement.lang=l;
  document.getElementById('ko').classList.toggle('on',l==='ko');
  document.getElementById('en').classList.toggle('on',l==='en');
  document.getElementById('ttl').textContent=T.title(T.cancer);
  document.getElementById('h1').textContent=T.title(T.cancer);
  document.getElementById('sub').textContent=T.sub(T.cancer);
  document.getElementById('lbl_sponsor').textContent=T.sponsor;
  document.getElementById('reset').textContent=T.reset;
  document.getElementById('cl_head').textContent=T.clusters;
  document.getElementById('wlView').textContent=T.wlview;
  document.getElementById('wlClear').textContent=T.wlclear;
  document.getElementById('exCsv').textContent=T.csv;
  document.getElementById('exXlsx').textContent=T.xlsx;
  syncWatchUI();rebuildSelect();renderClusters();rerender();}

document.getElementById('ko').onclick=()=>applyLang('ko');
document.getElementById('en').onclick=()=>applyLang('en');
document.getElementById('company').addEventListener('change',e=>showCompany(e.target.value));
document.getElementById('reset').addEventListener('click',()=>{document.getElementById('company').value='__ALL__';showCompany('__ALL__');});
document.getElementById('cards').addEventListener('click',e=>{const b=e.target.closest('.star');if(!b)return;
  e.preventDefault();toggleWatch(b.dataset.nct);});
document.getElementById('cards').addEventListener('contextmenu',e=>{const c=e.target.closest('.card');if(!c)return;
  e.preventDefault();toggleWatch(c.dataset.nct);});
document.getElementById('wlView').addEventListener('click',showWatch);
document.getElementById('wlClear').addEventListener('click',clearWatch);
document.getElementById('exCsv').addEventListener('click',exportCSV);
document.getElementById('exXlsx').addEventListener('click',exportXLSX);
if(gd&&gd.on){gd.on('plotly_click',ev=>{const nct=ev.points[0].customdata?ev.points[0].customdata[0]:null;if(!nct)return;
  const t=TRIALS.find(x=>x.nct===nct);if(!t)return;view={type:'trial',nct};
  document.getElementById('company').value=t.sponsor;highlight([[t.x,t.y]]);renderCards([t],`${t.nct} · ${t.sponsor}`);
  document.getElementById('cards').scrollTop=0;});}
applyLang(lang);
initWatchTrace();
</script>
</body></html>"""
