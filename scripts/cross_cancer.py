"""Cross-cancer meta-analysis: rank cancer types and (cancer x standard-of-care B)
data assets by BD value for the spatial-transcriptomics + AI TME service.

Aggregates every cancer's features + per-B value ranking into:
  - outputs/cross_cancer_summary.csv   (cancer-level opportunity ranking)
  - outputs/global_dataset_ranking.csv ((cancer x B) datasets, globally scored)
  - outputs/cross_cancer_report.md      (narrative report)
  - outputs/cross_cancer_opportunity.png (chart)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cancer_config import CONFIGS
from drug_normalize import soc_from_text, soc_classes

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "outputs"


def nz(s):
    s = s.astype(float)
    return (s - s.min()) / (s.max() - s.min() + 1e-9)


def main():
    rows, gparts = [], []
    for name, cfg in CONFIGS.items():
        slug = cfg["slug"]
        f = DATA / f"cohorts_{slug}_features.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        df = df[df.is_industry == 1].copy()
        soc = df.comparator_B.fillna("").apply(soc_from_text)
        df["has_soc"] = soc.apply(lambda s: bool(soc_classes(s)))
        df["B_soc"] = soc.apply(lambda s: ";".join(sorted(s)))
        d = df[df.has_soc].copy()
        if len(d) == 0:
            continue
        # cancer-level metrics
        rows.append({
            "cancer": name,
            "soc_cohorts": len(d),
            "phase3": int((d.phase_simple == "Phase 3").sum()),
            "bigpharma_trials": int(d.is_big_pharma.sum()),
            "sponsors": int(d.lead_sponsor.nunique()),
            "enroll_total": int(d.enrollment.fillna(0).sum()),
            "tme_rel%": round(100 * d.tme_relevant.mean()),
            "distinct_B": len({x for s in d.B_soc for x in s.split(";") if x}),
        })
        # (cancer x B) rows from the per-cancer value ranking
        vr = OUT / f"value_ranking_{slug}.csv"
        if vr.exists():
            v = pd.read_csv(vr)
            v.insert(0, "cancer", name)
            gparts.append(v)

    summ = pd.DataFrame(rows)
    # top standard-of-care B per cancer (by demand) for context
    glob = pd.concat(gparts, ignore_index=True) if gparts else pd.DataFrame()
    topB = (glob.sort_values("demand_breadth", ascending=False)
            .groupby("cancer").first()["B (standard therapy)"]
            if len(glob) else pd.Series(dtype=str))
    summ["top_B"] = summ.cancer.map(topB).fillna("-")
    summ["top_B_breadth"] = summ.cancer.map(
        glob.groupby("cancer")["demand_breadth"].max() if len(glob) else {}).fillna(0).astype(int)

    # cancer opportunity score: demand + urgency + partnerability + volume
    summ["opportunity"] = (
        0.30 * nz(summ.soc_cohorts) + 0.15 * nz(summ.top_B_breadth)
        + 0.20 * nz(summ.phase3)
        + 0.15 * nz(summ.bigpharma_trials) + 0.10 * nz(summ.sponsors)
        + 0.10 * nz(np.log1p(summ.enroll_total))
    ).round(3)
    summ = summ.sort_values("opportunity", ascending=False).reset_index(drop=True)
    summ.to_csv(OUT / "cross_cancer_summary.csv", index=False)

    # global (cancer x B) dataset ranking, re-normalized across ALL cancers
    if len(glob):
        g = glob.copy()
        g["global_value"] = (
            0.45 * nz(g.demand_breadth) + 0.15 * nz(g.sponsors) + 0.13 * nz(g.bigpharma)
            + 0.15 * nz(g.n_phase3) + 0.07 * nz(np.log1p(g.enroll_sum))
            + 0.05 * nz(g["tme_rel%"])
        ).round(3)
        g = g.sort_values("global_value", ascending=False).reset_index(drop=True)
        cols = ["cancer", "B (standard therapy)", "demand_breadth", "sponsors", "bigpharma",
                "n_phase3", "enroll_sum", "tme_rel%", "global_value"]
        g[cols].to_csv(OUT / "global_dataset_ranking.csv", index=False)

    # chart
    fig, ax = plt.subplots(figsize=(9, 7))
    s2 = summ.iloc[::-1]
    ax.barh(s2.cancer, s2.opportunity, color="#4C78A8")
    ax.set_title("Cancer-type opportunity score (spatial-TME data assets)")
    ax.set_xlabel("opportunity score")
    plt.tight_layout()
    plt.savefig(OUT / "cross_cancer_opportunity.png", dpi=140, bbox_inches="tight")

    # markdown report
    _report(summ, g if len(glob) else pd.DataFrame())
    print("wrote cross_cancer_summary.csv, global_dataset_ranking.csv, cross_cancer_report.md")
    print("\n=== TOP CANCERS ===")
    print(summ[["cancer", "soc_cohorts", "phase3", "bigpharma_trials", "top_B", "top_B_breadth", "opportunity"]].head(10).to_string(index=False))
    print("\n=== TOP 15 (cancer x B) DATASETS ===")
    print(g[cols].head(15).to_string(index=False))


def _report(summ, g):
    from datetime import date
    lines = [f"# 암종·데이터셋 가치 종합 분석 (cross-cancer)\n",
             f"_생성일 {date.today().isoformat()} · 데이터 출처 ClinicalTrials.gov · 19개 고형암 종합_\n",
             "## 배경\n",
             "공간전사체 + AI 종양미세환경(TME) 분석 서비스 관점에서, **표준치료(B)로 치료받은 환자의 TME 데이터**는 "
             "그 B를 대조군/backbone으로 쓰는 진행 중 임상이 많을수록(여러 제약사 수요) 가치가 높다. "
             "아래는 19개 고형암 전체를 통합해 (1) 고가치 암종, (2) 고가치 (암종×표준치료 B) 데이터셋을 우선순위화한 결과다.\n",
             "## 1. 고가치 암종 순위\n",
             "지표: SoC 코호트 수(수요) · Phase 3 수(시급성) · 빅파마 참여·스폰서 다양성(협업성) · 총 등록 규모(데이터량).\n"]
    t = summ.copy()
    t = t[["cancer", "soc_cohorts", "phase3", "bigpharma_trials", "sponsors", "enroll_total", "top_B", "top_B_breadth", "opportunity"]]
    lines.append(t.to_markdown(index=False))
    lines.append("\n## 2. 고가치 데이터셋 순위 (암종 × 표준치료 B) — 전역 스코어 상위 20\n")
    if len(g):
        gg = g[["cancer", "B (standard therapy)", "demand_breadth", "sponsors", "bigpharma", "n_phase3", "global_value"]].head(20)
        lines.append(gg.to_markdown(index=False))
    lines.append("\n## 3. 해석\n")
    top3 = ", ".join(summ.cancer.head(3))
    lines.append(f"- **최우선 암종**: {top3} — SoC 코호트가 많고 표준치료가 집중돼 하나의 TME 데이터셋이 다수 스폰서 수요를 충족.\n")
    if len(g):
        b0 = g.iloc[0]
        lines.append(f"- **최우선 데이터셋**: {b0['cancer']} × {b0['B (standard therapy)']} "
                     f"(진행 중 {int(b0.demand_breadth)}개 임상·스폰서 {int(b0.sponsors)}·빅파마 {int(b0.bigpharma)}) — "
                     "이 표준치료로 치료받은 환자의 공간전사체 데이터를 만들면 가장 많은 스폰서가 즉시 관심.\n")
    lines.append("- 상세 대시보드는 암종별 페이지 참조.\n")
    (OUT / "cross_cancer_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
