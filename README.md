# Clinical Cohort Value Finder

진행 중인 오픈 임상 코호트(ClinicalTrials.gov)를 분석해, **공간전사체 + AI 종양미세환경(TME)
분석** 관점에서 가치가 높은 데이터 자산을 우선순위화한다.

## 핵심 가치 모델
`A + B vs B`(또는 `A vs B`) 구조의 임상에서 **B는 표준치료(대조군/backbone)**이고, 스폰서들은
"B 치료에 따른 TME 변화·반응 바이오마커"를 공통적으로 원한다. 따라서 특정 암종에서 **표준치료 B로
치료받은 환자의 TME 데이터**는, 그 B를 쓰는 진행 중 임상이 많을수록(수요 폭) 가치가 높다.

## 파이프라인 (암종 파라미터화 — `CANCER` 환경변수)
```
scripts/cancer_config.py     암종별 설정(slug·바이오마커 패턴). CANCER 환경변수로 선택(기본 Colorectal)
scripts/build_scope.py       pipeline.json → <CANCER>·진행 중·Phase1-3 → data/scope_<slug>.json
scripts/enrich_ctgov.py      CT.gov API v2 전량 수집(캐시 data/raw_ctgov/) → data/cohorts_<slug>.csv
scripts/drug_normalize.py    약물 정규화·레지멘(FOLFOX/FLOT/SOX/FP/XP…)·빅파마·TME·SoC 계열
scripts/extract_comparator.py arm type 기반 대조군 B 추출 + feature → data/cohorts_<slug>_features.csv
scripts/_build_notebook.py   분석 노트북 생성 → notebooks/cohort_analysis_<slug>.ipynb
scripts/interactive_map.py   인터랙티브 HTML 대시보드 빌더(제약사 드롭다운·하이라이트·상세카드·관심목록·export)
scripts/build_nct_links.py   NCT→논문/학회초록 연결맵 생성 → data/nct_links.json (oncology R2 코퍼스 기반)
```

## 대시보드 기능
- **관심 목록(장바구니)**: 카드 우클릭(또는 ☆)으로 담기 → 맵에 금색 별 표시, localStorage 저장, CSV/Excel export
- **관련 논문·초록**: 각 임상 카드에 연결된 논문(PubMed)·학회 초록(DOI) 표시. `build_nct_links.py`로 `data/nct_links.json` 생성 후 `interactive_map.py`가 대시보드에 임베드

## 실행 (예: 위암)
```bash
export CANCER=Gastric   # 또는 Colorectal (기본)
python3 scripts/build_scope.py
python3 scripts/enrich_ctgov.py        # 네트워크(재실행 시 캐시 스킵)
python3 scripts/extract_comparator.py
python3 scripts/_build_notebook.py
CANCER=$CANCER jupyter nbconvert --to notebook --execute --inplace notebooks/cohort_analysis_${CANCER,,}*.ipynb
```
슬러그: Colorectal→`crc`, Gastric→`gastric`. 새 암종은 `cancer_config.py`의 `CONFIGS`에 항목만 추가.

## 산출물 (`outputs/`, `<slug>` = crc | gastric)
- `value_ranking_<slug>.csv` / `.png` — 표준치료 B 데이터 자산 가치 랭킹(수요 폭 중심)
- `cluster_representatives_<slug>.csv` — 클러스터별 대표 빅파마 코호트
- `cohort_map_<slug>.png` — 정적 2D 맵 / `cohort_map_<slug>.html` — 인터랙티브 대시보드
- `report_<slug>.html` — 실행된 노트북 리포트

## 데이터 한계
`pipeline.json`은 약물·임상시험 단위 집계라 개별 환자 데이터는 없다. 가치는 구조화 필드(arm 구조,
대조군, enrollment, eligibility 바이오마커) 기반 휴리스틱 스코어링이다.
