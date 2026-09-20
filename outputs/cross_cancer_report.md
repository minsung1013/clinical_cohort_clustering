# 암종·데이터셋 가치 종합 분석 (cross-cancer)

_생성일 2026-09-19 · 데이터 출처 ClinicalTrials.gov · 19개 고형암 종합_

## 배경

공간전사체 + AI 종양미세환경(TME) 분석 서비스 관점에서, **표준치료(B)로 치료받은 환자의 TME 데이터**는 그 B를 대조군/backbone으로 쓰는 진행 중 임상이 많을수록(여러 제약사 수요) 가치가 높다. 아래는 19개 고형암 전체를 통합해 (1) 고가치 암종, (2) 고가치 (암종×표준치료 B) 데이터셋을 우선순위화한 결과다.

## 1. 고가치 암종 순위

지표: SoC 코호트 수(수요) · Phase 3 수(시급성) · 빅파마 참여·스폰서 다양성(협업성) · 총 등록 규모(데이터량).

| cancer         |   soc_cohorts |   phase3 |   bigpharma_trials |   sponsors |   enroll_total | top_B         |   top_B_breadth |   opportunity |
|:---------------|--------------:|---------:|-------------------:|-----------:|---------------:|:--------------|----------------:|--------------:|
| Lung           |           394 |      224 |                188 |        149 |         143628 | carboplatin   |             135 |         1     |
| Breast         |           228 |      142 |                115 |         83 |         128979 | fluorouracil  |              59 |         0.604 |
| Colorectal     |           111 |       37 |                 46 |         80 |          38156 | bevacizumab   |              40 |         0.314 |
| Ovarian        |            90 |       44 |                 38 |         56 |          30912 | paclitaxel    |              41 |         0.277 |
| Prostate       |            84 |       42 |                 47 |         43 |          43549 | abiraterone   |              38 |         0.272 |
| Liver          |            95 |       36 |                 34 |         59 |          32772 | bevacizumab   |              18 |         0.248 |
| Pancreatic     |            70 |       24 |                 13 |         61 |          22271 | gemcitabine   |              42 |         0.221 |
| Head & Neck    |            75 |       24 |                 30 |         59 |          22888 | pembrolizumab |              24 |         0.217 |
| Melanoma       |            66 |       25 |                 30 |         53 |          19874 | pembrolizumab |              34 |         0.216 |
| Bladder        |            62 |       35 |                 37 |         36 |          22261 | gemcitabine   |              26 |         0.209 |
| Gastric        |            66 |       30 |                 25 |         44 |          23498 | CAPOX         |              17 |         0.194 |
| Endometrial    |            53 |       28 |                 21 |         36 |          20458 | paclitaxel    |              23 |         0.177 |
| Glioma         |            51 |        9 |                  9 |         47 |           6703 | temozolomide  |              27 |         0.138 |
| Renal          |            40 |       17 |                 21 |         26 |          20033 | pembrolizumab |              10 |         0.135 |
| Esophageal     |            37 |       14 |                 13 |         29 |          12033 | paclitaxel    |              11 |         0.116 |
| Neuroendocrine |            17 |       11 |                  6 |         13 |           5693 | octreotide    |               5 |         0.059 |
| Sarcoma        |            15 |        6 |                  4 |         13 |           3023 | sunitinib     |               3 |         0.036 |
| Mesothelioma   |             7 |        1 |                  6 |          6 |           3049 | pembrolizumab |               3 |         0.022 |
| Skin           |            10 |        2 |                  5 |          8 |           1142 | pembrolizumab |               4 |         0.007 |

## 2. 고가치 데이터셋 순위 (암종 × 표준치료 B) — 전역 스코어 상위 20

| cancer     | B (standard therapy)   |   demand_breadth |   sponsors |   bigpharma |   n_phase3 |   global_value |
|:-----------|:-----------------------|-----------------:|-----------:|------------:|-----------:|---------------:|
| Lung       | carboplatin            |              135 |         58 |          78 |         97 |          0.984 |
| Lung       | pemetrexed             |              104 |         46 |          64 |         77 |          0.79  |
| Lung       | pembrolizumab          |              100 |         49 |          63 |         65 |          0.774 |
| Lung       | cisplatin              |               86 |         35 |          57 |         70 |          0.677 |
| Lung       | paclitaxel             |               53 |         26 |          36 |         36 |          0.457 |
| Breast     | fluorouracil           |               59 |         33 |          29 |         33 |          0.438 |
| Lung       | docetaxel              |               44 |         39 |          17 |         32 |          0.411 |
| Breast     | fulvestrant            |               53 |         28 |          29 |         32 |          0.401 |
| Breast     | paclitaxel             |               44 |         25 |          21 |         34 |          0.4   |
| Breast     | capecitabine           |               43 |         26 |          20 |         36 |          0.399 |
| Colorectal | bevacizumab            |               40 |         37 |          16 |         19 |          0.378 |
| Breast     | trastuzumab            |               45 |         20 |          20 |         29 |          0.376 |
| Ovarian    | paclitaxel             |               41 |         26 |          12 |         25 |          0.351 |
| Pancreatic | gemcitabine            |               42 |         37 |           4 |         16 |          0.336 |
| Pancreatic | paclitaxel             |               40 |         35 |           4 |         14 |          0.322 |
| Prostate   | abiraterone            |               38 |         21 |          23 |         25 |          0.32  |
| Melanoma   | pembrolizumab          |               34 |         29 |          18 |         11 |          0.32  |
| Lung       | etoposide              |               29 |         20 |          15 |         20 |          0.295 |
| Lung       | osimertinib            |               38 |         19 |          15 |         16 |          0.29  |
| Prostate   | enzalutamide           |               32 |         23 |          20 |         20 |          0.289 |

## 3. 해석

- **최우선 암종**: Lung, Breast, Colorectal — SoC 코호트가 많고 표준치료가 집중돼 하나의 TME 데이터셋이 다수 스폰서 수요를 충족.

- **최우선 데이터셋**: Lung × carboplatin (진행 중 135개 임상·스폰서 58·빅파마 78) — 이 표준치료로 치료받은 환자의 공간전사체 데이터를 만들면 가장 많은 스폰서가 즉시 관심.

- 상세 대시보드는 암종별 페이지 참조.
