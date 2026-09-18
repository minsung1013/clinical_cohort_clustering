#!/bin/bash
# Batch-run the full pipeline for all remaining solid-tumor cancers.
cd /Users/minsung.choi/clinical_data_finder || exit 1

pairs=(
  "Bladder:bladder" "Breast:breast" "Endometrial:endometrial" "Esophageal:esoph"
  "Glioma:glioma" "Liver:liver" "Melanoma:melanoma" "Mesothelioma:meso"
  "Neuroendocrine:net" "Ovarian:ovarian" "Pancreatic:panc" "Prostate:prostate"
  "Renal:renal" "Sarcoma:sarcoma" "Skin:skin"
)

for p in "${pairs[@]}"; do
  name="${p%%:*}"; slug="${p##*:}"
  echo "===== $name ($slug)  $(date +%H:%M:%S) ====="
  export CANCER="$name"
  python3 scripts/build_scope.py 2>&1 | sed 's/^/  /'
  python3 scripts/enrich_ctgov.py >/dev/null 2>&1 && echo "  enrich done"
  python3 scripts/extract_comparator.py >/dev/null 2>&1 && echo "  extract done"
  python3 scripts/_build_notebook.py >/dev/null 2>&1
  if CANCER="$name" jupyter nbconvert --to notebook --execute --inplace \
       "notebooks/cohort_analysis_${slug}.ipynb" --ExecutePreprocessor.timeout=700 >/dev/null 2>&1; then
    echo "  notebook OK"
    CANCER="$name" jupyter nbconvert --to html --no-input \
      "notebooks/cohort_analysis_${slug}.ipynb" --output "outputs/report_${slug}.html" >/dev/null 2>&1
  else
    echo "  notebook FAIL"
  fi
done

echo "===== rebuild site  $(date +%H:%M:%S) ====="
python3 scripts/build_site.py 2>&1 | sed 's/^/  /'
echo "ALL DONE $(date +%H:%M:%S)"
