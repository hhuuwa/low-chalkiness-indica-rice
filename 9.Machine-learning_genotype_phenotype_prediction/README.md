# Machine-learning genotype-to-phenotype prediction

This directory contains scripts for predicting chalkiness-related phenotypes from genotype-derived features.

## Scripts

- `1-vcf2table.py`
  - extracts SNP genotypes from selected VCF files and converts them into tabular feature matrices.
- `2-get-inputdata.py`
  - merges genotype feature tables with phenotype values and exports model input files.
- `3-run_pred_result.py`
  - runs multiple regression models with 5-fold cross-validation and saves prediction performance.
- `4-plot-DEC.py`
  - visualizes prediction results for DEC-related traits.
- `5-plot-PGWC.py`
  - visualizes prediction results for PGWC-related traits.

## Notes

- The scripts currently use hard-coded input and output paths.
- Please update all file paths before running them in a new environment.
- Required Python packages include `numpy`, `pandas`, `scikit-learn`, `scipy`, `lightgbm`, `catboost`, `cyvcf2`, and `mpire`.

