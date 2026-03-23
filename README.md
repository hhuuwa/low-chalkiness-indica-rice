# Population Genomics Reveals Selection Signatures and Favorable Haplotypes for Low Chalkiness in Indica Rice

Code repository for the study **Population Genomics Reveals Selection Signatures and Favorable Haplotypes for Low Chalkiness in Indica Rice**.

This repository organizes the main analysis scripts used in the manuscript, from phenotype processing and variant discovery to population genomics, haplotype analysis, pedigree tracing, and genotype-to-phenotype prediction.


## Overview

Rice grain chalkiness is a key appearance and quality trait in indica rice. This repository collects the scripts used to:

- estimate BLUE values for chalkiness-related traits,
- call and annotate variants,
- infer population structure,
- quantify gene expression,
- identify favorable haplotypes and haplotype combinations,
- detect selection signatures,
- trace evolutionary origin and pedigree relationships, and
- evaluate machine-learning models for genotype-to-phenotype prediction.

## Repository structure

```text
.
|-- 1.Phenotyping of grain chalkiness
|-- 2.Variant discovery
|-- 3.Population structure analysis
|-- 4.Gene expression analysis
|-- 5.Haplotype analysis
|-- 6.Selection signals analysis
|-- 7.Evolutionary origin and pedigree tracing
|-- 8.Haplotype-combination analysis
`-- 9.Machine-learning_genotype_phenotype_prediction
```

## Analysis modules

### 1. Phenotyping of grain chalkiness

- `1.1.cala_BLUE.R`
  - fits a mixed-effects model and estimates BLUE values for traits with replicate measurements;
  - exports BLUE tables and histogram plots.
- `1.2.plot_pgwc_dec_boxplot.R`
  - compares PGWC and DEC between two subpopulations using boxplots and t-tests.

### 2. Variant discovery

- `run_variant_calling_snpeff.sh`
  - performs read mapping, BAM processing, per-chromosome GVCF calling, joint genotyping, SNP/INDEL filtering, and functional annotation with SnpEff.

### 3. Population structure analysis

- `3.1.run_popgen_tree_pca.sh`
  - filters SNPs, generates phylogeny input, builds a FastTree tree, and runs PLINK PCA.
- `3.2.plot_pca.R`
  - plots PCA results and writes the variance explained by each principal component.

### 4. Gene expression analysis

- `4.1.run_rnaseq_expression.sh`
  - runs `fastp`, `STAR`, and `featureCounts` for RNA-seq read processing and gene-level read counting.
- `4.2.calc_tpm.R`
  - calculates TPM values from the featureCounts output and GTF annotation.

### 5. Haplotype analysis

- `5.1.extract_promoter_cds_vcf.sh`
  - extracts promoter and CDS variants for candidate genes from a whole-genome VCF.
- `5.2.batch_chalk_haplotype.R`
  - performs gene-wise haplotype analysis, frequency comparison, phenotype comparison, and summary output.

### 6. Selection signals analysis

- `detect_highFST_genes.sh`
  - calculates windowed FST between Guangdong indica and international indica groups;
  - identifies top 10% high-FST regions and overlaps them with candidate genes.

### 7. Evolutionary origin and pedigree tracing

- `genehapr_network_temporal_pedigree.R`
  - builds haplotype networks;
  - summarizes temporal trends and maps pedigree-panel haplotypes to the expanded panel.

### 8. Haplotype-combination analysis

- `8.1.major_haplotype.R`
  - identifies major haplotypes and evaluates haplotype effects on chalkiness.
- `8.2.haplotype_combination.R`
  - combines major haplotypes across genes and compares phenotypic effects among combinations.

### 9. Machine-learning genotype-phenotype prediction

- `1-vcf2table.py`
  - converts VCF files into feature matrices.
- `2-get-inputdata.py`
  - merges genotype features with phenotype data to build model input tables.
- `3-run_pred_result.py`
  - evaluates multiple machine-learning regression models with 5-fold cross-validation.
- `4-plot-DEC.py`
  - plots prediction results for DEC.
- `5-plot-PGWC.py`
  - plots prediction results for PGWC.

## Recommended execution order

The analyses are modular, but a typical workflow is:

1. phenotype processing,
2. variant discovery,
3. population structure analysis,
4. expression analysis,
5. haplotype analysis,
6. selection signal analysis,
7. evolutionary/pedigree analysis,
8. haplotype-combination analysis,
9. machine-learning prediction.

Not every figure or table in the manuscript requires all modules to be run in one continuous pipeline.

## Software requirements

The repository contains scripts in `R`, `bash`, and `Python`. The exact runtime environment depends on the module, but the main tools used include:

### Command-line tools

- `bwa-mem2`
- `samtools`
- `GATK`
- `Picard`
- `SnpEff`
- `vcftools`
- `bcftools`
- `TASSEL5`
- `FastTree`
- `PLINK`
- `fastp`
- `STAR`
- `featureCounts`
- `bedtools`
- `bgzip`
- `tabix`

### R packages

- `argparser`
- `magrittr`
- `lsmeans`
- `lme4`
- `ggplot2`
- `ggpubr`
- `cowplot`
- `GenomicFeatures`
- `GenomicRanges`
- `geneHapR`
- `data.table`
- `dplyr`
- `tidyr`
- `ggsignif`
- `agricolae`

### Python packages

- `numpy`
- `pandas`
- `scikit-learn`
- `scipy`
- `lightgbm`
- `catboost`
- `cyvcf2`
- `mpire`
- `joblib`

## Important notes before running

### 1. Paths are currently environment-specific

Several scripts still contain hard-coded paths such as:

- `/project/data/...`
- `/data2/...`

Please update these paths to match your local or cluster environment before running the scripts.

### 2. Input files are not bundled here

This repository focuses on analysis code. Large raw datasets such as:

- FASTQ files,
- BAM files,
- genome references,
- full VCF files,
- expression matrices, and
- intermediate result files

are not included in this repository and should be prepared separately.

### 3. Some scripts are designed as reusable templates

Several scripts are written for specific datasets, sample lists, or candidate gene sets used in the manuscript. Before reuse, please check:

- input file names,
- sample metadata column names,
- subpopulation labels,
- output directories, and
- candidate gene lists.

### 4. Encoding

Some historical comments in shell scripts may show encoding artifacts. This does not affect the core logic, but cleaning comments and normalizing file encoding is recommended for long-term maintenance.

## Suggested data organization

For easier reuse, a project-specific directory structure like the following is recommended outside this repository:

```text
project/
|-- data/
|   |-- ref/
|   |-- fastq/
|   |-- vcf/
|   |-- phenotype/
|   `-- metadata/
|-- results/
|-- scripts/
`-- logs/
```

## Citation

If you use this code, please cite:

> Population Genomics Reveals Selection Signatures and Favorable Haplotypes for Low Chalkiness in Indica Rice

If the manuscript is published, please replace this section with the final journal citation, DOI, and year.

## Contact

For questions about the code or the manuscript, please contact the repository maintainer or corresponding author(s) associated with this study.

