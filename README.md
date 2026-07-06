# RNA-seq Diagnostic Plots in Python

A Python implementation of key bulk RNA-seq quality control diagnostic plots,
organized in workflow order from raw sequencing quality through outlier
detection. 
Python equivalent of the R version in
[rnaseq-diagnostic-plots](https://github.com/pdeosthale/rnaseq-diagnostic-plots).

## Biological Question

How do you know if your RNA-seq data is high quality before running
differential expression analysis? This repo provides a six-plot QC workflow
that catches the most common problems — failed sequencing runs, contamination,
outlier samples, and biotype composition issues — before they affect downstream
results.

## Disclaimer

This is an independent secondary analysis of publicly available data from GEO
(GSE161731) for educational and portfolio purposes. The original study was
designed and conducted by McClain et al. (2021).

## Dataset

**GSE161731** — Blood RNA-seq from patients with respiratory infections and
healthy controls. McClain et al. *Nature Communications* 12, 1079 (2021).
https://doi.org/10.1038/s41467-021-21289-y

Data downloads automatically from GEO on first run — no manual download
required.

## Plots Covered

### Stage 1 — Pre-alignment QC (raw FASTQ)

These plots require raw FASTQ or BAM files in practice. Simulated data illustrates typical patterns. Generate with FastQC and MultiQC.

| Plot | Data | What it shows | Tool (in practice) |
  |---|---|---|---|
  | 01. Per-base sequence quality | Simulated | Phred quality scores across 150bp read positions — good vs 3' degradation | FastQC + MultiQC |
| 02. Genomic feature assignment | Simulated | Exon/intron/intergenic/rRNA fractions — good vs four failure modes | featureCounts (Subread) |

### Stage 2 — Count-level QC (after alignment and counting)

| Plot | Data | What it shows |
|---|---|---|
| 03. Library size distribution | Real | Total mapped reads per sample — flags low-depth samples |
| 04. Gene detection per sample | Real | Number of genes detected per sample |
| 05. RNA biotype distribution | Real | Proportion of reads mapping to protein-coding, ncRNA, rRNA, etc. |

### Stage 3 — Sample-level QC (outlier detection)

| Plot | Data | What it shows |
|---|---|---|
| 06. PCA — outlier detection | Real | Identifies extreme samples more than 3 SD from mean on PC1 or PC2 |

## Key Interpretations

**Per-base sequence quality (Plot 01):** Good: Phred scores consistently above
Q30 with a slight acceptable decline at the 3' end. Bad: sharp drop below Q20
after position 100 — trim with Trimmomatic or Cutadapt before alignment.

**Genomic feature assignment (Plot 02):** Good: more than 60% exonic reads.
High intronic (more than 40%) suggests DNA contamination. High rRNA (more than
                                                                     20%) indicates failed rRNA depletion.

**Library size distribution (Plot 03):** Good: consistent depth, median 20-80M
reads. Flag samples below 20% of the median (red line). Cross-reference with
gene detection — low library size AND low gene detection = failed sample.

**Gene detection per sample (Plot 04):** Good: 12,000-18,000 genes detected per
sample, consistent within cohorts. Samples below 5,000 genes should be
investigated.

**RNA biotype distribution (Plot 05):** Good: more than 70% protein-coding
(dashed line), consistent across all samples. High rRNA suggests failed rRNA
depletion.

**PCA outlier detection (Plot 06):** Good: samples cluster loosely by biological
group. Red circles flag samples more than 3 SD from the mean on PC1 or PC2.
Cross-reference with library size and gene detection before removing.

## R vs Python Comparison

This repo intentionally mirrors the R version
([rnaseq-diagnostic-plots](https://github.com/pdeosthale/rnaseq-diagnostic-plots))
to demonstrate equivalent workflows in both languages.

| Step | R | Python |
  |---|---|---|
  | Data loading | `read.csv` | `pandas.read_csv` |
  | Normalization | `edgeR::cpm(log=TRUE)` | Custom `log2_cpm()` function |
  | PCA | `prcomp()` | `sklearn.decomposition.PCA` |
  | Biotype annotation | `org.Hs.eg.db` | `mygene` API |
  | Plotting | `ggplot2` + `pheatmap` | `matplotlib` + `seaborn` |
  
  **One difference:** Python PCA uses log2-CPM without TMM normalization, while
the R version uses TMM-normalized CPM from edgeR. This causes a minor
difference in outlier detection — the R version identifies 3 outliers
(X94478, X434482, X434741) while Python identifies 2 (X94478, X434741).
X434482 falls just inside the 3 SD threshold without TMM normalization.

## Requirements

```bash
pip install pandas numpy matplotlib seaborn scikit-learn mygene requests
```

Python 3.9 or later recommended.

## How to Reproduce

```bash
git clone https://github.com/pdeosthale/rnaseq-qc-python.git
cd rnaseq-qc-python
pip install pandas numpy matplotlib seaborn scikit-learn mygene requests
python3 scripts/01_rnaseq_diagnostic_plots.py
```

Data downloads automatically on first run (~500MB). All 6 figures saved to
`figures/`.

## Project Structure

```
scripts/
  └── 01_rnaseq_diagnostic_plots.py   # All 6 plots, self-contained

figures/
  ├── 01_per_base_quality_simulated.png
├── 02_genomic_feature_assignment_simulated.png
├── 03_library_size.png
├── 04_gene_detection.png
├── 05_biotype_distribution.png
└── 06_pca_outlier_detection.png

data/                               # Downloaded from GEO on first run (not tracked)
  ```

## Tools

pandas, numpy, matplotlib, seaborn, scikit-learn, mygene, urllib, ssl

**External tools (for pre-alignment QC, not in script):** FastQC, MultiQC,
featureCounts (Subread)

## Limitations

- Plots 01 and 02 use simulated data. In practice these require raw FASTQ
files (FastQC/MultiQC) and alignment BAM files (featureCounts/RSeQC).
- Biotype annotations use the mygene API which requires an internet connection.
The API returns broad gene type categories; more granular annotations are
available from Ensembl BioMart.
- PCA uses log2-CPM without TMM normalization — minor difference from the R
version which uses TMM-normalized CPM.
- QC thresholds (20% median library size, 5,000 gene minimum, 70%
                   protein-coding, 3 SD for PCA) are appropriate for blood RNA-seq and may
need adjustment for other tissue types.

## References

McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional
responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079
(2021). https://doi.org/10.1038/s41467-021-21289-y

Andrews, S. FastQC: A quality control tool for high throughput sequence data.
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  
  Liao, Y., Smyth, G.K. & Shi, W. featureCounts: an efficient general purpose
program for assigning sequence reads to genomic features. *Bioinformatics* 30,
923-930 (2014). https://doi.org/10.1093/bioinformatics/btt656