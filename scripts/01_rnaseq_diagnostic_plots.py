# ==============================================================================
# RNA-seq Diagnostic Plots: A Practical Guide (Python)
# Dataset: GSE161731 — Blood Transcriptomics in Respiratory Infections
#
# This script generates key diagnostic plots used in bulk RNA-seq analysis,
# organized in workflow order from raw sequencing quality through outlier
# detection. Python equivalent of the R version in rnaseq-diagnostic-plots.
#
# Plots 01-02 use simulated data (require raw FASTQ/BAM files in practice).
# Plots 03-06 use real data from GSE161731 downloaded directly from GEO.
#
# Plots covered:
#   01. Per-base sequence quality    (simulated — FastQC)
#   02. Genomic feature assignment   (simulated — featureCounts)
#   03. Library size distribution    (real data)
#   04. Gene detection per sample    (real data)
#   05. RNA biotype distribution     (real data)
#   06. PCA — outlier detection      (real data)
#
# Requirements:
#   pip install pandas numpy matplotlib seaborn scikit-learn mygene requests
# ==============================================================================

import os
import urllib.request
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import mygene
import warnings
warnings.filterwarnings("ignore")

## Set working directory to project root (one level up from scripts/)
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(project_root)

os.makedirs("figures", exist_ok=True)
os.makedirs("data",    exist_ok=True)

# Color palette — matches R version
PALETTE = {
    "healthy"  : "#2ecc71",
    "COVID19"  : "#e74c3c",
    "Influenza": "#3498db",
    "Bacterial": "#e67e22",
    "CoV_other": "#9b59b6"
}

print("Setup complete\n")


# ==============================================================================
# DOWNLOAD DATA FROM GEO
# ==============================================================================

COUNT_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731"
             "/suppl/GSE161731_counts.csv.gz")
KEY_URL   = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731"
             "/suppl/GSE161731_counts_key.csv.gz")

import ssl

def download_if_missing(url, path):
    if not os.path.exists(path):
        print(f"Downloading {os.path.basename(path)}...")
        ssl_context = ssl._create_unverified_context()
        with urllib.request.urlopen(url, context=ssl_context) as response:
            with open(path, "wb") as f:
                f.write(response.read())
        print("  Done")
    else:
        print(f"  {os.path.basename(path)} already exists, skipping download")

download_if_missing(COUNT_URL, "data/GSE161731_counts.csv.gz")
download_if_missing(KEY_URL,   "data/GSE161731_counts_key.csv.gz")

print("\nLoading data...")
counts_raw = pd.read_csv("data/GSE161731_counts.csv.gz", index_col=0)
key        = pd.read_csv("data/GSE161731_counts_key.csv.gz", dtype=str)

# Fix sample ID mismatch — R make.names() replaces - with . and adds X prefix
# Python equivalent: replaces - with . and prepend X if starts with digit
def make_names(s):
    s = str(s).replace("-", ".").replace(" ", ".")
    if s[0].isdigit():
        s = "X" + s
    return s

key["rna_id_fixed"] = key["rna_id"].apply(make_names)
key = key.set_index("rna_id_fixed")

counts_raw.columns = [make_names(c) for c in counts_raw.columns]
overlap = list(set(counts_raw.columns) & set(key.index))
counts  = counts_raw[overlap].copy()
meta    = key.loc[overlap].copy()
overlap = list(set(counts_raw.columns) & set(key.index))
print(f"Overlap: {len(overlap)} samples")
print(f"counts_raw columns sample: {list(counts_raw.columns[:3])}")
print(f"key index sample: {list(key.index[:3])}")
counts  = counts_raw[overlap].copy()
meta    = key.loc[overlap].copy()
# Clean cohort names
meta["cohort"] = meta["cohort"].replace({
    "COVID-19" : "COVID19",
    "CoV other": "CoV_other"
})
meta["age"]   = pd.to_numeric(meta["age"],   errors="coerce")
meta["batch"] = pd.to_numeric(meta["batch"], errors="coerce")

print(f"Counts: {counts.shape[0]:,} genes x {counts.shape[1]} samples")
print(f"Metadata: {meta.shape[0]} samples")

# Remove all-zero genes
counts_filtered = counts.loc[counts.sum(axis=1) > 0]
print(f"After removing zero-count genes: {counts_filtered.shape[0]:,} genes")


# ==============================================================================
# ==============================================================================
# QC — REMOVE OUTLIERS AND APPLY GENE FILTER
# ==============================================================================

def log2_cpm(counts_df, prior=1):
    """Compute log2(CPM + prior) normalization."""
    lib_sizes = counts_df.sum(axis=0)
    cpm = counts_df.divide(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + prior)

log_cpm        = log2_cpm(counts_filtered)
log_cpm_matrix = log_cpm.values.astype(float)

pca        = PCA(n_components=2)
pca_coords = pca.fit_transform(log_cpm_matrix.T)
pct_var    = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame({
    "PC1"   : pca_coords[:, 0],
    "PC2"   : pca_coords[:, 1],
    "sample": counts_filtered.columns.tolist(),
    "cohort": meta.loc[counts_filtered.columns, "cohort"].values
})

for pc in ["PC1", "PC2"]:
    mean, sd = pca_df[pc].mean(), pca_df[pc].std()
    pca_df[f"{pc}_outlier"] = (pca_df[pc] - mean).abs() > 3 * sd

pca_df["outlier"] = pca_df["PC1_outlier"] | pca_df["PC2_outlier"]
outlier_samples   = pca_df.loc[pca_df["outlier"], "sample"].tolist()
print(f"Outliers identified: {len(outlier_samples)}: {outlier_samples}")

# Step 1: Remove outliers
counts_qc = counts_filtered.drop(columns=outlier_samples)
meta_qc   = meta.loc[counts_qc.columns].copy()
print(f"After outlier removal: {meta_qc.shape[0]} samples")

# Step 2: Remove missing gender and age
meta_qc   = meta_qc.dropna(subset=["gender", "age"])
meta_qc   = meta_qc[meta_qc["gender"] != ""]
counts_qc = counts_qc[meta_qc.index]
print(f"After gender/age filter: {meta_qc.shape[0]} samples")

# Step 3: Keep earliest timepoint per subject
# Keep rna_id_fixed as a column before groupby so we can restore it as index
time_map          = {"early": 0, "middle": 1, "late": 2}
meta_qc = meta_qc.copy()
meta_qc["rna_id_fixed"] = meta_qc.index  # save index as column explicitly
meta_qc = meta_qc.reset_index(drop=True)
meta_qc["time_rank"] = meta_qc["time_since_onset"].map(time_map).fillna(99)
meta_qc = (meta_qc
           .sort_values("time_rank")
           .groupby("subject_id", sort=False)
           .first()
           .reset_index())
meta_qc = meta_qc.set_index("rna_id_fixed")
print(f"After deduplication: {meta_qc.shape[0]} samples")
print(f"meta_qc index sample: {list(meta_qc.index[:3])}")

# Step 4: Re-align counts to QC metadata
shared    = [s for s in meta_qc.index if s in counts_qc.columns]
meta_qc   = meta_qc.loc[shared]
counts_qc = counts_qc[shared]
print(f"Shared samples: {len(shared)}")

# Step 5: Gene filtering
min_samples = int(0.1 * counts_qc.shape[1])
counts_qc   = counts_qc.loc[(counts_qc >= 10).sum(axis=1) >= min_samples]

cohort_order          = ["healthy", "COVID19", "Influenza", "Bacterial", "CoV_other"]
meta_qc["cohort"]     = pd.Categorical(meta_qc["cohort"], categories=cohort_order)

print(f"After QC: {meta_qc.shape[0]} samples, {counts_qc.shape[0]:,} genes")
print(f"counts_qc columns sample: {list(counts_qc.columns[:3])}\n")

# ==============================================================================
# PLOT 01: PER-BASE SEQUENCE QUALITY (SIMULATED)
# ==============================================================================
#
# WHAT IT SHOWS:
#   Phred quality scores at each position along the sequencing read.
#   Q = -10 * log10(error probability). Q30 = 1 error per 1,000 bases.
#   Generated by FastQC; aggregated across samples with MultiQC.
#
#   Simulated here because this plot requires raw FASTQ files.
#
# WHAT GOOD LOOKS LIKE:
#   Scores above Q30 across all positions. Slight decline at 3' end is normal.
#
# WHAT BAD LOOKS LIKE:
#   Sharp quality drop after position 100 (3' degradation) → trim with
#   Trimmomatic or Cutadapt. Low quality throughout (<Q20) → run may have failed.

np.random.seed(42)
positions = np.arange(1, 151)

good_mean = 36 - (positions / 150) * 4
bad_mean  = np.where(
    positions <= 100,
    35 - (positions / 100) * 5,
    30 - ((positions - 100) / 50) * 18
)

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
fig.suptitle("01. Per-Base Sequence Quality (Simulated)", fontsize=13, fontweight="bold")

for ax, mean, q_label, color, title in zip(
    axes,
    [good_mean, bad_mean],
    ["Good quality", "Poor quality (3' degradation)"],
    ["#27ae60", "#c0392b"],
    ["Good quality", "Poor quality (3' degradation)"]
):
    # Quality zones
    ax.axhspan(28, 42, color="#2ecc71", alpha=0.10)
    ax.axhspan(20, 28, color="#f39c12", alpha=0.10)
    ax.axhspan(0,  20, color="#e74c3c", alpha=0.10)

    # IQR ribbon and mean line
    ax.fill_between(positions, mean - 2, mean + 1.5, color=color, alpha=0.25)
    ax.plot(positions, mean, color=color, linewidth=2)

    # Reference lines
    ax.axhline(30, linestyle="--", color="grey", linewidth=0.8)
    ax.axhline(20, linestyle="--", color="grey", linewidth=0.8)

    # Zone labels
    ax.text(148, 38, "Good (Q≥28)",        color="#27ae60", ha="right", fontsize=9, fontweight="bold")
    ax.text(148, 24, "Acceptable (Q20-28)", color="#e67e22", ha="right", fontsize=9, fontweight="bold")
    ax.text(148, 10, "Poor (Q<20)",         color="#c0392b", ha="right", fontsize=9, fontweight="bold")

    ax.set_title(title, fontweight="bold", fontsize=11)
    ax.set_xlabel("Position in Read (bp)")
    ax.set_xlim(0, 151)
    ax.set_ylim(0, 42)
    ax.set_yticks(range(0, 41, 10))
    ax.set_xticks(range(0, 151, 25))
    ax.grid(True, alpha=0.3)

axes[0].set_ylabel("Phred Quality Score")
fig.text(0.5, 0.01,
         "Simulated 150bp reads | Line = mean Phred score | Ribbon = IQR\n"
         "In practice: generate with FastQC, aggregate with MultiQC",
         ha="center", fontsize=9, color="grey")

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig("figures/01_per_base_quality_simulated.png", dpi=300, bbox_inches="tight")
plt.close()
print("Plot 01 saved: per-base sequence quality")


# ==============================================================================
# PLOT 02: GENOMIC FEATURE ASSIGNMENT (SIMULATED)
# ==============================================================================
#
# WHAT IT SHOWS:
#   Proportion of reads assigned to exons, introns, intergenic regions, rRNA.
#   Generated by featureCounts (Subread) or RSeQC — requires BAM files.
#
# WHAT GOOD LOOKS LIKE:
#   >60% exonic. <20% intronic. <10% intergenic. Low rRNA (<5%).
#
# WHAT BAD LOOKS LIKE:
#   High intronic (>40%) → DNA contamination.
#   High intergenic (>30%) → degraded RNA or wrong genome.
#   High rRNA (>20%) → rRNA depletion failed.

scenarios = {
    "Good\n(typical blood RNA-seq)":     {"Exonic":72, "Intronic":18, "Intergenic":5,  "rRNA":2,  "Unassigned":3},
    "High Intronic\n(DNA contamination)":{"Exonic":38, "Intronic":48, "Intergenic":7,  "rRNA":2,  "Unassigned":5},
    "High Intergenic\n(degraded RNA)":   {"Exonic":35, "Intronic":20, "Intergenic":38, "rRNA":3,  "Unassigned":4},
    "High rRNA\n(failed rRNA depletion)":{"Exonic":30, "Intronic":10, "Intergenic":5,  "rRNA":52, "Unassigned":3},
}

feature_colors = {
    "Exonic"    : "#2ecc71",
    "Intronic"  : "#3498db",
    "Intergenic": "#e67e22",
    "rRNA"      : "#e74c3c",
    "Unassigned": "#95a5a6"
}

features = ["Exonic", "Intronic", "Intergenic", "rRNA", "Unassigned"]
x        = np.arange(len(scenarios))

fig, ax = plt.subplots(figsize=(12, 6))
fig.suptitle("02. Genomic Feature Assignment (Simulated)",
             fontsize=13, fontweight="bold")

bottoms = np.zeros(len(scenarios))
for feature in features:
    vals = [scenarios[s][feature] for s in scenarios]
    bars = ax.bar(x, vals, bottom=bottoms,
                  color=feature_colors[feature], label=feature, width=0.6)
    # Label each segment
    for i, (v, b) in enumerate(zip(vals, bottoms)):
        if v > 3:
            ax.text(i, b + v/2, f"{v}%",
                    ha="center", va="center",
                    color="white", fontsize=9, fontweight="bold")
    bottoms += np.array(vals)

ax.axhline(60, linestyle="--", color="black", linewidth=1)
ax.text(-0.4, 62, "60% exonic threshold", fontsize=9)
ax.set_xticks(x)
ax.set_xticklabels(scenarios.keys(), fontsize=10, fontweight="bold")
ax.set_ylabel("Percentage of Reads (%)")
ax.set_ylim(0, 110)
ax.legend(loc="upper right", fontsize=9)
ax.set_title("Simulated scenarios showing common failure modes\n"
             "In practice: generate with featureCounts (Subread) or RSeQC",
             fontsize=9, color="grey")
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("figures/02_genomic_feature_assignment_simulated.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("Plot 02 saved: genomic feature assignment")


# ==============================================================================
# PLOT 03: LIBRARY SIZE DISTRIBUTION
# ==============================================================================
#
# WHAT IT SHOWS:
#   Total mapped reads per sample ordered from lowest to highest.
#
# WHAT GOOD LOOKS LIKE:
#   Consistent depth, median 20-80M reads. No samples below 20% of median.
#
# WHAT BAD LOOKS LIKE:
#   Samples below red threshold → failed sequencing or poor RNA quality.
#   Bimodal distribution → batch effects.

lib_sizes = counts_filtered.sum(axis=0) / 1e6
lib_df    = pd.DataFrame({
    "sample"  : lib_sizes.index,
    "lib_size": lib_sizes.values,
    "cohort"  : meta.loc[lib_sizes.index, "cohort"].values
}).sort_values("lib_size").reset_index(drop=True)

median_lib    = lib_df["lib_size"].median()
low_threshold = median_lib * 0.2

fig, ax = plt.subplots(figsize=(14, 5))
colors = [PALETTE.get(c, "#aaaaaa") for c in lib_df["cohort"]]
ax.bar(range(len(lib_df)), lib_df["lib_size"], color=colors, width=1.0)
ax.axhline(median_lib,    linestyle="--", color="black", linewidth=1.0,
           label=f"Median: {median_lib:.1f}M")
ax.axhline(low_threshold, linestyle=":",  color="red",   linewidth=1.0,
           label="Low threshold (20% of median)")
ax.text(5, median_lib + 1.5, f"Median: {median_lib:.1f}M", fontsize=9)
ax.text(5, low_threshold + 1.5, "Low threshold (20% of median)",
        fontsize=8, color="red")

# Legend for cohorts
handles = [mpatches.Patch(color=PALETTE[c], label=c) for c in cohort_order
           if c in lib_df["cohort"].values]
ax.legend(handles=handles, loc="upper left", fontsize=9, title="Cohort")

ax.set_title("03. Library Size Distribution", fontweight="bold", fontsize=13)
ax.set_xlabel(f"Sample (ordered by library size) | n={len(lib_df)}")
ax.set_ylabel("Library Size (millions of reads)")
ax.set_xlim(-1, len(lib_df))
ax.set_xticks([])
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("figures/03_library_size.png", dpi=300, bbox_inches="tight")
plt.close()
print("Plot 03 saved: library size distribution")


# ==============================================================================
# PLOT 04: GENE DETECTION PER SAMPLE
# ==============================================================================
#
# WHAT IT SHOWS:
#   Number of genes with count > 0 per sample, grouped by cohort.
#
# WHAT GOOD LOOKS LIKE:
#   12,000-18,000 genes detected. Consistent within cohorts.
#
# WHAT BAD LOOKS LIKE:
#   Samples below 5,000 genes → poor RNA quality or failed extraction.

n_genes = (counts_filtered > 0).sum(axis=0)
gene_df = pd.DataFrame({
    "sample" : n_genes.index,
    "n_genes": n_genes.values,
    "cohort" : meta.loc[n_genes.index, "cohort"].values
})
gene_df["cohort"] = pd.Categorical(gene_df["cohort"], categories=cohort_order)
gene_df = gene_df.sort_values("cohort")

fig, ax = plt.subplots(figsize=(9, 5))

for i, cohort in enumerate(cohort_order):
    subset = gene_df[gene_df["cohort"] == cohort]["n_genes"]
    if len(subset) == 0:
        continue
    x_pos = np.random.normal(i, 0.08, size=len(subset))
    ax.scatter(x_pos, subset, color=PALETTE[cohort],
               alpha=0.6, s=40, zorder=3)
    bp = ax.boxplot(subset, positions=[i], widths=0.4,
                    patch_artist=True, zorder=2,
                    medianprops=dict(color="black", linewidth=2))
    bp["boxes"][0].set_facecolor(PALETTE[cohort])
    bp["boxes"][0].set_alpha(0.3)

ax.axhline(5000, linestyle="--", color="red", linewidth=1)
ax.text(-0.4, 5300, "5,000 gene minimum", fontsize=8, color="red")
ax.set_xticks(range(len(cohort_order)))
ax.set_xticklabels(cohort_order, fontsize=10)
ax.set_title("04. Gene Detection per Sample", fontweight="bold", fontsize=13)
ax.set_xlabel("")
ax.set_ylabel("Number of Genes Detected (count > 0)")
ax.set_subtitle = ax.set_title
ax.text(0.5, 1.02,
        "GSE161731 | Samples with fewer than 5,000 genes detected should be investigated",
        transform=ax.transAxes, ha="center", fontsize=9, color="grey")
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("figures/04_gene_detection.png", dpi=300, bbox_inches="tight")
plt.close()
print("Plot 04 saved: gene detection per sample")


# ==============================================================================
# PLOT 05: RNA BIOTYPE DISTRIBUTION
# ==============================================================================
#
# WHAT IT SHOWS:
#   Proportion of reads mapping to protein-coding, ncRNA, pseudogene, rRNA etc.
#   Uses mygene to annotate Ensembl IDs with gene types.
#
# WHAT GOOD LOOKS LIKE:
#   >70% protein-coding. Consistent across samples. Low rRNA (<5%).
#
# WHAT BAD LOOKS LIKE:
#   High rRNA → rRNA depletion failed. Unusual sample profile → mislabeled.

print("\nFetching biotype annotations from mygene...")
mg          = mygene.MyGeneInfo()
ensembl_ids = counts_qc.index.tolist()

# Query in batches of 1000
batch_size = 1000
results    = []
for i in range(0, len(ensembl_ids), batch_size):
    batch = ensembl_ids[i:i + batch_size]
    res   = mg.querymany(batch, scopes="ensembl.gene",
                         fields="type_of_gene", species="human",
                         returnall=False, verbose=False)
    results.extend(res)

biotype_dict = {}
for r in results:
    if "query" in r and "type_of_gene" in r:
        biotype_dict[r["query"]] = r["type_of_gene"]

def simplify_biotype(b):
    if b == "protein-coding": return "Protein-coding"
    if b == "pseudo":         return "Pseudogene"
    if b == "ncRNA":          return "ncRNA"
    if b in ("snRNA","snoRNA"): return "snRNA/snoRNA"
    if b == "rRNA":           return "rRNA"
    if b == "tRNA":           return "tRNA"
    if b is None:             return "Unknown"
    return "Other"

biotype_series = pd.Series(
    {g: simplify_biotype(biotype_dict.get(g)) for g in ensembl_ids}
)
print(f"Biotype annotations retrieved: {len(biotype_dict)}")
print(f"Sample biotypes: {list(biotype_dict.items())[:5]}")
# Sum counts per biotype per sample
counts_qc_bt = counts_qc.copy()
counts_qc_bt["biotype"] = biotype_series
biotype_sums = counts_qc_bt.groupby("biotype").sum()
print("meta_qc index sample:", list(meta_qc.index[:3]))
print("biotype_sums columns sample:", list(biotype_sums.columns[:3]))
biotype_pct  = biotype_sums.divide(biotype_sums.sum(axis=0), axis=1) * 100

# Order samples by cohort
sample_order_qc = (meta_qc.sort_values("cohort").index.tolist())
sample_order_qc = [s for s in sample_order_qc if s in biotype_pct.columns]
biotype_pct     = biotype_pct[sample_order_qc]
print("biotype_pct shape:", biotype_pct.shape)
print("biotype_pct columns sample:", list(biotype_pct.columns[:3]))
print("sample_order_qc sample:", sample_order_qc[:3])
print("matching samples:", len([s for s in sample_order_qc if s in biotype_pct.columns]))
biotype_colors = {
    "Protein-coding": "#2ecc71",
    "ncRNA"         : "#3498db",
    "Pseudogene"    : "#e67e22",
    "snRNA/snoRNA"  : "#9b59b6",
    "rRNA"          : "#e74c3c",
    "tRNA"          : "#f39c12",
    "Other"         : "#95a5a6",
    "Unknown"       : "#bdc3c7"
}

biotype_order_plot = ["Protein-coding","ncRNA","Pseudogene",
                      "snRNA/snoRNA","rRNA","tRNA","Other","Unknown"]
biotype_order_plot = [b for b in biotype_order_plot if b in biotype_pct.index]

fig, ax = plt.subplots(figsize=(14, 6))
bottoms = np.zeros(biotype_pct.shape[1])

for biotype in biotype_order_plot:
    if biotype not in biotype_pct.index:
        continue
    vals = biotype_pct.loc[biotype].values
    ax.bar(range(len(sample_order_qc)), vals, bottom=bottoms,
           color=biotype_colors.get(biotype, "#aaaaaa"),
           label=biotype, width=1.0)
    bottoms += vals
print(f"  {biotype}: vals range {vals.min():.1f} - {vals.max():.1f}")
ax.axhline(70, linestyle="--", color="black", linewidth=1)
ax.text(3, 72, "70% protein-coding threshold", fontsize=8)
ax.set_xticks([])
ax.set_xlabel(f"Sample (ordered by cohort) | n={len(sample_order_qc)}")
ax.set_ylabel("Percentage of Reads (%)")
ax.set_title(
    f"05. RNA Biotype Distribution\n"
    f"GSE161731 | {len(sample_order_qc)} samples | "
    "More than 70% protein-coding expected for blood RNA-seq",
    fontsize=11, fontweight="bold"
)
ax.legend(loc="lower right", fontsize=8, ncol=2)
ax.set_ylim(0, 105)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("figures/05_biotype_distribution.png", dpi=300, bbox_inches="tight")
plt.close()
print("Plot 05 saved: biotype distribution")


# ==============================================================================
# PLOT 06: PCA — OUTLIER DETECTION
# ==============================================================================
#
# WHAT IT SHOWS:
#   PCA of all samples using log2(CPM+1), run BEFORE removing outliers.
#   Samples >3 SD from mean on PC1 or PC2 flagged in red.
#
# WHAT GOOD LOOKS LIKE:
#   Samples cluster loosely by cohort. No isolated outliers.
#
# WHAT BAD LOOKS LIKE:
#   Isolated samples far from the main cluster → investigate and likely remove.
#   All samples from one batch cluster separately → batch correction needed.

fig, ax = plt.subplots(figsize=(10, 7))

for cohort in cohort_order:
    subset = pca_df[pca_df["cohort"] == cohort]
    ax.scatter(subset["PC1"], subset["PC2"],
               color=PALETTE.get(cohort, "#aaaaaa"),
               s=60, alpha=0.8, label=cohort, zorder=2)

# Highlight outliers
outliers = pca_df[pca_df["outlier"]]
ax.scatter(outliers["PC1"], outliers["PC2"],
           facecolors="none", edgecolors="red",
           s=200, linewidths=2, zorder=3)

# Label outliers
for _, row in outliers.iterrows():
    ax.annotate(row["sample"],
                xy=(row["PC1"], row["PC2"]),
                xytext=(10, 5), textcoords="offset points",
                fontsize=8, color="red", fontweight="bold")

ax.set_xlabel(f"PC1 ({pct_var[0]:.0f}%)", fontsize=11)
ax.set_ylabel(f"PC2 ({pct_var[1]:.0f}%)", fontsize=11)
ax.set_title(
    "06. PCA — Outlier Detection\n"
    f"GSE161731 | log2(CPM+1) | PC1={pct_var[0]:.0f}%, PC2={pct_var[1]:.0f}% | "
    "Red circles = potential outliers (>3 SD from mean)",
    fontsize=11, fontweight="bold"
)
ax.legend(fontsize=9, title="Cohort")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figures/06_pca_outlier_detection.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Plot 06 saved: PCA outlier detection | Outliers flagged: {pca_df['outlier'].sum()}")


print("\n=== All 6 plots saved to figures/ ===")
