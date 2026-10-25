# Single-Cell RNA-Seq Analysis in Gastric Cancer (SPP1+ / C1QC+ Macrophages)

## Overview
This project analyzes single-cell RNA sequencing (scRNA-Seq) data from two gastric cancer samples
(**GSM5573504 / sample 39** and **GSM5573499 / sample 34**) obtained from the Gene Expression Omnibus (GEO).
The goal is to characterize the cellular composition of the tumor microenvironment, with a focus on
**macrophages** and **SPP1 signaling** toward T cells. The pipeline is implemented in **R** using the **Seurat** package.

---

## Key Features
- **Quality Control and Doublet Removal**: mitochondrial/feature filtering and doublet detection with scDblFinder.
- **Sample Integration**: Harmony-based integration of the two samples.
- **Cell Type Annotation**: automated reference-based annotation with SingleR (Human Primary Cell Atlas).
- **Differential Expression Analysis**: macrophages in sample CIN_39 vs sample CIN_34, with a volcano plot.
- **Functional Enrichment**: GO (BP, MF, CC) and KEGG pathway analysis with clusterProfiler.
- **Cell-Cell Communication**: CellChat analysis, including the SPP1 signaling pathway (macrophage to T cells).

---

## Project Goals
1. Characterize the cellular composition of gastric cancer samples.
2. Identify differentially expressed genes (DEGs) in macrophages between the two samples.
3. Explore enriched pathways and cell-cell communication, especially SPP1 signaling.

---

## Requirements

### Software and Libraries
- R (>= 4.0.0)
- Seurat, dplyr, ggplot2, ggrepel
- scDblFinder, SingleCellExperiment
- harmony
- SingleR, celldex
- miloR
- clusterProfiler, org.Hs.eg.db, enrichplot, msigdbr
- CellChat

### Data Source
The project uses two scRNA-Seq samples from GEO: **GSM5573504** (sample 39) and **GSM5573499** (sample 34).
Download them and set the data path at the top of the script.

---

## Pipeline Overview
1. **Data Loading**: read 10X-format data for each sample with `Read10X`.
2. **Quality Control**: filter cells (200 < nFeature_RNA < 6000, percent.mt < 10).
3. **Doublet Removal**: detect and remove doublets with scDblFinder.
4. **Ambient RNA**: CellBender workflow explained in the script (run outside R, not executed here).
5. **Normalization and Feature Selection**: LogNormalize, 2000 highly variable genes, scaling.
6. **Dimensionality Reduction and Clustering**: PCA, neighbor graph, clustering, UMAP.
7. **Annotation**: SingleR on each sample, then on the integrated object.
8. **Integration**: merge the samples and correct batch effects with Harmony; re-cluster at higher resolution.
9. **Marker Genes**: `FindAllMarkers`, DotPlots, and macrophage marker expression (CD68, CD14, CD163, C1QA, C1QB, ...).
10. **Differential Expression**: macrophage CIN_39 vs CIN_34 (adjusted p < 0.05, |log2FC| > 0.25) and volcano plot.
11. **Functional Analysis**: GO and KEGG enrichment for up- and down-regulated genes.
12. **Cell-Cell Communication**: CellChat, SPP1 pathway, macrophage to T cell interactions.

---

## Outputs
- `merged_CIN_39_34_integrated.rds`: integrated Seurat object
- `Macrophage_CIN39_vs_CIN34_DEG_results.csv`: full DE results
- `UPregulated_genes.csv`, `DOWNregulated_genes.csv`
- `GO_BP_*.csv`, `GO_MF_UPregulated.csv`, `GO_CC_UPregulated.csv`
- `KEGG_UPregulated.csv`, `KEGG_DOWNregulated.csv`

<!-- Add your figures here after uploading them to the figures/ folder, for example:
![UMAP](figures/umap.png)
![Volcano plot](figures/volcano.png)
-->

---

## Usage
1. Clone this repository:

```bash
git clone https://github.com/bassemtech2005-hub/Single-Cell-RNA-Seq-Analysis-in-Gastric-Cancer-SPP1-C1QC-.git
```

2. Download the two GEO samples and set the data paths in the script.
3. Install the required R packages.
4. Run `scrna_analysis.R`.

---

## Notes and Limitations
- The analysis uses **one sample per group**, so the differential expression results are exploratory
  (no biological replicates).
- The data and the large `.rds` object are not included in this repository.

---

## Author
Bassem
