# Paramo_SkyIslands_RootFungi

This repository contains the data, scripts, and analysis documents supporting the manuscript:

> **Parallel elevation filtering of Ericaceae root-associated fungal communities in Andean páramo ecosystems** — *[Authors]* — *[Journal, Year]*

## Supplementary HTML

- Analysis Appendix Page URL: https://ssangulo.github.io/Paramo_SkyIslands_RootFungi/appendix.html
- Repository copy of rendered HTML: `Appendix_HTML.html`
- Quarto source: `Analysis_pipeline.qmd`

## Repository contents

### Bioinformatic pipeline

Scripts to process Illumina ITS2 reads from *Gaultheria myrsinoides* roots sampled across four Colombian páramo regions. The workflow includes read concatenation, demultiplexing, adapter/primer trimming, quality filtering and denoising with DADA2, chimera removal, LULU curation, subtraction of OTUs from controls, replicate filtering, and taxonomic assignment with the UNITE database. Final outputs are curated OTU tables and phyloseq objects for downstream ecological analyses.

### Supplementary analyses and reproducibility

| File | Description |
| ---- | ----------- |
| `Analysis_pipeline.qmd` | Quarto source file for the supplementary HTML document. |
| `Appendix_HTML.html` | **Main reproducibility document.** Rendered supplementary HTML with all analyses, figures, and tables. Open in any browser. |
| `Appendix_PDF.pdf` | PDF code-free version of the supplementary material for journal submission. |

### Data and outputs

| Folder | Contents |
| ------ | -------- |
| `Data_files/` | Sample metadata and supporting data tables |
| `figures/` | All figures (PNG) referenced in the analysis documents |
| `tables/` | All result tables (CSV) referenced in the analysis documents |
| `objects/` | Saved R objects (phyloseq, phylogenetic tree, GLLVM heatmap) |
| `Scripts/` | R scripts for the full bioinformatic and analytical pipeline |

## Reproducibility

All analyses can be reproduced by cloning this repository and rendering `Analysis_pipeline.qmd` from the repository root. Computationally intensive steps (GLLVM fitting, iNEXT3D, NRI/NTI) were executed on a remote HPC; their outputs are archived in `figures/`, `tables/`, and `objects/` and embedded as static results.

Sequence data are deposited at ENA under accession **PRJEB107725**.
