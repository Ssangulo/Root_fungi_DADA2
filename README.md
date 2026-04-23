# Root_fungi_DADA2

This repository contains the data, scripts, and analysis documents supporting the manuscript:

> **Parallel elevation filtering of Ericaceae root-associated fungal communities in Andean páramo ecosystems** — *[Authors]* — *[Journal, Year]*

## Supplementary HTML

- GitHub Pages URL (after enabling Pages): `https://<username>.github.io/Root_fungi_DADA2/appendix.html`
- Repository copy of rendered HTML: `docs/appendix.html`
- Quarto source: `Analysis_pipeline_HTMLformat - original.qmd`

The HTML appendix is designed for readability and navigation (table of contents, section anchors, collapsible code, and copyable code blocks).

## Repository contents

### Bioinformatic pipeline

Scripts to process Illumina ITS2 reads from *Gaultheria myrsinoides* roots sampled across four Colombian páramo regions. The workflow includes read concatenation, demultiplexing, adapter/primer trimming, quality filtering and denoising with DADA2, chimera removal, LULU curation, subtraction of OTUs from controls, replicate filtering, and taxonomic assignment with the UNITE database. Final outputs are curated OTU tables and phyloseq objects for downstream ecological analyses.

### Supplementary analyses and reproducibility

| File | Description |
| ---- | ----------- |
| `Analysis_pipeline_HTMLformat - original.qmd` | Quarto source file for the supplementary HTML document. |
| `docs/appendix.html` | **Main reproducibility document.** Rendered supplementary HTML with all analyses, figures, and tables. Open in any browser. |
| `Analysis_pipeline_PDF.qmd` | Quarto source file for the PDF appendix. Renders a clean, code-free version of the supplementary material for journal submission. |

### Data and outputs

| Folder | Contents |
| ------ | -------- |
| `Data_files/` | Sample metadata and supporting data tables |
| `figures/` | All figures (PNG) referenced in the analysis documents |
| `tables/` | All result tables (CSV) referenced in the analysis documents |
| `objects/` | Saved R objects (phyloseq, phylogenetic tree, GLLVM heatmap) |
| `Scripts/` | R scripts for the full bioinformatic and analytical pipeline |

## Reproducibility

All analyses can be reproduced by cloning this repository and rendering `Analysis_pipeline_HTMLformat - original.qmd` from the repository root. Computationally intensive steps (GLLVM fitting, iNEXT3D, NRI/NTI) were executed on a remote Ubuntu server; their outputs are archived in `figures/`, `tables/`, and `objects/` and embedded as static results.

Sequence data are deposited at ENA under accession **PRJEB107725**.
