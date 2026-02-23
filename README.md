# Root_fungi_DADA2

This repository contains the data, scripts, and analysis documents supporting the manuscript:

> **[Manuscript title]** — *[Authors]* — *[Journal, Year]*

## Repository contents

### Bioinformatic pipeline
Scripts to process Illumina ITS2 reads from *Gaultheria myrsinoides* roots sampled across four Colombian páramo regions. The workflow includes read concatenation, demultiplexing, adapter/primer trimming, quality filtering and denoising with DADA2, chimera removal, LULU curation, subtraction of OTUs from controls, replicate filtering, and taxonomic assignment with the UNITE database. Final outputs are curated OTU tables and phyloseq objects for downstream ecological analyses.

### Supplementary analyses and reproducibility

| File | Description |
|------|-------------|
| `Analysis_pipeline_HTMLformat.html` | **Main reproducibility document.** Self-contained HTML file with all statistical analyses, figures, and tables, including collapsible code blocks for all methods executed on the remote server. Open in any browser — no additional software required. |
| `Analysis_pipeline_HTMLformat.qmd` | Quarto source file for the HTML document. Requires R and the packages listed in the setup chunk to render. |
| `Analysis_pipeline_PDF.qmd` | Quarto source file for the PDF appendix. Renders a clean, code-free version of the supplementary material for journal submission. |

### Data and outputs

| Folder | Contents |
|--------|----------|
| `Data_files/` | Sample metadata and supporting data tables |
| `figures/` | All figures (PNG) referenced in the analysis documents |
| `tables/` | All result tables (CSV) referenced in the analysis documents |
| `objects/` | Saved R objects (phyloseq, phylogenetic tree) |
| `Scripts/` | R scripts for the full bioinformatic and analytical pipeline |

## Reproducibility

All analyses can be reproduced by cloning this repository and rendering `Analysis_pipeline_HTMLformat.qmd` from the repository root. Computationally intensive steps (GLLVM fitting, iNEXT3D, NRI/NTI) were executed on a remote Ubuntu server; their outputs are archived in `figures/`, `tables/`, and `objects/` and embedded as static results.

Sequence data are deposited at ENA under accession **PRJEB107725**.
