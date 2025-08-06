# DMA of *D. magna* - Paper - Data Analysis Code

This repository provides the R code to reproduce the data analysis and figures for the Deep Metabolome Annotation (DMA) of *Daphnia magna* paper.

## Overview

The repository contains three main analysis workflows:

1. **Daphnia annotation summary** - Analysis of metabolite annotations from *D. magna* samples
2. **Metabolite reference standards analysis summary** - Analysis of metabolite standard mixture (MSM) data
3. **Phylo analysis** - Phylogenetic/metabolomics analysis across species

## Project Structure

```
├── input/
│   ├── input_for_summary_plots/          # Data for Daphnia and MSM analysis
│   │   ├── merged_annotations_all_classified.zip
│   │   ├── metabolite_standard_mixture_details.csv
│   │   └── pubchem_set.zip
│   └── input_for_phylometab_plot/        # Data for phylometab analysis
│       ├── chebi_with_inchikey_source_classyfire.csv
│       ├── phyloT_generated_tree_1734701763_newick.txt
│       ├── pubchem_kegg_hmdb_expanded.zip
│       └── Supplementary Data - *.csv
├── output/                               # Generated figures and summary tables
├── paper_summarise_daphnia.R            # Main Daphnia annotation analysis
├── paper_summarise_msm.R                # Metabolite standard mixture analysis
├── paper_phylometab.R                   # Phylometab metabolomics analysis
└── renv.lock                            # R package dependencies
```

## Requirements

- R (>= 4.4.3)
- RStudio (recommended)
- Required R packages are managed via `renv` (see Installation section)

## Installation

1. Clone this repository
2. Open the R project in RStudio: `dmagna-dma-paper.Rproj`
3. init the R environment using renv:

```r
renv::init()
```

This will install all required packages with their exact versions as specified in `renv.lock`.

## Usage

### 1. Daphnia Annotation Analysis

Run the main Daphnia annotation summarization:

```r
source("paper_summarise_daphnia.R")
```

**Generates:**
- Summary statistics and visualizations of metabolite annotations
- Classification analysis (superclass, class, subclass)
- Workflow comparison plots
- Venn diagrams for extraction methods, chromatography types, and polarity
- PCA analysis of annotations
- Tree maps and upset plots

**Key outputs:**
- `FIG_4a_tree_map.pdf` - Tree map visualization
- `FIG_4b_annotations_all_pca.pdf` - PCA plot of annotations
- `FIG_4c-e_*_bar.pdf` - Bar charts for chemical classifications
- `FIG_5a-e_*.pdf` - Workflow and method comparison plots
- `daphnia_annotation_summary.csv` - Summary statistics table

### 2. Metabolite Standard Mixture Analysis

Run the metabolite reference standards analysis:

```r
source("paper_summarise_msm.R")
```

**Generates:**
- Analysis of metabolite standard mixture (MSM) annotations
- Workflow-specific analysis for MSM data

**Key outputs:**
- `FIG_S10a_galaxy_msms_workflow_bar.pdf` - MSM workflow analysis
- `FIG_S10b_treemap_msm.pdf` - MSM tree map
- `FIG_S11_presence_absence_match_type_msm.pdf` - Match type analysis
- `msm_annotations_summary.csv` - MSM summary statistics

### 3. Phylo Analysis

Run the phylogenetic/ metabolomics analysis:

```r
source("paper_phylometab.R")
```

**Generates:**
- Phylogenetic tree with metabolite presence/absence data
- Cross-species metabolite comparison
- Database mapping analysis (KEGG, HMDB, MTox, ChEBI)

**Key output:**
- `FIG_6_phylomet.pdf` - Phylogenetic metabolomics plot

## Key Dependencies

The analysis relies on several R packages:

- **Data manipulation:** `dplyr`, `tidyr`, `data.table`, `stringr`
- **Visualization:** `ggplot2`, `cowplot`, `treemap`, `VennDiagram`, `UpSetR`
- **Chemical informatics:** `ChemmineR`
- **Phylogenetics:** `ape`, `ggtree`, `aplot`
- **Data import:** `openxlsx`, `jsonlite`

## Data Sources

The analysis uses several types of input data:

- **Metabolite annotations:** Classified metabolite annotations from DMA experiments
- **Reference standards:** Metabolite standard mixture composition and annotations
- **Database mappings:** Links to KEGG, HMDB, MTox, and ChEBI databases
- **Phylogenetic data:** Species tree in Newick format
- **Chemical classifications:** ClassyFire-based chemical taxonomy

## Figures Generated

The code reproduces the following figures from the paper:

**Main Figures:**
- Figure 4: Metabolite annotation overview (tree map, PCA, classification bars)
- Figure 5: Workflow and method comparisons (Venn diagrams, upset plots)
- Figure 6: Phylogenetic metabolomics analysis

**Supplementary Figures:**
- Figure S7-S11: Additional method comparisons and MSM analysis

## Output Files

All generated figures are saved as PDF files in the `output/` directory. Summary tables are saved as CSV files for further analysis or inclusion in manuscripts.

## License

See LICENSE file for details.

## Citation

If you use this code, please cite the corresponding *D. magna* DMA paper.