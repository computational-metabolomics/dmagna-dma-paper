# DMA of *D. magna* - Paper - Data Analysis Code

This repository provides a comprehensive set of R analysis code, supporting source data, and result files for the data-analysis summaries and figures in the Deep Metabolome Annotation (DMA) of *Daphnia magna* paper. Where complete figure-generation workflows are not included, the numerical input/source data are provided to support transparency and reuse.

## Overview

The repository provides analysis outputs and supporting files for the DMA of *D. magna* paper:

1. **Daphnia annotation summary** - R code and data for analysis of metabolite annotations from *D. magna* samples
2. **Metabolite reference standards analysis summary** - R code and data for analysis of metabolite standard mixture (MSM) data
3. **Phylo analysis** - R code and data for phylogenetic/metabolomics analysis across species
4. **Molecular network files** - Network/graph files used to generate molecular network figures
5. **IPA analysis** - IPA pathway results used to generate summary values and supplementary figure
6. **LC-MS optimisation analysis** - Numerical input/source data used to generate LC-MS optimisation analyis supplementary figures
7. **Example Feature Check (Galaxy workflow history access)** - Example for how the Galaxy workflow histories can be investigated

## Figures covered 

This repository provides the numerical data (and relevant code where highlighted) used to generate the following figures:

**Main Figures:**
- Figure 5 (`FIG_5a_tree_map.pdf`, `FIG_5b_annotations_all_pca.pdf`, `FIG_5c_annotations_all_superclass_bar.pdf`, `FIG_5d_annotations_all_class_bar.pdf`, `FIG_5e_annotations_all_subclass_bar.pdf`): Metabolite annotation overview (tree map, PCA, classification bars)
- Figure 6 (`FIG_6a_extraction_venn.pdf`, `FIG_6b_chromatography_type_venn.pdf`, `FIG_6c_polarity_venn.pdf`, `FIG_6d_annotations_all_workflow_bar.pdf`, `FIG_6e_annotations_all_upset.pdf`): Workflow and method comparisons (Venn diagrams, workflow bars, upset plots)
- Figure 7 (`FIG_7_phylomet.pdf`): Phylogenetic metabolomics analysis
- Figure 8: Molecular network analysis (negative ionisation; source files in `network_files/`)

**Supplementary Figures:**
- Figures S9-S26: LC-MS optimisation analysis (numerical input/source data only in `LC-MS_optimisation_analysis/`; figure-generation code is not included)
- Figure S27 (`FIG_S27_dims_lcms_venn.pdf`): DIMS vs LC-MS derived annotation Venn diagram
- Figure S28a (`FIG_S28a_annotation_type_venn.pdf`): Venn diagram of annotation computational approaches
- Figure S28b (`FIG_S28b_spectral_matching_venn.pdf`): Venn diagram of spectral matching approaches
- Figure S29 (`FIG_S29_annotations_all_mass_hist.pdf`): Histogram of the mass ranges of annotations observed
- Figure S30a (`FIG_S30a_galaxy_msms_workflow_bar.pdf`): MSM workflow analysis
- Figure S30b (`FIG_S30b_treemap_msm.pdf`): MSM tree map
- Figure S31 (`FIG_S31_presence_absence_match_type_msm.pdf`): MSM presence/absence and match-type analysis
- Figure S32 (`IPA_analysis/DMA_IPA_output_pathways.xls`): IPA pathway analysis source file
- Figure S33: Molecular network analysis (source files in `network_files/`)


## Project Structure

```
├── input/
│   ├── input_for_feature_check/          # Inputs for Galaxy workflow feature check
│   │   ├── galaxy_peaklist_references.csv
│   │   └── GalaxyNone-[samplelist_dma_daphnia_magna.tabular].tabular
│   ├── input_for_summary_plots/          # Data for Daphnia and MSM analysis
│   │   ├── merged_annotations_all_classified.zip
│   │   ├── metabolite_standard_mixture_details.csv
│   │   └── pubchem_set.zip
│   └── input_for_phylometab_plot/        # Data for phylometab analysis
│       ├── chebi_with_inchikey_source_classyfire.csv
│       ├── Daphnia_ChEBI.csv
│       ├── MTox.csv
│       ├── phyloT_generated_tree_1734701763_newick.txt
│       └── pubchem_kegg_hmdb_expanded.zip
├── IPA_analysis/                         # IPA source data for Supplemental Figure S32
│   └── DMA_IPA_output_pathways.xls
├── LC-MS_optimisation_analysis/          # Numerical input/source data for Supplemental Figures S9-S26 (no figure-generation code)
├── network_files/                        # Files used for molecular network figures
├── output/                               # Generated figures and summary tables
├── example_feature_check.R               # Galaxy workflow feature check example
├── paper_summarise_daphnia.R             # Main Daphnia annotation analysis
├── paper_summarise_msm.R                 # Metabolite standard mixture analysis
└── paper_phylometab.R                    # Phylometab metabolomics analysis
```

## Requirements for R analysis

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
- `FIG_5a_tree_map.pdf` - Tree map visualization
- `FIG_5b_annotations_all_pca.pdf` - PCA plot of annotations
- `FIG_5c-e_*_bar.pdf` - Bar charts for chemical classifications
- `FIG_6a-e_*.pdf` - Workflow and method comparison plots
- `daphnia_annotation_summary.csv` - Summary statistics table
- Figure S27 (`FIG_S27_dims_lcms_venn.pdf`) - DIMS vs LC-MS derived annotation Venn diagram
- Figure S28a (`FIG_S28a_annotation_type_venn.pdf`) - Venn diagram of annotation computational approaches
- Figure S28b (`FIG_S28b_spectral_matching_venn.pdf`) - Venn diagram of spectral matching approaches
- Figure S29 (`FIG_S29_annotations_all_mass_hist.pdf`) - Histogram of the mass ranges of annotations observed

### 2. Metabolite Standard Mixture Analysis

Run the metabolite reference standards analysis:

```r
source("paper_summarise_msm.R")
```

**Generates:**
- Analysis of metabolite standard mixture (MSM) annotations
- Workflow-specific analysis for MSM data

**Key outputs:**
- `FIG_S30a_galaxy_msms_workflow_bar.pdf` - MSM workflow analysis
- `FIG_S30b_treemap_msm.pdf` - MSM tree map
- `FIG_S31_presence_absence_match_type_msm.pdf` - Match type analysis
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
- `FIG_7_phylomet.pdf` - Phylogenetic metabolomics plot

### 4. Example Feature Check (Galaxy workflow history access)

Use the example feature check to show how readers can directly access files from Galaxy workflows and verify LC-MS feature details against blank-filtered XCMS features.

```r
source("example_feature_check.R")
```

**What it does:**
- Downloads XCMS peak lists and xcmsSet objects from Galaxy URLs
- Rebuilds RT windows and performs blank filtering
- Links the XCMS features from the Galaxy workflow to full annotation list


**Inputs:**
- Galaxy workflow file URLs in `input/input_for_feature_check/galaxy_peaklist_references.csv`
- Sample metadata in `input/input_for_feature_check/GalaxyNone-[samplelist_dma_daphnia_magna.tabular].tabular`
- Full annotation list in [input/input_for_summary_plots/merged_annotations_all_classified.zip](input/input_for_summary_plots/merged_annotations_all_classified.zip)

**Key outputs (per assay in `output/<assay_name>/`):**
- `*_DE_blank_filtered.RDS` and `*_blank_filtered_peak_matrix.csv`
- `*_xcms_passed_annos.csv`

## Key Dependencies

The analysis relies on several R packages:

- **Data manipulation:** `dplyr`, `tidyr`, `data.table`, `stringr`
- **Visualization:** `ggplot2`, `cowplot`, `treemap`, `VennDiagram`, `UpSetR`
- **Chemical informatics:** `ChemmineR`
- **Phylogenetics:** `ape`, `ggtree`, `aplot`
- **Data import:** `openxlsx`, `jsonlite`

## Output Files

All generated figures are saved as PDF files in the `output/` directory. Summary tables are saved as CSV files for further analysis or inclusion in manuscripts.

Also includes an updated metabolites file created for the MetaboLights study MTBLS2273.

## License

See LICENSE file for details.

## Citation

If you use this code or data within this repo please cite the corresponding *D. magna* DMA paper.
