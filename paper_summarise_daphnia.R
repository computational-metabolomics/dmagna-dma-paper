# =============================================================================
# Daphnia magna DMA Paper - Annotation Summary Analysis (Refactored)
# =============================================================================
# This script analyzes and summarizes metabolite annotations from Daphnia magna
# DMA experiments, generating figures and summary tables for the manuscript.
#
# Main outputs:
# - Figure 4: Metabolite annotation overview (tree map, PCA, classification)
# - Figure 5: Workflow and method comparisons  
# - Figure S7-S11: Supplementary analyses
# - Summary CSV files with annotation statistics
# =============================================================================

# Load required libraries ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(cowplot)
  library(RColorBrewer)
  
  # Visualization libraries
  library(treemap)
  library(UpSetR)
  library(VennDiagram)
  
  # Chemical informatics
  library(ChemmineR)
  
  # Data import
  library(openxlsx)
  library(jsonlite)
})

# Configuration ----
OUTPUT_DIR <- "output"
INPUT_DIR <- "input/input_for_summary_plots"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Define color palette ----
create_color_palette <- function() {
  colors <- c(
    "dodgerblue2", "green4", "#FB9A99", "#6A3D9A", "#FF7F00", 
    "grey70", "gold1", "skyblue2", "#E31A1C", "palegreen2",
    "#CAB2D6", "#FDBF6F", "grey15", "khaki2", "maroon", 
    "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", 
    "green1", "yellow4", "yellow3", "darkorange4", "brown", "grey30"
  )
  return(gplots::col2hex(colors))
}

c25 <- create_color_palette()

# Utility Functions ----

#' Standardize match type annotations based on hierarchical reliability
#' @param df Dataframe with match_type column
#' @return Dataframe with standardized match_type values
standardize_match_types <- function(df) {
  # Hierarchy: NMR > GC-MS > Spectral matching > Combined in silico > Individual in silico
  match_type_clean <- rep('', nrow(df))
  
  # Start with lowest priority and overwrite with higher priority
  match_type_clean[grepl("metfrag", df$match_type)] <- 'metfrag'
  match_type_clean[grepl("sirius", df$match_type)] <- 'sirius'
  match_type_clean[grepl("metfrag", df$match_type) & grepl("sirius", df$match_type)] <- 'metfrag & sirius'
  match_type_clean[grepl("mzcloud|gnps|sm", df$match_type)] <- 'sm'
  match_type_clean[grepl("GC-MS", df$match_type)] <- 'gcms'
  match_type_clean[grepl("NMR", df$match_type)] <- 'nmr'
  
  df$match_type <- match_type_clean
  return(df)
}

#' Clean and standardize classification data
#' @param df Dataframe with classification columns
#' @return Dataframe with cleaned classification columns
clean_classification_data <- function(df) {
  df$subclass[df$subclass==''] <- 'Undefined'
  df$subclass[is.na(df$subclass)] <- 'Undefined'
  df$subclass[df$subclass=='NA'] <- 'Undefined'
  
  return(df)
}

#' Save plot with consistent formatting
#' @param plot_obj ggplot object or grid object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot <- function(plot_obj, filename, width = 6, height = 6) {
  filepath <- file.path(OUTPUT_DIR, filename)
  ggsave(filename = filepath, plot = plot_obj, width = width, height = height, device = "pdf")
  cat("Saved:", filepath, "\n")
}

# Data Loading and Preprocessing ----

#' Load and preprocess annotation data
#' @return Preprocessed dataframe
load_annotation_data <- function() {
  cat("Loading annotation data...\n")
  
  # Load main annotation data
  df <- read.csv(unz(file.path(INPUT_DIR, "merged_annotations_all_classified.zip"),
                     "merged_annotations_all_classified.csv"))
  
  # Filter out metabolite standard mixture (MSM) data
  df <- df[!df$msm, ]
  
  # Remove unknown compounds
  df <- df[!grepl('UNKNOWN', df$inchikey), ]
  
  # Handle infinite values
  df[sapply(df, is.infinite)] <- NA
  
  cat("Loaded", nrow(df), "annotation records\n")
  return(df)
}

# Data Summarization Functions ----

#' Create main compound summary by InChI key
#' @param df Input annotation dataframe
#' @return Summary dataframe grouped by InChI key
create_compound_summary <- function(df) {
  cat("Creating compound summary by InChI key...\n")
  
  summary_df <- df %>%
    group_by(inchikey) %>%
    summarise(
      # Compound identifiers
      inchikey = toString(sort(unique(inchikey))),
      inchikey1 = toString(sort(unique(inchikey1))),
      smiles = toString(smiles[1]),
      molecular_formula = toString(sort(unique(molecular_formula))[1]),
      monoisotopic_exact_mass = toString(sort(unique(monoisotopic_exact_mass))[1]),
      compound_name = toString(sort(unique(compound_name))),
      
      # Database IDs
      pubchem_cids = toString(sort(unique(pubchem_cids))),
      pubchem_cid = toString(sort(unique(pubchem_cid))[1]),
      hmdb_ids = toString(sort(unique(cts_hmdb_ids))[1]),
      kegg_ids = toString(sort(unique(cts_kegg_ids))[1]),
      chebi_ids = toString(sort(unique(cts_chebi_ids))[1]),
      
      # Chemical classification
      kingdom = toString(sort(unique(kingdom))[1]),
      superclass = toString(sort(unique(superclass))[1]),
      class = toString(sort(unique(class))[1]),
      subclass = toString(sort(unique(subclass))[1]),
      direct_parent = toString(sort(unique(direct_parent))[1]),
      molecular_framework = toString(sort(unique(molecular_framework))[1]),
      predicted_lipidmaps_terms = toString(sort(unique(predicted_lipimaps_terms))[1]),
      
      # Experimental details
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac_number = toString(sort(unique(spe_frac))),
      spe_frac = toString(sort(unique(spe_frac_full))),
      chromatography = toString(sort(unique(chromatography))),
      measurement = toString(sort(unique(measurement))),
      polarity = toString(sort(unique(polarity))),
      
      # Annotation methods (boolean flags)
      lcmsdimsbool = max(lcmsdimsbool, na.rm = TRUE),
      nmrbool = max(nmrbool, na.rm = TRUE),
      gcmsbool = max(gcmsbool, na.rm = TRUE),
      smbool = max(smbool, na.rm = TRUE),
      mfbool = max(mfbool, na.rm = TRUE),
      siriusbool = max(siriusbool, na.rm = TRUE),
      mzcloudsmbool = max(mzcloudsmbool, na.rm = TRUE),
      galaxysmbool = max(galaxysmbool, na.rm = TRUE),
      galaxybool = max(galaxybool, na.rm = TRUE),
      gnpssmbool = max(gnpssmbool, na.rm = TRUE),
      
      # Annotation quality
      match_type = toString(unique(match_type)),
      msi_level = toString(sort(unique(msi_level))),
      .groups = 'drop'
    )
  
  # Clean and standardize data
  summary_df <- standardize_match_types(summary_df)
  summary_df <- clean_classification_data(summary_df)
  
  # Clean mass data
  summary_df$monoisotopic_exact_mass[summary_df$monoisotopic_exact_mass == 'NA'] <- ''
  
  cat("Created summary for", nrow(summary_df), "unique compounds\n")
  return(summary_df)
}

#' Create summary by compound and experimental factors
#' @param df Input annotation dataframe
#' @param group_vars Variables to group by (in addition to inchikey)
#' @return Summary dataframe
create_factor_summary <- function(df, group_vars) {
  group_vars_all <- c("inchikey", group_vars)
  
  summary_df <- df %>%
    group_by(across(all_of(group_vars_all))) %>%
    summarise(
      compoundname = toString(sort(unique(name))),
      exact_mass = toString(unique(exact_mass_temp)),
      match_type = toString(unique(match_type)),
      
      smbool = max(smbool),
      mfbool = max(mfbool),
      siriusbool = max(siriusbool),
      mzcloudsmbool =  max(mzcloudsmbool),
      galaxysmbool = max(galaxysmbool),
      nmrbool = max(nmrbool),
      gnpssmbool = max(gnpssmbool),
      galaxybool = max(galaxybool),
      match_type = toString(unique(match_type)),
      
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac = toString(sort(unique(spe_frac))),
      chromatography = toString(sort(unique(chromatography))),
      
      kingdom = toString(sort(unique(kingdom))[1]),
      superclass = toString(sort(unique(superclass))[1]),
      class = toString(sort(unique(class))[1]),
      subclass = toString(sort(unique(subclass))[1]),
      direct_parent = toString(sort(unique(direct_parent))[1]),
      molecular_framework = toString(sort(unique(molecular_framework))[1]),
      predicted_lipimaps_terms = toString(sort(unique(predicted_lipimaps_terms))[1]),
      
      measurement = toString(sort(unique(measurement))),
      polarity = toString(sort(unique(polarity))),
      msi_level = toString(sort(unique(msi_level))),
      .groups = 'drop'
    )
  
  

  
  return(clean_classification_data(summary_df))
}

# Analysis Functions ----

#' Generate annotation method statistics
#' @param summary_df Compound summary dataframe
#' @return List of statistics
calculate_annotation_stats <- function(summary_df) {
  cat("Calculating annotation statistics...\n")
  
  stats <- list()
  
  # Total annotation count
  stats$total_annotations <- nrow(summary_df)
  
  # ClassyFire coverage
  missing_classyfire <- summary_df[summary_df$superclass == 'Undefined' |summary_df$superclass == '', ]
  stats$with_classyfire <- nrow(summary_df) - nrow(missing_classyfire)
  stats$classyfire_percent <- (stats$with_classyfire / stats$total_annotations) * 100
  
  # Method counts
  stats$ms_based <- nrow(summary_df %>% filter(lcmsdimsbool == TRUE))
  
  stats$ms_based_only <- nrow(summary_df %>% filter(lcmsdimsbool == TRUE &
                                               gcmsbool == FALSE &
                                               nmrbool == FALSE))
  
  
  compound_summary_inchikey1 <- summary_df  %>% 
    group_by(inchikey1) %>%
    summarise(
      compoundname = toString(sort(unique(compound_name))),
    )
  
  stats$count_inchikey1 <- nrow(compound_summary_inchikey1)
  

  compound_summary_inchikey1_ms <- summary_df  %>% 
    filter(lcmsdimsbool == TRUE) %>%
    group_by(inchikey1) %>%
    summarise(
      compoundname = toString(sort(unique(compound_name))),
    )
  
  stats$count_inchikey1_ms <- nrow(compound_summary_inchikey1_ms)
  
  

  
  
  return(stats)
}

# Plotting Functions ----

#' Create bar plot for chemical classifications
#' @param summary_df Compound summary dataframe  
#' @param class_level Classification level ("superclass", "class", "subclass")
#' @param top_n Number of top categories to show
#' @return ggplot object
create_classification_barplot <- function(summary_df, class_level, top_n = 12) {
  
  # Create count data
  count_data <- summary_df %>%
    group_by(!!sym(class_level), match_type) %>%
    summarise(count = n_distinct(inchikey), .groups = 'drop') %>%
    arrange(desc(count))
  
  # Get top categories
  top_categories <- count_data %>%
    group_by(!!sym(class_level)) %>%
    summarise(total = sum(count), .groups = 'drop') %>%
    arrange(desc(total)) %>%
    slice_head(n = top_n) %>%
    pull(!!sym(class_level))
  
  # Filter and format data
  plot_data <- count_data %>%
    filter(!!sym(class_level) %in% top_categories) %>%
    mutate(
      !!class_level := factor(!!sym(class_level), levels = rev(top_categories)),
      match_type = factor(match_type, levels = c('sirius', 'metfrag', 'metfrag & sirius', 'sm', 'gcms', 'nmr'))
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = count, y = !!sym(class_level), fill = match_type)) +
    geom_bar(stat = 'identity', width = 0.7) +
    theme_bw() +
    scale_fill_manual("Annotation Method", values = unname(c25)) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Annotation Count", y = str_to_title(class_level))
  
  return(p)
}

#' Create tree map visualization
#' @param summary_df Compound summary dataframe
#' @return Nothing (saves PDF directly)
create_treemap <- function(summary_df) {
  cat("Creating tree map...\n")
  
  compound_superclass_count_only <- summary_df %>%
    group_by(superclass) %>%
    summarise(
      count = n_distinct(inchikey),
    ) %>%
    arrange(desc(count)) 
  
  compound_class_superclass_count <- summary_df %>%
    group_by(class) %>%
    summarise(
      superclass = unique(superclass)[1],
      count = n_distinct(inchikey),
    ) %>%
    arrange(desc(count)) 
  
  compound_class_superclass_count$superclass[compound_class_superclass_count$superclass==''] <- 'Undefined'
  compound_class_superclass_count$superclass[is.na(compound_class_superclass_count$superclass)] <- 'Undefined'
  
  compound_class_superclass_count$superclass <- factor(compound_class_superclass_count$superclass,
                                                       levels=unique(compound_superclass_count_only$superclass))
  pdf('output/FIG_4a_tree_map.pdf', width=15, height = 9)
  
  c25_treemap <- c(
    "dodgerblue2", 
    "green4",
    "#FB9A99",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "grey70",
    "#E31A1C",
    "skyblue2", 
    "gold", # red
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "grey15", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown", "grey30"
  )
  
  treemap(compound_class_superclass_count ,
          index=c("superclass","class"),
          vSize="count",
          type="index",
          palette =  gplots::col2hex(c25_treemap),
          
          bg.labels=244,
          align.labels=list(
            c("left", "top"), 
            c("right", "bottom")
            
          ),
          
          na.rm=F
  )
  dev.off()
}

#' Create PCA plot using PubChem fingerprints
#' @param summary_df Compound summary dataframe
#' @return Nothing (saves plots directly)
create_pca_plot <- function(summary_df) {
  cat("Creating PCA plot...\n")
  fpset <- read.SDFset(unzip("input/input_for_summary_plots/pubchem_set.zip", "pubchem_set.sdf"))
  cid(fpset) <- sdfid(fpset, tag=2)
  
  pubchem_fps <- fp2bit(fpset, type = 2, fptag = "PUBCHEM_CACTVS_SUBSKEYS")
  
  pca_res_scale_f <- prcomp(pubchem_fps , scale=FALSE)
  
  pca_res <- pca_res_scale_f
  pca_res_summary <- summary(pca_res)
  
  fp_meta <- merge(data.frame('cid'=rownames(pca_res$x)), compound_summary, by.x='cid', by.y='pubchem_cid', all.x = T,no.dups = T )
  
  compound_superclass_count_pca <- fp_meta %>%
    group_by(superclass) %>%
    summarise(
      count = n_distinct(inchikey),
    ) %>%
    arrange(desc(count)) 
  
  fp_meta$superclass[fp_meta$superclass == ''] = 'Undefined'
  fp_meta$superclass <- factor(fp_meta$superclass, 
                               levels=compound_superclass_count_pca$superclass[1:9])
  
  
  dtp <- data.frame('Superclass' =fp_meta$superclass, pca_res$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
  ggplot(data = dtp) + 
    geom_point(aes(x = PC1, y = PC2, col = Superclass)) + 
    scale_color_manual(values = c25)+
    xlab(paste('PC1 (', round(pca_res_summary$importance[,1][2]*100,1), '%)', sep='')) +
    ylab(paste('PC2 (', round(pca_res_summary$importance[,2][2]*100,1), '%)', sep='')) + 
    theme_minimal()  +
    theme(legend.position="")
  
  ggsave("output/FIG_4b_annotations_all_pca.pdf", width=6, height=6)
  
  ggplot(data = dtp) + 
    geom_point(aes(x = PC1, y = PC2, col = Superclass)) + 
    scale_color_manual(values = c25)+
    xlab(paste('PC1 (', round(pca_res_summary$importance[,1][2]*100,1), '%)', sep='')) +
    ylab(paste('PC2 (', round(pca_res_summary$importance[,2][2]*100,1), '%)', sep='')) + 
    theme_minimal()
  
  ggsave("output/FIG_4b_annotations_all_pca_legend.pdf", width=8, height=6)
  
}

#' Create Venn diagrams for method comparisons
#' @param summary_df Compound summary dataframe
#' @param df_merged_filtered Original annotation dataframe
#' @return Nothing (saves plots directly)
create_venn_diagrams <- function(summary_df, df_merged_filtered) {
  cat("Creating Venn diagrams...\n")
  
  # Color palette for Venn diagrams
  venn_colors <- brewer.pal(4, "Pastel2")
  
  # Helper function to create and save Venn diagrams
  create_venn <- function(sets, labels, filename, positions = c(-27, 26), distances = c(0.055, 0.055)) {
    vp <- venn.diagram(
      x = sets,
      category.names = labels,
      filename = NULL,
      output = TRUE,
      print.mode = c("raw", "percent"),
      width = 10, height = 10, units = "in",
      lwd = 2, lty = 'blank',
      fill = venn_colors[1:length(sets)],
      cex = 1.7, fontfamily = "sans",
      cat.cex = 2, cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = positions,
      cat.dist = distances,
      cat.fontfamily = "sans"
    )
    
    ggsave(vp, file = file.path(OUTPUT_DIR, filename), device = "pdf", width = 6, height = 6)
    cat("Saved:", file.path(OUTPUT_DIR, filename), "\n")
  }
  
  # Create factor summaries for Venn diagrams
  polarity_summary <- create_factor_summary(df_merged_filtered, "polarity")
  extraction_summary <- create_factor_summary(df_merged_filtered, "extraction")
  
  # 1. Polarity comparison
  pos_compounds <- unique(polarity_summary[polarity_summary$polarity == 'POSITIVE', ]$inchikey)
  neg_compounds <- unique(polarity_summary[polarity_summary$polarity == 'NEGATIVE', ]$inchikey)
  create_venn(list(pos_compounds, neg_compounds), 
              c("Positive", "Negative"), 
              "FIG_5c_polarity_venn.pdf")
  
  # 2. Extraction comparison  
  apolar_compounds <- unique(extraction_summary[extraction_summary$extraction == 'APOLAR', ]$inchikey)
  polar_compounds <- unique(extraction_summary[extraction_summary$extraction == 'POLAR', ]$inchikey)
  create_venn(list(apolar_compounds, polar_compounds), 
              c("Apolar", "Polar"), 
              "FIG_5a_extraction_venn.pdf")
  
  # 3. Annotation method comparison
  sm_compounds <- unique(summary_df[summary_df$smbool == 1, ]$inchikey)
  mf_compounds <- unique(summary_df[summary_df$mfbool == 1, ]$inchikey)
  sirius_compounds <- unique(summary_df[summary_df$siriusbool == 1, ]$inchikey)
  create_venn(sets=list(sm_compounds, mf_compounds, sirius_compounds),
              labels=c("Spectral \n matching", "MetFrag", "SIRIUS \n CSI:FingerID"),
              "FIG_S8a_annotation_type_venn.pdf",
              positions = c(-20, 26, 125),
              distances = c(0.055, 0.055, 0.055))
  
  # 4. Spectral matching tools comparison
  gnps_compounds <- unique(summary_df[summary_df$gnpssmbool == 1, ]$inchikey)
  mzcloud_compounds <- unique(summary_df[summary_df$mzcloudsmbool == 1, ]$inchikey)
  galaxy_compounds <- unique(summary_df[summary_df$galaxysmbool == 1, ]$inchikey)
  create_venn(list(gnps_compounds, mzcloud_compounds, galaxy_compounds),
              c("Spectral matching\n(GNPS workflow)", "Spectral matching\n(mzCloud)", "Spectral matching\n(Galaxy)"),
              "FIG_S8b_spectral_matching_venn.pdf",
              positions = c(-20, 26, 125),
              distances = c(0.055, 0.055, 0.055))
  
  # 5. Chromatography comparison
  c30_compounds <- unique(summary_df[grepl('C30', summary_df$chromatography), ]$inchikey)
  amd_compounds <- unique(summary_df[grepl('AMD', summary_df$chromatography), ]$inchikey)
  phe_compounds <- unique(summary_df[grepl('PHE', summary_df$chromatography), ]$inchikey)
  create_venn(list(c30_compounds, amd_compounds, phe_compounds),
              c("C30", "AMD", "PHE"),
              "FIG_5b_chromatography_type_venn.pdf",
              positions = c(-20, 26, 125),
              distances = c(0.055, 0.055, 0.055))
  
  # 6. Analytical platform comparison
  lcdims_compounds <- unique(summary_df[grepl('DIM|LC', summary_df$measurement), ]$inchikey)
  nmr_compounds <- unique(summary_df[grepl('NMR', summary_df$measurement), ]$inchikey)
  gcms_compounds <- unique(summary_df[grepl('GC', summary_df$measurement), ]$inchikey)
  create_venn(list(lcdims_compounds, nmr_compounds, gcms_compounds),
              c("LC-MS/MS, DIMSn", "NMR", "GC-MS"),
              "FIG_S7_dims_lcms_venn.pdf",
              positions = c(-20, 26, 125),
              distances = c(0.055, 0.055, 0.055))
}

# Workflow Analysis Functions ----

#' Create workflow summary and bar plots
#' @param df Input annotation dataframe
#' @return Nothing (saves plots directly)
create_workflow_analysis <- function(df) {
  cat("Creating workflow analysis...\n")
  
  df_merged_filtered_for_hist <- df[df$galaxysmbool |
                                    df$siriusbool |
                                    df$mfbool |
                                    df$mzcloudsmbool,
  ]
  
  workflow_summary1 <-  df_merged_filtered_for_hist %>%
    group_by(inchikey, extraction, polarity, spe, spe_frac, chromatography) %>%
    summarise(
      compoundname = toString(sort(unique(name))),
      hmdb_id = toString(sort(unique(hmdb_ids))),
      exact_mass = toString(unique(exact_mass_temp)),
      match_type = toString(unique(match_type)),
      
      inchikey = toString(sort(unique(inchikey))),
      
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac = toString(sort(unique(spe_frac))),
      chromatography = toString(sort(unique(chromatography))),
      
      measurement = toString(sort(unique(measurement))),
      polarity = toString(sort(unique(polarity))),
      count = n_distinct(inchikey1),
      
      kingdom = toString(sort(unique(kingdom))[1]),
      superclass = toString(sort(unique(superclass))[1]),
      class = toString(sort(unique(class))[1]),
      subclass = toString(sort(unique(subclass))[1]),
      direct_parent = toString(sort(unique(direct_parent))[1]),
      molecular_framework = toString(sort(unique(molecular_framework))[1]),
      predicted_lipimaps_terms = toString(sort(unique(predicted_lipimaps_terms))[1]),
      
      
      
    )
  
  workflow_summary2 <- df_merged_filtered_for_hist  %>%
    group_by(inchikey, extraction, spe, spe_frac, chromatography) %>%
    summarise(
      compoundname = toString(sort(unique(name))),
      hmdb_id = toString(sort(unique(hmdb_ids))),
      exact_mass = toString(unique(exact_mass_temp)),
      match_type = toString(unique(match_type)),
      
      inchikey = toString(sort(unique(inchikey))),
      
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac = toString(sort(unique(spe_frac))),
      chromatography = toString(sort(unique(chromatography))),
      measurement = toString(sort(unique(measurement))),
      polarity = toString(sort(unique(polarity))),
      polarity = 'All',
      count = n_distinct(inchikey1),
      
      kingdom = toString(sort(unique(kingdom))[1]),
      superclass = toString(sort(unique(superclass))[1]),
      class = toString(sort(unique(class))[1]),
      subclass = toString(sort(unique(subclass))[1]),
      direct_parent = toString(sort(unique(direct_parent))[1]),
      molecular_framework = toString(sort(unique(molecular_framework))[1]),
      predicted_lipimaps_terms = toString(sort(unique(predicted_lipimaps_terms))[1]),
      
    )
  
  workflow_summary3 <- rbind(workflow_summary1, workflow_summary2)
  
  workflow_summary3$extract_spe_lc <- paste(workflow_summary3$extraction,
                                            workflow_summary3$spe,
                                            workflow_summary3$spe_frac,
                                            workflow_summary3$chromatography)
  
  
  workflow_summary3$extract_spe_lc <- gsub('D-', '', workflow_summary3$extract_spe_lc)
  workflow_summary3$extract_spe_lc[workflow_summary3$extract_spe_lc=='APOLAR  1 C30'] <- 'APOLAR   C30'
  
  workflow_summary3$extract_spe_lc <- factor(workflow_summary3$extract_spe_lc, 
                                             levels=c("POLAR   AMD",
                                                      "POLAR   PHE",
                                                      "POLAR WCX 1 PHE",
                                                      "POLAR WCX 1 AMD",
                                                      "POLAR WCX 2 PHE",
                                                      "POLAR WCX 2 AMD",
                                                      "POLAR WCX 3 PHE",
                                                      "POLAR WCX 3 AMD",
                                                      "POLAR WCX 4 PHE",
                                                      "POLAR WCX 4 AMD",
                                                      "POLAR WAX 1 PHE",
                                                      "POLAR WAX 1 AMD",
                                                      "POLAR WAX 2 PHE",
                                                      "POLAR WAX 2 AMD",
                                                      "POLAR WAX 3 PHE",
                                                      "POLAR WAX 3 AMD",
                                                      "POLAR WAX 4 PHE",
                                                      "POLAR WAX 4 AMD",
                                                      "APOLAR C18 1 C30",
                                                      "APOLAR C18 2 C30",
                                                      "APOLAR C18 3 C30",
                                                      "APOLAR AMP 1 C30",
                                                      "APOLAR AMP 2 C30",
                                                      "APOLAR AMP 3 C30",
                                                      "APOLAR AMP 4 C30",
                                                      "APOLAR   C30"
                                             ))
  
  
  workflow_summary3$superclass[workflow_summary3$superclass == ''] = 'Undefined'
  workflow_summary3$class[workflow_summary3$class == ''] = 'Undefined'
  
  compound_superclass_count_for_bar <- workflow_summary3 %>%
    group_by(superclass) %>%
    summarise(
      count = n_distinct(inchikey),
      
      n = n()) %>%
    arrange(desc(n)) 
  
  # use this as the naming category
  names(c25) <- compound_superclass_count_for_bar$superclass[1:25]
  
  workflow_summary_for_bar <- workflow_summary3
  workflow_summary_for_bar$superclass <- factor(workflow_summary_for_bar$superclass, 
                                                levels=compound_superclass_count_for_bar$superclass[1:8])
  
  
  workflow_bar <- ggplot(workflow_summary_for_bar, aes(x=extract_spe_lc, y=count, fill=superclass, label = count))+
    geom_bar(stat='identity', width =0.7) +
    theme_bw()+
    scale_y_continuous(name="Unique metabolite annotation count", limits=c(0, 1600))+
    scale_fill_manual("legend", values = c25)+
    theme(legend.position="",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_discrete(guide = guide_axis(angle = 90), drop=FALSE) +
    facet_wrap(~polarity, ncol= 1)
  
  
  ggsave("output/FIG_5d_annotations_all_workflow_bar.pdf", width=10, height=10)
  
  workflow_bar_with_legend <- workflow_bar  +   theme(legend.position="bottom",
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank())
  
  ggsave("output/FIG_5d_annotations_all_workflow_bar_with_legend.pdf", width=10, height=10)
  
  # Create UpSet plot
  create_upset_plot(workflow_summary3)
}

#' Create UpSet plot for workflow intersections
#' @param workflow_summary Workflow summary dataframe
#' @return Nothing (saves plot directly)
create_upset_plot <- function(workflow_summary) {
  cat("Creating UpSet plot...\n")
  
  dl = plyr::dlply(workflow_summary, ~ extract_spe_lc , function(x){
    unique(x$inchikey)
  })
  
  upsetr <- upset(fromList(dl),mb.ratio = c(0.35, 0.65), 
                  
                  order.by = "freq", sets=as.character(unique(workflow_summary$extract_spe_lc)),
                  nintersects=140,
                  sets.bar.color='lightblue',
                  show.numbers=FALSE)
  print(upsetr)
  pdf('output/FIG_5e_annotations_all_upset.pdf', width=13, height=5)
  upsetr
  dev.off()
}

# Mass Distribution Analysis ----

#' Create mass distribution histogram
#' @param df Input annotation dataframe
#' @param summary_df Compound summary dataframe
#' @return Nothing (saves plot directly)
create_mass_histogram <- function(df, summary_df) {
  cat("Creating mass distribution histogram...\n")
  
  compound_superclass_count_only <- summary_df %>%
    group_by(superclass) %>%
    summarise(
      count = n_distinct(inchikey),
    ) %>%
    arrange(desc(count)) 
  
  df_mass_hist_filtered <- df[!is.na(df$polarity),]
  df_mass_hist_filtered <- merge(df_mass_hist_filtered, summary_df[,c('inchikey', 'monoisotopic_exact_mass')], all.x = T)
  mass_hist_summary1 <- df_mass_hist_filtered   %>%
    group_by(inchikey1) %>%
    summarise(
      monoisotopic_exact_mass = unique(monoisotopic_exact_mass),
      polarity = 'all',
      kingdom = toString(sort(unique(kingdom))[1]),
      superclass = toString(sort(unique(superclass))[1]),
      class = toString(sort(unique(class))[1]),
      subclass = toString(sort(unique(subclass))[1]),
    )
  
  mass_hist_summary2 <- df_mass_hist_filtered  %>%
    group_by(inchikey1, polarity) %>%
    summarise(
      monoisotopic_exact_mass = unique(monoisotopic_exact_mass),
      kingdom = toString(sort(unique(kingdom))[1]),
      superclass = toString(sort(unique(superclass))[1]),
      class = toString(sort(unique(class))[1]),
      subclass = toString(sort(unique(subclass))[1]),
    )
  mass_hist_summary3 <- rbind(mass_hist_summary1, mass_hist_summary2)
  
  mass_hist_summary3$superclass[mass_hist_summary3$superclass == ''] = 'Undefined'
  mass_hist_summary3$class[mass_hist_summary3$class == ''] = 'Undefined'
  
  mass_hist_summary3$superclass <- factor(mass_hist_summary3$superclass, 
                                          levels=compound_superclass_count_only$superclass[1:9])
  
  ggplot(mass_hist_summary3 , aes(as.numeric(monoisotopic_exact_mass), fill=superclass))+
    theme_bw()+
    geom_histogram(color = "black", bins=75, size=0.4) +
    scale_fill_manual("legend", values = c25)+
    theme(legend.position="bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    facet_wrap(~polarity, ncol = 1)
  
  ggsave("output/FIG_S9_annotations_all_mass_hist.pdf", width=20, height=10)
}

# Export Functions ----

#' Export summary tables
#' @param summary_df Compound summary dataframe
#' @return Nothing (saves CSV files)
export_summary_tables <- function(summary_df) {
  cat("Exporting summary tables...\n")
  
  # Format main summary for export
  export_columns <- c(
    'inchikey', 'inchikey1', 'smiles', 'molecular_formula', 'monoisotopic_exact_mass',
    'compound_name', 'pubchem_cids', 'hmdb_ids', 'kegg_ids', 'chebi_ids',
    'kingdom', 'superclass', 'class', 'subclass', 'direct_parent', 'molecular_framework',
    'predicted_lipidmaps_terms', 'assays', 'extraction', 'spe', 'spe_frac',
    'chromatography', 'measurement', 'polarity', 'lcmsdimsbool', 'gcmsbool',
    'nmrbool', 'smbool', 'mfbool', 'siriusbool', 'mzcloudsmbool', 'galaxysmbool',
    'gnpssmbool', 'msi_level'
  )
  
  summary_export <- summary_df[, export_columns]
  write.csv(summary_export, file.path(OUTPUT_DIR, 'daphnia_annotation_summary.csv'), na = "", row.names = FALSE)
  
  # Export classification counts
  export_classification_counts(summary_df, 'superclass', 'daphnia_annotations_all_superclass_count.csv')
  export_classification_counts(summary_df, 'class', 'daphnia_annotations_all_class_count.csv')
  export_classification_counts(summary_df, 'subclass', 'daphnia_annotations_all_subclass_count.csv')
  
  cat("Exported summary tables to output directory\n")
}

#' Export classification count tables
#' @param summary_df Compound summary dataframe
#' @param class_level Classification level
#' @param filename Output filename
export_classification_counts <- function(summary_df, class_level, filename) {
  counts <- summary_df %>%
    group_by(!!sym(class_level)) %>%
    summarise(count = n_distinct(inchikey), .groups = 'drop') %>%
    arrange(desc(count))
  
  write.csv(counts, file.path(OUTPUT_DIR, filename), row.names = FALSE)
}

# Main Analysis Pipeline ----

#' Run complete analysis pipeline
run_daphnia_analysis <- function() {
  cat("=== Starting Daphnia magna annotation analysis ===\n")
  start_time <- Sys.time()
  
  # 1. Load and preprocess data
  df_merged_filtered <<- load_annotation_data()
  
  # 2. Create main compound summary
  compound_summary <<- create_compound_summary(df_merged_filtered)
  
  # 3. Calculate and display statistics
  stats <- calculate_annotation_stats(compound_summary)
  cat("\n=== Analysis Statistics ===\n")
  cat("Total annotations:", stats$total_annotations, "\n")
  cat("With ClassyFire details:", stats$with_classyfire, "(", round(stats$classyfire_percent, 1), "%)\n")
  cat("Total inchikey1:", stats$count_inchikey1, "\n")
  cat("MS-based inchikey1:", stats$count_inchikey1_ms, "\n")
  cat("MS-based annotations:", stats$ms_based, "\n")
  cat("MS-based annotations (unique to MS):", stats$ms_based_only, "\n")
  
  
  # 4. Generate all figures
  cat("\n=== Generating Figures ===\n")
  
  # Figure 4: Main annotation overview
  create_treemap(compound_summary)
  create_pca_plot(compound_summary)
  
  # Classification bar plots
  for (level in c("superclass", "class", "subclass")) {
    plot <- create_classification_barplot(compound_summary, level)
    filename <- paste0("FIG_4", switch(level, "superclass" = "c", "class" = "d", "subclass" = "e"), 
                      "_annotations_all_", level, "_bar.pdf")
    save_plot(plot, filename, width = 4, height = 4)
  }
  
  # Save legend separately
  legend_plot <- create_classification_barplot(compound_summary, "superclass")
  legend <- cowplot::get_legend(legend_plot)
  pdf(file.path(OUTPUT_DIR, 'FIG_4_legend.pdf'))
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
  cat("Saved: output/FIG_4_legend.pdf\n")
  
  # Figure 5: Workflow and method comparisons
  create_workflow_analysis(df_merged_filtered)
  create_venn_diagrams(compound_summary, df_merged_filtered)
  
  # Supplementary figures
  create_mass_histogram(df_merged_filtered, compound_summary)
  
  # 5. Export summary tables
  export_summary_tables(compound_summary)
  
  # Analysis complete
  end_time <- Sys.time()
  cat("\n=== Analysis Complete ===\n")
  cat("Runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
  cat("All outputs saved to:", OUTPUT_DIR, "\n")
}

# Execute Analysis ----
# Run the analysis when script is sourced
run_daphnia_analysis()
