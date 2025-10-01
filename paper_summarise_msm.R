# =============================================================================
# Daphnia magna DMA Paper - Metabolite Standard Mixture (MSM) Analysis (Clean)
# =============================================================================
# This script analyzes metabolite annotations that match to the Metabolite 
# Standard Mixture (MSM) used in the DMA experiments.
#
# Main outputs:
# - Figure S10a: Workflow analysis for MSM matches
# - Figure S10b: Tree map of MSM annotations by chemical class
# - Figure S11: Presence/absence matrix of annotation methods for MSM
# - CSV summary of MSM annotation results
# =============================================================================

cat("=============================================================================\n")
cat("METABOLITE STANDARD MIXTURE (MSM) ANALYSIS\n") 
cat("=============================================================================\n")

# Load required libraries ----
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(VennDiagram)
  library(treemap)
})

# Configuration ----
setwd('.')
OUTPUT_DIR <- "output"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Define color palette ----
c25 <- c(
  "dodgerblue2", "green4", "#FB9A99", "#6A3D9A", "#FF7F00", 
  "grey70", "gold1", "skyblue2", "#E31A1C", "palegreen2",
  "#CAB2D6", "#FDBF6F", "black", "khaki2", "maroon", 
  "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", 
  "green1", "yellow4", "yellow3", "darkorange4", "brown", "grey30"
)
c25 <- gplots::col2hex(c25)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Load and preprocess MSM annotation data
#' @return List containing standardmix_details and dfmsm_merged
load_msm_data <- function() {
  cat("\n=== Loading Data ===\n")
  
  # Load metabolite standard mixture details
  cat("Loading metabolite standard mixture details...\n")
  standardmix_details <- read.csv("input/input_for_summary_plots/metabolite_standard_mixture_details.csv")
  cat("Standards in mixture:", nrow(standardmix_details), "\n")
  
  # Load merged annotation data
  cat("Loading merged annotation data...\n")
  df_merged_filtered <- read.csv(unz("input/input_for_summary_plots/merged_annotations_all_classified.zip",
                                     "merged_annotations_all_classified.csv"))
  
  # Filter for MSM data only
  cat("Filtering for MSM annotations...\n")
  dfmsm_merged <- df_merged_filtered[df_merged_filtered$msm,]
  cat("Total MSM annotation records:", nrow(dfmsm_merged), "\n")
  cat("Unique MSM compounds detected:", length(unique(dfmsm_merged$inchikey1)), "\n")
  
  return(list(
    standardmix_details = standardmix_details,
    dfmsm_merged = dfmsm_merged
  ))
}

#' Create main compound summary for MSM data
#' @param dfmsm_merged Filtered MSM annotation data
#' @param standardmix_details Standard mixture details
#' @return List containing compound summary and statistics
create_msm_summary <- function(dfmsm_merged, standardmix_details) {
  cat("\n=== Creating Main Annotation Summary ===\n")
  
  compound_summary1 <- dfmsm_merged  %>%
    group_by(inchikey1) %>%
    summarise(
      compoundname = toString(sort(unique(name))),
      match_type = toString(unique(match_type)),
      inchikey2d = toString(unique(inchikey1)),
      
      smbool = max(smbool, na.rm = T),
      mfbool = max(mfbool, na.rm = T),
      siriusbool = max(siriusbool, na.rm = T),
      mzcloudsmbool =  max(mzcloudsmbool, na.rm = T),
      galaxysmbool = max(galaxysmbool, na.rm = T),
      gnpssmbool = max(gnpssmbool, na.rm = T),
      match_type = toString(unique(match_type)),
      
      rank_min = min(rank, na.rm=T),
      wscore_max = max(is.finite(wscore),na.rm=T),
      metfrag_score_max = max(is.finite(metfrag_score), na.rm=T),
      sirius_score_max = max(is.finite(sirius_score), na.rm=T),
      sm_score_max = round(max(is.finite(sm_score), na.rm=T), 4),
      inchikey = toString(sort(unique(inchikey))),
      
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac = toString(sort(unique(spe_frac))),
      spe_frac = toString(sort(unique(spe_frac))),
      chromatography = toString(sort(unique(chromatography))),
      
      measurement = toString(sort(unique(measurement))),
      polarity = toString(sort(unique(polarity)))
    )
  
  # Merge with the standard mix (so only rows which match to the inchikey 2D are shown)
  compound_summary1 <- merge(compound_summary1, standardmix_details[,c('Compound', 'inchikey2d', 'superclass', 'class', 'smiles')],  by.x = 'inchikey1', by.y = 'inchikey2d')
  
  # Calculate match statistics
  total_count_matched <- nrow(compound_summary1)
  perc_matched_msm <- round(nrow(compound_summary1)/48*100, 1)
  
  cat("Standards successfully annotated:", total_count_matched, "out of 48\n")
  cat("Percentage matched:", perc_matched_msm, "%\n")
  
  # Export main summary
  cat("Exporting MSM annotation summary...\n")
  write.csv(compound_summary1, file.path(OUTPUT_DIR, 'msm_annotations_summary.csv'), row.names = FALSE)
  
  # Identify missing standards
  missing_mix <- standardmix_details[!standardmix_details$inchikey2d %in% unique(compound_summary1$inchikey2d),]
  missing_mix <- missing_mix[,-which(colnames(missing_mix) %in% c('IUPAC.name'))]
  cat("Standards not detected:", nrow(missing_mix), "\n")
  if(nrow(missing_mix) > 0) {
    cat("Missing compounds:", paste(missing_mix$Compound, collapse=", "), "\n")
  }
  
  return(list(
    compound_summary = compound_summary1,
    total_matched = total_count_matched,
    perc_matched = perc_matched_msm,
    missing_standards = missing_mix
  ))
}

#' Create workflow analysis summary
#' @param dfmsm_merged Filtered MSM annotation data
#' @param standardmix_details Standard mixture details
#' @return Workflow summary dataframe
create_workflow_summary <- function(dfmsm_merged, standardmix_details) {
  cat("\n=== Creating Workflow Analysis ===\n")
  
  # Summary by workflow components (with polarity)
  cat("Creating workflow summary with polarity...\n")
  compound_summary2 <- dfmsm_merged  %>%
    group_by(inchikey1, polarity, spe, spe_frac, chromatography) %>%
    summarise(
      compoundname = toString(sort(unique(name))),
      hmdb_id = toString(sort(unique(hmdb_ids))),
      match_type = toString(unique(match_type)),
      
      # Scoring information
      rank_min = min(rank, na.rm=T),
      wscore_max = max(wscore,na.rm=T),
      sirius_score_max = max(sirius_score, na.rm=T),
      metfrag_score_max = max(metfrag_score, na.rm=T),
      sm_score = round(max(sm_score, na.rm=T), 4),
      inchikey = toString(sort(unique(inchikey))),
      
      # Experimental details
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac = toString(sort(unique(spe_frac))),
      chromatography = toString(sort(unique(chromatography))),
      measurement = toString(sort(unique(measurement))),
      polarity = toString(sort(unique(polarity))),
      count = n_distinct(inchikey1),
      .groups = 'drop'
    )
  
  # Summary by workflow components (all polarities combined)
  cat("Creating workflow summary for all polarities...\n")
  compound_summary3 <- dfmsm_merged  %>%
    group_by(inchikey1, spe, spe_frac, chromatography) %>%
    summarise(
      compoundname = toString(sort(unique(name))),
      hmdb_id = toString(sort(unique(hmdb_ids))),
      match_type = toString(unique(match_type)),
      
      # Scoring information
      rank_min = min(rank, na.rm=T),
      wscore_max = max(wscore,na.rm=T),
      sirius_score_max = max(sirius_score, na.rm=T),
      metfrag_score_max = max(metfrag_score, na.rm=T),
      sm_score = round(max(sm_score, na.rm=T), 4),
      inchikey = toString(sort(unique(inchikey))),
      
      # Experimental details
      assays = toString(sort(unique(assay))),
      extraction = toString(sort(unique(extraction))),
      spe = toString(sort(unique(spe))),
      spe_frac = toString(sort(unique(spe_frac))),
      chromatography = toString(sort(unique(chromatography))),
      measurement = toString(sort(unique(measurement))),
      polarity = 'All',  # Combined polarity label
      count = n_distinct(inchikey1),
      .groups = 'drop'
    )
  
  # Combine workflow summaries
  compound_summary4 <- rbind(compound_summary2, compound_summary3)
  
  # Match to standard mixture for workflow analysis
  compound_summary4 <- merge(compound_summary4, 
                            standardmix_details[,c('Compound', 'inchikey2d', 'superclass', 'class')],  
                            by.x = 'inchikey1', by.y = 'inchikey2d')
  
  # Clean up workflow labels
  cat("Cleaning workflow labels...\n")
  compound_summary4$extract_spe_lc <- paste(compound_summary4$extraction, 
                                            compound_summary4$spe,
                                            compound_summary4$spe_frac, 
                                            compound_summary4$chromatography)
  compound_summary4$extract_spe_lc <- gsub('D-', '', compound_summary4$extract_spe_lc)
  compound_summary4$extract_spe_lc <- gsub('APOL ', 'APOLAR ', compound_summary4$extract_spe_lc)
  compound_summary4$extract_spe_lc <- gsub('POL ', 'POLAR ', compound_summary4$extract_spe_lc)
  compound_summary4$polarity <- gsub('neg', 'NEGATIVE', compound_summary4$polarity)
  compound_summary4$polarity <- gsub('pos', 'POSITIVE', compound_summary4$polarity)
  
  # Define workflow order
  workflow_levels <- c("POLAR   AMD", "POLAR   PHE", "POLAR WCX 1 PHE", "POLAR WCX 1 AMD",
                      "POLAR WCX 2 PHE", "POLAR WCX 2 AMD", "POLAR WCX 3 PHE", "POLAR WCX 3 AMD",
                      "POLAR WCX 4 PHE", "POLAR WCX 4 AMD", "POLAR WAX 1 PHE", "POLAR WAX 1 AMD",
                      "POLAR WAX 2 PHE", "POLAR WAX 2 AMD", "POLAR WAX 3 PHE", "POLAR WAX 3 AMD",
                      "POLAR WAX 4 PHE", "POLAR WAX 4 AMD", "APOLAR C18 1 C30", "APOLAR C18 2 C30",
                      "APOLAR C18 3 C30", "APOLAR AMP 1 C30", "APOLAR AMP 2 C30", "APOLAR AMP 3 C30",
                      "APOLAR AMP 4 C30", "APOLAR   C30")
  
  compound_summary4$extract_spe_lc <- factor(compound_summary4$extract_spe_lc, levels = workflow_levels)
  
  # Order superclasses by frequency
  compound_superclass_count <- compound_summary4 %>%
    group_by(superclass) %>%
    summarise(count = n_distinct(inchikey1), n = n(), .groups = 'drop') %>%
    arrange(desc(n)) 
  
  compound_summary4$superclass <- factor(compound_summary4$superclass, 
                                         levels = compound_superclass_count$superclass)
  
  return(compound_summary4)
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

#' Create workflow bar plot (Figure S10a)
#' @param compound_summary4 Workflow summary dataframe
#' @param c25 Color palette
create_workflow_barplot <- function(compound_summary4, c25) {
  cat("Creating Figure S10a: Workflow analysis bar plot...\n")
  
  galaxy_msm_bar <- ggplot(compound_summary4, aes(x = extract_spe_lc, y = count, fill = superclass)) +
    geom_col(width = 0.7) +
    theme_bw() +
    scale_fill_manual("Superclass", values = c25) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(~polarity, ncol = 1) +
    labs(x = "Workflow", y = "Unique compound count")
  
  ggsave(file.path(OUTPUT_DIR, "FIG_S10a_galaxy_msms_workflow_bar.pdf"), 
         plot = galaxy_msm_bar, width = 10, height = 10)
  cat("Saved: FIG_S10a_galaxy_msms_workflow_bar.pdf\n")
  
  return(galaxy_msm_bar)
}

#' Create chemical class tree map (Figure S10b)
#' @param compound_summary4 Workflow summary dataframe
#' @param c25 Color palette
create_msm_treemap <- function(compound_summary4, c25) {
  cat("Creating Figure S10b: Chemical class tree map...\n")
  
  compound_class_count <- compound_summary4 %>%
    group_by(class) %>%
    summarise(count = n_distinct(inchikey1),
             superclass = unique(superclass),
             n = n(),
             .groups = 'drop') %>%
    arrange(desc(n))
  
  pdf(file.path(OUTPUT_DIR, 'FIG_S10b_treemap_msm.pdf'))
  treemap(compound_class_count,
          index = c("superclass", "class"),
          vSize = "count",
          type = "index",
          palette = c25,
          bg.labels = 244,
          align.labels = list(c("left", "top"), c("right", "bottom")),
          na.rm = FALSE)
  dev.off()
  cat("Saved: FIG_S10b_treemap_msm.pdf\n")
}

#' Create presence/absence plot (Figure S11)
#' @param compound_summary1 Main compound summary
create_presence_absence_plot <- function(compound_summary1) {
  cat("Creating Figure S11: Annotation method presence/absence plot...\n")
  
  # Prepare data for presence/absence plot
  summary_bool <- compound_summary1[,c('Compound', 'smbool', 'siriusbool',
                       'mfbool', 'galaxysmbool', 'mzcloudsmbool', 'gnpssmbool')] 
  colnames(summary_bool) <- c('Compound', 'SM (all)', 'Sirius CSI:FingerID',
                              'MetFrag', 'SM (Galaxy)', 'SM (mzCloud)', 'SM (GNPS workflow)')
  
  # Convert to long format and clean compound names
  booldf <- summary_bool %>% pivot_longer(cols = -Compound) 
  booldf$Compound <- str_trunc(booldf$Compound, 20)
  
  # Set factor levels for consistent ordering
  booldf$Compound <- factor(booldf$Compound, levels = sort(unique(booldf$Compound)))
  booldf$name <- factor(booldf$name, levels = c('SM (Galaxy)', 'SM (mzCloud)', 'SM (GNPS workflow)',
                                              'SM (all)', 'MetFrag', 'Sirius CSI:FingerID'))
  
  # Calculate method summary statistics
  bool_summary <- booldf %>%
    group_by(name) %>%
    summarise(sum = sum(value), .groups = 'drop')
  
  cat("Annotation method coverage:\n")
  for(i in 1:nrow(bool_summary)) {
    cat(sprintf("  %s: %d compounds\n", bool_summary$name[i], bool_summary$sum[i]))
  }
  
  # Create presence/absence plot
  presence_plot <- ggplot(booldf, aes(x = name, y = Compound, shape = factor(value))) +
    theme_bw() +
    geom_point(size = 2) +
    scale_shape_manual(values = c(1, 16)) +  # 1 = open circle, 16 = filled circle
    xlab("Annotation method") + 
    ylab("Metabolite Reference Standard") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, "FIG_S11_presence_absence_match_type_msm.pdf"), 
         plot = presence_plot, width = 4, height = 6.5)
  cat("Saved: FIG_S11_presence_absence_match_type_msm.pdf\n")
  
  return(presence_plot)
}

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

#' Run complete MSM analysis pipeline
run_msm_analysis <- function() {
  cat("=== Starting MSM Analysis ===\n")
  start_time <- Sys.time()
  
  # 1. Load data
  data_list <- load_msm_data()
  standardmix_details <- data_list$standardmix_details
  dfmsm_merged <- data_list$dfmsm_merged
  
  # 2. Create main annotation summary
  summary_results <- create_msm_summary(dfmsm_merged, standardmix_details)
  compound_summary1 <- summary_results$compound_summary
  
  # 3. Create workflow analysis
  compound_summary4 <- create_workflow_summary(dfmsm_merged, standardmix_details)
  
  # 4. Generate all figures
  cat("\n=== Generating Figures ===\n")
  create_workflow_barplot(compound_summary4, c25)
  create_msm_treemap(compound_summary4, c25)
  create_presence_absence_plot(compound_summary1)
  
  # 5. Final summary
  end_time <- Sys.time()
  runtime <- round(difftime(end_time, start_time, units = "secs"), 2)
  
  cat("\n=== MSM Analysis Complete ===\n")
  cat("Key Results:\n")
  cat(sprintf("  Total standards in mixture: %d\n", 48))
  cat(sprintf("  Standards successfully annotated: %d (%.1f%%)\n", 
              summary_results$total_matched, summary_results$perc_matched))
  cat(sprintf("  Standards not detected: %d\n", nrow(summary_results$missing_standards)))
  cat("\nRuntime:", runtime, "seconds\n")
  cat("\nOutputs saved:\n")
  cat("  - msm_annotations_summary.csv\n")
  cat("  - FIG_S10a_galaxy_msms_workflow_bar.pdf\n")
  cat("  - FIG_S10b_treemap_msm.pdf\n")
  cat("  - FIG_S11_presence_absence_match_type_msm.pdf\n")
  cat("=============================================================================\n")
}

# Execute Analysis ----
# Run the analysis when script is sourced
run_msm_analysis()
