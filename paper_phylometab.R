# =============================================================================
# Daphnia magna DMA Paper - Phylogenetic Metabolomics Analysis (Clean)
# =============================================================================
# This script creates phylogenetic metabolomics visualizations by combining
# metabolite chemical classifications with species phylogenetic trees.
#
# Main outputs:
# - Figure 6: Combined phylogenetic-metabolomic heatmap with chemical tree,
#   species phylogenetic tree, and database mapping results
# =============================================================================

cat("=============================================================================\n")
cat("PHYLOGENETIC METABOLOMICS ANALYSIS\n") 
cat("=============================================================================\n")

# Load required libraries ----
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ape)
  library(aplot)
  library(ggtree)
  library(stringi)
  library(ggpubr)
  library(data.table)
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
  "grey70", "skyblue2", "gold1", "#E31A1C", "palegreen2",
  "#CAB2D6", "#FDBF6F", "grey15", "khaki2", "maroon", 
  "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", 
  "green1", "yellow4", "yellow3", "darkorange4", "brown", "grey30"
)

# =============================================================================
# DATA LOADING AND PREPROCESSING FUNCTIONS
# =============================================================================

#' Load and preprocess DMA annotation data for phylometab analysis
#' @return Preprocessed DMA dataframe with  scores
load_phylometab_data <- function() {
  cat("\n=== Loading Phylometab Data ===\n")
  
  # Load DMA annotation summary
  cat("Loading DMA annotation summary...\n")
  dma_df <- read.csv("input/input_for_phylometab_plot/Supplementary Data - 2 - Metabolite Annotation Summary.csv", 
                     stringsAsFactors = FALSE)
  
  # Filter out empty superclass entries
  dma_df <- dma_df[dma_df$superclass != '', ]
  cat("Total DMA annotations loaded:", nrow(dma_df), "\n")
  
  # Select relevant columns and calculate reliability score
  dma_df <- dma_df[, c('kingdom', 'superclass', 'class', 'subclass', 'inchikey', 'inchikey1', 
                       'gcmsbool', 'lcmsdimsbool', 'nmrbool', 'smbool', 'mfbool', 'siriusbool')]
  
  # Calculate reliability score as sum of detection methods
  # IMPORTANT - WE DO NOT USE RELIABILTIY SCORE FOR HEATMAP - JUST CONVERT TO PRESENCE ABSENCE 
  #           - initially used as alternative to count
  dma_df$reliability_score <- dma_df$gcmsbool + dma_df$nmrbool + dma_df$smbool + 
                              dma_df$mfbool + dma_df$siriusbool
  dma_df$reliability_score <- as.character(dma_df$reliability_score)
  
  # Make names compatible with tree construction
  dma_df <- make_classyfire_tree_compatible(dma_df)
  
  return(dma_df)
}

#' Make ClassyFire names compatible with phylogenetic tree construction
#' @param df Dataframe with ClassyFire classifications
#' @return Cleaned dataframe with tree-compatible names
make_classyfire_tree_compatible <- function(df) {
  cat("Making ClassyFire names tree-compatible...\n")
  
  # Store original full names before cleaning
  superclass_name_full <- df$superclass
  kingdom_full <- df$kingdom
  
  # Remove problematic characters for tree construction
  df <- as.data.frame(sapply(df, function(x) gsub(" ", "", x)))  
  df <- as.data.frame(sapply(df, function(x) gsub("\'", "", x)))
  df <- as.data.frame(sapply(df, function(x) gsub("->", "", x)))
  df <- as.data.frame(sapply(df, function(x) gsub("[()]", "", x)))
  df <- as.data.frame(sapply(df, function(x) gsub(",", "", x)))
  df <- as.data.frame(sapply(df, function(x) gsub("-", "", x)))
  
  # Convert to factors for tree construction
  df$kingdom <- as.factor(df$kingdom)
  df$superclass <- as.factor(df$superclass)
  df$class <- as.factor(df$class)
  df$subclass <- as.factor(df$subclass)
  df$inchikey <- as.factor(df$inchikey)
  
  # Restore original names for display
  df$superclass_name_full <- superclass_name_full
  df$kingdom_full <- kingdom_full
  
  return(df)
}

# =============================================================================
# DATABASE MAPPING FUNCTIONS
# =============================================================================

#' Map DMA annotations to external databases and species
#' @param dma_df DMA annotation dataframe
#' @param kegg_hmdb_pth Path to KEGG/HMDB mapping file
#' @param mtox_pth Path to MTox mapping file  
#' @param all_daphnia_pth Path to Daphnia ChEBI mapping file
#' @param all_species_pth Path to species mapping file
#' @return List containing mapped dataframe and species summary
map_dma_to_databases <- function(dma_df, kegg_hmdb_pth, mtox_pth, all_daphnia_pth, all_species_pth) {
  cat("\n=== Mapping to External Databases ===\n")
  
  # IMPORTANT - WE DO NOT USE RELIABILTIY SCORE FOR HEATMAP - JUST CONVERT TO PRESENCE ABSENCE 
  #           - initially used as alternative to count
  
  # Initialize database mapping columns
  dma_df$kegg <- '0'
  dma_df$hmdb <- '0'
  dma_df$mtox <- '0'
  dma_df$chebi_daphnia_all <- '0'
  
  # Map to KEGG and HMDB
  cat("Mapping to KEGG and HMDB databases...\n")
  mapping_hmdb_kegg <- read.csv(unz(kegg_hmdb_pth, paste(tools::file_path_sans_ext(basename(kegg_hmdb_pth)), '.csv', sep="")))
  mapping_hmdb_kegg$inchikey1 <- sapply(str_split(mapping_hmdb_kegg$inchikey, '-'), function(x) x[1])
  
  # KEGG mapping
  kegg_bool <- dma_df$inchikey1 %in% mapping_hmdb_kegg[!is.na(mapping_hmdb_kegg$kegg), ]$inchikey1
  dma_df[kegg_bool, ]$kegg <- dma_df[kegg_bool, ]$reliability_score
  
  # HMDB mapping
  hmdb_bool <- dma_df$inchikey1 %in% mapping_hmdb_kegg[!is.na(mapping_hmdb_kegg$hmdb), ]$inchikey1
  dma_df[hmdb_bool, ]$hmdb <- dma_df[hmdb_bool, ]$reliability_score
  
  
  # Map to MTox
  cat("Mapping to MTox database...\n")
  mapping_mtox <- read.csv(mtox_pth)
  mtox_bool <- dma_df$inchikey1 %in% mapping_mtox$InchiKey1
  dma_df[mtox_bool, ]$mtox <- dma_df[mtox_bool, ]$reliability_score
  
  # Map to Daphnia ChEBI
  cat("Mapping to Daphnia ChEBI database...\n")
  mapping_daphnia_all <- read.csv(all_daphnia_pth)
  daphnia_all_bool <- dma_df$inchikey1 %in% mapping_daphnia_all$INCHIKEY1
  dma_df[daphnia_all_bool, ]$chebi_daphnia_all <- dma_df[daphnia_all_bool, ]$reliability_score
  
  # Map to species data
  species_results <- map_to_species_data(dma_df, all_species_pth)
  dma_df_with_species <- species_results$dma_df
  species_summary <- species_results$species_summary
  
  return(list(dma_df_with_species, species_summary))
}

#' Map DMA data to species occurrence data
#' @param dma_df DMA annotation dataframe
#' @param all_species_pth Path to species data file
#' @return List with updated dataframe and species summary
map_to_species_data <- function(dma_df, all_species_pth) {
  cat("Mapping to species occurrence data...\n")
  
  # Load species data
  species_all <- fread(all_species_pth, sep = ',')
  
  # Define relevant species for analysis
  species_keep <- c('Homo sapiens', 'Caenorhabditis elegans', 'Escherichia coli',
                    'Mus musculus', 'Daphnia pulex', 'Daphnia galeata', 'Daphnia magna',
                    'Daphnia tenebrosa', 'Rattus norvegicus', 'Drosophila melanogaster',
                    'Danio rerio', 'Xenopus laevis', 'Xenopus', 'Mus musculus domesticus',
                    'Tabernaemontana elegans', 'Escherichia coli O25')
  
  # Filter species: keep top 50 by count or those in our relevant list
  species_summary <- species_all %>% 
    group_by(SPECIES_TEXT) %>% 
    filter(n() > 50 | SPECIES_TEXT %in% species_keep) %>% 
    summarize(count = n(), .groups = 'drop')
  
  # Remove problematic species for dendrogram construction
  remove_species <- c('Mycoplasma genitalium', "Homo Sapiens", "Cordyceps sinensis", 
                      "Alstonia spatulata", "Dysoxylum lenticellatum", "Juglans sinensis", 
                      "Panax japonicus var. major", "Piper boehmeriaefolium", "Sinocalamus affinis")
  species_summary <- species_summary %>% filter(!SPECIES_TEXT %in% remove_species)
  
  cat("Selected species for analysis:", nrow(species_summary), "\n")
  
  # Export species list for dendrogram creation - only required if a different species list needed
  # - For this script we have already expored and produced a dendrogram/tree in newick format for the species 
  #   of interest
  # write.csv(species_summary, file.path(OUTPUT_DIR, "species_to_get_dendrogram.csv"), row.names = FALSE)
  # cat("Exported species list for dendrogram creation\n")
  
  # Filter species data and prepare for mapping
  species_filt <- species_all %>% 
    group_by(SPECIES_TEXT) %>% 
    filter(SPECIES_TEXT %in% unique(species_summary$SPECIES_TEXT))
  
  species_filt$inchikey1 <- sapply(str_split(species_filt$inchikey, '-'), function(x) x[1])
  
  # Create species mapping function
  species_check <- function(x, check_df) {
    check_bool <- check_df$inchikey1 %in% x$inchikey1
    species_x <- rep(0, nrow(check_df))
    species_x[check_bool] <- check_df[check_bool, ]$reliability_score 
    return(tibble(species_x))
  }
  
  # Apply mapping to all species
  species_all_mapping <- species_filt %>%
    group_by(SPECIES_TEXT) %>%
    group_map(~species_check(.x, check_df = dma_df)) %>% 
    bind_cols()
  
  # Clean species names for column headers
  species_summary$SPECIES_TEXT_ <- gsub(" ", "_", species_summary$SPECIES_TEXT)
  colnames(species_all_mapping) <- species_summary$SPECIES_TEXT_
  
  # Combine with original dataframe
  dma_df_with_species <- cbind(dma_df, species_all_mapping)
  
  return(list(dma_df = dma_df_with_species, species_summary = species_summary))
}

# =============================================================================
# TREE CONSTRUCTION FUNCTIONS
# =============================================================================

#' Create phylogenetic tree from newick file
#' @param newick_pth Path to newick format tree file
#' @return ggplot tree object
create_taxa_tree <- function(newick_pth) {
  cat("\n=== Creating Taxa Tree ===\n")
  cat("Loading phylogenetic tree from:", newick_pth, "\n")
  
  phylogen_tree <- read.tree(newick_pth)
  p <- ggtree(phylogen_tree) + layout_dendrogram()
  
  cat("Taxa tree created successfully\n")
  return(p)
}

#' Create chemical classification tree from ClassyFire data
#' @param df Dataframe with ClassyFire classifications
#' @return List containing tree plot, node details, and colors
create_classyfire_tree <- function(df) {
  cat("\n=== Creating Chemical Classification Tree ===\n")
  
  # Handle special cases in classification
  df$superclass <- as.character(df$superclass)
  df[df$class == 'Organonitrogencompounds', ]$superclass <- 'Organonitrogencompounds'
  df[df$superclass == 'Hydrocarbons', ]$superclass <- 'Unsaturatedhydrocarbons'
  df$superclass <- as.factor(df$superclass)
  
  # Summarize superclasses by frequency
  superclass_sum <- df %>% 
    filter(kingdom == "Organiccompounds") %>% 
    group_by(superclass) %>%
    summarize(count = n(), full_name = unique(superclass_name_full), .groups = 'drop') %>% 
    arrange(desc(count))
  
  cat("Top chemical superclasses:\n")
  for(i in 1:min(5, nrow(superclass_sum))) {
    cat(sprintf("  %s: %d compounds\n", superclass_sum$superclass[i], superclass_sum$count[i]))
  }
  
  # Create phylogenetic tree from classification hierarchy
  tree <- as.phylo(~superclass/class/subclass/inchikey, data = df)
  
  # Highlight top superclasses
  max_superclass_num <- 8
  node_update <- data.frame(
    node = unlist(lapply(as.character(superclass_sum[1:max_superclass_num, ]$superclass), 
                        function(x){ nodeid(tree, x) })), 
    full_name = superclass_sum[1:max_superclass_num, ]$full_name,
    color = c25[1:max_superclass_num],
    lty = rep(0, max_superclass_num)
  )
  
  # Remove any nodes that couldn't be found
  node_update <- node_update[!is.na(node_update$node), ]
  hilight_colours <- node_update[order(node_update$full_name), ]$color
  
  # Prepare node details for plotting
  node_details <- data.frame(
    node = 1:treeio::getNodeNum(tree), 
    color = 'black', 
    label2 = '', 
    lty = 1
  )
  node_details[node_update$node, ]$color <- node_update$color
  node_details[node_update$node, ]$label2 <- as.character(node_update$full_name)
  node_details[node_update$node, ]$lty <- node_update$lty
  
  # Create tree plot with highlights
  p_1 <- ggtree(tree, layout = "slanted")
  p_2 <- p_1 + geom_hilight(data = node_update, aes(node = node, fill = full_name), 
                           align = "left", alpha = 0.80) + 
              scale_fill_manual(values = hilight_colours) + 
              theme(legend.position = "none")
  p_3 <- p_2 %<+% node_details + geom_text(aes(label = I(label2)), hjust = 1, size = 3.5, color = 'black')
  
  cat("Chemical classification tree created with", length(hilight_colours), "highlighted superclasses\n")
  
  return(list(p_3, node_details, hilight_colours))
}

# =============================================================================
# HEATMAP VISUALIZATION FUNCTIONS
# =============================================================================

#' Create comprehensive phylogenetic-metabolomic heatmap
#' @param dma_tree Chemical classification tree results
#' @param dma_df Complete DMA dataframe with mappings
#' @param species_summary Species summary statistics
#' @param newick_pth Path to phylogenetic tree file
#' @return Combined heatmap visualization
create_phylometab_heatmap <- function(dma_tree, dma_df, species_summary, newick_pth) {
  cat("\n=== Creating Phylogenetic-Metabolomic Heatmap ===\n")
  
  metabolite_tree <- dma_tree[[1]]
  
  # Create taxa tree and get species order
  p_tax <- create_taxa_tree(newick_pth)
  ordered_species <- rev(get_taxa_name(p_tax))
  ordered_species <- ordered_species[ordered_species %in% colnames(dma_df)]
  
  # Prepare species occurrence data
  species_long <- gather(dma_df[, c('inchikey', 'subclass', 'class', ordered_species)], 
                        category, value, ordered_species[1]:ordered_species[length(ordered_species)])
  species_long$value <- as.numeric(species_long$value)
  species_long$count <- 0
  # IMPORTANT - WE DO NOT USE RELIABILTIY SCORE FOR HEATMAP - JUST CONVERT TO PRESENCE ABSENCE
  species_long$count[species_long$value > 0] <- 1  # Convert to presence/absence
  
  # Calculate counts per subclass
  species_long <- species_long %>%
    group_by(subclass, category) %>%
    mutate(total_count = sum(count)) %>%
    ungroup() %>%
    arrange(category, class, subclass)
  
  # Group compounds for visualization (every 25 rows)
  n <- 25
  species_long$sum_n_rows <- rep(tapply(species_long$count, (seq_along(species_long$count) - 1) %/% n, sum), 
                                each = n)[1:nrow(species_long)]
  species_long$category <- factor(species_long$category, levels = ordered_species)
  
  # Create main species heatmap
  p_main <- ggplot(species_long, aes(x = category, y = inchikey, fill = sum_n_rows)) +
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "#6666FF") +  
    theme_minimal() + xlab(NULL) + ylab(NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Prepare database mapping data
  misc_long <- gather(dma_df[, c('inchikey', 'hmdb', 'kegg', 'mtox', 'class', 'subclass')], 
                     category, value, hmdb:mtox)
  misc_long$value <- as.numeric(misc_long$value)
  misc_long$count <- 0
  misc_long$count[misc_long$value > 0] <- 1
  
  # Calculate database mapping counts
  misc_long <- misc_long %>%
    group_by(subclass, category) %>%
    mutate(total_count = sum(count)) %>%
    ungroup() %>%
    arrange(category, class, subclass)
  
  misc_long$sum_n_rows <- rep(tapply(misc_long$count, (seq_along(misc_long$count) - 1) %/% n, sum), 
                             each = n)[1:nrow(misc_long)]
  
  # Create database mapping heatmap
  p_right <- ggplot(misc_long, aes(x = category, y = inchikey, fill = sum_n_rows)) +
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "#6666FF") +  
    theme_minimal() + xlab(NULL) + ylab(NULL) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Create species summary bar plots
  colnames(species_summary) <- c('category', 'count')
  species_summary$category <- gsub(" ", "_", species_summary$category)
  
  p_bottom1 <- ggplot(species_summary, aes(x = category, y = count)) + 
    ggtitle('Total count of metabolites published in ChEBI') +
    geom_bar(stat = "identity") + xlab(NULL) + ylab(NULL) + 
    theme_minimal() + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_y_reverse()
  
  species_long_numeric <- species_long
  species_long_numeric$value <- as.numeric(species_long$value)
  p_bottom2 <- ggplot(species_long_numeric, aes(x = category, y = value)) + 
    ggtitle('Count of matches to D. magna DMA annotations') +
    geom_bar(stat = 'identity') + xlab(NULL) + ylab(NULL) + 
    theme_minimal() + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_y_reverse()
  
  # Combine all plot elements
  cat("Combining plot elements...\n")
  heatmap_combined <- p_main %>% 
    insert_left(metabolite_tree, width = 0.6) %>% 
    insert_right(p_right, width = 0.2) %>% 
    insert_top(p_tax, height = 0.25) %>%
    insert_bottom(p_bottom1, height = 0.1) %>%
    insert_bottom(p_bottom2, height = 0.1)
  
  # Add legends
  heatmap_combined <- plot_list(
    heatmap_combined, 
    as_ggplot(get_legend(p_main + theme(legend.position = "bottom"))),
    as_ggplot(get_legend(p_bottom1 + theme(legend.position = "bottom"))), 
    ncol = 1, heights = c(0.9, 0.05, 0.01)
  )
  
  cat("Phylogenetic-metabolomic heatmap created successfully\n")
  return(heatmap_combined)
}

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

#' Run complete phylogenetic metabolomics analysis
run_phylometab_analysis <- function() {
  cat("=== Starting Phylogenetic Metabolomics Analysis ===\n")
  start_time <- Sys.time()
  
  # 1. Load and preprocess data
  dma_df <- load_phylometab_data()
  
  # 2. Map to external databases and species
  mapping_results <- map_dma_to_databases(
    dma_df,
    kegg_hmdb_pth = "input/input_for_phylometab_plot/pubchem_kegg_hmdb_expanded.zip",
    mtox_pth = "input/input_for_phylometab_plot/Supplementary Data - 5 - Matching to MTox.csv",
    all_daphnia_pth = "input/input_for_phylometab_plot/Supplementary Data - 4 - Matching to Daphnia ChEBI.csv",
    all_species_pth = "input/input_for_phylometab_plot/chebi_with_inchikey_source_classyfire.csv"
  )
  
  dma_df_mapped <- mapping_results[[1]]
  species_summary <- mapping_results[[2]]
  
  # 3. Create chemical classification tree
  dma_tree <- create_classyfire_tree(dma_df_mapped)
  
  # 4. Create combined phylogenetic-metabolomic heatmap
  cat("\n=== Generating Final Visualization ===\n")
  phylometab_heatmap <- create_phylometab_heatmap(
    dma_tree, 
    dma_df_mapped, 
    species_summary, 
    newick_pth = 'input/input_for_phylometab_plot/phyloT_generated_tree_1734701763_newick.txt'
  )
  
  # 5. Save final output
  cat("Saving phylogenetic metabolomics plot...\n")
  ggsave(file.path(OUTPUT_DIR, "FIG_6_phylomet.pdf"), phylometab_heatmap, width = 10, height = 10)
  
  # 6. Analysis summary
  end_time <- Sys.time()
  runtime <- round(difftime(end_time, start_time, units = "secs"), 2)
  
  cat("\n=== Phylogenetic Metabolomics Analysis Complete ===\n")
  cat("\nRuntime:", runtime, "seconds\n")
  cat("\nOutput saved:\n")
  cat("  - FIG_6_phylomet.pdf\n")
  cat("=============================================================================\n")
}

# Execute Analysis ----
# Run the analysis when script is sourced
run_phylometab_analysis()
