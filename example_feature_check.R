library(tidyverse)
library(structToolbox)
library(xcms)

feature_check <- function(galaxy_file_pl_url,
                                 galaxy_file_xdata_url, 
                                 sample_metadata_pth, 
                                 merged_annotation_list_pth,
                                 out_dir='.', 
                                 assay_name){
  
  dir.create(out_dir)
  
  sample_metadata <- read.table(sample_metadata_pth, sep = "\t", header = 1)
  blank_file_names <- sample_metadata[sample_metadata$blank=="yes",]$original_filename
  
  blank_file_names <- gsub("-", ".", blank_file_names)
  blank_file_names <- gsub(".raw", "", blank_file_names, ignore.case = T)
  blank_file_names <- gsub(".mzML", "", blank_file_names, ignore.case = T)
  blank_file_names <- unique(blank_file_names)
  
  
  
  # download relevant files
  peaklist_pth <- file.path(out_dir, paste(assay_name, "peaklist.tsv", sep="_"))
  download.file(galaxy_file_pl_url, destfile =peaklist_pth, method="curl")
  
  xdata_pth <- file.path(out_dir, paste(assay_name, "xdata.RData", sep="_"))
  download.file(galaxy_file_xdata_url, destfile =xdata_pth, method="curl")
  
  ##############################################################################
  # First get XCMS xset object and get full peak widths
  ##############################################################################
  load(xdata_pth)
  
  # Get the legacy xcmsSet object
  xdata <- as(xdata, "xcmsSet")
  
  xcms::sampclass(xdata) <- xdata@phenoData$sample_group
  
  rt_min <- xcms::groupval(xdata, method = "medret", value = "rtmin", intensity = "into")
  rt_max <-  xcms::groupval(xdata, method = "medret", value = "rtmax", intensity = "into")
  
  xcms_names <- xcms::groupnames(xdata)
  full_rtmin_rtmax <- data.frame(name=xcms_names)
  # warnings expected if no sample feature available to calculate rtmin and rtmax - convert them to NA
  safe_min <- function(x) if (all(is.na(x))) NA else min(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE)
  full_rtmin_rtmax$rtmin_full <- apply(rt_min[, !grepl("blank", colnames(rt_min), ignore.case = TRUE)], 1, safe_min)
  full_rtmin_rtmax$rtmax_full <- apply(rt_max[, !grepl("blank", colnames(rt_max), ignore.case = TRUE)], 1, safe_max)
  
  xcms_grp_info <- xcms::groups(xdata)[,c("mzmed", "rtmed")]
  colnames(xcms_grp_info) <- c("mzmed", "rtmed")
  full_rtmin_rtmax <- cbind(xcms_grp_info, full_rtmin_rtmax)
  
  ##############################################################################
  # Then load the xcms peak list generated from the XCMS W4M galaxy tool
  # And save as StructToolbox object for blank filtering
  ##############################################################################
  peak_matrix <- read.table(peaklist_pth, sep = "\t", header = 1, row.names = 1, check.names = F)
  
  peak_matrix[peak_matrix==0] <- NA
  
  VM = data.frame(name=rownames(peak_matrix ))
  
  VM <- merge(VM, full_rtmin_rtmax, all.x = T, sort = F)
  
  rownames(VM) <- VM$name
  
  SM = data.frame(id=colnames(peak_matrix ))
  SM$x <- "sample"
  # Use file name to check for blanks
  SM$x[grepl("blank", SM$id, ignore.case = T)] <- "blank"
  # AND use sample metadata (as some cases recorded as just B_). Include both
  # approaches just to ensure blanks are chosen
  SM$x[SM$id %in% blank_file_names] <- "blank" 
  
  SM$x <- factor(SM$x, levels = c("blank", "sample"))
  rownames(SM) <- SM$id
  
  DE <- 
    DatasetExperiment(
      data = as.data.frame(t(peak_matrix)), # samples in ROWS
      sample_meta = SM,
      variable_meta = VM
    )
  
  ##############################################################################
  # Perform blank filtering
  ##############################################################################
  M = blank_filter(
    fold_change = 10,
    blank_label = "blank",
    qc_label = "sample",
    factor_name = "x",
    fraction_in_blank=0)
  
  M = model_apply(M,DE)
  
  DE_filt <- struct::predicted(M)
  
  saveRDS(DE_filt, file.path(out_dir, paste(assay_name, "DE_blank_filtered.RDS", sep="_")))
  write.csv(DE_filt$data, file.path(out_dir, paste(assay_name, "blank_filtered_peak_matrix.csv", sep="_")))
  write.csv(DE_filt$variable_meta, file.path(out_dir, paste(assay_name, "blank_filtered_variable_metadata.csv", sep="_")))
  
  
  ##############################################################################
  # Find Galaxy annotations that can match to blank filtered features
  ##############################################################################
  # for specific assay
  all_annotations <-  read.csv(unzip(merged_annotation_list_pth, "merged_annotations_all_classified.csv"))
  all_annotations_assay_all <- all_annotations %>% filter(assay == assay_name) 
  all_annotations_assay <- all_annotations_assay_all %>% filter(measurement == "LCMSMS")
  
  # just xcms
  all_annotations_xcms <- all_annotations_assay %>% filter(grp_name != "")
  
  blank_annotations <- all_annotations_xcms[!all_annotations_xcms$grp_name %in% DE_filt$variable_meta$name,]
  
  # Save the annotations IDs that pass
  all_annos_xcms_passed <- all_annotations_xcms[all_annotations_xcms$grp_name %in% DE_filt$variable_meta$name,]
  xcms_pass_merged_filt_id <- all_annos_xcms_passed$merged_id
  write.csv(all_annos_xcms_passed, file.path(out_dir, paste(assay_name, "xcms_passed_annos.csv", sep="_")))
  
}


sample_metadata_pth <- "input/input_for_feature_check/GalaxyNone-[samplelist_dma_daphnia_magna.tabular].tabular"
galaxy_references <- read.csv('input/input_for_feature_check/galaxy_peaklist_references.csv')
merged_annotation_list_pth <- "input/input_for_summary_plots/merged_annotations_all_classified.zip"


for (i in seq_len(nrow(galaxy_references))) {
  
  assay_name <- galaxy_references$assay_name[i]
  print(paste(i, assay_name))
  if (assay_name == "") {
    break
  }

  galaxy_file_pl_url <- galaxy_references$galaxy_file_pl_url[i]
  galaxy_file_xdata_url <- galaxy_references$galaxy_file_xdata_url[i]
  out_dir <- file.path("output", assay_name)
  feature_check(
    galaxy_file_pl_url,
    galaxy_file_xdata_url,
    sample_metadata_pth,
    merged_annotation_list_pth,
    out_dir,
    assay_name
  )
  
}

