dir_list <- c(
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GastricCancerCellLines/SNU668/A00_IncucyteRawDataLiveDead_varyGlucose_250324/raw_images",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/MCF10A/A00_IncucyteRawDataLiveDead_varyGlucose_241015/raw_images/",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/MCF10A/A00b_IncucyteRawDataLiveDead_varyGlucose_241015",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/MCF10A/A00c_IncucyteRawDataLiveDead_varyGlucose_241015/raw_images",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/MDA-MB-231/B00_IncucyteRawDataLiveDead_varyGlucose_250213/raw_images",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/M00b_IncucyteRawDataLiveDead_varyGlucose/raw_images",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/C00_IncucyteRawDataLiveDead_varyGlucose/mCherry Orange calibrated all cells",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/C00_IncucyteRawDataLiveDead_varyGlucose/NIR calibrated dead cells",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/C00_IncucyteRawDataLiveDead_varyGlucose/phase_orig",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/I00_IncucyteRawDataLiveDead_varyGlucose/Red",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/I00_IncucyteRawDataLiveDead_varyGlucose/NIR",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/I00_IncucyteRawDataLiveDead_varyGlucose/phase"
)


summarize_files <- function(dir) {
  cat("\nDirectory:", dir, "\n")
  
  # 1) subdirectories
  subs <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  cat("Subdirs:", if(length(subs)) paste(subs, collapse = ", ") else "(none)", "\n")
  
  # 2) extension counts
  files <- list.files(dir, recursive = FALSE, all.files = FALSE)
  exts  <- tolower(tools::file_ext(files))
  cat("Extension counts:\n")
  print(table(exts))
  
  # 3) sample filenames
  cat("Sample files:\n")
  print(head(files, 10))
  
  # 4) split on '_' and extract 4th-from-last field
  frags <- strsplit(files, "_", fixed = TRUE)
  fourth_from_last <- sapply(frags, function(x) {
    if (length(x) >= 4) x[length(x) - 3] else NA_character_
  })
  uniques <- unique(na.omit(fourth_from_last))
  
  cat("Unique 4th-from-last fields:\n")
  if (length(uniques)) {
    cat(paste0("  - ", uniques), sep = "\n")
  } else {
    cat("  (none)\n")
  }
}

invisible(lapply(dir_list, summarize_files))


all_names <- do.call(c,lapply(dir_list, function(path){
  prefix <- paste(gsub("_","-",unlist(strsplit(path,split="/"))[8:9]),collapse="_")
  paste(prefix,list.files(path),sep="_")
}))



dups <- all_names[duplicated(all_names)]

## ------------------------------------------------------------------
## 0.  destination for the links
dest <- "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GlucoseStarvation/all_raw"
if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)

## ------------------------------------------------------------------
## 1.  channel → class dictionary  (add rows if a new experiment appears)
channel_map <- list(
  "SNU668_A00-IncucyteRawDataLiveDead-varyGlucose-250324" =
    c("Green-nucleus" = "alive",   "NIR-dead" = "dead",
      "phase" = "phase",           "Phase"    = "phase"),
  
  "MCF10A_A00-IncucyteRawDataLiveDead-varyGlucose-241015" =
    c("Red-Nuclei"  = "alive", "NIR-dead" = "dead",
      "Phase" = "phase"),
  "MCF10A_A00b-IncucyteRawDataLiveDead-varyGlucose-241015" =
    c("Red-Nuclei"  = "alive", "NIR-dead" = "dead",
      "Phase" = "phase"),  
  "MCF10A_A00c-IncucyteRawDataLiveDead-varyGlucose-241015" =
    c("Red-Nuclei"  = "alive", "NIR-dead" = "dead",
      "phase" = "phase"),
  
  "MDA-MB-231_B00-IncucyteRawDataLiveDead-varyGlucose-250213" =
    c("red"  = "aliveC1", "green" = "aliveC2",
      "NIR-dead" = "dead", "phase" = "phase", "Phase" = "phase"),
  
  "SUM-159_M00a-IncucyteRawDataLiveDead-varyGlucose" =
    c("Red-nuclei" = "alive", "Green-dead" = "dead", "phase" = "phase"),
  
  "SUM-159_M00b-IncucyteRawDataLiveDead-varyGlucose" =
    c("Red-nuclei" = "alive", "NIR-dead"   = "dead",
      "Phase" = "phase", "phase" = "phase"),
  
  "SUM-159_C00-IncucyteRawDataLiveDead-varyGlucose" =
    c("OrangeCalibrated" = "alive", "NIRCalibrated" = "dead",
      "phase" = "phase"),
  
  "SUM-159_I00-IncucyteRawDataLiveDead-varyGlucose" =
    c("Red" = "alive", "NIR" = "dead",
      "Phase" = "phase", "phase" = "phase")
)

## helper that extracts the same prefix you already use ----------------
prefix_from_dir <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  paste(gsub("_", "-", parts[8:9]), collapse = "_")
}

## helper that pulls the raw channel tag from a filename ----------------
raw_tag <- function(fname) {
  x <- strsplit(fname, "_", fixed = TRUE)[[1]]
  if (length(x) >= 4) x[length(x) - 3] else NA_character_
}

## ------------------------------------------------------------------
## 2.  walk through every directory and link its files
for (dir in dir_list) {
  
  prefix <- prefix_from_dir(dir)
  files  <- list.files(dir, full.names = TRUE)
  map    <- channel_map[[prefix]]
  
  
  
  if (is.null(map)) {
    message("No mapping for ", prefix, "; skipped.")
    next
  }
  
  for (src in files) {
    fname  <- basename(src)
    tag    <- raw_tag(fname)
    
    
    if (is.na(tag)) {
      message("Could not extract tag from: ", fname, "; skipped.")
      next
    }
    
    if (!tag %in% names(map)) {
      message("Tag ‘", tag, "’ not in map for ", prefix, "; skipped.")
      next
    }
    class  <- map[[tag]]
    new_fname <- sub(tag, class, fname, fixed = TRUE)
    ## keep your original prefix in front
    link_name <- file.path(dest, paste0(prefix, "_", new_fname))
    file.symlink(src, link_name)
  }
}