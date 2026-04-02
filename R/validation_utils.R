## Utilities for lightweight reproducibility manifests stored in a dedicated
## directory. Expected layout:
## - paths.txt: one-line JSON object or vector of canonical file paths to track
## - hashes.txt: sha256 hashes written once from those paths
##
## Typical use:
## 1. Create a reproducibility directory for one locked analysis bundle.
## 2. Write its canonical tracked files into paths.txt.
## 3. Call save_hashes() once to freeze the recorded hashes.
## 4. Call check_hashes() later to verify the tracked artifacts still match.

load_paths <- function(reproducibility_dir){
    # Load the canonical tracked-file list for one reproducibility bundle.
    paths <- readLines(file.path(reproducibility_dir,"paths.txt"))
    paths <- jsonlite::fromJSON(paths)
    paths
}

hash_paths <- function(paths){
  # Compute a sha256 hash for each tracked path in the declared order.
  hashes <- sapply(unlist(paths),digest::digest,algo="sha256",file=TRUE)
  as.character(hashes)
}

save_hashes <- function(reproducibility_dir){
  # Freeze the current hashes for the bundle. Refuse to overwrite an existing
  # manifest so the directory stays append-only unless edited intentionally.
  output_file <- file.path(reproducibility_dir,"hashes.txt")
  if(file.exists(output_file)) stop("hashes already exists!!")
  
  paths <- load_paths(reproducibility_dir)
  
  hashes <- hash_paths(paths)
  writeLines(hashes,output_file)
  
}

check_hashes <- function(reproducibility_dir){
  # Recompute hashes from the current files and compare them against the saved
  # manifest. Returns TRUE only when every tracked entry matches exactly.
  hash_file <- file.path(reproducibility_dir,"hashes.txt")
  stored_hashes <- readLines(hash_file)
  
  paths <- load_paths(reproducibility_dir)
  hashes <- hash_paths(paths)
  
  identical(hashes,stored_hashes)
}

get_spec_path <- function(reproducibility_dir){
  paths <- load_paths(reproducibility_dir)
  paths$spec_path
}
