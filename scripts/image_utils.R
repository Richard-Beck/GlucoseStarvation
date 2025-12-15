## convenience functions for displaying images in R
library(data.table)
library(tiff)
library(EBImage)
make_lookup <- function(img_dir  = "all_raw/",
                        channels = c("phase","alive","aliveC1","aliveC2","dead")){
  fts <- list.files(img_dir, pattern="\\.tif$", full.names=TRUE)
  fn  <- basename(fts)
  noe <- sub("\\.tif$", "", fn)
  
  # regex for any of the channels
  chan_pat <- paste(channels, collapse="|")
  
  # remove "_<channel>_" to recover base; also extract channel
  lookup <- data.table(
    full_path = fts,
    base      = gsub(paste0("_(", chan_pat, ")_"), "_", noe),
    channel   = sub(paste0(".*_(", chan_pat, ")_.*"), "\\1", noe)
  )
}



# helper to read only the first valid frame
safe_read <- function(fp){
  imgs <- suppressWarnings(tiff::readTIFF(fp, all=TRUE))
  # drop any 1×1 junk frames
  real <- Filter(function(x) all(dim(x)>1), imgs)
  real[[1]]
}

make_composite <- function(base_name,lookup, displayIm =FALSE,
                           channels=c("phase","alive","aliveC1","aliveC2","dead")){
  dt0 <- lookup[ base==base_name ]
  if (nrow(dt0)==0) stop("No images found for ", base_name)
  
  # read + normalize each channel into a named list
  img_list <- setNames(
    lapply(channels, function(ch){
      fp <- dt0[channel==ch, full_path]
      if (length(fp)==0) return(NULL)
      x <- safe_read(fp[1])
      #x <- x - min(x); if (max(x)>0) x <- x/max(x)
      x
    }),
    channels
  )
  
  rescale_im <- function(x){
    x <- x - min(x); if (max(x)>0) x <- x/max(x)
  }
  
  # phase background
  phase <- rescale_im(img_list[["phase"]])
  if (is.null(phase)) stop("No phase image for ", base_name)
  
  # merge all alive signals by simply adding
  alive_max <- rescale_im(Reduce(function(a,b){
    if (is.null(b)) a else a+b
  }, img_list[c("alive","aliveC1","aliveC2")], init = 0))
  
  dead <- rescale_im(img_list[["dead"]])
  if (is.null(dead)) dead <- array(0, dim(phase))
  
  # build RGB:
  #   R = phase + dead
  #   G = phase + alive
  #   B = phase
  R <- pmin(phase + dead, 1)
  G <- pmin(phase + alive_max, 1)
  B <- phase
  
  rgb <- EBImage::rgbImage(red=R, green=G, blue=B)
  if(displayIm) EBImage::display(rgb, method="raster")
  invisible(rgb)
}

library(data.table)
library(magick)
# library(EBImage) # Required for make_composite

process_timecourse <- function(data, 
                               sel_cellLine, 
                               sel_ploidy, 
                               sel_glucose, 
                               img_lookup_table,
                               selection_idx = c(1, 1), # c(Plate_Index, Position_Index)
                               output_folder = NULL,
                               fps = 1) { # Frames Per Second for the GIF
  
  # --- 1. Helper: Circular Indexing ---
  get_wrapped_item <- function(items, idx) {
    if(length(items) == 0) stop("List is empty, cannot select.")
    real_idx <- (idx - 1) %% length(items) + 1
    return(items[real_idx])
  }
  
  # --- 2. Filter Metadata ---
  subset_dt <- data[cellLine == sel_cellLine & 
                      ploidy == sel_ploidy & 
                      glucose == sel_glucose]
  
  if(nrow(subset_dt) == 0) stop("No data found for this condition.")
  
  # Map Experiments to Plates
  plate_map <- unique(subset_dt[, .(experiment, plateID)])
  unique_plates <- plate_map$plateID
  
  # --- 3. Select Plate (Modulo) ---
  target_plate <- get_wrapped_item(unique_plates, selection_idx[1])
  target_exp   <- plate_map[plateID == target_plate, experiment][1]
  
  plate_msg <- sprintf("Plate: %s [%d/%d]", target_plate, which(unique_plates==target_plate), length(unique_plates))
  
  # --- 4. Find Files & Select Position (Modulo) ---
  well_pattern <- paste0("_", target_plate, "_")
  matches <- img_lookup_table[grepl(target_exp, base, fixed=TRUE) & 
                                grepl(well_pattern, base, fixed=TRUE)]
  
  if(nrow(matches) == 0) stop("Metadata exists, but no images found on disk.")
  
  # Regex to find position number
  pos_pattern <- paste0(target_plate, "_([0-9]+)_")
  all_bases <- unique(matches$base)
  positions <- unique(sub(paste0(".*", pos_pattern, ".*"), "\\1", all_bases))
  
  if(all(positions == all_bases)) {
    chosen_pos <- "All"
    final_bases <- sort(all_bases)
  } else {
    positions <- positions[order(as.numeric(positions))]
    chosen_pos <- get_wrapped_item(positions, selection_idx[2])
    strict_pattern <- paste0("_", target_plate, "_", chosen_pos, "_")
    final_bases <- grep(strict_pattern, all_bases, value=TRUE)
  }
  
  final_bases <- sort(final_bases)
  pos_msg <- sprintf("Pos: %s [%d/%d]", chosen_pos, which(positions==chosen_pos), length(positions))
  
  message(paste("Selected:", plate_msg, "|", pos_msg, "| Frames:", length(final_bases)))
  
  # --- 5. EXECUTION BRANCHING ---
  
  if (is.null(output_folder)) {
    # === VIEW MODE (Immediate) ===
    for(b in final_bases) {
      # Assuming make_composite displays the image if displayIm=TRUE
      make_composite(base_name = b, lookup = img_lookup_table, displayIm = TRUE)
      # Sleep removed as requested
    }
    
  } else {
    # === SAVE MODE (GIF) ===
    if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    
    # Construct Filename
    safe_line <- gsub("[^A-Za-z0-9-]", "", sel_cellLine)
    fname <- sprintf("%s_%s_Glc%s_%s_Pos%s.gif", 
                     safe_line, sel_ploidy, sel_glucose, target_plate, chosen_pos)
    full_path <- file.path(output_folder, fname)
    
    message("Compiling GIF frames...")
    
    # Initialize an empty magick stack
    img_stack <- magick::image_blank(width=1, height=1) 
    # (We will drop the first junk frame later, or use list accumulation)
    
    # Better approach: Accumulate in a list
    magick_list <- list()
    
    for(i in seq_along(final_bases)) {
      if(i %% 5 == 0) cat(".")
      
      # Get raw array data [0..1]
      raw_img <- make_composite(base_name = final_bases[i], lookup = img_lookup_table, displayIm = FALSE)
      
      # Convert to Magick image
      # image_read handles [0,1] arrays automatically
      magick_list[[i]] <- magick::image_read(as.array(raw_img))
    }
    cat("\n")
    
    # Combine list into a single magick stack
    stack <- magick::image_join(magick_list)
    
    # Animate
    message("Animating at ", fps, " fps...")
    anim <- magick::image_animate(stack, fps = fps)
    
    # Save
    message("Writing GIF to: ", full_path)
    magick::image_write(anim, path = full_path, format = "gif")
    message("Done.")
  }
}

