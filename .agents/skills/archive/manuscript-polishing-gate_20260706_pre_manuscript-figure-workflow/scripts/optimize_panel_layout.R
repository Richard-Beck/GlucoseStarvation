#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  invisible(utils::globalVariables(c("figure", "panel")))
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  optimize_panel_layout.R --input subpanel_dimensions.csv --output-dir layout_dir [options]",
      "",
      "Required input columns:",
      "  figure, panel, and either width_in/height_in or width_px/height_px.",
      "  Optional columns: subpanel_png, order, dpi.",
      "",
      "Options:",
      "  --target-width <in>   Maximum/recommended figure width in inches (default: 7)",
      "  --max-height <in>     Maximum figure height in inches (default: 9.25)",
      "  --gap <in>            Gap between panel slots in inches (default: 0.08)",
      "  --sx-lo <num>         Lower x-scale bound per panel (default: 0.9)",
      "  --sx-hi <num>         Upper x-scale bound per panel (default: 1.05)",
      "  --sy-lo <num>         Lower y-scale bound per panel (default: 0.9)",
      "  --sy-hi <num>         Upper y-scale bound per panel (default: 1.05)",
      "  --top-k <int>         Number of unscaled trees to optimize (default: 25)",
      "  --max-panels <int>    Max panels per figure for exhaustive tree search (default: 8)",
      "  --dpi <num>           DPI used when only pixel dimensions are supplied (default: 300)",
      "  --layout-out <path>   Layout-plan CSV path (default: <output-dir>/layout_plan.csv)",
      "  --report-out <path>   Markdown report path (default: <output-dir>/layout_report.md)",
      "  --preview <true|false> Write simple layout preview PNGs (default: true)",
      "",
      "Coordinate convention:",
      "  x_in/y_in and x_npc/y_npc are lower-left panel origins.",
      "  Use y_npc directly in assembly; do not convert with 1 - y_npc.",
      "",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  out <- list(
    input = NULL,
    output_dir = NULL,
    target_width = 7,
    max_height = 9.25,
    gap = 0.08,
    sx_lo = 0.90,
    sx_hi = 1.05,
    sy_lo = 0.90,
    sy_hi = 1.05,
    top_k = 25,
    max_panels = 8,
    dpi = 300,
    layout_out = NULL,
    report_out = NULL,
    preview = TRUE,
    width_overflow = 1000,
    height_overflow = 1000,
    space_weight = 1,
    scale_weight = 0.02,
    width_fill_weight = 0.05
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- gsub("-", "_", sub("^--", "", key))
    if (name %in% c("help", "h")) {
      usage()
      quit(status = 0)
    }
    value <- TRUE
    if (i + 1 <= length(args) && !startsWith(args[[i + 1]], "--")) {
      value <- args[[i + 1]]
      i <- i + 1
    }
    out[[name]] <- value
    i <- i + 1
  }

  numeric_fields <- c(
    "target_width", "max_height", "gap", "sx_lo", "sx_hi", "sy_lo", "sy_hi",
    "top_k", "max_panels", "dpi", "width_overflow", "height_overflow",
    "space_weight", "scale_weight", "width_fill_weight"
  )
  for (field in numeric_fields) {
    out[[field]] <- as.numeric(out[[field]])
    if (is.na(out[[field]])) stop("Invalid numeric option --", gsub("_", "-", field), call. = FALSE)
  }
  out$top_k <- as.integer(out$top_k)
  out$max_panels <- as.integer(out$max_panels)
  out$preview <- tolower(as.character(out$preview)) %in% c("1", "true", "yes", "y")

  if (is.null(out$input)) stop("Missing required --input", call. = FALSE)
  if (is.null(out$output_dir)) stop("Missing required --output-dir", call. = FALSE)
  if (out$sx_lo > out$sx_hi || out$sy_lo > out$sy_hi) {
    stop("Scale lower bounds must not exceed upper bounds", call. = FALSE)
  }
  out
}

safe_name <- function(x) {
  gsub("_+$", "", gsub("[^A-Za-z0-9_-]+", "_", x))
}

read_panel_dimensions <- function(path, fallback_dpi) {
  panels <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(panels) <- gsub("-", "_", names(panels))
  if (!"figure" %in% names(panels)) panels$figure <- "Figure"
  if (!"panel" %in% names(panels)) panels$panel <- letters[seq_len(nrow(panels))]
  if (!"subpanel_png" %in% names(panels)) {
    for (candidate in c("subpanel_image", "path", "png", "file")) {
      if (candidate %in% names(panels)) {
        panels$subpanel_png <- panels[[candidate]]
        break
      }
    }
  }
  if (!"subpanel_png" %in% names(panels)) panels$subpanel_png <- ""

  has_inches <- all(c("width_in", "height_in") %in% names(panels))
  has_pixels <- all(c("width_px", "height_px") %in% names(panels))
  if (!has_inches && !has_pixels) {
    stop("Input must contain width_in/height_in or width_px/height_px columns", call. = FALSE)
  }

  if (has_inches) {
    panels$width_in <- as.numeric(panels$width_in)
    panels$height_in <- as.numeric(panels$height_in)
  } else {
    dpi <- if ("dpi" %in% names(panels)) as.numeric(panels$dpi) else rep(fallback_dpi, nrow(panels))
    dpi[is.na(dpi) | dpi <= 0] <- fallback_dpi
    panels$width_px <- as.numeric(panels$width_px)
    panels$height_px <- as.numeric(panels$height_px)
    panels$width_in <- panels$width_px / dpi
    panels$height_in <- panels$height_px / dpi
  }

  if (any(is.na(panels$width_in) | is.na(panels$height_in) | panels$width_in <= 0 | panels$height_in <= 0)) {
    stop("Panel dimensions must be positive numeric values", call. = FALSE)
  }
  panels$..input_order <- seq_len(nrow(panels))
  panels
}

mk <- function(dir, left, right) {
  list(dir = dir, left = left, right = right)
}

make_layout_enumerator <- function(panel_ids) {
  cache <- new.env(parent = emptyenv())

  layouts <- function(i, j) {
    key <- paste(i, j, sep = ":")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    if (i == j) {
      out <- list(i)
    } else {
      out <- list()
      for (k in i:(j - 1)) {
        for (left in layouts(i, k)) {
          for (right in layouts(k + 1, j)) {
            out <- c(out, list(mk("beside", left, right), mk("stack", left, right)))
          }
        }
      }
    }
    assign(key, out, envir = cache)
    out
  }

  tree_to_string <- function(x) {
    if (is.numeric(x)) return(panel_ids[[x]])
    if (x$dir == "beside") {
      paste0("[", tree_to_string(x$left), " | ", tree_to_string(x$right), "]")
    } else {
      paste0("[", tree_to_string(x$left), " / ", tree_to_string(x$right), "]")
    }
  }

  list(layouts = layouts, tree_to_string = tree_to_string)
}

unpack_scales <- function(p, n) {
  list(sx = p[seq_len(n)], sy = p[n + seq_len(n)])
}

dim_layout <- function(x, panels, gap, sx = rep(1, nrow(panels)), sy = rep(1, nrow(panels))) {
  if (is.numeric(x)) {
    width <- sx[[x]] * panels$width_in[[x]]
    height <- sy[[x]] * panels$height_in[[x]]
    return(list(W = width, H = height, A = width * height))
  }

  left <- dim_layout(x$left, panels, gap, sx, sy)
  right <- dim_layout(x$right, panels, gap, sx, sy)
  if (x$dir == "beside") {
    list(W = left$W + gap + right$W, H = max(left$H, right$H), A = left$A + right$A)
  } else {
    list(W = max(left$W, right$W), H = left$H + gap + right$H, A = left$A + right$A)
  }
}

score_layout <- function(x, panels, args, sx = rep(1, nrow(panels)), sy = rep(1, nrow(panels)),
                         penalize_scaling = FALSE) {
  d <- dim_layout(x, panels, args$gap, sx, sy)
  figure_area <- d$W * d$H
  wasted_frac <- if (figure_area > 0) max(0, figure_area - d$A) / figure_area else Inf
  score <- args$width_overflow * max(0, d$W / args$target_width - 1)^2 +
    args$height_overflow * max(0, d$H / args$max_height - 1)^2 +
    args$space_weight * wasted_frac +
    args$width_fill_weight * max(0, 1 - d$W / args$target_width)^2
  if (penalize_scaling) {
    score <- score + args$scale_weight * sum((sx - 1)^2 + (sy - 1)^2)
  }
  score
}

coords_layout <- function(x, panels, gap, sx = rep(1, nrow(panels)), sy = rep(1, nrow(panels)),
                          x0 = 0, y0 = 0, W = NULL, H = NULL) {
  d <- dim_layout(x, panels, gap, sx, sy)
  if (is.null(W)) W <- d$W
  if (is.null(H)) H <- d$H

  if (is.numeric(x)) {
    width <- sx[[x]] * panels$width_in[[x]]
    height <- sy[[x]] * panels$height_in[[x]]
    return(data.frame(
      panel = panels$panel[[x]],
      subpanel_png = panels$subpanel_png[[x]],
      x_in = x0,
      y_in = y0,
      width_in = width,
      height_in = height,
      sx = sx[[x]],
      sy = sy[[x]],
      stringsAsFactors = FALSE
    ))
  }

  left <- dim_layout(x$left, panels, gap, sx, sy)
  right <- dim_layout(x$right, panels, gap, sx, sy)
  if (x$dir == "beside") {
    rbind(
      coords_layout(x$left, panels, gap, sx, sy, x0, y0 + (H - left$H) / 2),
      coords_layout(x$right, panels, gap, sx, sy, x0 + left$W + gap, y0 + (H - right$H) / 2)
    )
  } else {
    rbind(
      coords_layout(x$left, panels, gap, sx, sy, x0 + (W - left$W) / 2, y0 + right$H + gap),
      coords_layout(x$right, panels, gap, sx, sy, x0 + (W - right$W) / 2, y0)
    )
  }
}

plot_preview <- function(path, figure_id, layout_df) {
  layout_width <- unique(layout_df$layout_width_in)
  layout_height <- unique(layout_df$layout_height_in)
  grDevices::png(path, width = 1800, height = 1400, res = 200, type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(1, 1, 3, 1), xaxs = "i", yaxs = "i")
  graphics::plot(
    NA,
    xlim = c(0, layout_width),
    ylim = c(0, layout_height),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = sprintf("%s layout: %.2f x %.2f in", figure_id, layout_width, layout_height)
  )
  graphics::rect(0, 0, layout_width, layout_height, border = "gray70", lty = 2)
  for (i in seq_len(nrow(layout_df))) {
    row <- layout_df[i, ]
    graphics::rect(
      row$x_in, row$y_in,
      row$x_in + row$width_in, row$y_in + row$height_in,
      col = "gray92", border = "black", lwd = 2
    )
    graphics::text(
      row$x_in + row$width_in / 2,
      row$y_in + row$height_in / 2,
      sprintf("%s\nsx %.2f\nsy %.2f", row$panel, row$sx, row$sy),
      cex = 0.8
    )
  }
}

optimize_figure <- function(figure_id, panels, args) {
  if ("order" %in% names(panels)) {
    panels <- panels[order(as.numeric(panels$order), panels$..input_order), , drop = FALSE]
  } else {
    panels <- panels[order(panels$..input_order), , drop = FALSE]
  }
  rownames(panels) <- NULL
  n <- nrow(panels)
  if (n > args$max_panels) {
    stop(
      sprintf(
        "Figure '%s' has %d panels; max-panels is %d. Split the figure or raise --max-panels knowingly.",
        figure_id, n, args$max_panels
      ),
      call. = FALSE
    )
  }

  enumerator <- make_layout_enumerator(panels$panel)
  all_trees <- enumerator$layouts(1, n)
  unscaled_scores <- vapply(all_trees, score_layout, numeric(1), panels = panels, args = args)
  top_idx <- order(unscaled_scores)[seq_len(min(args$top_k, length(unscaled_scores)))]

  best <- NULL
  lower <- c(rep(args$sx_lo, n), rep(args$sy_lo, n))
  upper <- c(rep(args$sx_hi, n), rep(args$sy_hi, n))
  p0 <- rep(1, 2 * n)
  for (idx in top_idx) {
    tree <- all_trees[[idx]]
    opt <- try(
      stats::optim(
        par = p0,
        fn = function(p) {
          scales <- unpack_scales(p, n)
          score_layout(tree, panels, args, scales$sx, scales$sy, penalize_scaling = TRUE)
        },
        method = "L-BFGS-B",
        lower = lower,
        upper = upper
      ),
      silent = TRUE
    )
    if (inherits(opt, "try-error")) next
    if (is.null(best) || opt$value < best$score) {
      scales <- unpack_scales(opt$par, n)
      best <- list(tree = tree, score = opt$value, sx = scales$sx, sy = scales$sy)
    }
  }

  if (is.null(best)) stop("No valid optimized layout found for ", figure_id, call. = FALSE)
  d <- dim_layout(best$tree, panels, args$gap, best$sx, best$sy)
  coords <- coords_layout(best$tree, panels, args$gap, best$sx, best$sy)
  figure_area <- d$W * d$H
  panel_area <- sum(coords$width_in * coords$height_in)
  wasted_frac <- if (figure_area > 0) max(0, figure_area - panel_area) / figure_area else NA_real_

  out <- coords
  out$figure <- figure_id
  out$input_width_in <- panels$width_in[match(out$panel, panels$panel)]
  out$input_height_in <- panels$height_in[match(out$panel, panels$panel)]
  out$x_npc <- out$x_in / d$W
  out$y_npc <- out$y_in / d$H
  out$width_npc <- out$width_in / d$W
  out$height_npc <- out$height_in / d$H
  out$layout_width_in <- d$W
  out$layout_height_in <- d$H
  out$target_width_in <- args$target_width
  out$max_height_in <- args$max_height
  out$layout_tree <- enumerator$tree_to_string(best$tree)
  out$layout_score <- best$score
  out$wasted_fraction <- wasted_frac
  out$over_width_in <- max(0, d$W - args$target_width)
  out$over_height_in <- max(0, d$H - args$max_height)
  out$scale_note <- ifelse(
    abs(out$sx - 1) > 0.01 | abs(out$sy - 1) > 0.01,
    "revise subpanel dimensions or accept optimizer scale",
    "near native size"
  )
  out[, c(
    "figure", "panel", "subpanel_png", "x_in", "y_in", "width_in", "height_in",
    "sx", "sy", "x_npc", "y_npc", "width_npc", "height_npc",
    "input_width_in", "input_height_in", "layout_width_in", "layout_height_in",
    "target_width_in", "max_height_in", "layout_tree", "layout_score",
    "wasted_fraction", "over_width_in", "over_height_in", "scale_note"
  )]
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(args$layout_out)) args$layout_out <- file.path(args$output_dir, "layout_plan.csv")
  if (is.null(args$report_out)) args$report_out <- file.path(args$output_dir, "layout_report.md")

  panels <- read_panel_dimensions(args$input, args$dpi)
  figures <- unique(panels$figure)
  layout_list <- vector("list", length(figures))
  names(layout_list) <- figures

  for (figure_id in figures) {
    layout_list[[figure_id]] <- optimize_figure(
      figure_id,
      panels[panels$figure == figure_id, , drop = FALSE],
      args
    )
  }

  layout_plan <- do.call(rbind, layout_list)
  rownames(layout_plan) <- NULL
  utils::write.csv(layout_plan, args$layout_out, row.names = FALSE)

  if (args$preview) {
    for (figure_id in figures) {
      preview_path <- file.path(args$output_dir, paste0(safe_name(figure_id), "_layout_preview.png"))
      plot_preview(preview_path, figure_id, layout_plan[layout_plan$figure == figure_id, , drop = FALSE])
    }
  }

  report <- c(
    "# Panel layout optimization report",
    "",
    sprintf("- Input dimensions: `%s`", args$input),
    sprintf("- Layout plan: `%s`", args$layout_out),
    sprintf("- Target width: %.2f in", args$target_width),
    sprintf("- Maximum height: %.2f in", args$max_height),
    sprintf("- Gap: %.2f in", args$gap),
    "- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.",
    "",
    "## Figures",
    ""
  )
  for (figure_id in figures) {
    one <- layout_plan[layout_plan$figure == figure_id, , drop = FALSE]
    report <- c(
      report,
      sprintf(
        "- `%s`: %.2f x %.2f in, wasted %.1f%%, tree `%s`",
        figure_id,
        unique(one$layout_width_in),
        unique(one$layout_height_in),
        100 * unique(one$wasted_fraction),
        unique(one$layout_tree)
      )
    )
  }
  report <- c(
    report,
    "",
    "## Scale recommendations",
    "",
    "Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated and this optimizer rerun before final assembly."
  )
  writeLines(report, args$report_out)

  cat("Wrote layout plan: ", args$layout_out, "\n", sep = "")
  cat("Wrote layout report: ", args$report_out, "\n", sep = "")
  for (figure_id in figures) {
    one <- layout_plan[layout_plan$figure == figure_id, , drop = FALSE]
    cat(
      sprintf(
        "%s: %.2f x %.2f in, %d panels, tree %s\n",
        figure_id,
        unique(one$layout_width_in),
        unique(one$layout_height_in),
        nrow(one),
        unique(one$layout_tree)
      )
    )
  }
}

tryCatch(
  main(),
  error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
    quit(status = 1)
  }
)
