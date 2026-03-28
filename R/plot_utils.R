library(ggplot2)
library(patchwork)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

make_plot_group_id <- function(df, group_cols = NULL) {
  if (is.null(group_cols) || !length(group_cols)) {
    return(rep("all", nrow(df)))
  }

  valid_cols <- group_cols[group_cols %in% names(df)]
  if (!length(valid_cols)) {
    return(rep("all", nrow(df)))
  }

  interaction(df[, valid_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
}

make_facet_label <- function(values, facet_by) {
  if (is.null(facet_by) || !length(facet_by)) {
    return(rep("'all'", length(values)))
  }

  if (facet_by == "G0") {
    vals <- as.character(values)
    levs <- as.character(sort(unique(values)))
    return(factor(
      paste0("G[0]==", vals),
      levels = paste0("G[0]==", levs)
    ))
  }

  if (facet_by == "ploidy") {
    vals <- as.character(values)
    levs <- as.character(sort(unique(values)))
    return(factor(
      paste0("ploidy==", vals),
      levels = paste0("ploidy==", levs)
    ))
  }

  vals <- gsub("'", "\\\\'", as.character(values), fixed = TRUE)
  levs <- gsub("'", "\\\\'", as.character(sort(unique(values))), fixed = TRUE)
  factor(
    paste0("'", vals, "'"),
    levels = paste0("'", levs, "'")
  )
}

get_facet_values <- function(df, facet_by) {
  if (is.null(facet_by) || !length(facet_by) || !(facet_by %in% names(df))) {
    return(rep("all", nrow(df)))
  }
  df[[facet_by]]
}

get_variable_title <- function(variable_name) {
  if (identical(variable_name, "N")) {
    return("N (alive cells)")
  }
  if (identical(variable_name, "R1")) {
    return("R1: Glucose (mM)")
  }
  variable_name
}

plot_fit_overlays <- function(
  sim_df,
  obs_df = NULL,
  color_by = NULL,
  linetype_by = NULL,
  ribbon_fill_by = NULL,
  color_values = NULL,
  linetype_values = NULL,
  ribbon_fill_values = NULL,
  color_label = NULL,
  linetype_label = NULL,
  ribbon_fill_label = NULL,
  facet_rows = "G0",
  facet_cols = "ploidy",
  title = NULL,
  subtitle = NULL,
  scales = "free",
  collect_guides = TRUE,
  legend_position = "right",
  legend_from_first_panel_only = TRUE,
  clip_sim_to_obs_max = TRUE,
  obs_max_multiplier = 1.5,
  line_width = 0.9,
  point_size = 1.1,
  line_alpha = 0.85,
  point_alpha = 0.75,
  ribbon_alpha = 0.18,
  print_plot = TRUE
) {
  stopifnot(all(c("time", "variable", "value") %in% names(sim_df)))

  sim_df <- sim_df
  color_levels <- NULL
  linetype_levels <- NULL
  ribbon_fill_levels <- NULL

  if (!is.null(color_by) && color_by %in% names(sim_df)) {
    color_levels <- unique(as.character(sim_df[[color_by]]))
    sim_df[[color_by]] <- factor(sim_df[[color_by]], levels = color_levels)
    if (!is.null(obs_df) && nrow(obs_df) && color_by %in% names(obs_df)) {
      obs_df[[color_by]] <- factor(obs_df[[color_by]], levels = color_levels)
    }
  }
  if (!is.null(linetype_by) && linetype_by %in% names(sim_df)) {
    linetype_levels <- unique(as.character(sim_df[[linetype_by]]))
    sim_df[[linetype_by]] <- factor(sim_df[[linetype_by]], levels = linetype_levels)
    if (!is.null(obs_df) && nrow(obs_df) && linetype_by %in% names(obs_df)) {
      obs_df[[linetype_by]] <- factor(obs_df[[linetype_by]], levels = linetype_levels)
    }
  }
  if (!is.null(ribbon_fill_by) && ribbon_fill_by %in% names(sim_df)) {
    ribbon_fill_levels <- unique(as.character(sim_df[[ribbon_fill_by]]))
    sim_df[[ribbon_fill_by]] <- factor(sim_df[[ribbon_fill_by]], levels = ribbon_fill_levels)
  }

  variable_y_max <- NULL
  if (isTRUE(clip_sim_to_obs_max) && !is.null(obs_df) && nrow(obs_df)) {
    variable_y_max <- tapply(obs_df$value, obs_df$variable, function(x) {
      if (!length(x) || all(is.na(x))) {
        return(NA_real_)
      }
      max(x, na.rm = TRUE) * obs_max_multiplier
    })

    for (var_name in names(variable_y_max)) {
      ymax <- variable_y_max[[var_name]]
      if (is.na(ymax)) {
        next
      }
      idx <- sim_df$variable == var_name & is.finite(sim_df$value) & sim_df$value > ymax
      sim_df$value[idx] <- NA_real_
      if ("value_lower" %in% names(sim_df)) {
        sim_df$value_lower[idx] <- NA_real_
      }
      if ("value_upper" %in% names(sim_df)) {
        sim_df$value_upper[idx] <- NA_real_
      }
    }
  }

  sim_df$.group_id <- make_plot_group_id(
    sim_df,
    unique(c("well_idx", "line_id", "ploidy", "G0", color_by, linetype_by, ribbon_fill_by, "group_1", "group_2"))
  )
  sim_df$.facet_row <- make_facet_label(get_facet_values(sim_df, facet_rows), facet_rows)
  sim_df$.facet_col <- make_facet_label(get_facet_values(sim_df, facet_cols), facet_cols)

  if (!is.null(obs_df) && nrow(obs_df)) {
    stopifnot(all(c("time", "variable", "value") %in% names(obs_df)))
    obs_df <- obs_df
    obs_df$.facet_row <- make_facet_label(get_facet_values(obs_df, facet_rows), facet_rows)
    obs_df$.facet_col <- make_facet_label(get_facet_values(obs_df, facet_cols), facet_cols)
  }

  has_interval <- all(c("value_lower", "value_upper") %in% names(sim_df))
  if (isTRUE(has_interval)) {
    has_interval <- any(!is.na(sim_df$value_lower) & !is.na(sim_df$value_upper))
  }

  vars_to_plot <- unique(sim_df$variable)

  plot_list <- lapply(seq_along(vars_to_plot), function(i) {
    var_name <- vars_to_plot[[i]]
    sub_sim <- sim_df[sim_df$variable == var_name, , drop = FALSE]
    sub_obs <- NULL
    if (!is.null(obs_df) && nrow(obs_df)) {
      sub_obs <- obs_df[obs_df$variable == var_name, , drop = FALSE]
    }
    show_panel_legend <- !isTRUE(collect_guides) || !isTRUE(legend_from_first_panel_only) || i == 1L

    p <- ggplot()

    if (isTRUE(has_interval)) {
      ribbon_map <- aes(
        x = .data$time,
        ymin = .data$value_lower,
        ymax = .data$value_upper,
        group = .data$.group_id
      )
      if (!is.null(ribbon_fill_by) && ribbon_fill_by %in% names(sub_sim)) {
        ribbon_map$fill <- as.name(ribbon_fill_by)
      }

      p <- p + geom_ribbon(
        data = sub_sim,
        mapping = ribbon_map,
        alpha = ribbon_alpha,
        show.legend = show_panel_legend,
        inherit.aes = FALSE
      )
    }

    line_map <- aes(
      x = .data$time,
      y = .data$value,
      group = .data$.group_id
    )
    if (!is.null(color_by) && color_by %in% names(sub_sim)) {
      line_map$colour <- as.name(color_by)
    }
    if (!is.null(linetype_by) && linetype_by %in% names(sub_sim)) {
      line_map$linetype <- as.name(linetype_by)
    }

    p <- p + geom_line(
      data = sub_sim,
      mapping = line_map,
      linewidth = line_width,
      alpha = line_alpha,
      show.legend = show_panel_legend
    )

    if (!is.null(sub_obs) && nrow(sub_obs)) {
      p <- p + geom_point(
        data = sub_obs,
        mapping = aes(x = .data$time, y = .data$value),
        color = "black",
        size = point_size,
        alpha = point_alpha
      )
    }

    p <- p +
      facet_grid(
        rows = vars(.facet_row),
        cols = vars(.facet_col),
        scales = scales,
        labeller = ggplot2::label_parsed
      ) +
      theme_minimal() +
      labs(
        title = get_variable_title(var_name),
        x = "time (h)",
        y = "",
        color = color_label %||% color_by,
        linetype = linetype_label %||% linetype_by,
        fill = ribbon_fill_label %||% ribbon_fill_by
      ) +
      theme(
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "grey80", fill = NA)
      )

    if (!is.null(color_values)) {
      p <- p + scale_color_manual(values = color_values, limits = color_levels, drop = FALSE)
    } else if (!is.null(color_levels)) {
      p <- p + scale_color_discrete(limits = color_levels, drop = FALSE)
    }
    if (!is.null(linetype_values)) {
      p <- p + scale_linetype_manual(values = linetype_values, limits = linetype_levels, drop = FALSE)
    } else if (!is.null(linetype_levels)) {
      p <- p + scale_linetype_discrete(limits = linetype_levels, drop = FALSE)
    }
    if (!is.null(ribbon_fill_values)) {
      p <- p + scale_fill_manual(values = ribbon_fill_values, limits = ribbon_fill_levels, drop = FALSE)
    } else if (!is.null(ribbon_fill_levels)) {
      p <- p + scale_fill_discrete(limits = ribbon_fill_levels, drop = FALSE)
    }

    p
  })

  combined_plot <- wrap_plots(plot_list, nrow = 1) +
    plot_layout(guides = if (isTRUE(collect_guides)) "collect" else "keep") +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  if (isTRUE(collect_guides)) {
    combined_plot <- combined_plot & theme(legend.position = legend_position)
  }

  if (isTRUE(print_plot)) {
    print(combined_plot)
  }

  invisible(combined_plot)
}
