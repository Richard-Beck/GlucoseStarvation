plot_slider_export <- function(plots, filepath) {
  suppressWarnings({
    library(ggplot2)
    library(plotly)
    library(htmlwidgets)
    
    # Convert ggplot list to plotly
    plotly_list <- lapply(plots, ggplotly)
    n <- length(plotly_list)
    
    # Set visibility and frame names
    for (i in seq_len(n)) {
      plotly_list[[i]]$x$layout$visible <- if (i == 1) TRUE else FALSE
      plotly_list[[i]]$x$layout$name <- paste("Frame", i)
    }
    
    # Use first plot as base
    base_plot <- plotly_list[[1]]
    
    # Add frames with layout for autoscaling axes
    base_plot$x$frames <- lapply(seq_len(n), function(i) {
      list(
        name = paste("Frame", i),
        data = plotly_list[[i]]$x$data,
        layout = plotly_list[[i]]$x$layout
      )
    })
    
    # Add slider and animation controls
    base_plot <- base_plot %>%
      layout(
        updatemenus = list(
          list(
            type = "buttons",
            showactive = FALSE,
            buttons = list(
              list(
                label = "Play",
                method = "animate",
                args = list(
                  NULL,
                  list(
                    mode = "immediate",
                    frame = list(duration = 0, redraw = TRUE),
                    transition = list(duration = 0)
                  )
                )
              ),
              list(
                label = "Pause",
                method = "animate",
                args = list(NULL, list(frame = list(duration = 0), mode = "immediate"))
              )
            )
          )
        ),
        sliders = list(
          list(
            active = 0,
            steps = lapply(seq_len(n), function(i) {
              list(
                label = paste("Plot", i),
                method = "animate",
                args = list(
                  list(paste("Frame", i)),
                  list(mode = "immediate", frame = list(duration = 0, redraw = TRUE))
                )
              )
            })
          )
        )
      )
    
    # Save as standalone HTML
    saveWidget(base_plot, filepath, selfcontained = TRUE)
  })
}
