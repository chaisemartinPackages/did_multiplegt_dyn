#' Internal function of did_multiplegt_dyn
#' @param data data
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @noRd
did_multiplegt_dyn_graph <- function(data) {
  grmat <- data$mat_res_XX
  l_XX <- data$l_XX 
  grmat[l_XX + 1, c(1,3,4)] <- 0
  grmat <- data.frame(grmat[, c(1,3,4,7)])
  did_multiplegt_dyn_plot <- ggplot(grmat,aes(x = .data$Time, y = .data$Estimate, group = 1)) + 
    geom_line(colour = "blue") + 
    geom_errorbar(data = ~dplyr::filter(.x, grmat$Estimate != 0),aes(ymin = .data$LB.CI, ymax = .data$UB.CI), position=position_dodge(0.05), width = 0.2, colour = "red")  + 
    geom_point(colour = "blue") + 
    ggtitle("DID, from last period before treatment changes (t=0) to t") + 
    xlab("Relative time to last period before treatment changes (t=0)") +
    theme(plot.title = element_text(hjust = 0.5))
  return(did_multiplegt_dyn_plot)
}

#' Internal function of did_multiplegt_dyn
#' @param obj A did_multiplegt_dyn object
#' @import ggplot2
#' @noRd
combine_plot <- function(obj) {
  color_set <- get_colors(length(obj$by_levels))
  base_plot <- obj$by_level_1$plot
  plot <-  ggplot(data = base_plot$data, aes(x = base_plot$data$Time, y = base_plot$data$Estimate))  +
            geom_point(colour = color_set[1]) + geom_line(aes(colour = paste0(obj$args$by," = ", obj$by_levels[1]))) +
            geom_errorbar(data = base_plot$data, aes(ymin = base_plot$data$LB.CI, ymax = base_plot$data$UB.CI), position=position_dodge(0.05), width = 0.2, colour = color_set[1])
  if (length(obj$by_levels) > 1) {
    for (j in 2:length(obj$by_levels)) {
      add_plot <-  obj[[paste0("by_level_",j)]][["plot"]]
      plot <- plot + 
          geom_point(data = add_plot$data, aes_(x =  add_plot$data$Time, y = add_plot$data$Estimate), colour = color_set[j]) +
          geom_errorbar(data = add_plot$data, aes_(ymin = add_plot$data$LB.CI, ymax = add_plot$data$UB.CI), position=position_dodge(0.05), width = 0.2, colour =  color_set[j]) +
          geom_line(data =  add_plot$data, aes_(x =  add_plot$data$Time, y =  add_plot$data$Estimate, 
              colour = paste0(obj$args$by, " = ", obj$by_levels[j])))
    }
  }
  plot <- plot + ylab("Estimate") + ggtitle("DID, from last period before treatment changes (t=0) to t") + 
    xlab("Relative time to last period before treatment changes (t=0)") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
    scale_colour_manual("", breaks = paste0(obj$args$by, " = ", obj$by_levels), values = color_set)
  return(plot)
}

#' Internal function to retrieve plot colors
#' @param N Number of colors to retrieve
#' @import ggplot2
#' @importFrom grDevices colors
#' @noRd
get_colors <- function(N) {
  must_color <- c(552, 26, 81, 68, 450, 640, 24, 498) 
  # Indices of the following colors in ggplot's colors()
  # Red, blue, green, cyan, magenta, violet, black, orange
  other_color <- 1:657
  other_color <- other_color[!(other_color %in% must_color)]
  if (N > length(must_color)) {
    index_colors <- c(must_color, sample(other_color, N - length(must_color)))
  } else {
    index_colors <- must_color[1:N]
  }
  return(colors()[index_colors])
}