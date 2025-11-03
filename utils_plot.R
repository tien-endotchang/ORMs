# Return the mean squared error
mse = function(e) return ( sum(e^2) / length(e) )

# xxx
test_relation <- function(x, y) {
  # Perform the Mann-Whitney test (wilcox.test in R)
  result <- wilcox.test(x, y, alternative = "two.sided", paired = TRUE)
  
  # Interpret the p-value and test statistic to classify the relationship
  if (result$p.value > 0.1 || is.na(result$p.value)) {
    return("tie")  # No significant difference
  } else {
    # Re-run the test with directional hypotheses
    greater_test <- wilcox.test(x, y, alternative = "greater", paired = TRUE)
    smaller_test <- wilcox.test(x, y, alternative = "less", paired = TRUE)
    
    if (greater_test$p.value <= 0.05) { # if x is greater than y
      return("greater")
    } else if (smaller_test$p.value <= 0.05) { # if x is smaller than y
      return("smaller")
    }
  }
}

# Plot the metrics (RMSE/Ktau) of the OTMs under Johnson's Z
plot_Metric_Z.J_All.reA = function(result, y.axis, log.scale = FALSE){
  y_pos = if(log.scale) log10(max(result$value)) else max(result$value)
  
  p = ggplot(data = result, aes(x = var, y = value)) + # x = variable
    geom_boxplot_pattern(
      aes(pattern = Method,
          pattern_angle = Method,
          pattern_density = Method),
      pattern_spacing = 0.015,
      pattern_colour  = 'black',
      alpha = 0.75,
      outlier.size = 1,
      outlier.alpha = 0.3, 
      position = position_dodge(width = 0.85)) +

    scale_pattern_manual(name="Method", values = c("none", "circle", "stripe", "crosshatch")) +
    scale_pattern_angle_manual(name="Method", values = c(0, 30, 30, -30)) +
    scale_pattern_density_manual(name="Method", values = c(0, 0.2, 0.02, 0.02)) +

    labs(x = "the number of predictor variables", y = y.axis) +
    theme(text = element_text(size = 15)) +
    
    stat_summary(
      fun = median, geom = "smooth", 
      aes(group = Method, colour = Method), 
      linewidth = 1.5, linetype = "solid", 
      position = position_dodge(width = 0.85),
      show.legend = FALSE) + 
    # 00000050
    scale_colour_manual(values = c("GDA" = "#00000050", "CorPA" = "#FF000050", 
                                   "RegPA" = "#00009950", "IdA" = "#FF000000")) +
    # "#F8766D" "#00BA38" "#619CFF"
    guides(pattern = guide_legend(override.aes = list(fill = "white", pattern_colour = "black"))) +
    theme_minimal(base_size = 18) +
    apa_theme
  
  weird <- scales::trans_new("weird_log",
                             transform = function(x) ifelse(x !=  log10(max(result$value)), log10(x), x),
                             inverse = function(x) ifelse(x != log10(max(result$value)), 10^x, x))
  
  if(log.scale){
    p = p + scale_y_continuous(breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x)),
                               transform = weird)
  }
  return( p )
}

# Plot the bias of the of all reallocation methods under Johnson's Z
plot_bias_Z.J_All.reA_wrap = function(result, y.axis){
  
  limits = c(-0.15, 0.15)
  p = ggplot(data = result, aes(x = variable, y = value)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#000000", linewidth = 1) +
    geom_boxplot_pattern(
      aes(pattern = Method,
          pattern_angle = Method,
          pattern_density = Method),
      pattern_colour  = 'black', 
      alpha = 0.75,
      outlier.size = 1,
      outlier.alpha = 0.3,
      position = position_dodge(width = 0.85)) +
    
    scale_pattern_manual(name="Method", values = c("none", "circle", "stripe", "crosshatch")) +
    scale_pattern_angle_manual(name="Method", values = c(0, 30, 30, -30)) + 
    scale_pattern_density_manual(name="Method", values = c(0, 0.2, 0.02, 0.02)) +
    scale_pattern_spacing_manual(name = "Method", values = c(0.025, 0.015, 0.02, 0.025)) +

    facet_wrap(~ num_p, scales = "free", nrow = 3, labeller = labeller(num_p = 
                                                                         c("4" = "p = 4",
                                                                           "7" = "p = 7",
                                                                           "10" = "p = 10"))) +
    
    labs(x = "the i-th important predictor variable", y = y.axis, fill = "OT") +
    theme(text = element_text(size = 15)) +

    stat_summary(
      fun = median, geom = "smooth",
      aes(group = Method, colour = Method),
      linewidth = 1.5, linetype = "solid", 
      position = position_dodge(width = 0.85),
      show.legend = FALSE) + 

    scale_colour_manual(values = c("GDA" = "#00000050", "CorPA" = "#FF000050", 
                                   "RegPA" = "#00009950", "IdA" = "#FF000000")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white", pattern_colour = "black"))) +
    theme_minimal(base_size = 18) +
    apa_theme
  
  return( p )
}

# Set APA theme
apa_theme <- theme(
  plot.margin = unit(c(1, 1, 1, 1), "cm"),
  plot.background = element_rect(fill = "white", color = NA),
  plot.title = element_text(size = 22, face = "bold",
                            hjust = 0.5,
                            margin = margin(b = 15)),
  axis.line = element_line(color = "black", linewidth = .5),
  axis.title = element_text(size = 15, color = "black",
                            face = "bold"),
  axis.text = element_text(size = 15, color = "black"),
  axis.text.x = element_text(margin = margin(t = 10)),
  axis.title.y = element_text(margin = margin(r = 10)),
  axis.ticks = element_line(size = .5),
  panel.grid = element_blank(),
  # legend.position = c(0.20, 0.8),
  legend.background = element_rect(color = "black"),
  legend.text = element_text(size = 12), ### 15
  legend.margin = margin(t = 5, l = 5, r = 5, b = 5),
  legend.key = element_rect(color = NA, fill = NA),
  legend.key.size = unit(1, 'cm'), 
  strip.text = element_text(size = 15, color = "black",
                            face = "bold")
)


apply_filters <- function(data_list, condition) {
  lapply(data_list, function(data) {
    data %>% filter(!!condition)
  })
}