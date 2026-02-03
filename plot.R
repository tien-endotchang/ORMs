rm(list = ls())

library(ggplot2)
library(ggpattern)
library(reshape2)
library(patchwork)
library(dplyr)
library(scales)
library(tidyverse)
library(latex2exp)
source("utils_plot.R")

p_min = 3
p_max = 10
p_set = c(3:10)
num_ev_set = c(rep(1000, 4), rep(2500, 4))
num_rep = 10
num_y = 100

res_folder_name = "res"
fig_folder_name = "fig"
dir.create(fig_folder_name)
scenario_set = c("1.1", "1.2", "2.1", "2.2")
ev_all = readRDS(paste0(res_folder_name, "//ev_all.rds"))
vif_all = readRDS(paste0(res_folder_name, "//vif_all.rds"))
tau.long_all_mean_raw = readRDS(paste0(res_folder_name, "//mean_Ktau_pM", p_max,".rds")) %>%
  mutate(
    variable = factor(variable, levels = 3:p_max),
    Z = factor(Z, levels = c("Johnson", "GS", "PC", "VM")),
    Method = factor(Method, levels = c("GDA", "CorPA", "RegPA", "IdA")),
  )

# w/o filtering
rmse.long_all_mean = data.frame()
tau.long_all_mean = data.frame()
for (i in 1:length(p_set)) {
  p = p_set[i]
  print(p)
  
  num_ev = num_ev_set[p-p_min+1]
  
  rmse.long_all_mean.p = readRDS(paste0(res_folder_name, "//mean_rmse_p", p, ".rds")) 
  tau.long_all_mean.p = tau.long_all_mean_raw %>% 
    rename(value = Ktau_mean, var = variable) %>% 
    filter(var == p)
  
  rmse.long_all_mean = rbind(rmse.long_all_mean, rmse.long_all_mean.p)
  tau.long_all_mean = rbind(tau.long_all_mean, tau.long_all_mean.p)
}
rmse.long_all_mean_table <- rmse.long_all_mean %>%
  group_by(Method, Z, var) %>%
  summarize(value = mean(rmse_mean), .groups = "drop") %>%
  pivot_wider(names_from = var, values_from = value) %>%
  as.data.frame()

tau.long_all_mean_table <- tau.long_all_mean %>%
  group_by(Method, Z, var) %>%
  summarize(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = var, values_from = value) %>%
  as.data.frame()

combined <- rmse.long_all_mean_table
combined[, 3:10] <- mapply(
  function(a, b) paste(a, "/", b),
  round(rmse.long_all_mean_table[, 3:10],2), round(tau.long_all_mean_table[, 3:10],2)
)
colnames(combined) = c("Method", "Z", paste0("p = ", 3:10))
write.csv(combined, paste0(fig_folder_name, "//Table 1.csv"))

# w/ filtering
data_count = matrix(0, 8, 4)
p_rmse_list = p_tau_list = list()
t_ev = 1.5
condition_map <- list(
  "1.1" = quote(vif1 <  4 * p & ev1.p <  t_ev / sqrt(p)),
  "1.2" = quote(vif1 <  4 * p & ev1.p >= t_ev / sqrt(p)),
  "2.1" = quote(vif1 >= 4 * p & ev1.p <  t_ev / sqrt(p)),
  "2.2" = quote(vif1 >= 4 * p & ev1.p >= t_ev / sqrt(p))
)

for(j in 1:length(scenario_set)){
  s = scenario_set[j]
  print(s)
  rmse.long_all_mean = data.frame()
  tau.long_all_mean = data.frame()
  bias.long.J_p3710 = data.frame()
  
  for (i in 1:length(p_set)) {
    p = p_set[i]
    print(p)
    num_ev = num_ev_set[p-p_min+1]
    
    rmse.long_all_mean.p = readRDS(paste0(res_folder_name, "//mean_rmse_p", p, ".rds")) 
    tau.long_all_mean.p = tau.long_all_mean_raw %>% 
      rename(value = Ktau_mean, var = variable) %>% 
      filter(var == p)
    bias.long_all_mean.p = readRDS(paste0(res_folder_name, "//mean_bias_p", p, ".rds")) 
    
    if (s %in% names(condition_map)) {
      condition <- condition_map[[s]]
      filtered_data <- apply_filters(list(rmse.long_all_mean.p, tau.long_all_mean.p, bias.long_all_mean.p), condition)
      rmse.long_all_mean.p <- filtered_data[[1]]
      tau.long_all_mean.p <- filtered_data[[2]]
      bias.long_all_mean.p <- filtered_data[[3]]
    }
  
    rmse.long_all_mean = rbind(rmse.long_all_mean, rmse.long_all_mean.p)
    tau.long_all_mean = rbind(tau.long_all_mean, tau.long_all_mean.p)
    data_count[i, j] = dim(rmse.long_all_mean.p)[1]
    
    if(p %in% c(4, 7, 10)){
      bias.long.J_p3710 = rbind(bias.long.J_p3710, bias.long_all_mean.p)
    }
  }
  
  rmse.long.J_mean = rmse.long_all_mean %>%
    rename(value = rmse_mean) %>%
    filter(Z == "Johnson") %>%
    mutate(var = factor(var))
  
  p_rmse_z.J_all.reA = plot_Metric_Z.J_All.reA(rmse.long.J_mean, "RMSE", T)

  tau.long.J_mean = tau.long_all_mean %>%
    filter(Z == "Johnson") %>%
    mutate(var = factor(var))
  
  p_tau_z.J_all.reA = plot_Metric_Z.J_All.reA(tau.long.J_mean, "Kendall's tau")

  p_bias_Z.J_All.reA_wrap_1 = plot_bias_Z.J_All.reA_wrap(bias.long.J_p3710, "bias")
  
  if(s=="1.2"){
    ggsave(paste0(fig_folder_name, "//fig_6.pdf"),
           plot = p_bias_Z.J_All.reA_wrap_1, width = 12, height = 8, dpi = 300)
  }
  p_rmse_list[[s]] = p_rmse_z.J_all.reA
  p_tau_list[[s]] = p_tau_z.J_all.reA
}


fig4 = (p_rmse_list$`1.1` + p_tau_list$`1.1`) +
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') 
ggsave(paste0(fig_folder_name, "//fig_4.pdf"),
       plot = fig4, width = 12, height = 5, dpi = 300)

fig5 = (p_rmse_list$`1.2` + p_tau_list$`1.2`) +
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') 
ggsave(paste0(fig_folder_name, "//fig_5.pdf"),
       plot = fig5, width = 12, height = 5, dpi = 300)

fig8 = (p_rmse_list$`2.1` + p_tau_list$`2.1`) +
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') 
ggsave(paste0(fig_folder_name, "//fig_8.pdf"),
       plot = fig8, width = 12, height = 5, dpi = 300)

fig9 = (p_rmse_list$`2.2` + p_tau_list$`2.2`) +
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') 
ggsave(paste0(fig_folder_name, "//fig_9.pdf"),
       plot = fig9, width = 12, height = 5, dpi = 300)



#####
# Figure 7
# win-loss analysis
# RMSE
all_plot_data <- list()
p_set = c(4,7,10)

for (p in p_set) {
  print(p)
  num_ev = num_ev_set[p-p_min+1]
  threshold = 4 * p
  rmse.long.J_sub_mean = readRDS(paste0(res_folder_name, "//mean_rmse_p", p, ".rds")) %>%
    filter(Z == "Johnson", vif1 < threshold) %>%
    mutate(
      ev1.sqp = ev1.p * sqrt(p)
    )
  
  res = rmse.long.J_sub_mean %>%
    filter(Method %in% c("CorPA", "RegPA")) %>%
    pivot_wider(names_from = Method, values_from = rmse_mean) %>%
    rename(mRMSE.RW = CorPA, mRMSE.GCD = RegPA) 
  
  res.WinLoss = res  %>%
    group_by(ev_index, ev1.p, ev1.sqp) %>%
    summarise(WinLoss = test_relation(mRMSE.RW, mRMSE.GCD), .groups = "drop") %>%
    mutate(
      ev1.p_bin = cut(ev1.p, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE),
      ev1.sqp_bin = cut(ev1.sqp, breaks = seq(0.6, 2.1, by = 0.3), include.lowest = TRUE),
      WinLoss = recode(WinLoss, "greater" = "GCD wins", "tie" = "tie", "smaller" = "RW wins"),
      WinLoss = factor(WinLoss, levels = c("GCD wins", "tie", "RW wins"))
    )
  
  summary_data = res.WinLoss %>%
    group_by(ev1.sqp_bin) %>%
    count(WinLoss) %>%
    mutate(
      proportion = n / sum(n),
      position = 1 - cumsum(proportion) + proportion / 2,
      p = as.factor(p) # Add p to data
    ) %>%
    ungroup()
  
  all_plot_data[[p]] = summary_data
}
all_plot_data_rmse = bind_rows(all_plot_data)
levels(all_plot_data_rmse$p) = paste0("p = ", p_set)


###################################################
### K-tau
all_plot_data <- list()
tau.long_all_mean = readRDS(paste0(res_folder_name, "//mean_Ktau_pM", p_max, ".rds"))
for (p in p_set) {
  print(p)
  num_ev = num_ev_set[p-p_min+1]
  threshold = 4 * p

  ktau.long.J_sub_mean <- tau.long_all_mean %>%
    filter(Z == "Johnson", variable == p, vif1 < threshold) %>%
    mutate(
      ev1.sqp = ev1.p * sqrt(p)
    )

  res = ktau.long.J_sub_mean %>%
    filter(Method %in% c("CorPA", "RegPA")) %>%
    pivot_wider(names_from = Method, values_from = Ktau_mean) %>%
    rename(mKtau.RW = CorPA, mKtau.GCD = RegPA) 
  
  res.WinLoss = res  %>%
    group_by(ev_index, ev1.p, ev1.sqp) %>%
    summarise(WinLoss = test_relation(mKtau.RW, mKtau.GCD), .groups = "drop") %>%
    mutate(
      ev1.p_bin = cut(ev1.p, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE),
      ev1.sqp_bin = cut(ev1.sqp, breaks = seq(0.6, 2.1, by = 0.3), include.lowest = TRUE),
      WinLoss = recode(WinLoss, "greater" = "RW wins", "tie" = "tie", "smaller" = "GCD wins"),
      WinLoss = factor(WinLoss, levels = c("GCD wins", "tie", "RW wins"))
    )
  
  summary_data = res.WinLoss %>%
    group_by(ev1.sqp_bin) %>%
    count(WinLoss) %>%
    mutate(
      proportion = n / sum(n),
      position = 1 - cumsum(proportion) + proportion / 2,
      p = as.factor(p) # Add p to data
    ) %>%
    ungroup()
  
  all_plot_data[[p]] = summary_data
}
all_plot_data_ktau = bind_rows(all_plot_data)
levels(all_plot_data_ktau$p) = paste0("p = ", p_set)


all_plot_data_rmse$metric = "rmse"
all_plot_data_ktau$metric = "ktau"
all_plot_data = rbind(all_plot_data_rmse, all_plot_data_ktau)
all_plot_data$metric = factor(all_plot_data$metric, levels = c("rmse", "ktau"))
all_plot_data = all_plot_data[!is.na(all_plot_data$ev1.sqp_bin), ]
levels(all_plot_data$metric) = c("RMSE", "Kendall's tau") 
# round(proportion * 100)
fig7 = ggplot(all_plot_data, aes(x = ev1.sqp_bin, y = proportion, fill = WinLoss)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(y = position, label = ifelse(proportion >= 0.05, round(proportion * 100), "")),
    color = "white",
    size = 3
  ) +
  labs(
    title = "",
    x = TeX("\\textbf{Interval for} $\\lambda_1/\\sqrt{p}$"),
    y = "Proportion",
    fill = "Result"
  ) +
  scale_fill_manual(
    values = c("GCD wins" = "#4C72B0", "tie" = "#999999", "RW wins" = "#C44E52"),
    limits = c("GCD wins", "tie", "RW wins")
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black") +  # Add the dashed vertical line
  facet_grid(rows = vars(p), cols = vars(metric)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(color = "black")
  )
fig7
ggsave(paste0(fig_folder_name, "//fig_7.pdf"),
       plot = fig7, width = 12, height = 8, dpi = 300)

# Figure 10
# win-loss analysis
# RMSE
all_plot_data <- list()
p_set = c(4,7,10)

for (p in p_set) {
  print(p)
  num_ev = num_ev_set[p-p_min+1]
  threshold = 4 * p
  rmse.long.J_sub_mean = readRDS(paste0(res_folder_name, "//mean_rmse_p", p, ".rds")) %>%
    filter(Z == "Johnson", vif1 >= threshold) %>%
    mutate(
      ev1.sqp = ev1.p * sqrt(p)
    )
  
  res = rmse.long.J_sub_mean %>%
    filter(Method %in% c("CorPA", "RegPA")) %>%
    pivot_wider(names_from = Method, values_from = rmse_mean) %>%
    rename(mRMSE.RW = CorPA, mRMSE.GCD = RegPA) 
  
  res.WinLoss = res  %>%
    group_by(ev_index, ev1.p, ev1.sqp) %>%
    summarise(WinLoss = test_relation(mRMSE.RW, mRMSE.GCD), .groups = "drop") %>%
    mutate(
      ev1.p_bin = cut(ev1.p, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE),
      ev1.sqp_bin = cut(ev1.sqp, breaks = seq(0.6, 2.1, by = 0.3), include.lowest = TRUE),
      WinLoss = recode(WinLoss, "greater" = "GCD wins", "tie" = "tie", "smaller" = "RW wins"),
      WinLoss = factor(WinLoss, levels = c("GCD wins", "tie", "RW wins"))
    )
  
  summary_data = res.WinLoss %>%
    group_by(ev1.sqp_bin) %>%
    count(WinLoss) %>%
    mutate(
      proportion = n / sum(n),
      position = 1 - cumsum(proportion) + proportion / 2,
      p = as.factor(p) # Add p to data
    ) %>%
    ungroup()
  
  all_plot_data[[p]] = summary_data
}
all_plot_data_rmse = bind_rows(all_plot_data)
levels(all_plot_data_rmse$p) = paste0("p = ", p_set)


###################################################
### K-tau
all_plot_data <- list()
tau.long_all_mean = readRDS(paste0(res_folder_name, "//mean_Ktau_pM", p_max, ".rds"))
for (p in p_set) {
  print(p)
  num_ev = num_ev_set[p-p_min+1]
  threshold = 4 * p
  
  ktau.long.J_sub_mean <- tau.long_all_mean %>%
    filter(Z == "Johnson", variable == p, vif1 >= threshold) %>%
    mutate(
      ev1.sqp = ev1.p * sqrt(p)
    )
  
  res = ktau.long.J_sub_mean %>%
    filter(Method %in% c("CorPA", "RegPA")) %>%
    pivot_wider(names_from = Method, values_from = Ktau_mean) %>%
    rename(mKtau.RW = CorPA, mKtau.GCD = RegPA) 
  
  res.WinLoss = res  %>%
    group_by(ev_index, ev1.p, ev1.sqp) %>%
    summarise(WinLoss = test_relation(mKtau.RW, mKtau.GCD), .groups = "drop") %>%
    mutate(
      ev1.p_bin = cut(ev1.p, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE),
      ev1.sqp_bin = cut(ev1.sqp, breaks = seq(0.6, 2.1, by = 0.3), include.lowest = TRUE),
      WinLoss = recode(WinLoss, "greater" = "RW wins", "tie" = "tie", "smaller" = "GCD wins"),
      WinLoss = factor(WinLoss, levels = c("GCD wins", "tie", "RW wins"))
    )
  
  summary_data = res.WinLoss %>%
    group_by(ev1.sqp_bin) %>%
    count(WinLoss) %>%
    mutate(
      proportion = n / sum(n),
      position = 1 - cumsum(proportion) + proportion / 2,
      p = as.factor(p) # Add p to data
    ) %>%
    ungroup()
  
  all_plot_data[[p]] = summary_data
}
all_plot_data_ktau = bind_rows(all_plot_data)
levels(all_plot_data_ktau$p) = paste0("p = ", p_set)


all_plot_data_rmse$metric = "rmse"
all_plot_data_ktau$metric = "ktau"
all_plot_data = rbind(all_plot_data_rmse, all_plot_data_ktau)
all_plot_data$metric = factor(all_plot_data$metric, levels = c("rmse", "ktau"))
all_plot_data = all_plot_data[!is.na(all_plot_data$ev1.sqp_bin), ]
levels(all_plot_data$metric) = c("RMSE", "Kendall's tau") 

fig10 = ggplot(all_plot_data, aes(x = ev1.sqp_bin, y = proportion, fill = WinLoss)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(y = position, label = ifelse(proportion >= 0.05, round(proportion * 100), "")),
    color = "white",
    size = 3
  ) +
  labs(
    title = "",
    x = TeX("\\textbf{Interval for} $\\lambda_1/\\sqrt{p}$"),
    y = "Proportion",
    fill = "Result"
  ) +
  scale_fill_manual(
    values = c("GCD wins" = "#4C72B0", "tie" = "#999999", "RW wins" = "#C44E52"),
    limits = c("GCD wins", "tie", "RW wins")
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black") +  # Add the dashed vertical line
  facet_grid(rows = vars(p), cols = vars(metric)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(color = "black")
  )
fig10
ggsave(paste0(fig_folder_name, "//fig_10.pdf"),
       plot = fig10, width = 12, height = 8, dpi = 300)
