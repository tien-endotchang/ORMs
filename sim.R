rm(list = ls())
library(svd)
library(MASS)
library(relaimpo)
library(GPArotation)
library(reshape2)
library(fungible) 
library(dplyr)
source("utils_sim.R")


sim.master = function(p_min, p_max, num_ev, num_rep = 10, num_y = 100){
  p_set = p_min:p_max
  
  df.time = data.frame(matrix(0, nrow = length(p_set) * num_ev * num_rep, ncol = 4))
  colnames(df.time) = c("p", "ev_index", "seed", "time")

  num_methods = 4  # Number of methods
  tau_columns = p_max - p_min + 1 + 3
  num_iterations = num_y * num_rep * num_ev
  
  res.tau_GDA = array(NA, dim = c(num_iterations, tau_columns, num_methods))
  res.tau_CorPA = array(NA, dim = c(num_iterations, tau_columns, num_methods))
  res.tau_RegPA = array(NA, dim = c(num_iterations, tau_columns, num_methods))
  res.tau_IdA = array(NA, dim = c(num_iterations, tau_columns, num_methods))
  
  ev_all = array(NA, dim = c(num_ev, p_max+1, length(p_set)))
  vif_all = array(NA, dim = c(num_ev * num_rep, p_max+2, length(p_set)))
  
  start_total = Sys.time()
  time_count = 1
  
  for(k in 1:length(p_set)){
    count = 1
    start_p = Sys.time()
    p = p_set[k]
    num_comparisons = p + 3
    
    res.ordered_bias_GDA = array(NA, dim = c(num_iterations, num_comparisons, num_methods))
    res.ordered_bias_CorPA = array(NA, dim = c(num_iterations, num_comparisons, num_methods))
    res.ordered_bias_RegPA = array(NA, dim = c(num_iterations, num_comparisons, num_methods))
    res.ordered_bias_IdA = array(NA, dim = c(num_iterations, num_comparisons, num_methods))
    
    set.seed(2) ### for reproducibility
    ev_pre = runif_simplex(num_ev, p, p)
    res.ev = data.frame(matrix(0, nrow = num_ev, ncol = p + 1))
    res.vif = data.frame(matrix(0, nrow = num_rep * num_ev, ncol = p + 2))
    colnames(res.ev) = c(paste0("ev", 1:p), "ev_index")
    colnames(res.vif) = c(paste0("vif", 1:p), "ev_index", "seed")
    
    for(j in 1:num_ev){
      start_ev = Sys.time()
      
      for(c in 1:num_rep){
        start = Sys.time()
        
        cat(paste0("number of predictors: ", p, ": the ", j, " th correlation matrix", ": seed: ", c, "\n"))
        correlation = rMAP(ev_pre[j, ], Seed = c)$R
        
        ### COR ONLY ###
        corxx.SVD = svd(correlation)
        V = corxx.SVD$v
        sqrtD = sqrt(corxx.SVD$d)
        invSqrtD = 1 / sqrtD
        lambda = (V %*% diag(sqrtD) %*% t(V))
        
        Q_pc = V
        Q_G = V %*% diag(sqrtD) %*% t(V) %*% solve(chol(correlation)) # cholesky ~~
        Q_vm = varimax(lambda)$rotmat
        
        GD_ZX = GenerateGDZX(correlation, NULL, corxx.SVD)
        GD_ZpcX = GenerateGDZX(correlation, Q_pc, corxx.SVD)
        GD_ZGX = GenerateGDZX(correlation, Q_G, corxx.SVD)
        GD_ZvmX = GenerateGDZX(correlation, Q_vm, corxx.SVD)
        
        CorPA = GenerateCorPA(corxx=correlation, Q=NULL, corxx.SVD=corxx.SVD)
        CorPA_Zpc = GenerateCorPA(corxx=correlation, Q=Q_pc, corxx.SVD=corxx.SVD)
        CorPA_ZG = GenerateCorPA(corxx=correlation, Q=Q_G, corxx.SVD=corxx.SVD)
        CorPA_Zvm = GenerateCorPA(corxx=correlation, Q=Q_vm, corxx.SVD=corxx.SVD)

        RegPA = GenerateRegPA(corxx=correlation, Q=NULL, corxx.SVD=corxx.SVD)
        RegPA_Zpc = GenerateRegPA(corxx=correlation, Q=Q_pc, corxx.SVD=corxx.SVD)
        RegPA_ZG = GenerateRegPA(corxx=correlation, Q=Q_G, corxx.SVD=corxx.SVD)
        RegPA_Zvm = GenerateRegPA(corxx=correlation, Q=Q_vm, corxx.SVD=corxx.SVD)
        ### COR ONLY ###
        
        ev = eigen(correlation)$values
        vif = sort(diag(solve(correlation)), decreasing = TRUE)
        
        sphere_points = rsphere_n(n = num_y, dim = p) # different y
        for(i in 1:num_y){
          beta = sphere_points[i, ]
          
          ### COR ONLY ###
          ryx = (corxx.SVD$v %*% diag(sqrt(corxx.SVD$d)) %*% t(corxx.SVD$v)) %*% beta / sqrt(1 + 0.5^2)
          coryx = diag(1, p+1)
          coryx[2:(p+1), 1] = coryx[1, 2:(p+1)] = ryx
          coryx[2:(p+1), 2:(p+1)] = correlation
          GD = unname(calc.relimp(coryx)@lmg)
          
          beta_Z = as.numeric( (V %*% diag(invSqrtD) %*% t(V)) %*% ryx )
          beta_Z_pc = as.numeric( t(Q_pc) %*% (V %*% diag(invSqrtD) %*% t(V)) %*% ryx )
          beta_Z_G = as.numeric( t(Q_G) %*% (V %*% diag(invSqrtD) %*% t(V)) %*% ryx )
          beta_Z_vm = as.numeric( t(Q_vm) %*% (V %*% diag(invSqrtD) %*% t(V)) %*% ryx )
          
          epsilon = as.numeric((CorPA) %*% beta_Z^2)
          epsilon_pc = as.numeric((CorPA_Zpc) %*% beta_Z_pc^2) 
          epsilon_G = as.numeric((CorPA_ZG) %*% beta_Z_G^2)
          epsilon_vm = as.numeric((CorPA_Zvm) %*% beta_Z_vm^2)
          
          GDA = as.numeric((GD_ZX) %*% beta_Z^2)
          GDA_pc = as.numeric((GD_ZpcX) %*% beta_Z_pc^2)
          GDA_G = as.numeric((GD_ZGX) %*% beta_Z_G^2)
          GDA_vm = as.numeric((GD_ZvmX) %*% beta_Z_vm^2)
          
          delta = as.numeric(RegPA %*% beta_Z^2)
          delta_pc = as.numeric(RegPA_Zpc %*% beta_Z_pc^2)
          delta_G = as.numeric(RegPA_ZG %*% beta_Z_G^2)
          delta_vm = as.numeric(RegPA_Zvm %*% beta_Z_vm^2)
          ### COR ONLY ###
          
          # estimating bias
          res.ordered_bias_GDA[count, , 1] = c((GDA - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_GDA[count, , 2] = c((GDA_pc - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_GDA[count, , 3] = c((GDA_G - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_GDA[count, , 4] = c((GDA_vm - GD)[order(GD, decreasing = TRUE)], i, c, j)
          
          res.ordered_bias_CorPA[count, , 1] = c((epsilon - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_CorPA[count, , 2] = c((epsilon_pc - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_CorPA[count, , 3] = c((epsilon_G - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_CorPA[count, , 4] = c((epsilon_vm - GD)[order(GD, decreasing = TRUE)], i, c, j)
          
          res.ordered_bias_RegPA[count, , 1] = c((delta - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_RegPA[count, , 2] = c((delta_pc - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_RegPA[count, , 3] = c((delta_G - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_RegPA[count, , 4] = c((delta_vm - GD)[order(GD, decreasing = TRUE)], i, c, j)
          
          res.ordered_bias_IdA[count, , 1] = c((beta_Z^2 - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_IdA[count, , 2] = c((beta_Z_pc^2 - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_IdA[count, , 3] = c((beta_Z_G^2 - GD)[order(GD, decreasing = TRUE)], i, c, j)
          res.ordered_bias_IdA[count, , 4] = c((beta_Z_vm^2 - GD)[order(GD, decreasing = TRUE)], i, c, j)
          
          # Kendall's tau
          res.tau_GDA[count, p - p_min + 1, 1] = cor(GD, GDA, method = "kendall")
          res.tau_GDA[count, p - p_min + 1, 2] = cor(GD, GDA_pc, method = "kendall")
          res.tau_GDA[count, p - p_min + 1, 3] = cor(GD, GDA_G, method = "kendall")
          res.tau_GDA[count, p - p_min + 1, 4] = cor(GD, GDA_vm, method = "kendall")
          
          res.tau_CorPA[count, p - p_min + 1, 1] = cor(GD, epsilon, method = "kendall")
          res.tau_CorPA[count, p - p_min + 1, 2] = cor(GD, epsilon_pc, method = "kendall")
          res.tau_CorPA[count, p - p_min + 1, 3] = cor(GD, epsilon_G, method = "kendall")
          res.tau_CorPA[count, p - p_min + 1, 4] = cor(GD, epsilon_vm, method = "kendall")
          
          res.tau_RegPA[count, p - p_min + 1, 1] = cor(GD, delta, method = "kendall")
          res.tau_RegPA[count, p - p_min + 1, 2] = cor(GD, delta_pc, method = "kendall")
          res.tau_RegPA[count, p - p_min + 1, 3] = cor(GD, delta_G, method = "kendall")
          res.tau_RegPA[count, p - p_min + 1, 4] = cor(GD, delta_vm, method = "kendall")
          
          res.tau_IdA[count, p - p_min + 1, 1] = cor(GD, beta_Z^2, method = "kendall")
          res.tau_IdA[count, p - p_min + 1, 2] = cor(GD, beta_Z_pc^2, method = "kendall")
          res.tau_IdA[count, p - p_min + 1, 3] = cor(GD, beta_Z_G^2, method = "kendall")
          res.tau_IdA[count, p - p_min + 1, 4] = cor(GD, beta_Z_vm^2, method = "kendall")
          
          count = count + 1
        }
        
        res.vif[c+(j-1)*num_rep, 1:(p+2)] = c(vif, j, c)
        
        end = Sys.time()
        df.time[time_count, 1:4] = c(p, j, c, end - start)
        time_count = time_count + 1
      }
      
      res.ev[j, 1:(p+1)] = c(ev, j)
      end_ev = Sys.time()
      cat(paste0("Finish repetition: "))
      print(end_ev - start_ev)
    }
    ev_all[1:num_ev, 1:(p+1), k] = as.matrix(res.ev)
    vif_all[1:(num_ev*num_rep), 1:(p+2), k] = as.matrix(res.vif)
    
    bias.long_GDA = resArr2Long(res.ordered_bias_GDA, p, NULL,"bias", "GDA")
    bias.long_CorPA = resArr2Long(res.ordered_bias_CorPA, p, NULL, "bias", "CorPA")
    bias.long_RegPA = resArr2Long(res.ordered_bias_RegPA, p, NULL, "bias", "RegPA")
    bias.long_IdA = resArr2Long(res.ordered_bias_IdA, p, NULL, "bias", "IdA")
    bias.long_all = rbind(bias.long_GDA, bias.long_IdA, bias.long_RegPA, bias.long_CorPA)
    
    bias.long_all = bias.long_all %>%
      mutate(
        variable = factor(variable, levels = 1:p),
        Z = factor(Z, levels = c("Johnson", "GS", "PC", "VM")),
        Method = factor(Method, levels = c("GDA", "CorPA", "RegPA", "IdA")),
        ev1.p = ev_all[ev_index, 1, p - p_min + 1] / p,
        vif1 = vif_all[(ev_index - 1) * num_rep + seed, 1, p - p_min + 1],
        kappa = ev_all[ev_index, 1, p - p_min + 1] / ev_all[ev_index, p, p - p_min + 1]
      )
    
    rmse.long_all_mean = bias.long_all %>%
      group_by(idy, ev_index, seed, Z, Method, ev1.p, vif1, kappa) %>%
      summarize(rmse = sqrt(mse(value)), .groups = "drop") %>%
      group_by(ev_index, seed, Z, Method, ev1.p, vif1, kappa) %>%
      summarize(rmse_mean = mean(rmse), .groups = "drop") %>%
      mutate(var = factor(p, levels = 3:10)) %>%
      as.data.frame()

    bias.long.J_mean = bias.long_all %>%
      filter(Z == "Johnson", Method %in% c("GDA", "CorPA", "RegPA", "IdA")) %>%
      group_by(ev_index, seed, variable, Z, Method, ev1.p, vif1, kappa) %>%
      summarize(value = mean(value), .groups = "drop") %>%
      mutate(num_p = p)%>%
      as.data.frame()

    saveRDS(rmse.long_all_mean, paste0("mean_rmse_p", p, ".rds"), compress = FALSE)
    saveRDS(bias.long.J_mean, paste0("mean_bias_p", p, ".rds"), compress = FALSE)

    rm(bias.long_GDA, bias.long_CorPA, bias.long_RegPA, bias.long_IdA, bias.long_all)
    
    end_p = Sys.time()
    cat(paste0("Finish p = ", p, " : "))
    print(end_p - start_p)
  }
  
  for(i in 1:4){
    res.tau_GDA[, p - p_min + 2, i] = res.ordered_bias_GDA[, p+1, 1]
    res.tau_GDA[, p - p_min + 3, i] = res.ordered_bias_GDA[, p+2, 1]
    res.tau_GDA[, p - p_min + 4, i] = res.ordered_bias_GDA[, p+3, 1]
    
    res.tau_CorPA[, p - p_min + 2, i] = res.ordered_bias_GDA[, p+1, 1]
    res.tau_CorPA[, p - p_min + 3, i] = res.ordered_bias_GDA[, p+2, 1]
    res.tau_CorPA[, p - p_min + 4, i] = res.ordered_bias_GDA[, p+3, 1]
    
    res.tau_RegPA[, p - p_min + 2, i] = res.ordered_bias_GDA[, p+1, 1]
    res.tau_RegPA[, p - p_min + 3, i] = res.ordered_bias_GDA[, p+2, 1]
    res.tau_RegPA[, p - p_min + 4, i] = res.ordered_bias_GDA[, p+3, 1]
    
    res.tau_IdA[, p - p_min + 2, i] = res.ordered_bias_GDA[, p+1, 1]
    res.tau_IdA[, p - p_min + 3, i] = res.ordered_bias_GDA[, p+2, 1]
    res.tau_IdA[, p - p_min + 4, i] = res.ordered_bias_GDA[, p+3, 1]
  }
  
  tau.long_GDA = resArr2Long(res.tau_GDA, p, p_min, "tau", "GDA")
  tau.long_CorPA = resArr2Long(res.tau_CorPA, p, p_min, "tau", "CorPA")
  tau.long_RegPA = resArr2Long(res.tau_RegPA, p, p_min, "tau", "RegPA")
  tau.long_IdA = resArr2Long(res.tau_IdA, p, p_min, "tau", "IdA")
  tau.long_all = rbind(tau.long_GDA, tau.long_IdA, tau.long_RegPA, tau.long_CorPA)
  
  saveRDS(tau.long_all, paste0("Ktau_pm", p_min, "_pM", p_max, ".rds"), compress = TRUE)
  saveRDS(ev_all, paste0("ev_all_", p_min, "_", p_max, ".rds"), compress = FALSE)
  saveRDS(vif_all, paste0("vif_all_", p_min, "_", p_max, ".rds"), compress = FALSE)

  end_total = Sys.time()
  print(end_total - start_total) # 2.2 days
  
}

################################
folder_name = "res"
dir.create("res")
setwd(paste0("./", folder_name))
# folder_name = "Rtmp"
# dir.create("C:/Rtmp")
# setwd(paste0("C:/Rtmp"))
num_rep = 10
num_y = 100

# small
p_min = 3
p_max = 6
num_ev_s = 1000
sim.master(p_min=p_min, p_max=p_max, num_ev=num_ev_s, num_rep=num_rep, num_y=num_y)

# large - 1
p_min = 7
p_max = 9
num_ev_l = 2500
sim.master(p_min=p_min, p_max=p_max, num_ev=num_ev_l, num_rep=num_rep, num_y=num_y)

# large - 2
p_min = 10
p_max = 10
num_ev_l = 2500
sim.master(p_min=p_min, p_max=p_max, num_ev=num_ev_l, num_rep=num_rep, num_y=num_y)

tau.long_all_3_6 = readRDS(paste0("Ktau_pm", 3, "_pM", 6, ".rds"))
tau.long_all_7_9 = readRDS(paste0("Ktau_pm", 7, "_pM", 9, ".rds"))
tau.long_all_10_10 = readRDS(paste0("Ktau_pm", 10, "_pM", 10, ".rds"))
tau.long_all = rbind(tau.long_all_3_6, tau.long_all_7_9, tau.long_all_10_10)
saveRDS(tau.long_all, paste0("Ktau_pM", 10, ".rds"))

ev_all = array(NA, dim = c(num_ev_l, 11, 8))
ev_all_3_6 = readRDS(paste0("ev_all_", 3, "_", 6, ".rds"))
ev_all_7_9 = readRDS(paste0("ev_all_", 7, "_", 9, ".rds"))
ev_all_10_10 = readRDS(paste0("ev_all_", 10, "_", 10, ".rds"))
ev_all[1:num_ev_s, 1:7, 1:4] = ev_all_3_6
ev_all[1:num_ev_l, 1:10, 5:7] = ev_all_7_9
ev_all[1:num_ev_l, 1:11, 8:8] = ev_all_10_10
saveRDS(ev_all, paste0("ev_all.rds"))

vif_all = array(NA, dim = c(num_ev_l * num_rep, 12, 8))
vif_all_3_6 = readRDS(paste0("vif_all_", 3, "_", 6, ".rds"))
vif_all_7_9 = readRDS(paste0("vif_all_", 7, "_", 9, ".rds"))
vif_all_10_10 = readRDS(paste0("vif_all_", 10, "_", 10, ".rds"))
vif_all[1:(num_ev_s * num_rep), 1:8, 1:4] = vif_all_3_6
vif_all[1:(num_ev_l * num_rep), 1:11, 5:7] = vif_all_7_9
vif_all[1:(num_ev_l * num_rep), 1:12, 8:8] = vif_all_10_10
saveRDS(vif_all, paste0("vif_all.rds"))

################################
ev_all = readRDS(paste0("ev_all.rds"))
vif_all = readRDS(paste0("vif_all.rds"))
p_min = 3
p_max = 10

tau.long_all = readRDS(paste0("Ktau_pM", 10, ".rds"))
tau.long_all_mean.temp <- tau.long_all %>%
  group_by(variable, Z, Method, ev_index, seed) %>%
  summarize(Ktau_mean = mean(value), .groups = "drop") 
tau.long_all_mean = data.frame()

for(i in 1:length(p_min:p_max)){
  p = (p_min:p_max)[i]
  print(p)
  
  tau.long_all_mean.p <- tau.long_all_mean.temp %>%
    filter(variable == as.character(p)) %>%
    mutate(
      ev1.p = ev_all[ev_index, 1, p - p_min + 1] / p,
      vif1 = vif_all[(ev_index - 1) * num_rep + seed, 1, p - p_min + 1],
      kappa = ev_all[ev_index, 1, p - p_min + 1] / ev_all[ev_index, p, p - p_min + 1]
    ) %>%
    as.data.frame()
  
  tau.long_all_mean = rbind(tau.long_all_mean, tau.long_all_mean.p)
}

saveRDS(tau.long_all_mean, paste0("mean_Ktau_pM", p_max, ".rds"))