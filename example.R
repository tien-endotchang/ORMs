rm(list = ls())
library(relaimpo)
source("utils_sim.R")


# Suh et al. (1998)
# used in Azen (2003)
Suh_yx = matrix(c(1.0000, 0.2346, 0.3637, 0.4875, 0.3162, 0.6208,
                  0.2346, 1.0000, 0.1775, 0.1125, 0.1510, 0.3425,
                  0.3637, 0.1775, 1.0000, 0.2614, 0.2384, 0.2899,
                  0.4875, 0.1125, 0.2614, 1.0000, 0.2490, 0.3818,
                  0.3162, 0.1510, 0.2384, 0.2490, 1.0000, 0.1918,
                  0.6208, 0.3425, 0.2899, 0.3818, 0.1918, 1.0000), nrow = 6)
names = c("Y", "Health", "Finance", "Famliy", "Housing", "Self")
row.names(Suh_yx) = colnames(Suh_yx) = names
Suh_yx.obj = generate_realdata_result(Suh_yx)
print(Suh_yx.obj)
write.csv(data.frame(Suh_yx.obj$table), paste0("fig//Suh1998.csv"))
# In Johnson (2001), 
corxx = matrix(c(1.00, 0.72, 0.53, 0.49, 0.42, 0.56, 0.48,
                 0.72, 1.00, 0.66, 0.67, 0.51, 0.68, 0.60,
                 0.53, 0.66, 1.00, 0.59, 0.44, 0.52, 0.51,
                 0.49, 0.67, 0.59, 1.00, 0.58, 0.60, 0.61,
                 0.42, 0.51, 0.44, 0.58, 1.00, 0.53, 0.56,
                 0.56, 0.68, 0.52, 0.60, 0.53, 1.00, 0.58,
                 0.48, 0.60, 0.51, 0.61, 0.56, 0.58, 1.00), nrow = 7)

beta_yx4 = matrix(c(0.200, 0.251, 0.148, 0.135, 0.116, 0.122, 0.107), ncol = 1)
rho_yx4 = corxx %*% beta_yx4
Johnson_yx4 = diag(8)
Johnson_yx4[2:8, 2:8] = corxx
Johnson_yx4[1, 2:8] = rho_yx4
Johnson_yx4[2:8, 1] = rho_yx4
colnames(Johnson_yx4) = row.names(Johnson_yx4) = c("operator", "jstp", "njstp", "woctp", "icp", "ocp", "jtc", "hws")
Johnson_yx.obj = generate_realdata_result(Johnson_yx4)
print(Johnson_yx.obj)
write.csv(data.frame(Johnson_yx.obj$table), paste0("fig//Johnson2001.csv"))

# rw.obj = RelativeWeight(Johnson_yx4)
# rw.obj$RW / sum(rw.obj$RW) * 100

