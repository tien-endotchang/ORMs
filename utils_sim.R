runif_simplex = function(n, d, K){
  out = matrix(0, n, d)
  for(i in 1:n){
    s = runif(d - 1, 0, K)
    s = sort(c(0, s, K))
    out[i, 1:d] = diff(s)
  }
  return( t(apply(out, 1, sort, decreasing = TRUE)) )
}

rsphere_n = function(n, dim, r = 1.0) {
  v = matrix(0, nrow = n, ncol = dim)
  for(j in 1:n){
    for(i in 1:dim){v[j, i] = rnorm(1)}
    m = sqrt(sum(v[j, ]^2))
    v[j, ] = v[j, ] / m
  }
  v
}

GenerateGDZX = function(corxx, Q = NULL, corxx.SVD = NULL){
  p = dim(corxx)[1]
  if(is.null(corxx.SVD)) corxx.SVD = svd(corxx)
  if(is.null(Q)) Q = diag(p)
  V = corxx.SVD$v
  sqrtD = sqrt(corxx.SVD$d)
  corxz = (V %*% diag(sqrtD) %*% t(V)) %*% Q
  
  GD_ZX = matrix(0, p, p)
  for(i in 1:p){
    corzix = diag(1, p+1)
    corzix[2:(p+1), 1] = corzix[1, 2:(p+1)] = corxz[, i]
    corzix[2:(p+1), 2:(p+1)] = corxx    
    GD_zi = unname(calc.relimp(corzix + (1E-12)*diag(p+1))@lmg) # compromise
    GD_ZX[, i] = GD_zi
  }
  return(GD_ZX)
}

RelativeWeight = function(coryx, Q = NULL, corxx.SVD = NULL){
  p = dim(coryx)[1] - 1
  ryx = coryx[1, 2:(p+1)]
  corxx = coryx[2:(p+1), 2:(p+1)]
  
  if(is.null(corxx.SVD)) corxx.SVD = svd(corxx)
  V = corxx.SVD$v
  sqrtD = sqrt(corxx.SVD$d)
  invSqrtD = 1 / sqrtD
  
  if(is.null(Q)) Q = diag(p)
  ryz = t(Q) %*% (V %*% diag(invSqrtD) %*% t(V)) %*% ryx
  corxz = (V %*% diag(sqrtD) %*% t(V)) %*% Q
  corPA.raw = corxz^2
  corPA = corPA.raw / matrix(colSums(corPA.raw), p, p, byrow = TRUE)
  RW = as.numeric(corPA %*% ryz^2)
  return( list(ReA = corPA, RW = RW) )
}

GenerateCorPA = function(corxx, Q = NULL, corxx.SVD = NULL){
  p = dim(corxx)[1]
  
  if(is.null(corxx.SVD)) corxx.SVD = svd(corxx)
  V = corxx.SVD$v
  sqrtD = sqrt(corxx.SVD$d)
  invSqrtD = 1 / sqrtD
  
  if(is.null(Q)) Q = diag(p)
  corxz = (V %*% diag(sqrtD) %*% t(V)) %*% Q
  CorPA.raw = corxz^2
  CorPA = CorPA.raw / matrix(colSums(CorPA.raw), p, p, byrow = TRUE)
  return( CorPA )
}

GenerateGreen = function(coryx, Q = NULL, corxx.SVD = NULL){
  p = dim(coryx)[1] - 1
  ryx = coryx[1, 2:(p+1)]
  corxx = coryx[2:(p+1), 2:(p+1)]
  
  if(is.null(corxx.SVD)) corxx.SVD = svd(corxx)
  V = corxx.SVD$v
  sqrtD = sqrt(corxx.SVD$d)
  invSqrtD = 1 / sqrtD
  
  if(is.null(Q)) Q = diag(p)
  ryz = t(Q) %*% (V %*% diag(invSqrtD) %*% t(V)) %*% ryx
  
  regxz = (V %*% diag(invSqrtD) %*% t(V)) %*% Q
  regPA.raw = regxz^2
  regPA = regPA.raw / matrix(colSums(regPA.raw), p, p, byrow = TRUE)
  GCD = as.numeric(regPA %*% ryz^2)
  return( list(ReA = regPA, GCD = GCD) )
}

GenerateRegPA = function(corxx, Q = NULL, corxx.SVD = NULL){
  p = dim(corxx)[1]
  
  if(is.null(corxx.SVD)) corxx.SVD = svd(corxx)
  V = corxx.SVD$v
  sqrtD = sqrt(corxx.SVD$d)
  invSqrtD = 1 / sqrtD
  
  if(is.null(Q)) Q = diag(p)
  regxz = (V %*% diag(invSqrtD) %*% t(V)) %*% Q
  regPA.raw = regxz^2
  regPA = regPA.raw / matrix(colSums(regPA.raw), p, p, byrow = TRUE)
  return( regPA )
}

resArr2Long = function(resArr, p, p_min=NULL, metric = c("bias", "tau"),
                       method = c("GDA", "IdA", "RegPA", "CorPA"),
                       Z_list = c("Johnson", "GS", "PC", "VM")){
  res.long = data.frame()
  if(metric == "bias"){
    colname = c(1:p, "idy", "seed", "ev_index")
  }else{
    colname = c(p_min:p, "idy", "seed", "ev_index")
  }
  for(i in 1:length(Z_list)){
    temp = data.frame(resArr[, , i])
    colnames(temp) = colname
    long.temp = melt(temp, id.vars = c("idy", "ev_index", "seed"))
    if(i == 1){
      long.temp$Z = "Johnson"
    }
    if(i == 2){
      long.temp$Z = "GS"
    }
    if(i == 3){
      long.temp$Z = "PC"
    }
    if(i == 4){
      long.temp$Z = "VM"
    }
    res.long = rbind(res.long, long.temp)
  }
  res.long$variable = factor(res.long$variable)
  res.long$Z = factor(res.long$Z, levels = Z_list)
  res.long$Method = method
  res.long$Method = factor(res.long$Method, levels = c("GDA", "IdA", "RegPA", "CorPA"))
  return( res.long )
}

# Return the mean squared error
mse = function(e) return ( sum(e^2) / length(e) )

generate_realdata_result = function(coryx){
  p = dim(coryx)[1] - 1
  GD = calc.relimp(coryx)@lmg
  RW.obj = RelativeWeight(coryx, Q = NULL, corxx.SVD = NULL)
  RW = RW.obj$RW
  GCD.obj = GenerateGreen(coryx, Q = NULL, corxx.SVD = NULL)
  GCD = GCD.obj$GCD

  regPA = GCD.obj$ReA
  corxx = coryx[2:(p+1), 2:(p+1)]
  
  vif = diag(solve(corxx))
  ev = eigen(corxx)$values
  
  df.table = data.frame(matrix(0, nrow = 7, ncol = p))
  row.names(df.table) = c("GD", "GCD", "RW",
                          "row-sum of regPA",
                          "VIF",
                          "GCD-GD", "RW-GD")
  colnames(df.table) = colnames(coryx)[2:(p+1)][order(GD, decreasing = T)]
  df.table[1:7, 1:p] = matrix(c(round(GD[order(GD, decreasing = T)] / sum(GD) * 100, 5),
                                round(GCD[order(GD, decreasing = T)] / sum(GD) * 100, 5),
                                round(RW[order(GD, decreasing = T)] / sum(GD) * 100, 5),
                                round(apply(regPA, 1, sum)[order(GD, decreasing = T)], 5),
                                round(vif[order(GD, decreasing = T)], 5),
                                round((GCD - GD)[order(GD, decreasing = T)] / sum(GD) * 100, 5),
                                round((RW - GD)[order(GD, decreasing = T)] / sum(GD) * 100, 5)),
                              nrow = 7, byrow = TRUE)
  df.summary = c(p,
                 sum(GD),
                 round(ev[1]/sqrt(p), 5),
                 round(max(vif)/p, 5),
                 round(sqrt(mean((GD - GCD)^2)), 5),
                 round(sqrt(mean((GD - RW)^2)), 5),
                 round(cor(GD, GCD, method = "kendall"), 5),
                 round(cor(GD, RW, method = "kendall"), 5))
  names(df.summary) = c("p", "R2", "ev1/sqrt(p)", "vif1/p", "RMSE.GCD", "RMSE.RW", "Ktau.GCD", "Ktau.RW")
  
  
  return( list(table = df.table,
               summary = df.summary) )
}

# Generate optimal Z by optimization through Wolfram Mathematica
# generate_lambda_optimal.Z_corPA = function(correlation){
#   scriptPath = "MinSRT_Wrap.wl"
#   command <- paste0("wolframscript -file", 
#                     " ", 
#                     scriptPath, 
#                     "  (", 
#                     correlation[1, 2], ") (", 
#                     correlation[1, 3], ") (", 
#                     correlation[2, 3], ")")
#   result = system(command, intern = TRUE)
#   
#   if(length(result) != 1) return( "NG" ) # Constraints not satisfied
#   if(result == "Not Good") return( "NG" ) # Bad local optima
#   
#   parts = unlist(strsplit(result, "\\}\\{"))
#   parts = gsub("[{}]", "", parts)
#   numbers = as.numeric(unlist(strsplit(parts, ", ")))
#   lambda_BcorPA = matrix(numbers, nrow = 3, byrow = T)
#   return( lambda_BcorPA )
# }