rm(list=ls())

source("functions.txt")

library(multiColl) # install.packages("multiColl")
library(lrmest) # install.packages("lrmest")

#########################################################################################################
# data

data(pcd) # data in 'lrmtest' package
y = pcd[,1]
n <- length(y)
cte = rep(1, n)
x2 = pcd[,2]
x3 = pcd[,3]
x4 = pcd[,4]
x5 = pcd[,5]
X = cbind(cte, x2, x3, x4, x5)
p = ncol(X)

#########################################################################################################
# Multicollinearity detection: Table 1

  Z <- cbind(y,X[,-1]) 
  round(cor(Z),digits=4) # Matrix of correlation

  round(CVs(X),digits=4) # Coefficient of variation (function in 'multiColl' package)
  round(VIF(X), digits=4) # VIF (function in 'multiColl' package) 
  round(CN(X), digits=4) # Condition number (function in 'multiColl' package)

#########################################################################################################
# Estimation: Table 2

# OLS: Table 2, column 1

  reg = lm(y ~ X[,-1])
  summary(reg)
  beta = as.matrix(reg$coefficients)
  sigma2 = as.numeric(crossprod(reg$residuals)/reg$df.residual)
  
  MSE = sigma2 * sum(diag(solve(crossprod(X)))) # MSE for MCO
  MSE
  
# Ridge with k proposed by Hoerld and Kennard: Table 3, column 2
  
  # k calculation
  
    LD = eigen(crossprod(X))
    L = LD[[1]]
    D = LD[[2]]
  
    X_star = X%*%D
    initial_alpha = BetaR(y, X_star, 0) # Canonical estimators (function in 'functions.txt')
    alpha = max(initial_alpha)
  
    khk = round(p*sigma2/(alpha^2), digits=4)
    khk
    
  # estimation 
  
    rid(y~x2+x3+x4+x5, k=khk, data=pcd) # function in 'lrmtest'
  
# Ridge with Value of k that minimizes MSE: Table 2, column 3
  
  # k calculation
    
    salto = 0.0001
    tope = 0.002
    k_values = seq(0.0001,tope,salto)
    best_k = k_calculation(y, X[,-1], k_values, intercept=1, DATA=pcd)
    best_k
  
  # estimation
  
    rid(y~x2+x3+x4+x5, best_k, data=pcd)

## Ridge with STANDARDIZED VARIABLES and value that k that minimizes MSE: Table 2, column 4 
  
  # data standardization
    
    x2_est=estandarizar(x2) # function in 'functions.txt'
    x3_est=estandarizar(x3)
    x4_est=estandarizar(x4)
    x5_est=estandarizar(x5)
    X_est=cbind(x2_est,x3_est,x4_est,x5_est)

  # k calculation
    
    salto = 0.1
    tope = 21
    k_values = seq(20,tope,salto)
    best_k = k_calculation(y, X_est, k_values, intercept=0, DATA=pcd)
    best_k
  
  # estimation
  
    rid(y~0+X_est, best_k, data=pcd)

#########################################################################################################
## Limits for FIV/CN, landa min and MSE when raising each variable: Table 3
    
# Limits for FIV
    
  round(VIF(X[,-2]),digits=4) # function in 'multiColl' package
  round(VIF(X[,-3]),digits=4) 
  round(VIF(X[,-4]),digits=4) 
  round(VIF(X[,-5]),digits=4) 

# Limits for FIV
  
  round(CN(X[,-2]), digits=4) # function in 'multiColl' package
  round(CN(X[,-3]), digits=4)
  round(CN(X[,-4]), digits=4)
  round(CN(X[,-5]), digits=4)
  
# landa min and MSE
  
  # raising each variable
  
    landa_min = numeric() 
    MSE_raise = numeric() 
    j = 0  
    for (var_alzada in 2:p){
      j = j + 1
      a_r = auxiliary_regression(X, var_alzada, sigma2, beta) # function in 'functions.txt'
      alfa = a_r[[1]]
      landa_min[j] = a_r[[2]]
      MSE_raise[j] = MSEraise(X, sigma2, beta, landa_min[j], alfa, var_alzada) # function in 'functions.txt'
    }
    round(landa_min,digits=4)
    MSE_raise
    
#########################################################################################################
# Estimation raising each variable: Table 4
    
  vifs = matrix(, p-1, p-1)
  cns = numeric()
  cvs = numeric()
  j = 0  
  for (var_alzada in 2:p){
    j = j + 1
    Xalzada = X_raised(X, var_alzada, landa=landa_min[j]) # function in 'functions.txt'
    vifs[,j] = VIF(Xalzada)
    cns[j] = CN(Xalzada)
    cvs[j] = CV(Xalzada[,var_alzada])
    regresion = lm(y~Xalzada[,-1])
    print(summary(regresion))
  }
  round(vifs,digits=4)
  round(cns,digits=4)
  round(cvs,digits=4)
  
#########################################################################################################  
# landa than makes all VIFs lower than 10: Table 5
  
  # landa calculation 
  
  inicio = 4
  salto = 0.005
  tope = 5
  
  vif_raised3 = AlzadoFIV(X[,-1], alzar=2, inicio, salto, tope) # function in 'functions.txt'
  index = VIF_less10(vif_raised3)
  landa3 = vif_raised3[index,1]
  landa3
  
  vif_raised5 = AlzadoFIV(X[,-1], alzar=4, inicio, salto, tope)
  index = VIF_less10(vif_raised5)
  landa5 = vif_raised5[index,1]
  landa5
  
  # estimation
  
    vifs = matrix(, p-1, 2)
    cns = numeric()
    cvs = numeric()
    mses = numeric()
    landas_min = c(landa3, landa5)
    j = 0  
    for (var_alzada in c(3,5)){
      j = j + 1
      Xalzada = X_raised(X, var_alzada, landa=landas_min[j]) # function in 'functions.txt'
      vifs[,j] = VIF(Xalzada)
      cns[j] = CN(Xalzada)
      cvs[j] = CV(Xalzada[,var_alzada])
      a_r = auxiliary_regression(X, var_alzada, sigma2, beta) # function in 'functions.txt'
      alfa = a_r[[1]]
      mses[j] = MSEraise(X, sigma2, beta, landas_min[j], alfa, var_alzada) # function in 'functions.txt'
      regresion = lm(y~Xalzada[,-1])
      print(summary(regresion))
    }
    round(vifs,digits=4)
    round(cns,digits=4)
          round(cvs,digits=4)
                round(mses,digits=4)
  
#########################################################################################################  
# landa that makes CN<20 for variable 3: Table 6
  
  # landa calculation 
  
  inicio = 15
  salto = 0.01
  tope = 20
  landa_values = seq(inicio,tope,salto)
  
  best_landa = numeric()
  NC_opt = numeric()
  j = 0
  for (var_alzada in c(3,5)){
    j = j + 1
    for (landa in landa_values) {
      Xraised = X_raised(X, var_alzada, landa)
      NC <- CN(Xraised)
      if (NC < 20) {
        best_landa[j] <- landa
        NC_opt[j] = NC
        break
      }
    }  
  }
  best_landa
  NC_opt
  
  # estimation
  
    vifs = matrix(, p-1, 2)
    cns = numeric()
    cvs = numeric()
    mses = numeric()
    j = 0  
    for (var_alzada in c(3,5)){
      j = j + 1
      Xalzada = X_raised(X, var_alzada, landa=best_landa[j]) # function in 'functions.txt'
      vifs[,j] = VIF(Xalzada)
      cns[j] = CN(Xalzada)
      cvs[j] = CV(Xalzada[,var_alzada])
      a_r = auxiliary_regression(X, var_alzada, sigma2, beta) # function in 'functions.txt'
      alfa = a_r[[1]]
      mses[j] = MSEraise(X, sigma2, beta, best_landa[j], alfa, var_alzada) # function in 'functions.txt'
      regresion = lm(y~Xalzada[,-1])
      print(summary(regresion))
    }
    round(vifs,digits=4)
    round(cns,digits=4)
          round(cvs,digits=4)
                round(mses,digits=4)
  
#########################################################################################################  
# Residualization: Table 7
  
  vifs = matrix(, p-1, 2)
  cns = numeric()
  cvs = numeric()
  mses = numeric()
  j = 0   
  for (var_ort in c(3,5)){
    j = j + 1
    X_ort = X
    aux_orth <- lm(X[,var_ort] ~ X[,-var_ort])
    e_auxorth = aux_orth$residuals
    X_ort[,var_ort] = e_auxorth
    vifs[,j] = VIF(X_ort)
    cns[j] = CN(X_ort)
    cvs[j] = CV(X_ort[,var_ort])
    mod_orth <- lm(y ~ X_ort[,-1])
    print(summary(mod_orth))
    mses[j] = sigma2 * sum(diag(solve(crossprod(X_ort))))
  }
  round(vifs,digits=4)
  round(cns,digits=4)
  round(cvs,digits=4)
  round(mses,digits=4)