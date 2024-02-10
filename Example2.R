rm(list=ls())

library(latex2exp) # install.packages("latex2exp")
library(multiColl) # install.packages("multiColl")

source("functions.txt")

##################################################################
# data

Y<- c(3.626,3.899,3.447,3.011,2.543,2.234,2.201,2.364,2.140,2.290,2.350,
      2.320,2.320,2.190,2.200,2.630,2.950,3.310,3.620,3.600,4.088,4.381,
      4.650,4.677,4.476,5.054,5.368,4.339,2.218,1.673,1.337,1.239,1.224,
      1.253,1.405,1.522,1.736,2.128,2.114,2.056,1.671,1.281,0.897,0.600,
      0.572,0.5062,0.537)
x2<- c(92.9,93.9,93.9,94.4,95.1,95.7,95.9,96.4,96.8,98.0,98.1,98.7,98.8,
       100.0,100.3,101.0,101.1,102.4,102.5,102.8,103.0,104.4,104.5,105.8,
       106.4,108.2,108.5,108.2,107.5,108.4,108.1,108.7,108.7,110.1,110.0,
       110.9,111.4,113.2,112.9,114.1,114.4,115.9,115.8,116.8,116.5,117.6,117.3)
x3<- c(17211.0,2724.0,17232.0,9577.0,4117.0,-2134.0,6117.0,10949.0,18360.0,
       13646.0,8424.0,14319.0,3885.0,4493.0,-320.0,-2736.0,-6909.0,-4848.0,
       -4255.0,1347.0,8781.0,8723.0,3662.0,-17548.0,-37041.0,-27624.0,-37723.0,
       -43584.0,-16070.0,-5029.0,7294.0,85.0,-4399.0,-2431.0,2137.0,-4345.0,
       -12643.0,-2272.0,-3592.0,8071.0,12202.0,35619.0,42161.0,43880.0,52483.0,
       56376.0,48981.0)
x4<- c(-51384.0,-49567.1,-52128.4,-53593.3,-65480.0,-50343.8,-75646.4,-59120.8,
       -69246.3,-60313.8,-56782.9,-55313.1,-67034.4,-61942.8,-46258.4,-43761.4,
       -37562.6,-35609.6,-27064.0,-32497.2,-18389.0,-9923.5,-9727.0,-23729.9,
       -28909.3,-46527.0,-49654.0,-81729.7,-121227.5,-142580.9,-164699.2,-152269.2,
       -162477.4,-128366.4,-169848.0,-129290.2,-104646.7,-103143.8,-102621.8,
       -104240.4,-82309.3,-91620.9,-85054.4,-99998.2,-81287.1,-77738.8,-73003.3)

n <- length(Y)
cte = rep(1,n)
X = cbind(cte,x2,x3,x4)
p = ncol(X)

##################################################################
# Matrix of correlations: Table 8

RdetR(X) # function in 'multiColl' package
round(VIF(X),digits=4) # function in 'multiColl' package
round(CN(X),digits=4) # function in 'multiColl' package
round(CVs(X),digits=4) # function in 'multiColl' package

##################################################################
# Angles: Table 9

round(angle(cte,x2),digits=4) # function in 'functions.txt'
round(angle(cte,x3),digits=4)
round(angle(cte,x4),digits=4)
round(angle(x2,x3),digits=4)
round(angle(x2,x4),digits=4)
round(angle(x3,x4),digits=4)

#################################################################
# Table 10

var_alzada = 2

# OLS estimation: first column

  regresion <- lm(Y ~ X[,-1])
  summary(regresion)
  beta = coef(lm(Y ~ X[,-1]))
  
  sigma2 = as.numeric(summary(regresion)[6])^2
  MSE = sigma2 * sum(diag(solve(crossprod(X)))) # MSE for OLS
  MSE

  angle(cte,x2)

# landa calculation
  
  inicio = 0
  salto = 0.01
  tope = 40
  landa_values = seq(inicio,tope,salto)
  
  best_landa_cv <- numeric()
  best_landa_nc <- numeric()
  best_landa_norma <- numeric()
  
  normas <- numeric()
  diferencias_normas <- numeric()

  j1 = 0
  j2 = 0
  j3 = 0
  j4 = 0
  for (landa in landa_values) {
    X_alzada = X_raised(X, var_alzada, landa)
    #
    x_alzada = X_alzada[, var_alzada]
    cv = CV(x_alzada)
    if (cv > 0.1002506) {
      j1 = j1 + 1
      best_landa_cv[j1] <- landa # landa that makes cv<0.1002506: second column (RE-CV)
    }
    #
    nc <- CN(X_alzada)
    if (nc < 20) {
      j2 = j2 + 1
      best_landa_nc[j2] <- landa # landa that makes CN<20: third column (RE-CN)
    }
    #
    regresion_norma <- lm(Y ~ X_alzada[,-1])
    beta_norma <- as.matrix(regresion_norma$coef)
    normaF <- norm(beta_norma, "F")^2
    j3 = j3 + 1
    if (j3 == 1) normas[j3] = normaF
    if (j3 > 1) {
      normas[j3] <- normaF
      diferencias_normas[j3-1] <- normas[j3] - normas[j3-1]
      if (abs(diferencias_normas[j3-1])/0.1 < 0.001) {
        j4 = j4 + 1
        best_landa_norma[j4] <- landa # landa that stabilize norm: fourth column (RE-STABLE)
      }
    }  
  }
  
  # landa that minimizes MSE: fith column (RE-MSE)
  a_r = auxiliary_regression(X, var_alzada, sigma2, beta) # function in 'functions.txt'
  landa_minMSE = a_r[[2]]
  
  landas = c(0, min(best_landa_cv), min(best_landa_nc), min(best_landa_norma), landa_minMSE)
  landas
  
# Figures
  
  # Figure 1
    plot(landa_values, normas, xlab=TeX('$\\lambda$'), ylab=TeX('$||\\hat{\\beta}(\\lambda)||$'))
    
  # Figure 1
    plot(landa_values[-1], diferencias_normas/0.1, xlab=TeX('$\\lambda$'), ylab=TeX('$(||\\hat{\\beta}(\\lambda)||-||\\hat{\\beta}(\\lambda-0.1)||)/0.1$'))

# estimation

  vifs = matrix(, p-1, length(landas))
  cns = numeric()
  cvs = numeric()
  normas = numeric()
  angles = numeric()
  mses = numeric()
  j = 0
  for (best_landa in landas){
    j = j + 1
    X_alzada = X_raised(X, var_alzada, best_landa)
    regresion = lm(Y~X_alzada[,-1])
    print(summary(regresion))
    #
    beta_alzado = as.matrix(regresion$coef)
    normas[j] = norm(beta_alzado,"F")^2
    #
    vifs[,j] = VIF(X_alzada) # function in 'multiColl' package
    cns[j] = CN(X_alzada) # function in 'multiColl' package
    x_alzada = X_alzada[, var_alzada]
    cvs[j] = CV(x_alzada) # function in 'multiColl' package
    #
    angles[j] = angle(cte,x_alzada) 
    #
    a_r = auxiliary_regression(X, var_alzada, sigma2, beta) # function in 'functions.txt'
    alfa = a_r[[1]]
    mses[j] = MSEraise(X, sigma2, beta, best_landa, alfa, var_alzada) # function in 'functions.txt'
  }
round(landas,digits=4)
round(vifs,digits=4)
round(cns,digits=4)
round(mses,digits=4)
round(normas,digits=4)
round(cvs,digits=4)
round(angles,digits=4)
  
  