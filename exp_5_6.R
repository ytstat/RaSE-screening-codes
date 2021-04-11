# This is the R script used to compare the performance of various screening methods in Examples 5 and 6 in [TF21]
# 'v.list' records the order of p features provided by different methods.
# 'v.list[[1]]' and 'v.list[[2]]' represent results of Examples 5 and 6 in [TF21], respectively
# 'timing' records the computational times of different methods in Examples 5 and 6
# [TF21] sets seed = 1:200, then calculate the median MMS to obtain the results in Table 4.
#
# Reference:
# [TF21] Tian, Y. and Feng, Y., 2021. RaSE: A Variable Screening Framework via Random Subspace Ensembles. arXiv preprint arXiv:2102.03892.
#


library(RaSEn)
library(screening)
library(caret)
library(VariableScreening)
library(e1071)
library(dplyr)
library(SIS)
library(energy)
library(doParallel)
library(EDMeasure)




set.seed(seed, kind = "L'Ecuyer-CMRG")
cores <- detectCores()

# -----------------------------------------------------

CIN_KDE = function(x){
  x=as.vector(x)
  n=length(x)
  h=1.06*sd(x)/(n^(0.2))
  phi=D=matrix(0,n,n)
  for(j in 1:(n-1)){
    for(i in (j+1):n){
      a=x[i]-x[j]
      phi[i,j]=exp(-0.5*(a/h)^2)
      D[i,j]=a
    }
  }
  phi=phi+t(phi)+diag(n)
  D=D-t(D)
  Dphi=D*phi
  S=sum((colSums(Dphi)/colSums(phi))^2)
  C=n*(h^4)
  return(S/C)
}

CINscreenDiscrete = function(x1,y1){
  x1=as.vector(x1)
  y1=as.vector(y1)

  x=x1[!is.na(x1)]
  y=y1[!is.na(x1)]

  n=length(y)
  infxy=0
  uy = unique(y)

  for(j in 1:length(uy)){
    xs=x[(y == uy[j])]
    w = length(xs)/n
    infxy=infxy+w*CIN_KDE(xs)
  }

  infx=CIN_KDE(x)
  covxinf=infxy-infx
  covxinf.scaled=infxy/infx #Omitting the addtive constant (-1) since this won't affect the ranking
  output=list(infxy,infx,covxinf,covxinf.scaled,length(uy),n)
  names(output)=c("infxy","infx","covxinf","covxinf.scaled","numslice","n")

  return(output)
}

mat <- rep(list(NULL), 17)
names(mat) <- c("sis", "isis", "holp", "dc_sis", "sirs", "mv_sis", "ip", "mdc", "cin", "knn", "ebic", "bic", "svm", "knn1", "ebic1", "bic1", "svm1")
v.list <- rep(list(mat), 2)

# -------------------------------------------------------
# Example 5:   notes: when seed = 10, mv_sis will pop up feature 549 twice --- error!
# -------------------------------------------------------
n <- 200
p <- 2000
train.data <- RaModel("screening", 5, n = 200, p = 2000)
xtrain <- train.data$x
ytrain <- train.data$y
xtrain <- scale(xtrain)
D <- floor(sqrt(n))
registerDoParallel(cores)



fit.sis <- screening(xtrain, ytrain, method = 'sis', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
fit.isis <- SIS(xtrain, ytrain, family = "binomial", iter.max=2, penalty = "lasso")$ix_list[[2]]
fit.holp <- screening(xtrain, ytrain, method = 'holp', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
fit.dcsis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
fit.mvsis <- screenIID(xtrain, ytrain, method = "MV-SIS")$rank
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
fit.mdc <- order(apply(xtrain, 2, function(x){mdd(ytrain, x)^2}), decreasing = T)
fit.cin <- order(apply(xtrain, 2, function(x){CINscreenDiscrete(x,ytrain)$covxinf.scaled}), decreasing = T)
fit.ebic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "logistic", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "logistic", criterion = "bic", cores = cores, iteration =1)
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "loo", cores = cores, iteration =1, k =5)
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial",cores = cores, iteration =1, cv = 5)


v.list[[1]][['sis']] <- fit.sis
v.list[[1]][['isis']] <- fit.isis
v.list[[1]][['holp']] <- fit.holp
v.list[[1]][['dc_sis']] <- fit.dcsis
v.list[[1]][['sirs']] <- fit.sirs
v.list[[1]][['mv_sis']] <- unique(fit.mvsis)
v.list[[1]][['ip']] <- fit.ip
v.list[[1]][['mdc']] <- fit.mdc
v.list[[1]][['cin']] <- fit.cin
v.list[[1]][['ebic']] <- RaRank(fit.ebic, selected.num = "p", iteration = 0)
v.list[[1]][['bic']] <- RaRank(fit.bic, selected.num = "p", iteration = 0)
v.list[[1]][['knn']] <- RaRank(fit.knn, selected.num = "p", iteration = 0)
v.list[[1]][['svm']] <- RaRank(fit.svm, selected.num = "p", iteration = 0)
v.list[[1]][['ebic1']] <- RaRank(fit.ebic, selected.num = "p", iteration = 1)
v.list[[1]][['bic1']] <- RaRank(fit.bic, selected.num = "p", iteration = 1)
v.list[[1]][['knn1']] <- RaRank(fit.knn, selected.num = "p", iteration = 1)
v.list[[1]][['svm1']] <- RaRank(fit.svm, selected.num = "p", iteration = 1)


# -------------------------------------------------------
# Example 6
# -------------------------------------------------------
n <- 200
p <- 2000
train.data <- RaModel("screening", 6, n = 200, p = 2000)
xtrain <- train.data$x
ytrain <- train.data$y
D <- floor(sqrt(n))
xtrain <- scale(xtrain)

SIS.list <- SIS(x = xtrain, y = ytrain, family = 'multinom', penalty = 'lasso', iter.max = 2)$ix_list
fit.sis <- SIS.list[[1]]
fit.isis <- SIS.list[[2]]
# fit.holp <- NULL
fit.dcsis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
fit.mvsis <- screenIID(xtrain, ytrain, method = "MV-SIS")$rank
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
fit.mdc <- order(apply(xtrain, 2, function(x){mdd(ytrain, x)^2}), decreasing = T)
fit.cin <- order(apply(xtrain, 2, function(x){CINscreenDiscrete(x,ytrain)$covxinf.scaled}), decreasing = T)
fit.ebic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "logistic", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "logistic", criterion = "bic", cores = cores, iteration =1)
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "loo", cores = cores, iteration =1, k =5)
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial",cores = cores, iteration =1, cv = 5)

v.list[[2]][['sis']] <- fit.sis
v.list[[2]][['isis']] <- fit.isis
v.list[[2]][['holp']] <- NULL
v.list[[2]][['dc_sis']] <- fit.dcsis
v.list[[2]][['sirs']] <- fit.sirs
v.list[[2]][['mv_sis']] <- fit.mvsis
v.list[[2]][['ip']] <- fit.ip
v.list[[2]][['mdc']] <- fit.mdc
v.list[[2]][['cin']] <- fit.cin
v.list[[2]][['ebic']] <- RaRank(fit.ebic, selected.num = "p", iteration = 0)
v.list[[2]][['bic']] <- RaRank(fit.bic, selected.num = "p", iteration = 0)
v.list[[2]][['knn']] <- RaRank(fit.knn, selected.num = "p", iteration = 0)
v.list[[2]][['svm']] <- RaRank(fit.svm, selected.num = "p", iteration = 0)
v.list[[2]][['ebic1']] <- RaRank(fit.ebic, selected.num = "p", iteration = 1)
v.list[[2]][['bic1']] <- RaRank(fit.bic, selected.num = "p", iteration = 1)
v.list[[2]][['knn1']] <- RaRank(fit.knn, selected.num = "p", iteration = 1)
v.list[[2]][['svm1']] <- RaRank(fit.svm, selected.num = "p", iteration = 1)

