# This is the R script used to compare the performance of various screening methods in Examples 1-3 in [TF21]
# 'v.list' records the order of p features provided by different methods.
# 'v.list[[1]]', 'v.list[[2]]' and 'v.list[[3]]' represent results of Examples 1, 2 and 3 in [TF21], respectively
# 'timing' records the computational times of different methods in Example 1
# [TF21] sets seed = 1:200, then calculate the median MMS to obtain the results in Tables 1 and 3 and the average running time to generate Table 2.
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

CINscreen = function(x1,y1,numslice=5){
  x1=as.vector(x1)
  y1=as.vector(y1)
  # if(length(x1)!=length(y1)){
  #   print("X and Y are of different length!")
  # }
  #else{
  #xNA=which(is.na(x))
  #yNA=which(is.na(y))
  #xyNA=union(xNA,yNA)
  x=x1[!is.na(x1)]
  y=y1[!is.na(x1)]
  n=length(y)
  ind=order(y)
  numrep=round(n/numslice)
  infxy=0
  for(i in 1:(numslice-1)){
    xs=x[ind[((i-1)*numrep+1):(i*numrep)]]
    infxy=infxy+(numrep/n)*CIN_KDE(xs)
  }
  xs=x[ind[((numslice-1)*numrep+1):n]]
  infxy=infxy+((n-(numslice-1)*numrep)/n)*CIN_KDE(xs)
  infx=CIN_KDE(x)
  covxinf=infxy-infx
  covxinf.scaled=infxy/infx #Omitting the addtive constant (-1) since this won't affect the ranking
  output=list(infxy,infx,covxinf,covxinf.scaled,numslice,n)
  names(output)=c("infxy","infx","covxinf","covxinf.scaled","numslice","n")
  return(output)
  #}
}


mat <- rep(list(NULL), 17)
names(mat) <- c("sis", "isis", "holp", "dc_sis", "sirs", "mv_sis", "ip", "mdc", "cis", "knn", "ebic", "bic", "svm", "knn1", "ebic1", "bic1", "svm1")
v.list <- rep(list(mat), 3)


# -------------------------------------------------------
# Example 1: continuous response, simulation 2 in SIS (2008)
# -------------------------------------------------------
timing <- rep(NA, 16)
names(timing) <- c("sis", "isis", "holp", "dc_sis", "sirs", "ip", "mdc", "cis", "knn", "ebic", "bic", "svm", "knn1", "ebic1", "bic1", "svm1")

n <- 100
p <- 1000
train.data <- RaModel(model.no = 1, n = 100, p = 1000, model.type = "screening")
xtrain <- train.data$x
ytrain <- train.data$y
xtrain <- scale(xtrain)
D <- floor(sqrt(n))


registerDoParallel(cores)

t1 <- proc.time()
fit.sis <- screening(xtrain, ytrain, method = 'sis', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
timing['sis'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.isis <- SIS(xtrain, ytrain, iter.max=2, penalty = "lasso")$ix_list[[2]]
timing['isis'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.holp <- screening(xtrain, ytrain, method = 'holp', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
timing['holp'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
timing['sirs'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.dcsis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
timing['dc_sis'] <- (proc.time() - t1)['elapsed']


t1 <- proc.time()
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
timing['ip'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.mdc <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {mdd(ytrain, xtrain[, i])^2}, decreasing = T)
timing['mdc'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.cis <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {CINscreen(xtrain[, i], ytrain)$covxinf.scaled}, decreasing = T)
timing['cis'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.ebic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
timing['ebic1'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "bic", cores = cores, iteration =1)
timing['bic1'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "cv", cores = cores, cv = 5,iteration =1, k =5)
timing['knn1'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial",cores = cores, iteration =1, cv = 5)
timing['svm1'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "ebic", cores = cores, iteration =0, gam = 0.5)
timing['ebic'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "bic", cores = cores, iteration =0)
timing['bic'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "cv", cores = cores, cv = 5,iteration =0, k =5)
timing['knn'] <- (proc.time() - t1)['elapsed']

t1 <- proc.time()
fit <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial",cores = cores, iteration =0, cv = 5)
timing['svm'] <- (proc.time() - t1)['elapsed']


v.list[[1]][['sis']] <- fit.sis
v.list[[1]][['isis']] <- fit.isis
v.list[[1]][['holp']] <- fit.holp
v.list[[1]][['dc_sis']] <- fit.dcsis
v.list[[1]][['sirs']] <- fit.sirs
v.list[[1]][['mv_sis']] <- NULL
v.list[[1]][['ip']] <- fit.ip
v.list[[1]][['mdc']] <- fit.mdc
v.list[[1]][['cis']] <- fit.cis
v.list[[1]][['ebic']] <- RaRank(fit.ebic, selected.num="p", iteration=0)
v.list[[1]][['bic']] <- RaRank(fit.bic, selected.num="p", iteration=0)
v.list[[1]][['knn']] <- RaRank(fit.knn, selected.num="p", iteration=0)
v.list[[1]][['svm']] <- RaRank(fit.svm, selected.num="p", iteration=0)
v.list[[1]][['ebic1']] <- RaRank(fit.ebic, selected.num="p", iteration=1)
v.list[[1]][['bic1']] <- RaRank(fit.bic, selected.num="p", iteration=1)
v.list[[1]][['knn1']] <- RaRank(fit.knn, selected.num="p", iteration=1)
v.list[[1]][['svm1']] <- RaRank(fit.svm, selected.num="p", iteration=1)



# -------------------------------------------------------
# Example 2: continuous response, knn example
# -------------------------------------------------------
n <- 200
p <- 2000
train.data <- RaModel(model.type = "screening", model.no = 2, n = 200, p = 2000)
xtrain <- train.data$x
ytrain <- train.data$y
xtrain <- scale(xtrain)
D <- floor(sqrt(n))


fit.sis <- screening(xtrain, ytrain, method = 'sis', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
fit.isis <- SIS(xtrain, ytrain, iter.max=2, penalty = "lasso")$ix_list[[2]]
fit.holp <- screening(xtrain, ytrain, method = 'holp', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
fit.dcsis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
fit.mdc <- order(apply(xtrain, 2, function(x){mdd(ytrain, x)^2}), decreasing = T)
fit.cis <- order(apply(xtrain, 2, function(x){CINscreen(x,ytrain)$covxinf.scaled}), decreasing = T)
fit.ebic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "bic", cores = cores, iteration =1)
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "cv", cores = cores, cv = 5,iteration =1, k =5)
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial",cores = cores, iteration =1, cv = 5)


v.list[[2]][['sis']] <- fit.sis
v.list[[2]][['isis']] <- fit.isis
v.list[[2]][['holp']] <- fit.holp
v.list[[2]][['dc_sis']] <- fit.dcsis
v.list[[2]][['sirs']] <- fit.sirs
v.list[[2]][['mv_sis']] <- NULL
v.list[[2]][['ip']] <- fit.ip
v.list[[2]][['mdc']] <- fit.mdc
v.list[[2]][['cis']] <- fit.cis
v.list[[2]][['ebic']] <- RaRank(fit.ebic, selected.num="p", iteration=0)
v.list[[2]][['bic']] <- RaRank(fit.bic, selected.num="p", iteration=0)
v.list[[2]][['knn']] <- RaRank(fit.knn, selected.num="p", iteration=0)
v.list[[2]][['svm']] <- RaRank(fit.svm, selected.num="p", iteration=0)
v.list[[2]][['ebic1']] <- RaRank(fit.ebic, selected.num="p", iteration=1)
v.list[[2]][['bic1']] <- RaRank(fit.bic, selected.num="p", iteration=1)
v.list[[2]][['knn1']] <- RaRank(fit.knn, selected.num="p", iteration=1)
v.list[[2]][['svm1']] <- RaRank(fit.svm, selected.num="p", iteration=1)


# -------------------------------------------------------
# Example 3: continuous response, example 1.c in DC-SIS (2012)
# -------------------------------------------------------
n <- 200
p <- 2000
train.data <- RaModel(model.type = "screening", model.no = 3, n = 200, p = 2000)
xtrain <- train.data$x
ytrain <- train.data$y
xtrain <- scale(xtrain)
D <- floor(sqrt(n))

fit.sis <- screening(xtrain, ytrain, method = 'sis', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
fit.isis <- SIS(xtrain, ytrain, iter.max=2, penalty = "lasso")$ix_list[[2]]
fit.holp <- screening(xtrain, ytrain, method = 'holp', ebic.gamma = 0.5, num.select = ncol(xtrain))$screen
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
fit.dcsis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
fit.mdc <- order(apply(xtrain, 2, function(x){mdd(ytrain, x)^2}), decreasing = T)
fit.cis <- order(apply(xtrain, 2, function(x){CINscreen(x,ytrain)$covxinf.scaled}), decreasing = T)
fit.ebic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "bic", cores = cores, iteration =1)
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "cv", cores = cores, cv = 5, iteration =1, k =5)
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial",cores = cores, iteration =1, cv = 5)


v.list[[3]][['sis']] <- fit.sis
v.list[[3]][['isis']] <- fit.isis
v.list[[3]][['holp']] <- fit.holp
v.list[[3]][['dc_sis']] <- fit.dcsis
v.list[[3]][['sirs']] <- fit.sirs
v.list[[3]][['mv_sis']] <- NULL
v.list[[3]][['ip']] <- fit.ip
v.list[[3]][['mdc']] <- fit.mdc
v.list[[3]][['cis']] <- fit.cis
v.list[[3]][['ebic']] <- RaRank(fit.ebic, selected.num="p", iteration=0)
v.list[[3]][['bic']] <- RaRank(fit.bic, selected.num="p", iteration=0)
v.list[[3]][['knn']] <- RaRank(fit.knn, selected.num="p", iteration=0)
v.list[[3]][['svm']] <- RaRank(fit.svm, selected.num="p", iteration=0)
v.list[[3]][['ebic1']] <- RaRank(fit.ebic, selected.num="p", iteration=1)
v.list[[3]][['bic1']] <- RaRank(fit.bic, selected.num="p", iteration=1)
v.list[[3]][['knn1']] <- RaRank(fit.knn, selected.num="p", iteration=1)
v.list[[3]][['svm1']] <- RaRank(fit.svm, selected.num="p", iteration=1)


