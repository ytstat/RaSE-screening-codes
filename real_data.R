# This is the R script used to compare the performance of various screening methods in two real studies in [TF21]
# 'v.list' records the order of p features provided by different methods.
# 'v.list[[1]]' and 'v.list[[2]]' represent results of two real studies in [TF21], respectively
# [TF21] sets seed = 1:200, then calculate the median MMS to obtain the results in Table 5.
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
library(class)
library(datamicroarray)
library(glmnet)
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



# -----------------------------------------------------
s <- matrix(nrow = 2, ncol = 20)
colnames(s) <- c("lasso", "sis", "isis", "holp", "dc_sis", "mv_sis", "sirs", "ip", "mdc", "cin",
                 "ebic", "ebic1", "bic", "bic1",
                 "knn.model", "knn", "knn1",
                 "svm.model", "svm", "svm1")


mat <- rep(list(NULL), 17)
names(mat) <- c("sis", "isis", "holp", "dc_sis", "sirs", "mv_sis", "ip",  "mdc", "cin", "knn", "ebic", "bic", "svm", "knn1", "ebic1", "bic1", "svm1")
v.list <- rep(list(mat), 2)

# -----------------------------------------------------
# Example 1: alon
# -----------------------------------------------------
data('alon', package = 'datamicroarray')
x <- alon$x
y <- as.numeric(alon$y)-1


test_ind <- sample(1:length(y), length(y)/10)
xtrain <- x[-test_ind, ]
ytrain <- y[-test_ind]
xtest <- x[test_ind, ]
ytest <- y[test_ind]
scale0 <- scale(xtrain)
xtrain <- scale0
xtest <- scale(xtest, center = attributes(scale0)$`scaled:center`, scale = attributes(scale0)$`scaled:scale`)
num_va <- floor(length(ytrain)/log(length(ytrain)))
n <- length(ytrain)
p <- ncol(xtrain)
D <- floor(sqrt(n))

er <- numeric(20)
names(er) <- c("lasso", "sis", "isis", "holp", "dc_sis", "mv_sis", "sirs", "ip", "mdc", "cin",
               "ebic", "ebic1", "bic", "bic1",
               "knn.model", "knn", "knn1",
               "svm.model", "svm", "svm1")

registerDoParallel(cores)

fit.sis <- screening(xtrain, ytrain, method = "sis", family = "binomial", num.select = ncol(xtrain))$screen
fit.isis <- SIS(xtrain, ytrain, family = "binomial", iter.max=2, penalty = "lasso")$ix_list[[2]]
fit.holp <- screening(xtrain, ytrain, method = "holp", family = "binomial", num.select = ncol(xtrain))$screen
fit.dc_sis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
fit.mv_sis <- screenIID(xtrain, ytrain, method = "MV-SIS")$rank
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
fit.mdc <- order(apply(xtrain, 2, function(x){mdd(ytrain, x)^2}), decreasing = T)
fit.cin <- order(apply(xtrain, 2, function(x){CINscreenDiscrete(x,ytrain)$covxinf.scaled}), decreasing = T)
fit.ebic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "logistic", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "logistic", criterion = "bic", cores = cores, iteration =1)
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "loo", cores = cores, iteration =1, k =5)
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial", cores = cores, iteration =1, cv = 5)


lambda.cv <- cv.glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial")
lasso <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.sis <- fit.sis[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.sis]), y = ytrain, family = "binomial")
lasso.sis <- glmnet(x = as.matrix(xtrain[, va.sis]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)

va.isis <- fit.isis[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.isis], ytrain, family = "binomial")
lasso.isis <- glmnet(as.matrix(xtrain)[, va.isis], ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.holp <- fit.holp[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.holp]), y = ytrain, family = "binomial")
lasso.holp <- glmnet(x = as.matrix(xtrain[, va.holp]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.mv_sis <- fit.mv_sis[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.mv_sis]), y = ytrain, family = "binomial")
lasso.mv_sis <- glmnet(x = as.matrix(xtrain[, va.mv_sis]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)

va.dc_sis <- fit.dc_sis[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.dc_sis]), y = ytrain, family = "binomial")
lasso.dc_sis <- glmnet(x = as.matrix(xtrain[, va.dc_sis]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.sirs<- fit.sirs[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.sirs]), y = ytrain, family = "binomial")
lasso.sirs <- glmnet(x = as.matrix(xtrain[, va.sirs]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)

va.ip<- fit.ip[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.ip]), y = ytrain, family = "binomial")
lasso.ip <- glmnet(x = as.matrix(xtrain[, va.ip]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.mdc<- fit.mdc[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.mdc]), y = ytrain, family = "binomial")
lasso.mdc <- glmnet(x = as.matrix(xtrain[, va.mdc]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.cin<- fit.cin[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.cin]), y = ytrain, family = "binomial")
lasso.cin <- glmnet(x = as.matrix(xtrain[, va.cin]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.ebic <- RaRank(fit.ebic, selected.num="p", iteration=0)[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.ebic]), y = ytrain, family = "binomial")
lasso.ebic <- glmnet(x = as.matrix(xtrain[, va.ebic]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)

va.ebic1 <- va.ebic <- RaRank(fit.ebic, selected.num="p", iteration=1)[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.ebic1]), y = ytrain, family = "binomial")
lasso.ebic1 <- glmnet(x = as.matrix(xtrain[, va.ebic1]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.bic <- RaRank(fit.bic, selected.num="p", iteration=0)[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.bic]), y = ytrain, family = "binomial")
lasso.bic <- glmnet(x = as.matrix(xtrain[, va.bic]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


va.bic1 <- RaRank(fit.bic, selected.num="p", iteration=1)[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.bic1]), y = ytrain, family = "binomial")
lasso.bic1 <- glmnet(x = as.matrix(xtrain[, va.bic1]), y = ytrain, family = "binomial", lambda = lambda.cv$lambda.min)


knn.model <- knn(train = xtrain, test = xtest, cl = ytrain, k = 5)

va.knn <- RaRank(fit.knn, selected.num="p", iteration=0)[1:num_va]
knn.knn <- knn(train = xtrain[, va.knn], test = xtest[, va.knn], cl = ytrain, k = 5)

va.knn1 <- RaRank(fit.knn, selected.num="p", iteration=1)[1:num_va]
knn.knn1 <- knn(train = xtrain[, va.knn1], test = xtest[, va.knn1], cl = ytrain, k = 5)

svm.model <- svm(x = xtrain, y = ytrain, type = "C-classification", kernel = "radial")

va.svm <- RaRank(fit.svm, selected.num="p", iteration=0)[1:num_va]
svm.svm <- svm(x = xtrain[, va.svm], y = ytrain, type = "C-classification", kernel = "radial")

va.svm1 <- RaRank(fit.svm, selected.num="p", iteration=1)[1:num_va]
svm.svm1 <- svm(x = xtrain[, va.svm1], y = ytrain, type = "C-classification", kernel = "radial")

er['lasso'] <- mean(as.numeric(predict(lasso, as.matrix(xtest), type = "class")) != ytest)


er['sis'] <- mean(as.numeric(predict(lasso.sis, as.matrix(xtest[, va.sis]), type = "class")) != ytest)
er['isis'] <- mean(as.numeric(predict(lasso.isis, as.matrix(xtest[, va.isis]), type = "class")) != ytest)
er['holp'] <- mean(as.numeric(predict(lasso.holp, as.matrix(xtest[, va.holp]), type = "class")) != ytest)
er['dc_sis'] <- mean(as.numeric(predict(lasso.dc_sis, as.matrix(xtest[, va.dc_sis]), type = "class")) != ytest)
er['sirs'] <- mean(as.numeric(predict(lasso.sirs, as.matrix(xtest[, va.sirs]), type = "class")) != ytest)
er['mv_sis'] <- mean(as.numeric(predict(lasso.mv_sis, as.matrix(xtest[, va.mv_sis]), type = "class")) != ytest)
er['ip'] <- mean(as.numeric(predict(lasso.ip, as.matrix(xtest[, va.ip]), type = "class")) != ytest)
er['mdc'] <- mean(as.numeric(predict(lasso.mdc, as.matrix(xtest[, va.mdc]), type = "class")) != ytest)
er['cin'] <- mean(as.numeric(predict(lasso.cin, as.matrix(xtest[, va.cin]), type = "class")) != ytest)


er['ebic'] <- mean(as.numeric(predict(lasso.ebic, as.matrix(xtest[, va.ebic]), type = "class")) != ytest)
er['ebic1'] <- mean(as.numeric(predict(lasso.ebic1, as.matrix(xtest[, va.ebic1]), type = "class")) != ytest)
er['bic'] <- mean(as.numeric(predict(lasso.bic, as.matrix(xtest[, va.bic]), type = "class")) != ytest)
er['bic1'] <- mean(as.numeric(predict(lasso.bic1, as.matrix(xtest[, va.bic1]), type = "class")) != ytest)

er['knn.model'] <- mean(knn.model != ytest)
er['knn'] <- mean(knn.knn != ytest)
er['knn1'] <- mean(knn.knn1 != ytest)

er['svm.model'] <- mean(predict(svm.model, xtest) != ytest)
er['svm'] <- mean(predict(svm.svm, xtest[, va.svm]) != ytest)
er['svm1'] <- mean(predict(svm.svm1, xtest[, va.svm1]) != ytest)

s[1, ] <- er




v.list[[1]][['sis']] <- fit.sis
v.list[[1]][['isis']] <- fit.isis
v.list[[1]][['holp']] <- fit.holp
v.list[[1]][['dc_sis']] <- fit.dc_sis
v.list[[1]][['sirs']] <- fit.sirs
v.list[[1]][['mv_sis']] <- fit.mv_sis
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


# -----------------------------------------------------
# Example 2: expression
# -----------------------------------------------------
# D <- read.table("/Users/yetian/Desktop/Dropbox/Columbia/Research/Project/RaSE Screening/dataset/RatEyeExpression.txt", sep ="\t")
D <- read.table("/home/yt2170/work/screening/datasets/RatEyeExpression.txt", sep ="\t")
D <- t(D)
colnames(D) <- D[1,]
D <- D[-1,]
y <- as.numeric(D[,'1389163_at'])
D <- D[, which(colnames(D) != "1389163_at")]
x <- matrix(as.numeric(D), nrow = nrow(D))

x_var <- apply(x, 2, var)
va <- order(x_var, decreasing = T)[1:5000]
x <- x[, va]



test_ind <- sample(1:length(y), length(y)/10)
xtrain <- x[-test_ind, ]
ytrain <- y[-test_ind]
xtest <- x[test_ind, ]
ytest <- y[test_ind]
scale0 <- scale(xtrain)
xtrain <- scale0
xtest <- scale(xtest, center = attributes(scale0)$`scaled:center`, scale = attributes(scale0)$`scaled:scale`)
num_va <- floor(length(ytrain)/log(length(ytrain)))
n <- length(ytrain)
p <- ncol(xtrain)
D <- sqrt(n)

er <- numeric(20)
names(er) <- c("lasso", "sis", "isis", "holp", "dc_sis", "mv_sis", "sirs", "ip", "mdc", "cin",
               "ebic", "ebic1", "bic", "bic1",
               "knn.model", "knn", "knn1",
               "svm.model", "svm", "svm1")

fit.sis <- screening(xtrain, ytrain, method = "sis", family = "gaussian", num.select = ncol(xtrain))$screen
fit.isis <- SIS(xtrain, ytrain, family = "gaussian", iter.max = 2, penalty = "lasso")$ix_list[[2]]
fit.holp <- screening(xtrain, ytrain, method = "holp", family = "gaussian", num.select = ncol(xtrain))$screen
fit.dc_sis <- screenIID(xtrain, ytrain, method = "DC-SIS")$rank
fit.mv_sis <- NULL
fit.sirs <- screenIID(xtrain, ytrain, method = "SIRS")$rank
fit.ip <- order(foreach(i = 1:ncol(xtrain), .combine = "c") %dopar% {dcor(xtrain[,i]^2, ytrain^2)}, decreasing = T)
fit.mdc <- order(apply(xtrain, 2, function(x){mdd(ytrain, x)^2}), decreasing = T)
fit.cin <- order(apply(xtrain, 2, function(x){CINscreen(x,ytrain)$covxinf.scaled}), decreasing = T)
fit.ebic <- RaScreen(xtrain,ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "ebic", cores = cores, iteration =1, gam = 0.5)
fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "lm", criterion = "bic", cores = cores, iteration =1)
fit.knn <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "knn", criterion = "cv", cores = cores, cv = 5, iteration =1, k =5)
fit.svm <- RaScreen(xtrain, ytrain, B1 = 200, B2 = 20*floor(p/D), model = "svm", criterion = "cv", kernel = "radial", cores = cores, iteration =1, cv = 5)



lambda.cv <- cv.glmnet(as.matrix(xtrain), ytrain)
lasso <- glmnet(as.matrix(xtrain), ytrain, lambda = lambda.cv$lambda.min)

va.sis <- fit.sis[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.sis], ytrain)
lasso.sis <- glmnet(as.matrix(xtrain)[, va.sis], ytrain, lambda = lambda.cv$lambda.min)


va.isis <- fit.isis[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.isis], ytrain)
lasso.isis <- glmnet(as.matrix(xtrain)[, va.isis], ytrain, lambda = lambda.cv$lambda.min)

va.holp <- fit.holp[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.holp], ytrain)
lasso.holp <- glmnet(as.matrix(xtrain)[, va.holp], ytrain, lambda = lambda.cv$lambda.min)

va.sirs <- fit.sirs[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.sirs], ytrain)
lasso.sirs <- glmnet(as.matrix(xtrain)[, va.sirs], ytrain, lambda = lambda.cv$lambda.min)

va.dc_sis <- fit.dc_sis[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.dc_sis], ytrain)
lasso.dc_sis <- glmnet(as.matrix(xtrain)[, va.dc_sis], ytrain, lambda = lambda.cv$lambda.min)


va.ip<- fit.ip[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.ip], ytrain)
lasso.ip <- glmnet(as.matrix(xtrain)[, va.ip], ytrain, lambda = lambda.cv$lambda.min)

va.mdc<- fit.mdc[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.mdc]), y = ytrain)
lasso.mdc <- glmnet(x = as.matrix(xtrain[, va.mdc]), y = ytrain, lambda = lambda.cv$lambda.min)


va.cin<- fit.cin[1:num_va]
lambda.cv <- cv.glmnet(x = as.matrix(xtrain[, va.cin]), y = ytrain)
lasso.cin <- glmnet(x = as.matrix(xtrain[, va.cin]), y = ytrain, lambda = lambda.cv$lambda.min)


va.ebic <- RaRank(fit.ebic, selected.num="p", iteration=0)[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.ebic], ytrain)
lasso.ebic <- glmnet(as.matrix(xtrain)[, va.ebic], ytrain, lambda = lambda.cv$lambda.min)

va.ebic1 <- RaRank(fit.ebic, selected.num="p", iteration=1)[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.ebic1], ytrain)
lasso.ebic1 <- glmnet(as.matrix(xtrain)[, va.ebic1], ytrain, lambda = lambda.cv$lambda.min)

va.bic <- RaRank(fit.bic, selected.num="p", iteration=0)[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.bic], ytrain)
lasso.bic <- glmnet(as.matrix(xtrain)[, va.bic], ytrain, lambda = lambda.cv$lambda.min)

va.bic1 <- RaRank(fit.bic, selected.num="p", iteration=1)[1:num_va]
lambda.cv <- cv.glmnet(as.matrix(xtrain)[, va.bic1], ytrain)
lasso.bic1 <- glmnet(as.matrix(xtrain)[, va.bic1], ytrain, lambda = lambda.cv$lambda.min)


knn.model <- knnreg(x = xtrain, y = ytrain, k = 5)

va.knn <- RaRank(fit.knn, selected.num="p", iteration=0)[1:num_va]
knn.knn <- knnreg(x = xtrain[, va.knn], y = ytrain, k = 5)

va.knn1 <- RaRank(fit.knn, selected.num="p", iteration=1)[1:num_va]
knn.knn1 <- knnreg(x = xtrain[, va.knn1], y = ytrain, k = 5)

svm.model <- svm(x = xtrain, y = ytrain, type = "eps-regression", kernel = "radial")

va.svm <- RaRank(fit.svm, selected.num="p", iteration=0)[1:num_va]
svm.svm <- svm(x = xtrain[, va.svm], y = ytrain, type = "eps-regression", kernel = "radial")

va.svm1 <- RaRank(fit.svm, selected.num="p", iteration=1)[1:num_va]
svm.svm1 <- svm(x = xtrain[, va.svm1], y = ytrain, type = "eps-regression", kernel = "radial")



er['sis'] <- mean((predict(lasso.sis, as.matrix(xtest)[, va.sis]) - ytest)^2)
er['isis'] <- mean((predict(lasso.isis, as.matrix(xtest)[, va.isis]) - ytest)^2)
er['holp'] <- mean((predict(lasso.holp, as.matrix(xtest)[, va.holp]) - ytest)^2)
er['sirs'] <- mean((predict(lasso.sirs, as.matrix(xtest)[, va.sirs]) - ytest)^2)
er['dc_sis'] <- mean((predict(lasso.dc_sis, as.matrix(xtest)[, va.dc_sis]) - ytest)^2)
er['ip'] <- mean((predict(lasso.ip, as.matrix(xtest)[, va.ip]) - ytest)^2)

er['mdc'] <- mean((predict(lasso.mdc, as.matrix(xtest)[, va.mdc]) - ytest)^2)
er['cin'] <- mean((predict(lasso.cin, as.matrix(xtest)[, va.cin]) - ytest)^2)


er['mv_sis'] <- NA

er['lasso'] <- mean((predict(lasso, as.matrix(xtest)) - ytest)^2)
er['ebic'] <- mean((predict(lasso.ebic, as.matrix(xtest)[, va.ebic]) - ytest)^2)
er['ebic1'] <- mean((predict(lasso.ebic1, as.matrix(xtest)[, va.ebic1]) - ytest)^2)
er['bic'] <- mean((predict(lasso.bic, as.matrix(xtest)[, va.bic]) - ytest)^2)
er['bic1'] <- mean((predict(lasso.bic1, as.matrix(xtest)[, va.bic1]) - ytest)^2)

er['knn.model'] <- mean((predict(knn.model, xtest) - ytest)^2)
er['knn'] <- mean((predict(knn.knn, xtest[, va.knn]) - ytest)^2)
er['knn1'] <- mean((predict(knn.knn1, xtest[, va.knn1]) - ytest)^2)

er['svm.model'] <- mean((predict(svm.model, xtest) - ytest)^2)
er['svm'] <- mean((predict(svm.svm, xtest[, va.svm]) - ytest)^2)
er['svm1'] <- mean((predict(svm.svm1, xtest[, va.svm1]) - ytest)^2)

s[2, ] <- er




v.list[[2]][['sis']] <- fit.sis
v.list[[2]][['isis']] <- fit.isis
v.list[[2]][['holp']] <- fit.holp
v.list[[2]][['dc_sis']] <- fit.dc_sis
v.list[[2]][['sirs']] <- fit.sirs
v.list[[2]][['mv_sis']] <- fit.mv_sis
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







