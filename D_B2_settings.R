# This is the R script used to study the impact of (D, B2) (Figure 3 in [TF21])
# 'mat' and 'mat1' record the MMS of RaSE-BIC and RaSE_1-BIC under different (D, B2) settings, respectively
# [TF21] sets seed = 1:200, then take the median in each (D, B2) settings to plot Figure 3
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
n <- 100
p <- 1000

mat <- matrix(rep(list(NULL), 17*20), nrow = 17)
rownames(mat) <- seq(200, 5000, 300)
colnames(mat) <- seq(2, 40, 2)

mat1 <- matrix(rep(list(NULL), 17*20), nrow = 17)
rownames(mat1) <- seq(200, 5000, 300)
colnames(mat1) <- seq(2, 40, 2)


# -------------------------------------------------------
# Model 1: continuous response, simulation 2 in SIS (2008)
# -------------------------------------------------------

train.data <- RaModel(model.type = "screening", model.no = 1, n = 100, p = 1000)
xtrain <- train.data$x
ytrain <- train.data$y
xtrain <- scale(xtrain)



B2 <- seq(200, 5000, 300)
D <- seq(2, 40, 2)
for (b2 in 1:length(B2)) {
  for (d in 1:length(D)) {
      print(c(b2, d))
      fit.bic <- RaScreen(xtrain, ytrain, B1 = 200, B2 = B2[b2], D = D[d], model = "lm", criterion = "bic", cores = cores, iteration = 1)
      mat[[b2, d]] <- RaRank(fit.bic, selected.num="p", iteration=0)
      mat1[[b2, d]] <- RaRank(fit.bic, selected.num="p", iteration=1)
  }
}

