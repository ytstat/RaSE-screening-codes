# This is the R script used to study the impact of (B1, B2) (Figure 1 in [TF21])
# 'mat' and 'mat1' record the MMS of RaSE-BIC and RaSE_1-BIC under different (B1, B2) settings, respectively
# [TF21] sets seed = 1:200, then take the median in each (B1, B2) settings to plot Figure 1
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

set.seed(seed, kind = "L'Ecuyer-CMRG")
cores <- detectCores()

# -----------------------------------------------------

mat <- matrix(rep(list(NULL), 10*17), nrow = 10)
rownames(mat) <- seq(100, 1000, 100)
colnames(mat) <- seq(1000, 100000, 6000)

mat1 <- matrix(rep(list(NULL), 10*17), nrow = 10)
rownames(mat1) <- seq(100, 1000, 100)
colnames(mat1) <- seq(1000, 100000, 6000)



# -------------------------------------------------------
# Model 1: continuous response, simulation 2 in SIS (2008)
# -------------------------------------------------------
n <- 100
p <- 1000
train.data <- RaModel(model.type = "screening", model.no = 1, n = 100, p = 1000)
xtrain <- train.data$x
ytrain <- train.data$y
xtrain <- scale(xtrain)



B1 <- seq(100, 1000, 100)
B2 <- seq(1000, 100000, 6000)
for (b1 in 1:10) {
  for (b2 in 1:17) {
      print(c(b1, b2))
      fit.bic <- RaScreen(xtrain, ytrain, B1 = B1[b1], B2 = B2[b2], model = "lm", criterion = "bic", cores = cores, iteration = 1)
      mat[[b1, b2]] <- RaRank(fit.bic, selected.num="p", iteration=0)
      mat1[[b1, b2]] <- RaRank(fit.bic, selected.num="p", iteration=1)
  }
}

