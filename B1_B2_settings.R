
library(RaSEn)
library(screening)
library(caret)
library(VariableScreening)
library(e1071)
library(dplyr)
library(SIS)


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

