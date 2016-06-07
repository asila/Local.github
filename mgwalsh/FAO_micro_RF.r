#' Random Forest regression bioavailability predictions of wheat plant micronutrient contents from soil data
#' data courtesy FAO (1982) & ICRAF (2016)
#' M. Walsh, May 2016

# Required packages
# install.packages(c("devtools","caret","doParallel","randomForest")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(doParallel)
  require(randomForest)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Target elements: Wheat plant micro-nutrient concentrations (ppm)
pB  <- fao_cal$pB  ## Boron
pCu <- fao_cal$pCu ## Copper
pMn <- fao_cal$pMn ## Manganese
pMo <- fao_cal$pMo ## Molybdenum
pZn <- fao_cal$pZn ## Zinc
pFe <- fao_cal$pFe ## Iron

# Covariates
wetc <- fao_cal[c(4:24)] ## Wet chemistry calibration data
wetv <- fao_val[c(4:24)] ## Wet chemistry validation data from 8 randomly selected countries
mirc <- fao_cal[c(32:1795)] ## MIR calibration data
mirv <- fao_val[c(32:1795)] ## MIR validation data for 8 randomly selected countries

# RF models ---------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "oob", allowParallel = TRUE)
tg <- expand.grid(mtry=seq(20, 200, by=10))

# Plant Boron content (ppm) -----------------------------------------------
# Wet chemistry covariates
pB_wet.rfo <- train(wetc, pB,
                    preProc = c("center", "scale"),
                    method = "rf",
                    ntree = 501,
                    tuneGrid = tg,
                    trControl = tc)
print(pB_wet.rfo)
pB_wet.imp <- varImp(pB_wet.rfo, useModel = FALSE)
plot(pB_wet.imp, top=20)

# MIR covariates
pB_mir.rfo <- train(mirc, pB,
                    preProc = c("center", "scale"),
                    method = "rf",
                    ntree = 501,
                    tuneGrid = tg,
                    trControl = tc)
print(pB_mir.rfo)
pB_mir.imp <- varImp(pB_mir.rfo, useModel = FALSE)
plot(pB_mir.imp, top=20)

# Plant Copper content (ppm) ----------------------------------------------
# Wet chemistry covariates
pCu_wet.rfo <- train(wetc, pCu,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pCu_wet.rfo)
pCu_wet.imp <- varImp(pCu_wet.rfo, useModel = FALSE)
plot(pCu_wet.imp, top=20)

# MIR covariates
pCu_mir.rfo <- train(mirc, pCu,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pCu_mir.rfo)
pCu_mir.imp <- varImp(pCu_mir.rfo, useModel = FALSE)
plot(pCu_mir.imp, top=20)

# Plant Manganese content (ppm) -------------------------------------------
# Wet chemistry covariates
pMn_wet.rfo <- train(wetc, pMn,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pMn_wet.rfo)
pMn_wet.imp <- varImp(pMn_wet.rfo, useModel = FALSE)
plot(pMn_wet.imp, top=20)

# MIR covariates
pMn_mir.rfo <- train(mirc, pMn,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pMn_mir.rfo)
pMn_mir.imp <- varImp(pMn_mir.rfo, useModel = FALSE)
plot(pMn_mir.imp, top=20)

# Plant Molybdenum content (ppm) ------------------------------------------
# Wet chemistry covariates
pMo_wet.rfo <- train(wetc, pMo,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pMo_wet.rfo)
pMo_wet.imp <- varImp(pMo_wet.rfo, useModel = FALSE)
plot(pMo_wet.imp, top=20)

# MIR covariates
pMo_mir.rfo <- train(mirc, pMo,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pMo_mir.rfo)
pMo_mir.imp <- varImp(pMo_mir.rfo, useModel = FALSE)
plot(pMo_mir.imp, top=20)

# Plant Zinc content (ppm) ------------------------------------------------
# Wet chemistry covariates
pZn_wet.rfo <- train(wetc, pZn,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pZn_wet.rfo)
pZn_wet.imp <- varImp(pZn_wet.rfo, useModel = FALSE)
plot(pZn_wet.imp, top=20)

# MIR covariates
pZn_mir.rfo <- train(mirc, pZn,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pZn_mir.rfo)
pZn_mir.imp <- varImp(pZn_mir.rfo, useModel = FALSE)
plot(pZn_mir.imp, top=20)

# Plant Iron contents (ppm) -----------------------------------------------
# Wet chemistry covariates
pFe_wet.rfo <- train(wetc, pFe,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pFe_wet.rfo)
pFe_wet.imp <- varImp(pFe_wet.rfo, useModel = FALSE)
plot(pZn_wet.imp, top=20)

# MIR covariates
pFe_mir.rfo <- train(mirc, pFe,
                     preProc = c("center", "scale"),
                     method = "rf",
                     ntree = 501,
                     tuneGrid = tg,
                     trControl = tc)
print(pFe_mir.rfo)
pFe_mir.imp <- varImp(pFe_mir.rfo, useModel = FALSE)
plot(pFe_mir.imp, top=20)
