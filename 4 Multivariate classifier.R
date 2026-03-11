#### Multivariate classifier ####
# written by Yinhan Chen and supervised by Qiang Luo
# Email: yinhanchen23@m.fudan.edu.cn
# released on 23 Sep 2025
# please cite: XXXXXX

library(caret)
library(pROC)
library(tidyverse)
library(Boruta)
library(psych)
library(randomForest)
library(e1071)  # bayesian
library(xgboost)  # XGBoost

###############
# load data
###############

rm(list=ls())
gc()

PKUdata = read.table(file="D:/BSC_and_MDD_v1/data/PKU/PKUdata.csv",header = T,sep = ",")
XYdata = read.table(file="D:/BSC_and_MDD_v1/data/XY/XYdata.csv",header = T,sep = ",")


FC_PKU_age_adjust = read.table(file="D:/BSC_and_MDD_v1/data/PKU/PKU_FC_bl_and_fu_age_adjusted.csv",sep = ",")

PKU_FC_bl = FC_PKU_age_adjust[c(1:168),]
PKU_FC_fu = FC_PKU_age_adjust[c(169:258),]


FC_XY_age_adjust = read.table(file="D:/BSC_and_MDD_v1/data/XY/XY_FC_bl_and_fu_age_adjusted.csv",sep = ",")

XY_FC_bl = FC_XY_age_adjust[c(1:411),]
XY_FC_fu = FC_XY_age_adjust[c(411:528),]

PKUdata$meanFC = rowMeans(PKU_FC_bl)
XYdata$meanFC = rowMeans(XY_FC_bl)

####################################
# classification models
####################################

ref_formula = diag~age+sex
BSC_formula = diag~age+sex+BSC_0w
FC_formula = diag~age+sex+meanFC
BSC_FC_formula = diag~age+sex+meanFC+BSC_0w

train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 100,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = TRUE
)

##################################
# rf
##################################
# 1. train in XY---------------------

## basic model
model_rf_XY = train(ref_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_XY.rds")

auc_rf_XY = model_rf_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_XY$AUC)
model_rf_XY$pred$sex <- XYdata$sex[model_rf_XY$pred$rowIndex]
auc_rf_XY_female <- model_rf_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_XY_female$AUC)
auc_rf_XY_male <- model_rf_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_XY_male$AUC)


## meanFC model
model_rf_FC_XY = train(FC_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_FC_XY.rds")
model_rf_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_rf_FC_XY.rds")

auc_rf_FC_XY = model_rf_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_FC_XY$AUC)
model_rf_FC_XY$pred$sex <- XYdata$sex[model_rf_FC_XY$pred$rowIndex]
auc_rf_FC_XY_female <- model_rf_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_XY_female$AUC)
auc_rf_FC_XY_male <- model_rf_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_XY_male$AUC)

## BSC model
model_rf_BSC_XY = train(BSC_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_XY.rds")
model_rf_BSC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_XY.rds")

auc_rf_BSC_XY = model_rf_BSC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_XY$AUC)
model_rf_BSC_XY$pred$sex <- XYdata$sex[model_rf_BSC_XY$pred$rowIndex]
auc_rf_BSC_XY_female <- model_rf_BSC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_XY_female$AUC)
auc_rf_BSC_XY_male <- model_rf_BSC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_XY_male$AUC)

## meanFC+BSC model
model_rf_BSC_FC_XY = train(BSC_FC_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_FC_XY.rds")
model_rf_BSC_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_FC_XY.rds")

model_rf_BSC_FC_XY_female = train(diag~age+meanFC+BSC_0w, data = XYdata[XYdata$sex=="female",], method = "rf", trControl = train_control)
auc_model_rf_BSC_FC_XY_female_trained =model_rf_BSC_FC_XY_female$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))

auc_rf_BSC_FC_XY = model_rf_BSC_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_FC_XY$AUC)
model_rf_BSC_FC_XY$pred$sex <- XYdata$sex[model_rf_BSC_FC_XY$pred$rowIndex]
auc_rf_BSC_FC_XY_female <- model_rf_BSC_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_XY_female$AUC)
auc_rf_BSC_FC_XY_male <- model_rf_BSC_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_XY_male$AUC)

t.test(auc_rf_BSC_FC_XY$AUC,auc_rf_XY$AUC)
t.test(auc_rf_BSC_FC_XY_female$AUC,auc_rf_XY_female$AUC)
t.test(auc_rf_BSC_FC_XY_male$AUC,auc_rf_XY_male$AUC)

t.test(auc_rf_BSC_FC_XY$AUC,auc_rf_FC_XY$AUC)
t.test(auc_rf_BSC_FC_XY_female$AUC,auc_rf_FC_XY_female$AUC)
t.test(auc_rf_BSC_FC_XY_male$AUC,auc_rf_FC_XY_male$AUC)

## visualization

# all
roc_fc = roc(model_rf_FC_XY$pred$obs, model_rf_FC_XY$pred$MDD,
             levels = c("HC", "MDD"), direction = "<")
roc_bsc_fc = roc(model_rf_BSC_FC_XY$pred$obs, model_rf_BSC_FC_XY$pred$MDD,
                 levels = c("HC", "MDD"), direction = "<")
tiff("D:/BSC_and_MDD_v1/picture/AUC_of_rf_train_on_XY.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,     
     font.lab = 2 )
lines(roc_fc, col = "#7ABF98")     
legend("bottomright", 
       title = "XY Both Sexes",
       legend = c("BSC Model (AUC = 0.72)",
                  "Basic Model (AUC = 0.72)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2,cex = 0.9,bty = "n")
dev.off()

# female
pred_bsc_fc_female <- subset(model_rf_BSC_FC_XY$pred, sex == "female")
pred_fc_female  <- subset(model_rf_FC_XY$pred, sex == "female")

roc_bsc_fc <- roc(pred_bsc_fc_female$obs, pred_bsc_fc_female$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_female$obs, pred_fc_female$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC_and_MDD_v1/picture/AUC_of_rf_train_on_XY_female.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "XY Female",
       legend = c("BSC Model (AUC = 0.75)",
                  "Basic Model (AUC = 0.75)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()

# male
pred_bsc_fc_male <- subset(model_rf_BSC_FC_XY$pred, sex == "male")
pred_fc_male  <- subset(model_rf_FC_XY$pred, sex == "male")

roc_bsc_fc <- roc(pred_bsc_fc_male$obs, pred_bsc_fc_male$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_male$obs, pred_fc_male$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC_and_MDD_v1/picture/AUC_of_rf_train_on_XY_male.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "XY male",
       legend = c("BSC Model (AUC = 0.64)",
                  "Basic Model (AUC = 0.64)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()


# 2. train in PKU-----------------------------

## basic model
model_rf_PKU = train(ref_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_PKU.rds")

auc_rf_PKU = model_rf_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_PKU$AUC)
model_rf_PKU$pred$sex <- PKUdata$sex[model_rf_PKU$pred$rowIndex]
auc_rf_PKU_female <- model_rf_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_PKU_female$AUC)
auc_rf_PKU_male <- model_rf_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_PKU_male$AUC)

## meanFC model
model_rf_FC_PKU = train(FC_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_FC_PKU.rds")
#model_rf_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_rf_FC_PKU.rds")

auc_rf_FC_PKU = model_rf_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_FC_PKU$AUC)
model_rf_FC_PKU$pred$sex <- PKUdata$sex[model_rf_FC_PKU$pred$rowIndex]
auc_rf_FC_PKU_female <- model_rf_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_PKU_female$AUC)
auc_rf_FC_PKU_male <- model_rf_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_PKU_male$AUC)

## BSC model
model_rf_BSC_PKU = train(BSC_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_PKU.rds")
model_rf_BSC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_PKU.rds")

auc_rf_BSC_PKU = model_rf_BSC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_PKU$AUC)
model_rf_BSC_PKU$pred$sex <- PKUdata$sex[model_rf_BSC_PKU$pred$rowIndex]
auc_rf_BSC_PKU_female <- model_rf_BSC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_PKU_female$AUC)
auc_rf_BSC_PKU_male <- model_rf_BSC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_PKU_male$AUC)

## meanFC+BSC model
model_rf_BSC_FC_PKU = train(BSC_FC_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_FC_PKU.rds")
#model_rf_BSC_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_FC_PKU.rds")

auc_rf_BSC_FC_PKU = model_rf_BSC_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_FC_PKU$AUC)
model_rf_BSC_FC_PKU$pred$sex <- PKUdata$sex[model_rf_BSC_FC_PKU$pred$rowIndex]
auc_rf_BSC_FC_PKU_female <- model_rf_BSC_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_PKU_female$AUC)
auc_rf_BSC_FC_PKU_male <- model_rf_BSC_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_PKU_male$AUC)

t.test(auc_rf_BSC_FC_PKU$AUC,auc_rf_FC_PKU$AUC)

# visualization
roc_fc = roc(model_rf_FC_PKU$pred$obs, model_rf_FC_PKU$pred$MDD,
             levels = c("HC", "MDD"), direction = "<")
roc_bsc_fc = roc(model_rf_BSC_FC_PKU$pred$obs, model_rf_BSC_FC_PKU$pred$MDD,
                 levels = c("HC", "MDD"), direction = "<")
tiff("D:/BSC and MDD/picture/bandpass01/AUC_of_rf_train_on_PKU.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,     
     font.lab = 2 )
lines(roc_fc, col = "#7ABF98")     
legend("bottomright", 
       title = "PKU Both Sexes",
       legend = c("BSC Model (AUC = 0.62)",
                  "Basic Model (AUC = 0.50)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2,cex = 0.9,bty = "n")
dev.off()

# female
pred_bsc_fc_female <- subset(model_rf_BSC_FC_PKU$pred, sex == "female")
pred_fc_female  <- subset(model_rf_FC_PKU$pred, sex == "female")

roc_bsc_fc <- roc(pred_bsc_fc_female$obs, pred_bsc_fc_female$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_female$obs, pred_fc_female$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC and MDD/picture/bandpass01/AUC_of_rf_train_on_PKU_female.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "PKU Female",
       legend = c("BSC Model (AUC = 0.60)",
                  "Basic Model (AUC = 0.51)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()

# male
pred_bsc_fc_male <- subset(model_rf_BSC_FC_PKU$pred, sex == "male")
pred_fc_male  <- subset(model_rf_FC_PKU$pred, sex == "male")

roc_bsc_fc <- roc(pred_bsc_fc_male$obs, pred_bsc_fc_male$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_male$obs, pred_fc_male$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC and MDD/picture/bandpass01/AUC_of_rf_train_on_PKU_male.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "PKU male",
       legend = c("BSC Model (AUC = 0.67)",
                  "Basic Model (AUC = 0.54)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()

# 3. test model in PKU 8week---------------------------

model_rf_FC_PKU = readRDS("D:/BSC and MDD/data/model_rf_FC_PKU.rds")
model_rf_BSC_FC_PKU = readRDS("D:/BSC and MDD/data/model_rf_BSC_FC_PKU.rds")

PKUdata_test = PKUdata[!is.na(PKUdata$BSC_8w),c(10,3,4,27,39)]
colnames(PKUdata_test) = c('diag','age','sex','BSC_0w','meanFC')

prediction_FC = predict(model_rf_FC_PKU,PKUdata_test, type = "prob")[2]
prediction_BSC_FC = predict(model_rf_BSC_FC_PKU,PKUdata_test, type = "prob")[2]

auc(roc(PKUdata_test$diag, prediction_FC$MDD,levels = c("HC", "MDD"), direction = "<"))
auc(roc(PKUdata_test$diag, prediction_BSC_FC$MDD,levels = c("HC", "MDD"), direction = "<"))


##################################
# rf
##################################
# 1. train in XY---------------------

## basic model
model_rf_XY = train(ref_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_XY.rds")

auc_rf_XY = model_rf_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_XY$AUC)
model_rf_XY$pred$sex <- XYdata$sex[model_rf_XY$pred$rowIndex]
auc_rf_XY_female <- model_rf_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_XY_female$AUC)
auc_rf_XY_male <- model_rf_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_XY_male$AUC)


## meanFC model
model_rf_FC_XY = train(FC_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_FC_XY.rds")
model_rf_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_rf_FC_XY.rds")

auc_rf_FC_XY = model_rf_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_FC_XY$AUC)
model_rf_FC_XY$pred$sex <- XYdata$sex[model_rf_FC_XY$pred$rowIndex]
auc_rf_FC_XY_female <- model_rf_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_XY_female$AUC)
auc_rf_FC_XY_male <- model_rf_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_XY_male$AUC)

## BSC model
model_rf_BSC_XY = train(BSC_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_XY.rds")
model_rf_BSC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_XY.rds")

auc_rf_BSC_XY = model_rf_BSC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_XY$AUC)
model_rf_BSC_XY$pred$sex <- XYdata$sex[model_rf_BSC_XY$pred$rowIndex]
auc_rf_BSC_XY_female <- model_rf_BSC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_XY_female$AUC)
auc_rf_BSC_XY_male <- model_rf_BSC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_XY_male$AUC)

## meanFC+BSC model
model_rf_BSC_FC_XY = train(BSC_FC_formula, data = XYdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_FC_XY.rds")
model_rf_BSC_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_rf_BSC_FC_XY.rds")

auc_rf_BSC_FC_XY = model_rf_BSC_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_FC_XY$AUC)
model_rf_BSC_FC_XY$pred$sex <- XYdata$sex[model_rf_BSC_FC_XY$pred$rowIndex]
auc_rf_BSC_FC_XY_female <- model_rf_BSC_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_XY_female$AUC)
auc_rf_BSC_FC_XY_male <- model_rf_BSC_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_XY_male$AUC)

t.test(auc_rf_BSC_FC_XY$AUC,auc_rf_XY$AUC)
t.test(auc_rf_BSC_FC_XY_female$AUC,auc_rf_XY_female$AUC)
t.test(auc_rf_BSC_FC_XY_male$AUC,auc_rf_XY_male$AUC)

t.test(auc_rf_BSC_FC_XY$AUC,auc_rf_FC_XY$AUC)
t.test(auc_rf_BSC_FC_XY_female$AUC,auc_rf_FC_XY_female$AUC)
t.test(auc_rf_BSC_FC_XY_male$AUC,auc_rf_FC_XY_male$AUC)


# 2. train in PKU-----------------------------

## basic model
model_rf_PKU = train(ref_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_PKU.rds")

auc_rf_PKU = model_rf_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_PKU$AUC)
model_rf_PKU$pred$sex <- PKUdata$sex[model_rf_PKU$pred$rowIndex]
auc_rf_PKU_female <- model_rf_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_PKU_female$AUC)
auc_rf_PKU_male <- model_rf_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_PKU_male$AUC)

## meanFC model
model_rf_FC_PKU = train(FC_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_FC_PKU.rds")
#model_rf_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_rf_FC_PKU.rds")

auc_rf_FC_PKU = model_rf_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_FC_PKU$AUC)
model_rf_FC_PKU$pred$sex <- PKUdata$sex[model_rf_FC_PKU$pred$rowIndex]
auc_rf_FC_PKU_female <- model_rf_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_PKU_female$AUC)
auc_rf_FC_PKU_male <- model_rf_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_FC_PKU_male$AUC)

## BSC model
model_rf_BSC_PKU = train(BSC_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_PKU.rds")
model_rf_BSC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_PKU.rds")

auc_rf_BSC_PKU = model_rf_BSC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_PKU$AUC)
model_rf_BSC_PKU$pred$sex <- PKUdata$sex[model_rf_BSC_PKU$pred$rowIndex]
auc_rf_BSC_PKU_female <- model_rf_BSC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_PKU_female$AUC)
auc_rf_BSC_PKU_male <- model_rf_BSC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_PKU_male$AUC)

## meanFC+BSC model
model_rf_BSC_FC_PKU = train(BSC_FC_formula, data = PKUdata, method = "rf", trControl = train_control)
saveRDS(model_rf_BSC_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_FC_PKU.rds")
#model_rf_BSC_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_rf_BSC_FC_PKU.rds")

auc_rf_BSC_FC_PKU = model_rf_BSC_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_rf_BSC_FC_PKU$AUC)
model_rf_BSC_FC_PKU$pred$sex <- PKUdata$sex[model_rf_BSC_FC_PKU$pred$rowIndex]
auc_rf_BSC_FC_PKU_female <- model_rf_BSC_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_PKU_female$AUC)
auc_rf_BSC_FC_PKU_male <- model_rf_BSC_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_rf_BSC_FC_PKU_male$AUC)

t.test(auc_rf_BSC_FC_PKU$AUC,auc_rf_FC_PKU$AUC)

##################################
# glm
##################################
# 1. train in XY---------------------

## basic model
model_glm_XY = train(ref_formula, data = XYdata, method = "glm", trControl = train_control)
saveRDS(model_glm_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_glm_XY.rds")

auc_glm_XY = model_glm_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_XY$AUC)
model_glm_XY$pred$sex <- XYdata$sex[model_glm_XY$pred$rowIndex]
auc_glm_XY_female <- model_glm_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_XY_female$AUC)
auc_glm_XY_male <- model_glm_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_XY_male$AUC)


## meanFC model
model_glm_FC_XY = train(FC_formula, data = XYdata, method = "glm", trControl = train_control)
saveRDS(model_glm_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_glm_FC_XY.rds")
model_glm_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_glm_FC_XY.rds")

auc_glm_FC_XY = model_glm_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_FC_XY$AUC)
model_glm_FC_XY$pred$sex <- XYdata$sex[model_glm_FC_XY$pred$rowIndex]
auc_glm_FC_XY_female <- model_glm_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_FC_XY_female$AUC)
auc_glm_FC_XY_male <- model_glm_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_FC_XY_male$AUC)

## BSC model
model_glm_BSC_XY = train(BSC_formula, data = XYdata, method = "glm", trControl = train_control)
saveRDS(model_glm_BSC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_glm_BSC_XY.rds")
model_glm_BSC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_glm_BSC_XY.rds")

auc_glm_BSC_XY = model_glm_BSC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_BSC_XY$AUC)
model_glm_BSC_XY$pred$sex <- XYdata$sex[model_glm_BSC_XY$pred$rowIndex]
auc_glm_BSC_XY_female <- model_glm_BSC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_XY_female$AUC)
auc_glm_BSC_XY_male <- model_glm_BSC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_XY_male$AUC)

## meanFC+BSC model
model_glm_BSC_FC_XY = train(BSC_FC_formula, data = XYdata, method = "glm", trControl = train_control)
saveRDS(model_glm_BSC_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_glm_BSC_FC_XY.rds")
model_glm_BSC_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_glm_BSC_FC_XY.rds")

auc_glm_BSC_FC_XY = model_glm_BSC_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_BSC_FC_XY$AUC)
model_glm_BSC_FC_XY$pred$sex <- XYdata$sex[model_glm_BSC_FC_XY$pred$rowIndex]
auc_glm_BSC_FC_XY_female <- model_glm_BSC_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_FC_XY_female$AUC)
auc_glm_BSC_FC_XY_male <- model_glm_BSC_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_FC_XY_male$AUC)

t.test(auc_glm_BSC_FC_XY$AUC,auc_glm_XY$AUC)
t.test(auc_glm_BSC_FC_XY_female$AUC,auc_glm_XY_female$AUC)
t.test(auc_glm_BSC_FC_XY_male$AUC,auc_glm_XY_male$AUC)

t.test(auc_glm_BSC_FC_XY$AUC,auc_glm_FC_XY$AUC)
t.test(auc_glm_BSC_FC_XY_female$AUC,auc_glm_FC_XY_female$AUC)
t.test(auc_glm_BSC_FC_XY_male$AUC,auc_glm_FC_XY_male$AUC)

## visualization

# all
roc_fc = roc(model_glm_FC_XY$pred$obs, model_glm_FC_XY$pred$MDD,
             levels = c("HC", "MDD"), direction = "<")
roc_bsc_fc = roc(model_glm_BSC_FC_XY$pred$obs, model_glm_BSC_FC_XY$pred$MDD,
                 levels = c("HC", "MDD"), direction = "<")
tiff("D:/BSC_and_MDD_v1/picture/AUC_of_glm_train_on_XY.tiff", width = 1200, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,     
     font.lab = 2,
     lwd = 3,
     lty = 1)
lines(roc_fc, col = "#7ABF98",lty = 2,lwd = 3)     
title(main = "XYH Both Sexes",
      font.main = 2,
      cex.main = 0.9)
legend("bottomright", 
       inset = 0.02,
       legend = c("BSC Model (AUC = 0.60)",
                  "Basic Model (AUC = 0.58)"),
       lty = c(1,2),   
       col = c("#BF355E", "#7ABF98"), lwd = 3,cex = 0.9,bty = "n")
dev.off()

# female
pred_bsc_fc_female <- subset(model_glm_BSC_FC_XY$pred, sex == "female")
pred_fc_female  <- subset(model_glm_FC_XY$pred, sex == "female")

roc_bsc_fc <- roc(pred_bsc_fc_female$obs, pred_bsc_fc_female$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_female$obs, pred_fc_female$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC_and_MDD_v1/picture/AUC_of_glm_train_on_XY_female.tiff", width = 1200, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     lwd = 3,
     lty = 1,
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98",lty = 2,lwd = 3)  
title(main = "XYH Female",
      font.main = 2,
      cex.main = 0.9)
legend("bottomright", 
       inset = 0.02,
       lty = c(1,2),  
       legend = c("BSC Model (AUC = 0.61)",
                  "Basic Model (AUC = 0.59)"),
       col = c("#BF355E", "#7ABF98"), lwd = 3, cex = 0.9, bty = "n")
dev.off()

# male
pred_bsc_fc_male <- subset(model_glm_BSC_FC_XY$pred, sex == "male")
pred_fc_male  <- subset(model_glm_FC_XY$pred, sex == "male")

roc_bsc_fc <- roc(pred_bsc_fc_male$obs, pred_bsc_fc_male$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_male$obs, pred_fc_male$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC_and_MDD_v1/picture/AUC_of_glm_train_on_XY_male.tiff", width = 1200, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     lwd = 3,
     lty = 1,
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98",lty = 2,lwd = 3)   
title(main = "XYH Male",
      font.main = 2,
      cex.main = 0.9)
legend("bottomright", 
       inset = 0.02,
       lty = c(1,2),  
       legend = c("BSC Model (AUC = 0.59)",
                  "Basic Model (AUC = 0.58)"),
       col = c("#BF355E", "#7ABF98"), lwd = 3, cex = 0.9, bty = "n")
dev.off()



# 2. train in PKU-----------------------------

## basic model
model_glm_PKU = train(ref_formula, data = PKUdata, method = "glm", trControl = train_control)
saveRDS(model_glm_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_glm_PKU.rds")

auc_glm_PKU = model_glm_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_PKU$AUC)
model_glm_PKU$pred$sex <- PKUdata$sex[model_glm_PKU$pred$rowIndex]
auc_glm_PKU_female <- model_glm_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_PKU_female$AUC)
auc_glm_PKU_male <- model_glm_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_PKU_male$AUC)

## meanFC model
model_glm_FC_PKU = train(FC_formula, data = PKUdata, method = "glm", trControl = train_control)
saveRDS(model_glm_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_glm_FC_PKU.rds")
model_glm_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_glm_FC_PKU.rds")

auc_glm_FC_PKU = model_glm_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_FC_PKU$AUC)
model_glm_FC_PKU$pred$sex <- PKUdata$sex[model_glm_FC_PKU$pred$rowIndex]
auc_glm_FC_PKU_female <- model_glm_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_FC_PKU_female$AUC)
auc_glm_FC_PKU_male <- model_glm_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_FC_PKU_male$AUC)

## BSC model
model_glm_BSC_PKU = train(BSC_formula, data = PKUdata, method = "glm", trControl = train_control)
saveRDS(model_glm_BSC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_glm_BSC_PKU.rds")
model_glm_BSC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_glm_BSC_PKU.rds")

auc_glm_BSC_PKU = model_glm_BSC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_BSC_PKU$AUC)
model_glm_BSC_PKU$pred$sex <- PKUdata$sex[model_glm_BSC_PKU$pred$rowIndex]
auc_glm_BSC_PKU_female <- model_glm_BSC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_PKU_female$AUC)
auc_glm_BSC_PKU_male <- model_glm_BSC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_PKU_male$AUC)

## meanFC+BSC model
model_glm_BSC_FC_PKU = train(BSC_FC_formula, data = PKUdata, method = "glm", trControl = train_control)
saveRDS(model_glm_BSC_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_glm_BSC_FC_PKU.rds")
model_glm_BSC_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_glm_BSC_FC_PKU.rds")

auc_glm_BSC_FC_PKU = model_glm_BSC_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_glm_BSC_FC_PKU$AUC)
model_glm_BSC_FC_PKU$pred$sex <- PKUdata$sex[model_glm_BSC_FC_PKU$pred$rowIndex]
auc_glm_BSC_FC_PKU_female <- model_glm_BSC_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_FC_PKU_female$AUC)
auc_glm_BSC_FC_PKU_male <- model_glm_BSC_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_glm_BSC_FC_PKU_male$AUC)

t.test(auc_glm_BSC_FC_PKU$AUC,auc_glm_PKU$AUC)
t.test(auc_glm_BSC_FC_PKU_female$AUC,auc_glm_PKU_female$AUC)
t.test(auc_glm_BSC_FC_PKU_male$AUC,auc_glm_PKU_male$AUC)

t.test(auc_glm_BSC_FC_PKU$AUC,auc_glm_FC_PKU$AUC)
t.test(auc_glm_BSC_FC_PKU_female$AUC,auc_glm_FC_PKU_female$AUC)
t.test(auc_glm_BSC_FC_PKU_male$AUC,auc_glm_FC_PKU_male$AUC)

## visualization

# all
roc_fc = roc(model_glm_FC_PKU$pred$obs, model_glm_FC_PKU$pred$MDD,
             levels = c("HC", "MDD"), direction = "<")
roc_bsc_fc = roc(model_glm_BSC_FC_PKU$pred$obs, model_glm_BSC_FC_PKU$pred$MDD,
                 levels = c("HC", "MDD"), direction = "<")
tiff("D:/BSC_and_MDD_v1/picture/AUC_of_glm_train_on_PKU.tiff", width = 1200, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,   
     lwd = 3,
     lty = 1,
     font.lab = 2 )
lines(roc_fc, col = "#7ABF98",lty = 2,lwd = 3)     
title(main = "PKU Both Sexes",
      font.main = 2,
      cex.main = 0.9)
legend("bottomright", 
       inset = 0.02,
       legend = c("BSC Model (AUC = 0.45)",
                  "Basic Model (AUC = 0.62)"),
       lty = c(1,2), 
       col = c("#BF355E", "#7ABF98"), lwd = 3,cex = 0.9,bty = "n")
dev.off()

# female
pred_bsc_fc_female <- subset(model_glm_BSC_FC_PKU$pred, sex == "female")
pred_fc_female  <- subset(model_glm_FC_PKU$pred, sex == "female")

roc_bsc_fc <- roc(pred_bsc_fc_female$obs, pred_bsc_fc_female$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_female$obs, pred_fc_female$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC_and_MDD_v1/picture/AUC_of_glm_train_on_PKU_female.tiff", width = 1200, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     lwd = 3,
     lty = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98",lty = 2,lwd = 3)  
title(main = "PKU Female",
      font.main = 2,
      cex.main = 0.9)
legend("bottomright", 
       inset = 0.02,
       legend = c("BSC Model (AUC = 0.45)",
                  "Basic Model (AUC = 0.61)"),
       lty = c(1,2), 
       col = c("#BF355E", "#7ABF98"), lwd = 3, cex = 0.9, bty = "n")
dev.off()

# male
pred_bsc_fc_male <- subset(model_glm_BSC_FC_PKU$pred, sex == "male")
pred_fc_male  <- subset(model_glm_FC_PKU$pred, sex == "male")

roc_bsc_fc <- roc(pred_bsc_fc_male$obs, pred_bsc_fc_male$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_male$obs, pred_fc_male$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC_and_MDD_v1/picture/AUC_of_glm_train_on_PKU_male.tiff", width = 1200, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     lwd = 3,
     lty = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98",lty = 2,lwd = 3)   
title(main = "PKU Male",
      font.main = 2,
      cex.main = 0.9)
legend("bottomright", 
       inset = 0.02,
       legend = c("BSC Model (AUC = 0.49)",
                  "Basic Model (AUC = 0.65)"),
       lty = c(1,2), 
       col = c("#BF355E", "#7ABF98"), lwd = 3, cex = 0.9, bty = "n")
dev.off()

##################################
# xgbTree
##################################
# 1. train in XY---------------------

## basic model
model_xgbTree_XY = train(ref_formula, data = XYdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_xgbTree_XY.rds")

auc_xgbTree_XY = model_xgbTree_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_XY$AUC)
model_xgbTree_XY$pred$sex <- XYdata$sex[model_xgbTree_XY$pred$rowIndex]
auc_xgbTree_XY_female <- model_xgbTree_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_XY_female$AUC)
auc_xgbTree_XY_male <- model_xgbTree_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_XY_male$AUC)


## meanFC model
model_xgbTree_FC_XY = train(FC_formula, data = XYdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_xgbTree_FC_XY.rds")
model_xgbTree_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_xgbTree_FC_XY.rds")

auc_xgbTree_FC_XY = model_xgbTree_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_FC_XY$AUC)
model_xgbTree_FC_XY$pred$sex <- XYdata$sex[model_xgbTree_FC_XY$pred$rowIndex]
auc_xgbTree_FC_XY_female <- model_xgbTree_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_FC_XY_female$AUC)
auc_xgbTree_FC_XY_male <- model_xgbTree_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_FC_XY_male$AUC)

## BSC model
model_xgbTree_BSC_XY = train(BSC_formula, data = XYdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_BSC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_xgbTree_BSC_XY.rds")
model_xgbTree_BSC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_xgbTree_BSC_XY.rds")

auc_xgbTree_BSC_XY = model_xgbTree_BSC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_BSC_XY$AUC)
model_xgbTree_BSC_XY$pred$sex <- XYdata$sex[model_xgbTree_BSC_XY$pred$rowIndex]
auc_xgbTree_BSC_XY_female <- model_xgbTree_BSC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_XY_female$AUC)
auc_xgbTree_BSC_XY_male <- model_xgbTree_BSC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_XY_male$AUC)

## meanFC+BSC model
model_xgbTree_BSC_FC_XY = train(BSC_FC_formula, data = XYdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_BSC_FC_XY, file = "D:/BSC_and_MDD_v1/data/XY/model_xgbTree_BSC_FC_XY.rds")
model_xgbTree_BSC_FC_XY = readRDS("D:/BSC_and_MDD_v1/data/XY/model_xgbTree_BSC_FC_XY.rds")

auc_xgbTree_BSC_FC_XY = model_xgbTree_BSC_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_BSC_FC_XY$AUC)
model_xgbTree_BSC_FC_XY$pred$sex <- XYdata$sex[model_xgbTree_BSC_FC_XY$pred$rowIndex]
auc_xgbTree_BSC_FC_XY_female <- model_xgbTree_BSC_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_FC_XY_female$AUC)
auc_xgbTree_BSC_FC_XY_male <- model_xgbTree_BSC_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_FC_XY_male$AUC)

t.test(auc_xgbTree_BSC_FC_XY$AUC,auc_xgbTree_XY$AUC)
t.test(auc_xgbTree_BSC_FC_XY_female$AUC,auc_xgbTree_XY_female$AUC)
t.test(auc_xgbTree_BSC_FC_XY_male$AUC,auc_xgbTree_XY_male$AUC)

t.test(auc_xgbTree_BSC_FC_XY$AUC,auc_xgbTree_FC_XY$AUC)
t.test(auc_xgbTree_BSC_FC_XY_female$AUC,auc_xgbTree_FC_XY_female$AUC)
t.test(auc_xgbTree_BSC_FC_XY_male$AUC,auc_xgbTree_FC_XY_male$AUC)


# 2. train in PKU-----------------------------

## basic model
model_xgbTree_PKU = train(ref_formula, data = PKUdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_PKU.rds")

auc_xgbTree_PKU = model_xgbTree_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_PKU$AUC)
model_xgbTree_PKU$pred$sex <- PKUdata$sex[model_xgbTree_PKU$pred$rowIndex]
auc_xgbTree_PKU_female <- model_xgbTree_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_PKU_female$AUC)
auc_xgbTree_PKU_male <- model_xgbTree_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_PKU_male$AUC)

## meanFC model
model_xgbTree_FC_PKU = train(FC_formula, data = PKUdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_FC_PKU.rds")
#model_xgbTree_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_FC_PKU.rds")

auc_xgbTree_FC_PKU = model_xgbTree_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_FC_PKU$AUC)
model_xgbTree_FC_PKU$pred$sex <- PKUdata$sex[model_xgbTree_FC_PKU$pred$rowIndex]
auc_xgbTree_FC_PKU_female <- model_xgbTree_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_FC_PKU_female$AUC)
auc_xgbTree_FC_PKU_male <- model_xgbTree_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_FC_PKU_male$AUC)

## BSC model
model_xgbTree_BSC_PKU = train(BSC_formula, data = PKUdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_BSC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_BSC_PKU.rds")
model_xgbTree_BSC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_BSC_PKU.rds")

auc_xgbTree_BSC_PKU = model_xgbTree_BSC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_BSC_PKU$AUC)
model_xgbTree_BSC_PKU$pred$sex <- PKUdata$sex[model_xgbTree_BSC_PKU$pred$rowIndex]
auc_xgbTree_BSC_PKU_female <- model_xgbTree_BSC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_PKU_female$AUC)
auc_xgbTree_BSC_PKU_male <- model_xgbTree_BSC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_PKU_male$AUC)

## meanFC+BSC model
model_xgbTree_BSC_FC_PKU = train(BSC_FC_formula, data = PKUdata, method = "xgbTree", trControl = train_control)
saveRDS(model_xgbTree_BSC_FC_PKU, file = "D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_BSC_FC_PKU.rds")
#model_xgbTree_BSC_FC_PKU = readRDS("D:/BSC_and_MDD_v1/data/PKU/model_xgbTree_BSC_FC_PKU.rds")

auc_xgbTree_BSC_FC_PKU = model_xgbTree_BSC_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_xgbTree_BSC_FC_PKU$AUC)
model_xgbTree_BSC_FC_PKU$pred$sex <- PKUdata$sex[model_xgbTree_BSC_FC_PKU$pred$rowIndex]
auc_xgbTree_BSC_FC_PKU_female <- model_xgbTree_BSC_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_FC_PKU_female$AUC)
auc_xgbTree_BSC_FC_PKU_male <- model_xgbTree_BSC_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_xgbTree_BSC_FC_PKU_male$AUC)

t.test(auc_xgbTree_BSC_FC_PKU$AUC,auc_xgbTree_PKU$AUC)
t.test(auc_xgbTree_BSC_FC_PKU_female$AUC,auc_xgbTree_PKU_female$AUC)
t.test(auc_xgbTree_BSC_FC_PKU_male$AUC,auc_xgbTree_PKU_male$AUC)

t.test(auc_xgbTree_BSC_FC_PKU$AUC,auc_xgbTree_FC_PKU$AUC)
t.test(auc_xgbTree_BSC_FC_PKU_female$AUC,auc_xgbTree_FC_PKU_female$AUC)
t.test(auc_xgbTree_BSC_FC_PKU_male$AUC,auc_xgbTree_FC_PKU_male$AUC)
