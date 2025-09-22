#### Multivariate classifier ####
# written by Yinhan Chen and supervised by Qiang Luo
# Email: yinhanchen23@m.fudan.edu.cn
# released on 23 Sep 2025
# please cite: XXXXXX

library(R.matlab)
library(caret)
library(pROC)
library(tidyverse)
library(Boruta)

rm(list=ls())
gc()

#############################################
# load data
#############################################

# BSC data
PKUdata = read.table(file="PKUdata.csv",header = T,sep = ",")
XYdata = read.table(file="XYdata.csv",header = T,sep = ",")

XYdata$BSC_group_0w  = factor(XYdata$BSC_group_0w, levels = c("female-like", "androgynous", "male-like"))
XYdata$BSC_group_8w  = factor(XYdata$BSC_group_8w, levels = c("female-like", "androgynous", "male-like"))
XYdata$remission = factor(XYdata$remission,levels=c("nonremission","remission"))
XYdata$response = factor(XYdata$response,levels=c("nonresponse","response"))

PKUdata$BSC_group_0w  = factor(PKUdata$BSC_group_0w, levels = c("female-like", "androgynous", "male-like"))
PKUdata$BSC_group_8w  = factor(PKUdata$BSC_group_8w, levels = c("female-like", "androgynous", "male-like"))
PKUdata$remission = factor(PKUdata$remission,levels=c("nonremission","remission"))
PKUdata$response = factor(PKUdata$response,levels=c("nonresponse","response"))

# FC data

FC_data = readMat("alldata_for_combat_analysis.mat",dat)

FC = FC_data$dat
batch = as.vector(FC_data$batch)
mod = FC_data$mod

PKU_FC_bl = FC[c(1:168),]
PKU_FC_fu = FC[c(169:258),]
PKU_FC = FC[c(1:258),]
XY_FC_bl = FC[c(259:672),]
XY_FC_fu = FC[c(673:793),]
XY_FC = FC[c(259:793),]

# calculation of mean FC 

PKUdata$meanFC = rowMeans(PKU_FC_bl)
PKUdata$meanFC_8w[!is.na(PKUdata$BSC_8w)] = rowMeans(PKU_FC_fu)
XYdata$meanFC = rowMeans(XY_FC_bl)

PKUdata$sex = as.factor(PKUdata$sex)
XYdata$sex = as.factor(XYdata$sex)


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

# 1. train in XY---------------------

## meanFC model
model_LMT_FC_XY = train(FC_formula, data = XYdata, method = "LMT", trControl = train_control)
saveRDS(model_LMT_FC_XY, file = "D:/BSC and MDD/data/model_LMT_FC_XY.rds")
model_LMT_FC_XY = readRDS("D:/BSC and MDD/data/model_LMT_FC_XY.rds")

auc_LMT_FC_XY = model_LMT_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_LMT_FC_XY$AUC)
model_LMT_FC_XY$pred$sex <- XYdata$sex[model_LMT_FC_XY$pred$rowIndex]
auc_LMT_FC_XY_female <- model_LMT_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_FC_XY_female$AUC)
auc_LMT_FC_XY_male <- model_LMT_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_FC_XY_male$AUC)

## meanFC+BSC model
model_LMT_BSC_FC_XY = train(BSC_FC_formula, data = XYdata, method = "LMT", trControl = train_control)
saveRDS(model_LMT_BSC_FC_XY, file = "D:/BSC and MDD/data/model_LMT_BSC_FC_XY.rds")
model_LMT_BSC_FC_XY = readRDS("D:/BSC and MDD/data/model_LMT_BSC_FC_XY.rds")

auc_LMT_BSC_FC_XY = model_LMT_BSC_FC_XY$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_LMT_BSC_FC_XY$AUC)
model_LMT_BSC_FC_XY$pred$sex <- XYdata$sex[model_LMT_BSC_FC_XY$pred$rowIndex]
auc_LMT_BSC_FC_XY_female <- model_LMT_BSC_FC_XY$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_BSC_FC_XY_female$AUC)
auc_LMT_BSC_FC_XY_male <- model_LMT_BSC_FC_XY$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_BSC_FC_XY_male$AUC)

t.test(auc_LMT_BSC_FC_XY$AUC,auc_LMT_FC_XY$AUC)

## visualization

# all
roc_fc = roc(model_LMT_FC_XY$pred$obs, model_LMT_FC_XY$pred$MDD,
             levels = c("HC", "MDD"), direction = "<")
roc_bsc_fc = roc(model_LMT_BSC_FC_XY$pred$obs, model_LMT_BSC_FC_XY$pred$MDD,
                 levels = c("HC", "MDD"), direction = "<")
tiff("D:/BSC and MDD/picture/AUC_of_LMT_train_on_XY.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,     
     font.lab = 2 )
lines(roc_fc, col = "#7ABF98")     
legend("bottomright", 
       title = "XY Both Sexes",
       legend = c("BSC Model (AUC = 0.72)",
                  "Basic Model (AUC = 0.71)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2,cex = 0.9,bty = "n")
dev.off()

# female
pred_bsc_fc_female <- subset(model_LMT_BSC_FC_XY$pred, sex == "female")
pred_fc_female  <- subset(model_LMT_FC_XY$pred, sex == "female")

roc_bsc_fc <- roc(pred_bsc_fc_female$obs, pred_bsc_fc_female$MDD,
               levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_female$obs, pred_fc_female$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC and MDD/picture/AUC_of_LMT_train_on_XY_female.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "XY Female",
       legend = c("BSC Model (AUC = 0.75)",
                  "Basic Model (AUC = 0.74)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()

# male
pred_bsc_fc_male <- subset(model_LMT_BSC_FC_XY$pred, sex == "male")
pred_fc_male  <- subset(model_LMT_FC_XY$pred, sex == "male")

roc_bsc_fc <- roc(pred_bsc_fc_male$obs, pred_bsc_fc_male$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_male$obs, pred_fc_male$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC and MDD/picture/AUC_of_LMT_train_on_XY_male.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "XY male",
       legend = c("BSC Model (AUC = 0.64)",
                  "Basic Model (AUC = 0.63)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()


# 2. train in PKU-----------------------------

## meanFC model
model_LMT_FC_PKU = train(FC_formula, data = PKUdata, method = "LMT", trControl = train_control)
saveRDS(model_LMT_FC_PKU, file = "D:/BSC and MDD/data/model_LMT_FC_PKU.rds")
model_LMT_FC_PKU = readRDS("D:/BSC and MDD/data/model_LMT_FC_PKU.rds")

auc_LMT_FC_PKU = model_LMT_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_LMT_FC_PKU$AUC)
model_LMT_FC_PKU$pred$sex <- PKUdata$sex[model_LMT_FC_PKU$pred$rowIndex]
auc_LMT_FC_PKU_female <- model_LMT_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_FC_PKU_female$AUC)
auc_LMT_FC_PKU_male <- model_LMT_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_FC_PKU_male$AUC)

## meanFC+BSC model
model_LMT_BSC_FC_PKU = train(BSC_FC_formula, data = PKUdata, method = "LMT", trControl = train_control)
saveRDS(model_LMT_BSC_FC_PKU, file = "D:/BSC and MDD/data/model_LMT_BSC_FC_PKU.rds")
model_LMT_BSC_FC_PKU = readRDS("D:/BSC and MDD/data/model_LMT_BSC_FC_PKU.rds")

auc_LMT_BSC_FC_PKU = model_LMT_BSC_FC_PKU$pred %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD,levels = c("HC", "MDD"),direction = "<")$auc))
mean(auc_LMT_BSC_FC_PKU$AUC)
model_LMT_BSC_FC_PKU$pred$sex <- PKUdata$sex[model_LMT_BSC_FC_PKU$pred$rowIndex]
auc_LMT_BSC_FC_PKU_female <- model_LMT_BSC_FC_PKU$pred %>%
  filter(sex == "female") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_BSC_FC_PKU_female$AUC)
auc_LMT_BSC_FC_PKU_male <- model_LMT_BSC_FC_PKU$pred %>%
  filter(sex == "male") %>%
  group_by(Resample) %>%
  summarise(AUC = as.numeric(roc(obs, MDD, levels = c("HC", "MDD"), direction = "<")$auc))
mean(auc_LMT_BSC_FC_PKU_male$AUC)

t.test(auc_LMT_BSC_FC_PKU$AUC,auc_LMT_FC_PKU$AUC)

# visualization
roc_fc = roc(model_LMT_FC_PKU$pred$obs, model_LMT_FC_PKU$pred$MDD,
             levels = c("HC", "MDD"), direction = "<")
roc_bsc_fc = roc(model_LMT_BSC_FC_PKU$pred$obs, model_LMT_BSC_FC_PKU$pred$MDD,
                 levels = c("HC", "MDD"), direction = "<")
tiff("D:/BSC and MDD/picture/AUC_of_LMT_train_on_PKU.tiff", width = 1200, height = 1200, res = 300)
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
pred_bsc_fc_female <- subset(model_LMT_BSC_FC_PKU$pred, sex == "female")
pred_fc_female  <- subset(model_LMT_FC_PKU$pred, sex == "female")

roc_bsc_fc <- roc(pred_bsc_fc_female$obs, pred_bsc_fc_female$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_female$obs, pred_fc_female$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC and MDD/picture/AUC_of_LMT_train_on_PKU_female.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "PKU Female",
       legend = c("BSC Model (AUC = 0.59)",
                  "Basic Model (AUC = 0.52)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()

# male
pred_bsc_fc_male <- subset(model_LMT_BSC_FC_PKU$pred, sex == "male")
pred_fc_male  <- subset(model_LMT_FC_PKU$pred, sex == "male")

roc_bsc_fc <- roc(pred_bsc_fc_male$obs, pred_bsc_fc_male$MDD,
                  levels = c("HC", "MDD"), direction = "<")
roc_fc <- roc(pred_fc_male$obs, pred_fc_male$MDD,
              levels = c("HC", "MDD"), direction = "<")

tiff("D:/BSC and MDD/picture/AUC_of_LMT_train_on_PKU_male.tiff", width = 1200, height = 1200, res = 300)
plot(roc_bsc_fc, col = "#BF355E", 
     xlab = "1-Specificity", 
     ylab = "Sensitivity",
     cex.lab = 1,
     font.lab = 2)
lines(roc_fc, col = "#7ABF98")   
legend("bottomright", 
       title = "PKU male",
       legend = c("BSC Model (AUC = 0.67)",
                  "Basic Model (AUC = 0.53)"),
       col = c("#BF355E", "#7ABF98"), lwd = 2, cex = 0.9, bty = "n")
dev.off()

# 3. test model in PKU 8week---------------------------

model_LMT_FC_PKU = readRDS("D:/BSC and MDD/data/model_LMT_FC_PKU.rds")
model_LMT_BSC_FC_PKU = readRDS("D:/BSC and MDD/data/model_LMT_BSC_FC_PKU.rds")

PKUdata_test = PKUdata[!is.na(PKUdata$BSC_8w),c(10,3,4,27,39)]
colnames(PKUdata_test) = c('diag','age','sex','BSC_0w','meanFC')

prediction_FC = predict(model_LMT_FC_PKU,PKUdata_test, type = "prob")[2]
prediction_BSC_FC = predict(model_LMT_BSC_FC_PKU,PKUdata_test, type = "prob")[2]

auc(roc(PKUdata_test$diag, prediction_FC$MDD,levels = c("HC", "MDD"), direction = "<"))
auc(roc(PKUdata_test$diag, prediction_BSC_FC$MDD,levels = c("HC", "MDD"), direction = "<"))

