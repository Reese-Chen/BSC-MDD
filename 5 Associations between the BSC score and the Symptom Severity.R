#### Associations between the BSC score and the Symptom Severity ####
# written by Yinhan Chen and supervised by Qiang Luo
# Email: yinhanchen23@m.fudan.edu.cn
# released on 23 Sep 2025
# please cite: XXXXXX

library(dplyr)
library(lmerTest)
library(ggplot2)
library(ggsci)           
library(ggExtra)         
library(R.matlab)
library(lavaan)
library(semPlot)
library(readxl)

#############################################
# load data
#############################################

rm(list=ls())
gc()

# BSC data
PKUdata = read.table(file="PKUdata.csv",header = T,sep = ",")
XYdata = read.table(file="XYdata.csv",header = T,sep = ",")

# split HAMD for PKU
PKU_HAMD_data_baseline = read.table("PKUdata_MDD_HAMD_detail.csv",header = T,sep = ",")
PKU_HAMD_data_8week = read.table("PKUdata_MDD_HAMD_8w_detail.csv",header = T,sep = ",")

# merge data
PKUdata$ID = gsub("_0w$", "", PKUdata$ID)

PKUdata = merge(PKUdata,PKU_HAMD_data_baseline,by.x = "ID",by.y = "subID",all.x = T)
PKUdata <- PKUdata %>%
  rename_with(
    ~ paste0(gsub("_", "", .x), "_0w"),  
    matches("^HAMD_\\d+$")               
  )
colnames(PKUdata)[40:42] = c('HAMD_core_0w','HAMD_anxiety_0w','HAMD_sleep_0w')

PKUdata = merge(PKUdata,PKU_HAMD_data_8week,by.x = "ID",by.y = "subID",all.x = T)
colnames(PKUdata)[61:77] = paste0(colnames(PKUdata)[61:77],'_8w')

# HAMD core scores for PKU

PKUdata$HAMD_core_0w = with(PKUdata,HAMD1_0w + HAMD2_0w + HAMD13_0w + HAMD7_0w + HAMD8_0w + HAMD10_0w)
PKUdata$HAMD_anxiety_0w = with(PKUdata,HAMD10_0w + HAMD11_0w+ HAMD12_0w+ HAMD13_0w+ HAMD15_0w+ HAMD17_0w)
PKUdata$HAMD_sleep_0w = with(PKUdata,HAMD4_0w + HAMD5_0w + HAMD6_0w)

PKUdata$HAMD_core_8w = with(PKUdata,HAMD1_8w + HAMD2_8w + HAMD13_8w + HAMD7_8w + HAMD8_8w + HAMD10_8w)
PKUdata$HAMD_anxiety_8w = with(PKUdata,HAMD10_8w + HAMD11_8w+ HAMD12_8w+ HAMD13_8w+ HAMD15_8w+ HAMD17_8w)
PKUdata$HAMD_sleep_8w = with(PKUdata,HAMD4_8w + HAMD5_8w + HAMD6_8w)

PKUdata$HAMD_core_reduction = PKUdata$HAMD_core_0w-PKUdata$HAMD_core_8w
PKUdata$HAMD_anxiety_reduction = PKUdata$HAMD_anxiety_0w-PKUdata$HAMD_anxiety_8w
PKUdata$HAMD_sleep_reduction = PKUdata$HAMD_sleep_0w-PKUdata$HAMD_sleep_8w

# HAMD core scores for XY

XYdata$HAMD_core_0w = with(XYdata,HAMD1_0w + HAMD2_0w + HAMD13_0w + HAMD7_0w + HAMD8_0w + HAMD10_0w)
XYdata$HAMD_anxiety_0w = with(XYdata,HAMD10_0w + HAMD11_0w+ HAMD12_0w+ HAMD13_0w+ HAMD15_0w+ HAMD17_0w)
XYdata$HAMD_sleep_0w = with(XYdata,HAMD4_0w + HAMD5_0w + HAMD6_0w)

XYdata$HAMD_core_8w = with(XYdata,HAMD1_8w + HAMD2_8w + HAMD13_8w + HAMD7_8w + HAMD8_8w + HAMD10_8w)
XYdata$HAMD_anxiety_8w = with(XYdata,HAMD10_8w + HAMD11_8w+ HAMD12_8w+ HAMD13_8w+ HAMD15_8w+ HAMD17_8w)
XYdata$HAMD_sleep_8w = with(XYdata,HAMD4_8w + HAMD5_8w + HAMD6_8w)

XYdata$HAMD_core_reduction = XYdata$HAMD_core_0w-XYdata$HAMD_core_8w
XYdata$HAMD_anxiety_reduction = XYdata$HAMD_anxiety_0w-XYdata$HAMD_anxiety_8w
XYdata$HAMD_sleep_reduction = XYdata$HAMD_sleep_0w-XYdata$HAMD_sleep_8w

###############################################
# descriptive analysis of split HAMD
###############################################

mean(PKUdata$HAMD_core_0w[PKUdata$diag=="MDD"],na.rm=T)
sd(PKUdata$HAMD_core_0w[PKUdata$diag=="MDD"],na.rm=T)
mean(PKUdata$HAMD_sleep_0w[PKUdata$diag=="MDD"],na.rm=T)
sd(PKUdata$HAMD_sleep_0w[PKUdata$diag=="MDD"],na.rm=T)
mean(PKUdata$HAMD_anxiety_0w[PKUdata$diag=="MDD"],na.rm=T)
sd(PKUdata$HAMD_anxiety_0w[PKUdata$diag=="MDD"],na.rm=T)

mean(PKUdata$HAMD_core_8w[PKUdata$diag=="MDD"],na.rm=T)
sd(PKUdata$HAMD_core_8w[PKUdata$diag=="MDD"],na.rm=T)
mean(PKUdata$HAMD_sleep_8w[PKUdata$diag=="MDD"],na.rm=T)
sd(PKUdata$HAMD_sleep_8w[PKUdata$diag=="MDD"],na.rm=T)
mean(PKUdata$HAMD_anxiety_8w[PKUdata$diag=="MDD"],na.rm=T)
sd(PKUdata$HAMD_anxiety_8w[PKUdata$diag=="MDD"],na.rm=T)

mean(XYdata$HAMD_core_0w[XYdata$diag=="MDD"],na.rm=T)
sd(XYdata$HAMD_core_0w[XYdata$diag=="MDD"],na.rm=T)
mean(XYdata$HAMD_sleep_0w[XYdata$diag=="MDD"],na.rm=T)
sd(XYdata$HAMD_sleep_0w[XYdata$diag=="MDD"],na.rm=T)
mean(XYdata$HAMD_anxiety_0w[XYdata$diag=="MDD"],na.rm=T)
sd(XYdata$HAMD_anxiety_0w[XYdata$diag=="MDD"],na.rm=T)

mean(XYdata$HAMD_core_8w[XYdata$diag=="MDD"],na.rm=T)
sd(XYdata$HAMD_core_8w[XYdata$diag=="MDD"],na.rm=T)
mean(XYdata$HAMD_sleep_8w[XYdata$diag=="MDD"],na.rm=T)
sd(XYdata$HAMD_sleep_8w[XYdata$diag=="MDD"],na.rm=T)
mean(XYdata$HAMD_anxiety_8w[XYdata$diag=="MDD"],na.rm=T)
sd(XYdata$HAMD_anxiety_8w[XYdata$diag=="MDD"],na.rm=T)

###############################################
# relationship between HAMD and BSC at baseline
###############################################

HAMD_names = c(colnames(XYdata)[17:34],'HAMD_core_0w','HAMD_anxiety_0w','HAMD_sleep_0w')

# 1. relationship in XY------------------------------------

results = matrix(0,nrow=21,ncol=4)
colnames(results) = c('estimate','std.error','t_value','p')
rownames(results) = HAMD_names

for (i in 1:21){
  fit = summary(lm(as.formula(paste0(HAMD_names[i],'~BSC_0w+age+sex+edu+meanFD')),data = XYdata[which(XYdata$diag=="MDD"),]))
  results[i,1:4] = fit$coefficients[2,]
}
results = as.data.frame(results)
results$p_adjust = NA
results$p_adjust[1:17] = p.adjust(results$p[1:17],'fdr')
write.csv(results,"relationship between baseline BSC and splited HAMD in XY.csv")


# 2. relationship in PKU------------------------------------

results = matrix(0,nrow=21,ncol=4)
colnames(results) = c('estimate','std.error','t_value','p')
rownames(results) = HAMD_names

for (i in 1:21){
  fit = summary(lm(as.formula(paste0(HAMD_names[i],'~BSC_0w+age+sex+edu+meanFD')),data = PKUdata[which(PKUdata$diag=="MDD"),]))
  results[i,1:4] = fit$coefficients[2,]
}
results = as.data.frame(results)
results$p_adjust = NA
results$p_adjust[1:17] = p.adjust(results$p[1:17],'fdr')
write.csv(results,"relationship between baseline BSC and splited HAMD in PKU.csv")

######################################################
# relationship between BSC and HAMD_sleep
######################################################

# 1. XY--------------------------------

t.test(XYdata$HAMD_sleep_0w[XYdata$sex=="female" & XYdata$diag=="MDD"],
       XYdata$HAMD_sleep_0w[XYdata$sex=="male" & XYdata$diag=="MDD"],)

model = lm(HAMD_sleep_0w~BSC_0w+age+sex+edu+meanFD,XYdata[XYdata$diag=="MDD",])
summary(model)

model = lm(HAMD_sleep_0w~BSC_0w+age+edu+meanFD,XYdata[XYdata$diag=="MDD" & XYdata$sex=="female",])
summary(model)

model = lm(HAMD_sleep_0w~BSC_0w+age+edu+meanFD,XYdata[XYdata$diag=="MDD" & XYdata$sex=="male",])
summary(model)


model = lm(HAMD_sleep_8w~BSC_8w+age+sex+edu+meanFD,XYdata[XYdata$diag=="MDD",])
summary(model)
model = lm(HAMD_sleep_reduction~BSC_0w+age+sex+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+age+sex+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+age+sex+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+age+sex+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD" & XYdata$BSC_group_0w=="female-like",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+age+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD" & XYdata$sex=="female",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+age+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD" & XYdata$sex=="male",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+sex+age+edu+meanFD+HAMD_sleep_0w,XYdata[XYdata$diag=="MDD",])
summary(model)

# 2. PKU--------------------------------

model = lm(HAMD_sleep_0w~BSC_0w+age+sex+edu+meanFD,PKUdata[PKUdata$diag=="MDD",])
summary(model)

t.test(PKUdata$HAMD_sleep_0w[PKUdata$sex=="female" & PKUdata$diag=="MDD"],
       PKUdata$HAMD_sleep_0w[PKUdata$sex=="male" & PKUdata$diag=="MDD"],)

model = lm(HAMD_sleep_0w~BSC_0w+age+edu+meanFD,PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="female",])
summary(model)
model = lm(HAMD_sleep_0w~BSC_0w+age+edu+meanFD,PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="male",])
summary(model)

model = lm(HAMD_sleep_8w~BSC_8w+age+sex+edu+meanFD,PKUdata[PKUdata$diag=="MDD",])
summary(model)
model = lm(HAMD_sleep_reduction~BSC_0w+age+sex+edu+meanFD+HAMD_sleep_0w,PKUdata[PKUdata$diag=="MDD",])
summary(model)
model = lm(HAMD_sleep_reduction~BSC_reduction+age+sex+edu+meanFD+HAMD_sleep_0w,PKUdata[PKUdata$diag=="MDD",])
summary(model)

model = lm(HAMD_sleep_reduction~BSC_reduction+age+edu+meanFD+HAMD_sleep_0w,PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="male",])
summary(model)



######################################################
# change of HAMD_sleep at 8week and baseline
######################################################

# 1. XY-----------------------------
XYdata_bl = XYdata[,c(2:8,34,66:68)]
colnames(XYdata_bl) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','HAMD_core','HAMD_anxiety','HAMD_sleep')
XYdata_bl$time = "baseline"
XYdata_fu = XYdata[,c(2:7,53,52,69:71)]
colnames(XYdata_fu) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','HAMD_core','HAMD_anxiety','HAMD_sleep')
XYdata_fu$time = "8week"
XYdata_long = rbind(XYdata_bl[!is.na(XYdata_bl$BSC) & !is.na(XYdata_bl$HAMD_sleep),],
                    XYdata_fu[!is.na(XYdata_fu$BSC) & !is.na(XYdata_fu$HAMD_sleep),])
XYdata_long = XYdata_long[which(XYdata_long$diag=="MDD"),]
XYdata_long$time = as.factor(XYdata_long$time)
XYdata_long$time = relevel(XYdata_long$time, ref = "baseline")

model = lmer(HAMD~time+age+sex+edu+meanFD+(1|ID),XYdata_long[XYdata_long$diag =="MDD",])
summary(model)

model = lmer(HAMD_sleep~time+age+sex+edu+meanFD+(1|ID),XYdata_long[XYdata_long$diag =="MDD",])
summary(model)

model = lmer(HAMD_core~time+age+sex+edu+meanFD+(1|ID),XYdata_long[XYdata_long$diag =="MDD",])
summary(model)

model = lmer(HAMD_anxiety~time+age+sex+edu+meanFD+(1|ID),XYdata_long[XYdata_long$diag =="MDD",])
summary(model)

model = lmer(BSC~time+age+sex+edu+meanFD+(1+time|ID),XYdata_long[XYdata_long$diag =="MDD",])
summary(model)

model = lmer(BSC~time+age+edu+meanFD+(1|ID),XYdata_long[XYdata_long$diag =="MDD" & XYdata_long$sex=="female",])
summary(model)

model = lmer(BSC~time+age+edu+meanFD+(1+time|ID),XYdata_long[XYdata_long$diag =="MDD" & XYdata_long$sex=="male",])
summary(model)



# 2. PKU-----------------------------
PKUdata_bl = PKUdata[,c(1,10,4,3,6,5,12,7,40:42)]
colnames(PKUdata_bl) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','HAMD_core','HAMD_anxiety','HAMD_sleep')
PKUdata_bl$time = "baseline"
PKUdata_fu = PKUdata[,c(1,10,4,3,6,5,27,8,78:80)]
colnames(PKUdata_fu) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','HAMD_core','HAMD_anxiety','HAMD_sleep')
PKUdata_fu$time = "8week"
PKUdata_long = rbind(PKUdata_bl[!is.na(PKUdata_bl$BSC) & !is.na(PKUdata_bl$HAMD_sleep),],
                    PKUdata_fu[!is.na(PKUdata_fu$BSC) & !is.na(PKUdata_fu$HAMD_sleep),])
PKUdata_long = PKUdata_long[which(PKUdata_long$diag=="MDD"),]
PKUdata_long$time = as.factor(PKUdata_long$time)
PKUdata_long$time = relevel(PKUdata_long$time, ref = "baseline")

model = lmer(HAMD~time+age+sex+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD",])
summary(model)

model = lmer(HAMD_sleep~time+age+sex+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD",])
summary(model)

model = lmer(HAMD_core~time+age+sex+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD",])
summary(model)

model = lmer(HAMD_anxiety~time+age+sex+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD",])
summary(model)

model = lmer(BSC~time+age+sex+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD",])
summary(model)

model = lmer(BSC~time+age+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD" & PKUdata_long$sex=="female",])
summary(model)

model = lmer(BSC~time+age+edu+meanFD+(1|ID),PKUdata_long[PKUdata_long$diag =="MDD" & PKUdata_long$sex=="male",])
summary(model)


##############################################
# other sleep-related index in XY
##############################################

convertedIDs = read_excel("converted_IDs_Depression.xlsx",sheet=1)
XYdata = merge(XYdata,convertedIDs,by.x = 'ID',by.y = '匿名化ID',all.x = T)
XYdata$原ID = sub("(\\d{2})(\\d{3})", "ZJDEP-\\1-\\2",XYdata$原ID)

other_index1 = read_excel("baseline_MDD.xlsx",sheet=3)
XYdata_merged = merge(XYdata,other_index1[,c(3,9:11)],by.x = '原ID',by.y = '被试编号',all.x = T)

other_index2 = read_excel("baseline_MDD.xlsx",sheet=4)
XYdata_merged = merge(XYdata_merged,other_index2[,c(3,7)],by.x = '原ID',by.y = '被试编号',all.x = T)

model = lm(表二4您的睡眠如何~BSC_0w+sex+age+edu+meanFD,XYdata_merged[XYdata_merged$diag=="MDD",])
summary(model)

model = lm(是否有失眠~BSC_0w+sex+age+edu+meanFD,XYdata_merged[XYdata_merged$diag=="MDD",])
summary(model)

model = lm(是否睡眠过多~BSC_0w+sex+age+edu+meanFD,XYdata_merged[XYdata_merged$diag=="MDD",])
summary(model)

model = lm(HAMA4失眠~BSC_0w+sex+age+edu+meanFD,XYdata_merged[XYdata_merged$diag=="MDD",])
summary(model)

############################################
# relationship between FC and HAMD sleep
############################################

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
UKB_FC = FC[c(794:8193),]

PKUdata$meanFC = rowMeans(PKU_FC_bl)
XYdata$meanFC = rowMeans(XY_FC_bl)

# 1. relationship in XY--------------------

XY_FC_bl_MDD = XY_FC_bl[which(XYdata$diag=="MDD"),]
p_FC_HAMD_sleep_XY = matrix(0,nrow=dim(XY_FC_bl_MDD)[2],ncol=1)
for (j in 1:dim(XY_FC_bl_MDD)[2]){
  model = lm(HAMD_sleep_0w~XY_FC_bl_MDD[,j]+age+sex+meanFD+edu,XYdata[XYdata$diag=="MDD",])
  p_FC_HAMD_sleep_XY[j,1] = summary(model)$coefficient[2,4]
}
sum(p_FC_HAMD_sleep_XY<0.05)
p_FC_HAMD_sleep_XY_adjust = p.adjust(p_FC_HAMD_sleep_XY,method="fdr")
sum(p_FC_HAMD_sleep_XY_adjust<0.05)

model = lm(HAMD_sleep_0w~meanFC+age+sex+meanFD+edu,XYdata[XYdata$diag=="MDD",])
summary(model)

# 2. relationship in PKU---------------------

PKU_FC_bl_MDD = PKU_FC_bl[which(PKUdata$diag=="MDD"),]
p_FC_HAMD_sleep_PKU = matrix(0,nrow=dim(PKU_FC_bl_MDD)[2],ncol=1)
for (j in 1:dim(PKU_FC_bl_MDD)[2]){
  model = lm(HAMD_sleep_0w~PKU_FC_bl_MDD[,j]+age+sex+meanFD+edu,PKUdata[PKUdata$diag=="MDD",])
  p_FC_HAMD_sleep_PKU[j,1] = summary(model)$coefficient[2,4]
}
sum(p_FC_HAMD_sleep_PKU<0.05)
p_FC_HAMD_sleep_PKU_adjust = p.adjust(p_FC_HAMD_sleep_PKU,method="fdr")
sum(p_FC_HAMD_sleep_PKU_adjust<0.05)

model = lm(HAMD_sleep_0w~meanFC+age+sex+meanFD+edu,PKUdata[PKUdata$diag=="MDD",])
summary(model)


# 3. compare in two cohorts---------------

sum(p_FC_HAMD_sleep_XY<0.05 & p_FC_HAMD_sleep_PKU<0.05)



##################################################
# visualization
##################################################

XYdata$Cohort = 'XY'
PKUdata$Cohort = 'PKU'

plot_data = rbind(XYdata[,c(1,3,4,8,64,68,71,74,77)],
                  PKUdata[,c(1,10,4,12,36,42,80,83,85)])
plot_data=plot_data[plot_data$diag=="MDD",]

# 1. HAMD sleep~BSC---------------------

p1 = ggplot(plot_data,aes(x=BSC_0w,y=HAMD_sleep_0w,group=Cohort,color=Cohort))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic() +
  theme(
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      size = 1.5  
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    #legend.position = "none", 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )+
  theme(legend.position = "none")+
  labs(x = "BSC",
       y = "HAMD-sleep")
p = ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.4)
ggsave("scatterplot_of_BSC_HAMD_sleep.tiff",p, width = 4, height = 4, dpi = 300)


# 2. HAMD sleep reduction~BSC reduction---------------------

p1 = ggplot(plot_data,aes(x=BSC_reduction,y=HAMD_sleep_reduction,group=Cohort,color=Cohort))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic() +
  theme(
    panel.border = element_rect(
      colour = "black", 
      fill = NA, 
      size = 1.5  
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    #legend.position = "none", 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )+
  #theme(legend.position = "none")+
  labs(x = "BSC reduction",
       y = "HAMD-sleep reduction")
p = ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.4)
ggsave("scatterplot_of_BSC_reduction_HAMD_sleep_reduction_label.tiff",p, width = 4, height = 4, dpi = 300)
