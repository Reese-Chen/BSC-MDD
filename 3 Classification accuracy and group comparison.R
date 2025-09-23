#### Classification accuracy and group comparision####
# written by Yinhan Chen and supervised by Qiang Luo
# Email: yinhanchen23@m.fudan.edu.cn
# released on 23 Sep 2025
# please cite: XXXXXX

library(ggplot2)
library(ggsignif)
library(pheatmap)
library(reshape2)
library(boot)
library(dplyr)

rm(list=ls())
gc()

#############################################
# load data
#############################################

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

sum((PKUdata$diag=="HC") & (PKUdata$BSC_0w>=0.35) & (PKUdata$BSC_0w<=0.65))/sum(PKUdata$diag=="HC")
sum((XYdata$diag=="HC") & (XYdata$BSC_0w>=0.35) & (XYdata$BSC_0w<=0.65))/sum(XYdata$diag=="HC")

model = lm(BSC_0w~diag+age+edu+meanFD,data = PKUdata[which(PKUdata$sex=="female"),])
summary(model)
model = lm(BSC_0w~diag+age+edu+meanFD,data = XYdata[which(XYdata$sex=="female"),])
summary(model)
model = lm(BSC_8w~diag+age+sex+edu+meanFD,data = PKUdata)
summary(model)

############################################
# classification accuracy
############################################

accuracy_function = function(data, indices) {
  sample_data = data[indices, ]
  mean(sample_data$sex == sample_data$predictedsex)
}

set.seed(123)  

# 1. PKU------------------------------------------------------------------------

PKUdata$predictedsex = "female"
PKUdata$predictedsex[which(PKUdata$BSC_0w>0.5)] = "male"
PKUdata$predictedsex = as.factor(PKUdata$predictedsex)

class_acc_PKU = matrix(0,nrow=4,ncol=3)
count =0
for (sex in c('female','male')){
  for (diag in c('HC','MDD')){
    count  = count +1
    boot_result = boot(data=PKUdata[which(PKUdata$diag==diag & PKUdata$sex==sex),],
                       statistic = accuracy_function,
                       R = 1000)
    x = boot.ci(boot_result, type = "perc")
    class_acc_PKU[count,] = c((x$percent[4]+x$percent[5])/2,x$percent[4],x$percent[5])
  }
}

# 2. XY------------------------------------------------------------------------

XYdata$predictedsex = "female"
XYdata$predictedsex[which(XYdata$BSC_0w>0.5)] = "male"
XYdata$predictedsex = as.factor(XYdata$predictedsex)

class_acc_XY = matrix(0,nrow=4,ncol=3)
count =0
for (sex in c('female','male')){
  for (diag in c('HC','MDD')){
    count  = count +1
    boot_result = boot(data=XYdata[which(XYdata$diag==diag & XYdata$sex==sex),],
                       statistic = accuracy_function,
                       R = 1000)
    x = boot.ci(boot_result, type = "perc")
    class_acc_XY[count,] = c((x$percent[4]+x$percent[5])/2,x$percent[4],x$percent[5])
  }
}

############################################
# compare the BSC in HC and MDD
############################################

# 1. XY-------------------------------------------------------------------------

model = lm(BSC_0w~diag+age+sex+edu+meanFD,data = XYdata)
summary(model)

model1 = lm(BSC_0w~diag+edu+age+meanFD,data = XYdata[which(XYdata$sex=='female'),])
summary(model1)

model1 = lm(BSC_0w~diag+edu+age+meanFD,data = XYdata[which(XYdata$sex=='male'),])
summary(model1)

t.test(XYdata$BSC_0w[which(XYdata$sex=="female" & XYdata$diag=="HC")],XYdata$BSC_0w[which(XYdata$sex=="female" & XYdata$diag=="MDD")],paired = FALSE)
t.test(XYdata$BSC_0w[which(XYdata$sex=="male" & XYdata$diag=="HC")],XYdata$BSC_0w[which(XYdata$sex=="male" & XYdata$diag=="MDD")],paired = FALSE)

model = lm(BSC_0w~diag*sex+age+edu+meanFD,data = XYdata)
summary(model)

# 2. PKU------------------------------------------------------------------------

model = lm(BSC_0w~diag+age+sex+edu+meanFD,data = PKUdata)
summary(model)
confint(model)

model1 = lm(BSC_0w~diag+edu+age+meanFD,data = PKUdata[which(PKUdata$sex=='female'),])
summary(model1)

model1 = lm(BSC_0w~diag+edu+age+meanFD,data = PKUdata[which(PKUdata$sex=='male'),])
summary(model1)

t.test(PKUdata$BSC_0w[which(PKUdata$sex=="female" & PKUdata$diag=="HC")],PKUdata$BSC_0w[which(PKUdata$sex=="female" & PKUdata$diag=="MDD")],paired = FALSE)
t.test(PKUdata$BSC_0w[which(PKUdata$sex=="male" & PKUdata$diag=="HC")],PKUdata$BSC_0w[which(PKUdata$sex=="male" & PKUdata$diag=="MDD")],paired = FALSE)

model = lm(BSC_0w~diag*sex+age+edu+meanFD,data = PKUdata)
summary(model)

############################################
# classification proportion
############################################

# 1. PKU------------------------------------------------------------------------

result1 = matrix(0,nrow = 12, ncol = 5)
rownames(result1) = rep(c('HC','MDD'),each=6)
colnames(result1) = c('sex','BSC_group','proportion','lower','upper')
result1[,1] = rep(c('female','male'),2,each = 3)
result1[,2] = rep(c('female-like','androgynous','male-like'),4)
count = 0
for (diagnosis in c('HC','MDD')){
  for (sex in c('female','male')){
    for (group in c('female-like','androgynous','male-like')){
      count = count+1
      datax = PKUdata$BSC_group_0w[which(PKUdata$diag==diagnosis & PKUdata$sex==sex)]
      confid = binom.test(sum(datax == group), length(datax), conf.level = 0.95)
      result1[count,3:5] = c(mean(datax==group),confid$conf.int[1],confid$conf.int[2]) 
    }
  }
}
result1

table_BSC_group = table(PKUdata$BSC_group_0w[which(PKUdata$diag=='MDD')],PKUdata$sex[which(PKUdata$diag=='MDD')])
chi2_test = chisq.test(t(table_BSC_group))
chi2_test

# 2. XY------------------------------------------------------------------------

result1 = matrix(0,nrow = 12, ncol = 5)
rownames(result1) = rep(c('HC','MDD'),each=6)
colnames(result1) = c('sex','BSC_group','proportion','lower','upper')
result1[,1] = rep(c('female','male'),2,each = 3)
result1[,2] = rep(c('female-like','androgynous','male-like'),4)
count = 0
for (diagnosis in c('HC','MDD')){
  for (sex in c('female','male')){
    for (group in c('female-like','androgynous','male-like')){
      count = count+1
      datax = XYdata$BSC_group_0w[which(XYdata$diag==diagnosis & XYdata$sex==sex)]
      confid = binom.test(sum(datax == group), length(datax), conf.level = 0.95)
      result1[count,3:5] = c(mean(datax==group),confid$conf.int[1],confid$conf.int[2]) 
    }
  }
}
result1

table_BSC_group = table(XYdata$BSC_group_0w[which(XYdata$diag=='MDD')],XYdata$sex[which(XYdata$diag=='MDD')])
chi2_test = chisq.test(t(table_BSC_group))
chi2_test

################################################################################
# visualization of distribution
################################################################################

# 1. histgram and density plot of BSC in PKU----------------------------
p_PKU_female = 
ggplot(data = PKUdata[which(PKUdata$sex=="female"),], aes(x = BSC_0w,fill = diag,color=diag)) +
  geom_density(aes(y = ..density..*40), size = 1.2, alpha = 0.2) +
  geom_histogram(aes(y = (..count.. / sum(..count..)) * 100), binwidth = 0.15, alpha = 0.2, position = "identity",) +
  labs(
    x = "BSC",
    y = "Percentage (%)",
    title = "PKU Female"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0) 
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0) 
  ) +
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
  )

ggsave("Distribution of BSC in PKU Female.png", p_PKU_female, width = 4, height = 4, dpi = 300)


# 2. histgram and density plot of BSC in PKU----------------------------
p_PKU_male = 
  ggplot(data = PKUdata[which(PKUdata$sex=="male"),], aes(x = BSC_0w,fill = diag,color=diag)) +
  geom_density(aes(y = ..density..*40), size = 1.2, alpha = 0.2) +
  geom_histogram(aes(y = (..count.. / sum(..count..)) * 100), binwidth = 0.15, alpha = 0.2, position = "identity",) +
  labs(
    x = "BSC",
    y = "Value",
    title = "PKU Male"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0) 
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0) 
  ) +
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
  )

ggsave("Distribution of BSC in PKU Male.png", p_PKU_male, width = 4, height = 4, dpi = 300)

# 3. histgram and density plot of BSC in XY----------------------------
p_XY_female = 
  ggplot(data = XYdata[which(XYdata$sex=="female"),], aes(x = BSC_0w,fill = diag,color=diag)) +
  geom_density(aes(y = ..density..*40), size = 1.2, alpha = 0.2) +
  geom_histogram(aes(y = (..count.. / sum(..count..)) * 100), binwidth = 0.15, alpha = 0.2, position = "identity",) +
  labs(
    x = "BSC",
    y = "Value",
    title = "XY Female"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0) 
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0) 
  ) +
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
  )

ggsave("Distribution of BSC in XY Female.png", p_XY_female, width = 4, height = 4, dpi = 300)


# 4. histgram and density plot of BSC in XY----------------------------
p_XY_male = 
  ggplot(data = XYdata[which(XYdata$sex=="male"),], aes(x = BSC_0w,fill = diag,color=diag)) +
  geom_density(aes(y = ..density..*40), size = 1.2, alpha = 0.2) +
  geom_histogram(aes(y = (..count.. / sum(..count..)) * 100), binwidth = 0.15, alpha = 0.2, position = "identity",) +
  labs(
    x = "BSC",
    y = "Value",
    title = "XY Male"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0) 
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0) 
  ) +
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
  )

ggsave("Distribution of BSC in XY Male.png", p_XY_male, width = 4, height = 4, dpi = 300)

################################################################################
# visualization of change of BSC
################################################################################

# long data form
PKUdata_bl = PKUdata[,c(1,10,4,3,6,5,12,7,13)]
colnames(PKUdata_bl) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','BSC_group')
PKUdata_bl$time = "Baseline"
PKUdata_fu = PKUdata[,c(1,10,4,3,6,5,27,8,28)]
colnames(PKUdata_fu) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','BSC_group')
PKUdata_fu$time = "8-week"
PKUdata_long = rbind(PKUdata_bl[!is.na(PKUdata_bl$BSC),],
                     PKUdata_fu[!is.na(PKUdata_fu$BSC),])
PKUdata_long = PKUdata_long[which(PKUdata_long$diag=="MDD"),]
PKUdata_long$time = as.factor(PKUdata_long$time)
PKUdata_long$time = relevel(PKUdata_long$time, ref = "Baseline")

XYdata_bl = XYdata[,c(2:8,34,9)]
colnames(XYdata_bl) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','BSC_group')
XYdata_bl$time = "Baseline"
XYdata_fu = XYdata[,c(2:7,53,52,54)]
colnames(XYdata_fu) = c('ID','diag','sex','age','edu','meanFD','BSC','HAMD','BSC_group')
XYdata_fu$time = "8-week"
XYdata_long = rbind(XYdata_bl[!is.na(XYdata_bl$BSC),],
                    XYdata_fu[!is.na(XYdata_fu$BSC),])
XYdata_long = XYdata_long[which(XYdata_long$diag=="MDD"),]
XYdata_long$time = as.factor(XYdata_long$time)
XYdata_long$time = relevel(XYdata_long$time, ref = "Baseline")


# 1. PKU--------------------------

p1 = ggplot(PKUdata_long, aes(x = time, y = BSC, group = ID, color = as.factor(BSC_group))) +
  geom_point() +
  geom_line() +
  labs(x = "Time", y = "BSC",title='PKU') +
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
    legend.position = "none", 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )
ggsave("change of BSC at baseline and 8week PKU.tiff",p1, width = 5, height = 5, dpi = 300)


# 2. XY---------------------------

p2 = ggplot(XYdata_long, aes(x = time, y = BSC, group = ID, color = as.factor(BSC_group))) +
  geom_point() +
  geom_line() +
  labs(x = "Time", y = "BSC",title='XY') +
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
    legend.position = "none", 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )
ggsave("change of BSC at baseline and 8week XY.tiff",p2, width = 5, height = 5, dpi = 300)

############################################
# analyes for split BSC
############################################

# 1. compare between MDD and HC--------------------

# in XY

# baseline
result = matrix(0,nrow = length(BSC_names),ncol=4)
colnames(result) = c("Estimate","Std. Error", "t value", "Pr(>|t|)" )
for (i in 1:length(BSC_names)){
  model = lm(as.formula(paste0(BSC_names[i],'~diag+age+sex+edu+meanFD')),data = XYdata_bl)
  fit = summary(model)
  result[i,] = fit$coefficients[2,]
}

# in PKU

# baseline
result = matrix(0,nrow = length(BSC_names),ncol=4)
colnames(result) = c("Estimate","Std. Error", "t value", "Pr(>|t|)" )
for (i in 1:length(BSC_names)){
  model = lm(as.formula(paste0(BSC_names[i],'~diag+age+sex+edu+meanFD')),data = PKUdata_bl)
  fit = summary(model)
  result[i,] = fit$coefficients[2,]
}

# 8-week
result = matrix(0,nrow = length(BSC_names),ncol=4)
colnames(result) = c("Estimate","Std. Error", "t value", "Pr(>|t|)" )
for (i in 1:length(BSC_names)){
  model = lm(as.formula(paste0(BSC_names[i],'~diag+age+sex+edu+meanFD')),data = PKUdata_fu)
  fit = summary(model)
  result[i,] = fit$coefficients[2,]
}

# 2. visualization-----------------------

# in XY male/female
BSC_means = rbind(
  colMeans(XYdata[XYdata$diag=="HC" & XYdata$sex=="female",c(8,10:16)]),
  colMeans(XYdata[XYdata$diag=="MDD" & XYdata$sex=="female",c(8,10:16)]),
  colMeans(XYdata[XYdata$diag=="MDD" & XYdata$sex=="female",c(53,55:61)],na.rm=T)
)
BSC_means = as.data.frame(BSC_means)
BSC_means$diag = c('HC','MDD_baseline','MDD_8week')
BSC_means$diag = factor(BSC_means$diag, levels = c('HC', 'MDD_baseline', 'MDD_8week'))
BSC_long = melt(BSC_means, id.vars = "diag", variable.name = "Network", value.name = "Mean")
p_XY = ggplot(BSC_long, aes(x = Network, y = diag, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "red",limits = c(0.4, 0.80)) +
  labs(y = "ZIB-MDD female")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in XY baseline and 8week female.png", p_XY, width = 10, height = 3, dpi = 300)

# XY MDD-HC female

BSC_means = rbind(
  colMeans(XYdata[XYdata$diag=="HC" & XYdata$sex=="female",c(8,10:16)]),
  colMeans(XYdata[XYdata$diag=="MDD" & XYdata$sex=="female",c(8,10:16)]),
  colMeans(XYdata[XYdata$diag=="MDD" & XYdata$sex=="female",c(53,55:61)],na.rm=T)
)
BSC_means = as.data.frame(BSC_means)
BSC_means$diag = c('HC','MDD_baseline','MDD_8week')
BSC_means$diag = factor(BSC_means$diag, levels = c('HC', 'MDD_baseline', 'MDD_8week'))

BSC_MDD_HC = BSC_means
BSC_MDD_HC[2,] = BSC_MDD_HC[2,]-BSC_MDD_HC[1,]
BSC_MDD_HC[3,] = BSC_MDD_HC[3,]-BSC_MDD_HC[1,]
BSC_MDD_HC = BSC_MDD_HC[c(2:3),]
BSC_MDD_HC$diag = c('Baseline','8-week')
BSC_long = melt(BSC_MDD_HC, id.vars = "diag", variable.name = "Network", value.name = "Mean")
p_XY = ggplot(BSC_long, aes(x = Network, y = diag, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient2(
    low = "blue",      
    mid = "white",    
    high = "red",      
    midpoint = 0,     
    limits = c(-0.1, 0.2)  
  )+
  labs(y = "XY Female")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in XY baseline and 8week female MDD minus HC.png", p_XY, width = 10, height = 3, dpi = 300)

# XY MDD-HC male

BSC_means = rbind(
  colMeans(XYdata[XYdata$diag=="HC" & XYdata$sex=="male",c(8,10:16)]),
  colMeans(XYdata[XYdata$diag=="MDD" & XYdata$sex=="male",c(8,10:16)]),
  colMeans(XYdata[XYdata$diag=="MDD" & XYdata$sex=="male",c(53,55:61)],na.rm=T)
)
BSC_means = as.data.frame(BSC_means)
BSC_means$diag = c('HC','MDD_baseline','MDD_8week')
BSC_means$diag = factor(BSC_means$diag, levels = c('HC', 'MDD_baseline', 'MDD_8week'))

BSC_MDD_HC = BSC_means
BSC_MDD_HC[2,] = BSC_MDD_HC[2,]-BSC_MDD_HC[1,]
BSC_MDD_HC[3,] = BSC_MDD_HC[3,]-BSC_MDD_HC[1,]
BSC_MDD_HC = BSC_MDD_HC[c(2:3),]
BSC_MDD_HC$diag = c('Baseline','8-week')
BSC_long = melt(BSC_MDD_HC, id.vars = "diag", variable.name = "Network", value.name = "Mean")
p_XY = ggplot(BSC_long, aes(x = Network, y = diag, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient2(
    low = "blue",      
    mid = "white",    
    high = "red",      
    midpoint = 0,     
    limits = c(-0.1, 0.2)  
  )+
  labs(y = "XY Male")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in XY baseline and 8week male MDD minus HC.png", p_XY, width = 10, height = 3, dpi = 300)


# in PKU female/male
BSC_means = rbind(
  colMeans(PKUdata[PKUdata$diag=="HC" & PKUdata$sex=="female",c(12,14:20)]),
  colMeans(PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="female",c(12,14:20)]),
  colMeans(PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="female",c(27,29:35)],na.rm=T)
)
BSC_means = as.data.frame(BSC_means)
BSC_means$diag = c('HC','MDD_baseline','MDD_8week')
BSC_means$diag = factor(BSC_means$diag, levels = c('HC', 'MDD_baseline', 'MDD_8week'))
BSC_long = melt(BSC_means, id.vars = "diag", variable.name = "Network", value.name = "Mean")
p_PKU = ggplot(BSC_long, aes(x = Network, y = diag, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "red",limits = c(0.40, 0.80)) +
  labs(y = "PKU female")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in PKU baseline and 8week female.png", p_PKU, width = 10, height = 3, dpi = 300)


# PKU MDD-HC female

BSC_means = rbind(
  colMeans(PKUdata[PKUdata$diag=="HC" & PKUdata$sex=="female",c(12,14:20)]),
  colMeans(PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="female",c(12,14:20)]),
  colMeans(PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="female",c(27,29:35)],na.rm=T)
)
BSC_means = as.data.frame(BSC_means)
BSC_means$diag = c('HC','MDD_baseline','MDD_8week')
BSC_means$diag = factor(BSC_means$diag, levels = c('HC', 'MDD_baseline', 'MDD_8week'))

BSC_MDD_HC = BSC_means
BSC_MDD_HC[2,] = BSC_MDD_HC[2,]-BSC_MDD_HC[1,]
BSC_MDD_HC[3,] = BSC_MDD_HC[3,]-BSC_MDD_HC[1,]
BSC_MDD_HC = BSC_MDD_HC[c(2:3),]
BSC_MDD_HC$diag = c('Baseline','8-week')
BSC_long = melt(BSC_MDD_HC, id.vars = "diag", variable.name = "Network", value.name = "Mean")
p_XY = ggplot(BSC_long, aes(x = Network, y = diag, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient2(
    low = "blue",      
    mid = "white",    
    high = "red",      
    midpoint = 0,      
    limits = c(-0.1, 0.2)  
  )+
  labs(y = "PKU Female")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in PKU baseline and 8week female MDD minus HC.png", p_XY, width = 10, height = 3, dpi = 300)

# PKU MDD-HC male

BSC_means = rbind(
  colMeans(PKUdata[PKUdata$diag=="HC" & PKUdata$sex=="male",c(12,14:20)]),
  colMeans(PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="male",c(12,14:20)]),
  colMeans(PKUdata[PKUdata$diag=="MDD" & PKUdata$sex=="male",c(27,29:35)],na.rm=T)
)
BSC_means = as.data.frame(BSC_means)
BSC_means$diag = c('HC', 'MDD_baseline', 'MDD_8week')
BSC_means$diag = factor(BSC_means$diag, levels = c('HC', 'MDD_baseline', 'MDD_8week'))

BSC_MDD_HC = BSC_means
BSC_MDD_HC[2,] = BSC_MDD_HC[2,]-BSC_MDD_HC[1,]
BSC_MDD_HC[3,] = BSC_MDD_HC[3,]-BSC_MDD_HC[1,]
BSC_MDD_HC = BSC_MDD_HC[c(2:3),]
BSC_MDD_HC$diag = c('Baseline','8-week')
BSC_long = melt(BSC_MDD_HC, id.vars = "diag", variable.name = "Network", value.name = "Mean")
p_XY = ggplot(BSC_long, aes(x = Network, y = diag, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient2(
    low = "blue",      
    mid = "white",     
    high = "red",      
    midpoint = 0,      
    limits = c(-0.1, 0.2)  
  )+
  labs(y = "PKU Male")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in PKU baseline and 8week male MDD minus HC.png", p_XY, width = 10, height = 3, dpi = 300)

############################################
# splitted BSC at baseline for remission and nonremission
############################################

# in XY
XY_BSC = XYdata[which(!is.na(XYdata$remission)), c(8,10:16,63)]
BSC_means = aggregate(. ~ remission, data = XY_BSC, FUN = mean, na.rm = TRUE)
BSC_long = melt(BSC_means, id.vars = "remission", variable.name = "Network", value.name = "Mean")
p_XY = ggplot(BSC_long, aes(x = Network, y = remission, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "red",limits = c(0.45, 0.68)) +
  labs(y = "ZIB-MDD")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in XY for remission and nonremission.png", p_XY, width = 10, height = 3, dpi = 300)


# in PKU
PKU_BSC = PKUdata[which(!is.na(PKUdata$remission)), c(12,14:21)]
BSC_means = aggregate(. ~ remission, data = PKU_BSC, FUN = mean, na.rm = TRUE)
BSC_long = melt(BSC_means, id.vars = "remission", variable.name = "Network", value.name = "Mean")
p_PKU = ggplot(BSC_long, aes(x = Network, y = remission, fill = Mean)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradient(low = "white", high = "red",limits = c(0.45, 0.68)) +
  labs(y = "PKU")+ 
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c(
    "BSC_0w" = "BSC",
    "BSC_between_network_0w" = "BSC between network",
    "BSC_Default_mode_0w" = "BSC Default mode",
    "BSC_Dorsal_attention_0w" = "BSC Dorsal attention",
    "BSC_Language_0w" = "BSC Language",
    "BSC_Salience_0w" = "BSC Salience",
    "BSC_Sensorimotor_0w" = "BSC Sensorimotor",
    "BSC_Visual_0w" = "BSC Visual"
  ))
ggsave("Heatmap of split BSC in PKU remission and nonremission.png", p_PKU, width = 10, height = 3, dpi = 300)




