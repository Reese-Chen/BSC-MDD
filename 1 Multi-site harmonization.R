#### Multi-site harmonization ####
# written by Yinhan Chen and supervised by Qiang Luo
# Email: yinhanchen23@m.fudan.edu.cn
# released on 23 Sep 2025
# please cite: XXXXXX

library(R.matlab)
library(sva)
library(umap)
library(ggplot2)
library(forcats)

rm(list=ls())
gc()

#############################################
# load data
#############################################

FC_data = readMat("alldata_for_combat_analysis.mat",dat)

FC = FC_data$dat
batch = as.vector(FC_data$batch)
mod = FC_data$mod

############################################
# UMAP for original rsFC data
############################################

batch = as.factor(batch)
set.seed(123)
umap_result = umap(FC)

umap_data = as.data.frame(umap_result$layout)
colnames(umap_data) = c("UMAP1", "UMAP2")
umap_data$Dataset = batch
umap_data$Dataset = as.factor(umap_data$Dataset)
umap_data$Dataset = fct_recode(umap_data$Dataset,PKU = "1",XY = "2",UKB="3")

dataset_colors = c("PKU" = "#A6170A", "XY" = "#0079C9", "UKB" = "#A3A3A3")
umap_data$shape = ifelse(umap_data$Dataset == "UKB", 21, 16)


p1 =  ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(color = Dataset, shape = shape), 
                 alpha = 0.7, size = 2, stroke = 1) +
      scale_color_manual(values = dataset_colors) +
      scale_shape_identity() +  
      labs(x = "UMAP Dimension 1",
           y = "UMAP Dimension 2") + 
      theme_classic(base_size = 14) +
      theme(
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24, face = "bold")
      )

ggsave('UMAP_for_FC_before_combat.png',p1,width=10,height = 8)

#############################################
# harmonization and visualization
#############################################

# 1. combat-----------------------------
combat_mdata=ComBat(dat=t(FC),batch=batch,mod = mod,par.prior=T,prior.plot=F,ref.batch=3)

umap_result_after = umap(t(combat_mdata))

umap_data_after = as.data.frame(umap_result_after$layout)
colnames(umap_data_after) = c("UMAP1", "UMAP2")
umap_data_after$Dataset = batch
umap_data_after$Dataset = as.factor(umap_data_after$Dataset)
umap_data_after$Dataset = fct_recode(umap_data_after$Dataset,PKU = "1",XY = "2",UKB="3")

# 2. visualization with umap-----------------------
dataset_colors = c("PKU" = "#A6170A", "XY" = "#0079C9", "UKB" = "#A3A3A3")
umap_data_after$shape = ifelse(umap_data_after$Dataset == "UKB", 21, 16)


p2 =  ggplot(umap_data_after, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Dataset, shape = shape), 
             alpha = 0.7, size = 2, stroke = 1) +
  scale_color_manual(values = dataset_colors) +
  scale_shape_identity() +  
  labs(x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") + 
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24, face = "bold")
  )

ggsave('UMAP_for_FC_after_combat.png',p2,width=10,height = 8)

# 3. save results-------------------------
writeMat("PKU_FC_bl_after_SVA_combat.mat", PKU_FC_bl = combat_mdata[,c(1:168)])
writeMat("PKU_FC_fu_after_SVA_combat.mat", PKU_FC_fu = combat_mdata[,c(169:258)])
writeMat("XY_FC_bl_after_SVA_combat.mat", XY_FC_bl = combat_mdata[,c(259:672)])
writeMat("XY_FC_fu_after_SVA_combat.mat", XY_FC_fu = combat_mdata[,c(673:793)])

