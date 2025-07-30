rm(list = ls())
setwd('/media/ubuntu/XYZ/mine/newlgg')

library(org.Hs.eg.db)
library(MetaNeighbor)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)

load(file = 'GSE89567/data/GSE89567.tam.Rdata')
load(file = 'GSE152277/data/GSE152277.tam.Rdata')
load(file = 'GSE202096/data/GSE202096.tam.Rdata')
load(file = 'GSE222520/data/GSE222520.tam.Rdata')
load(file = 'GSE270109/data/GSE270109.tam.Rdata')
load(file = 'GSE227718/data/GSE227718.tam.Rdata')

seurat.datalist <- list(GSE89567.tam,GSE152277.tam,GSE202096.tam,GSE222520.tam,
                 GSE227718.tam,GSE270109.tam)
sce.datalist <- list()
names(seurat.datalist) <- c('GSE89567','GSE152277','GSE202096','GSE222520',
                 'GSE227718','GSE270109')
datas_name <- names(seurat.datalist)
for(i in 1:length(seurat.datalist)){
  meta.data <- seurat.datalist[[i]]@meta.data %>% dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,DF_hi.lo,sample,clusters,celltype)
  seurat.datalist[[i]]@meta.data <- meta.data
  sce.datalist[[i]] <- as.SingleCellExperiment(seurat.datalist[[i]])
  head(rownames(sce.datalist[[i]]))
  colnames(colData(sce.datalist[[i]]))
}
names(sce.datalist) <- datas_name
fused_data <- mergeSCE(sce.datalist)
dim(fused_data)
head(colData(fused_data))
table(fused_data$clusters,fused_data$study_id)
all_global_hvgs <- variableGenes(dat=fused_data,exp_labels=fused_data$study_id)
head(all_global_hvgs)
length(all_global_hvgs)

GSE.aurocs <- MetaNeighborUS(var_genes = all_global_hvgs,
                                   dat = fused_data,
                                   study_id = fused_data$study_id,
                                   cell_type = fused_data$clusters,
                                   fast_version=TRUE)
write.csv(GSE.aurocs,file = 'compare/all.aurocs.csv')


####筛选出5个公共数据中有相似簇的cluster，各个cluster相似性可视化
##cluster1:GSE202096|BMDM.c1	GSE222520|BMDM.c1	GSE227718|BMDM.c1	GSE152277|BMDM.c1	
cluster1 <- GSE.aurocs[rownames(GSE.aurocs) %in% c('GSE152277|BMDM.c1','GSE202096|BMDM.c1',
                                                   'GSE222520|BMDM.c1','GSE227718|BMDM.c1'),
                       colnames(GSE.aurocs) %in% c('GSE152277|BMDM.c1','GSE202096|BMDM.c1',
                                                   'GSE222520|BMDM.c1','GSE227718|BMDM.c1')]

df_long <- melt(cluster1)
colnames(df_long) <- c('Cluster1','Cluster2','Similarity')
df_half <- df_long[c(1:4,6:8,11:12,16), ]

# 绘制热图
ggplot(df_half, aes(x = Cluster1, y = Cluster2, fill = Similarity)) +
  geom_tile() + # 绘制色块
  geom_text(aes(label = round(Similarity, 2)), color = "white") + # 在色块上添加标签
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) + # 设置颜色梯度从0到1
  theme_minimal() + # 使用简约主题
  labs(x = NULL, y = NULL, fill = "Similarity") + # 添加标签和图例标题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 横坐标标签倾斜45°
    axis.text.y = element_text(hjust = 1), # 纵坐标标签靠右对齐
    axis.title.x = element_blank(), # 移除x轴标题
    axis.title.y = element_blank(), # 移除y轴标题
    panel.grid = element_blank(), # 移除网格线
    panel.border = element_blank() # 移除面板边框
  ) +
  coord_fixed() # 使坐标轴比例固定，以确保热图方形
ggsave("compare/heatmap/cluster1_heatmap.png", width = 8, height = 6, dpi = 300, bg = 'white')

##cluster2:GSE89567|MG.c1	GSE152277|MG.c1 GSE270109|MG.c1	
cluster2 <- GSE.aurocs[rownames(GSE.aurocs) %in% c('GSE89567|MG.c1','GSE152277|MG.c1','GSE270109|MG.c1'),
                       colnames(GSE.aurocs) %in% c('GSE89567|MG.c1','GSE152277|MG.c1','GSE270109|MG.c1')]

df_long <- melt(cluster2)
colnames(df_long) <- c('Cluster1','Cluster2','Similarity')
df_half <- df_long[c(1:3,5:6,9), ]

# 绘制热图
ggplot(df_half, aes(x = Cluster1, y = Cluster2, fill = Similarity)) +
  geom_tile() + # 绘制色块
  geom_text(aes(label = round(Similarity, 2)), color = "white") + # 在色块上添加标签
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) + # 设置颜色梯度从0到1
  theme_minimal() + # 使用简约主题
  labs(x = NULL, y = NULL, fill = "Similarity") + # 添加标签和图例标题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 横坐标标签倾斜45°
    axis.text.y = element_text(hjust = 1), # 纵坐标标签靠右对齐
    axis.title.x = element_blank(), # 移除x轴标题
    axis.title.y = element_blank(), # 移除y轴标题
    panel.grid = element_blank(), # 移除网格线
    panel.border = element_blank() # 移除面板边框
  ) +
  coord_fixed() # 使坐标轴比例固定，以确保热图方形
ggsave("compare/heatmap/cluster2_heatmap.png", width = 8, height = 6, dpi = 300, bg = 'white')

##cluster3:GSE89567|MG.c2	GSE222520|MG.c3	GSE270109|MG.c3
cluster3 <- GSE.aurocs[rownames(GSE.aurocs) %in% c('GSE89567|MG.c2','GSE222520|MG.c3','GSE270109|MG.c3'),
                       colnames(GSE.aurocs) %in% c('GSE89567|MG.c2','GSE222520|MG.c3','GSE270109|MG.c3')]

df_long <- melt(cluster3)
colnames(df_long) <- c('Cluster1','Cluster2','Similarity')
df_half <- df_long[c(1:3,5:6,9), ]

# 绘制热图
ggplot(df_half, aes(x = Cluster1, y = Cluster2, fill = Similarity)) +
  geom_tile() + # 绘制色块
  geom_text(aes(label = round(Similarity, 2)), color = "white") + # 在色块上添加标签
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) + # 设置颜色梯度从0到1
  theme_minimal() + # 使用简约主题
  labs(x = NULL, y = NULL, fill = "Similarity") + # 添加标签和图例标题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 横坐标标签倾斜45°
    axis.text.y = element_text(hjust = 1), # 纵坐标标签靠右对齐
    axis.title.x = element_blank(), # 移除x轴标题
    axis.title.y = element_blank(), # 移除y轴标题
    panel.grid = element_blank(), # 移除网格线
    panel.border = element_blank() # 移除面板边框
  ) +
  coord_fixed() # 使坐标轴比例固定，以确保热图方形
ggsave("compare/heatmap/cluster3_heatmap.png", width = 8, height = 6, dpi = 300, bg = 'white')

##cluster4:GSE152277|MG.c3	GSE202096|MG.c1
cluster4 <- GSE.aurocs[rownames(GSE.aurocs) %in% c('GSE152277|MG.c3','GSE202096|MG.c1'),
                       colnames(GSE.aurocs) %in% c('GSE152277|MG.c3','GSE202096|MG.c1')]

df_long <- melt(cluster4)
colnames(df_long) <- c('Cluster1','Cluster2','Similarity')
df_half <- df_long[c(1:2,4), ]

# 绘制热图
ggplot(df_half, aes(x = Cluster1, y = Cluster2, fill = Similarity)) +
  geom_tile() + # 绘制色块
  geom_text(aes(label = round(Similarity, 2)), color = "white") + # 在色块上添加标签
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) + # 设置颜色梯度从0到1
  theme_minimal() + # 使用简约主题
  labs(x = NULL, y = NULL, fill = "Similarity") + # 添加标签和图例标题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 横坐标标签倾斜45°
    axis.text.y = element_text(hjust = 1), # 纵坐标标签靠右对齐
    axis.title.x = element_blank(), # 移除x轴标题
    axis.title.y = element_blank(), # 移除y轴标题
    panel.grid = element_blank(), # 移除网格线
    panel.border = element_blank() # 移除面板边框
  ) +
  coord_fixed() # 使坐标轴比例固定，以确保热图方形
ggsave("compare/heatmap/cluster4_heatmap.png", width = 8, height = 6, dpi = 300, bg = 'white')

##cluster5:GSE152277|MG.c2	GSE222520|MG.c2
cluster5 <- GSE.aurocs[rownames(GSE.aurocs) %in% c('GSE152277|MG.c2','GSE222520|MG.c2'),
                       colnames(GSE.aurocs) %in% c('GSE152277|MG.c2','GSE222520|MG.c2')]

df_long <- melt(cluster5)
colnames(df_long) <- c('Cluster1','Cluster2','Similarity')
df_half <- df_long[c(1:2,4), ]

# 绘制热图
ggplot(df_half, aes(x = Cluster1, y = Cluster2, fill = Similarity)) +
  geom_tile() + # 绘制色块
  geom_text(aes(label = round(Similarity, 2)), color = "white") + # 在色块上添加标签
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) + # 设置颜色梯度从0到1
  theme_minimal() + # 使用简约主题
  labs(x = NULL, y = NULL, fill = "Similarity") + # 添加标签和图例标题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 横坐标标签倾斜45°
    axis.text.y = element_text(hjust = 1), # 纵坐标标签靠右对齐
    axis.title.x = element_blank(), # 移除x轴标题
    axis.title.y = element_blank(), # 移除y轴标题
    panel.grid = element_blank(), # 移除网格线
    panel.border = element_blank() # 移除面板边框
  ) +
  coord_fixed() # 使坐标轴比例固定，以确保热图方形
ggsave("compare/heatmap/cluster5_heatmap.png", width = 8, height = 6, dpi = 300, bg = 'white')




