setwd('E:/mine/newlgg')
setwd('/media/ubuntu/047E-B974//mine/newlgg')


library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(enrichplot)

load(file = 'GSE89567/data/GSE89567.tam.Rdata')
load(file = 'GSE152277/data/GSE152277.tam.Rdata')
load(file = 'GSE202096/data/GSE202096.tam.Rdata')
load(file = 'GSE222520/data/GSE222520.tam.Rdata')
load(file = 'GSE270109/data/GSE270109.tam.Rdata')
load(file = 'GSE227718/data/GSE227718.tam.Rdata')


######提取差异基因的函数
##提取tam中unique cluster
get_cluster_TAM <- function(seurat_obj, cell_column = 'celltype', cluster_column = 'clusters'){
  meta.data <- seurat_obj@meta.data
  unique_clusters <- unique(meta.data[[cluster_column]])
  return(unique_clusters)
}

##提取差异基因列表
extract_DEGs <- function(seurat_obj, cluster_id, tam_clusters){
  comparison_clusters <- setdiff(tam_clusters,cluster_id)
  #查找差异基因
  Idents(seurat_obj) <- 'clusters'
  degs <- FindMarkers(seurat_obj, ident.1 = cluster_id, ident.2 = comparison_clusters, min.pct = 0.25, logfc.threshold = 0.25)
  return(degs)
}

tam_clusters_GSE89567 <- get_cluster_TAM(GSE89567.tam)
tam_clusters_GSE152277 <- get_cluster_TAM(GSE152277.tam)
tam_clusters_GSE202096 <- get_cluster_TAM(GSE202096.tam)
tam_clusters_GSE222520 <- get_cluster_TAM(GSE222520.tam)
tam_clusters_GSE227718 <- get_cluster_TAM(GSE227718.tam)
tam_clusters_GSE270109 <- get_cluster_TAM(GSE270109.tam)

########提取每个类群的各cluster的差异基因
#########cluster1 GSE152277_BMDM.c1 GSE202096_BMDM.c1 GSE222520_BMDM.c1 GSE227718_BMDM.c1
deg_GSE152277_BMDM.c1 <-extract_DEGs(GSE152277.tam,'BMDM.c1',tam_clusters_GSE152277)
deg_GSE202096_BMDM.c1 <-extract_DEGs(GSE202096.tam,'BMDM.c1',tam_clusters_GSE202096)
deg_GSE222520_BMDM.c1 <-extract_DEGs(GSE222520.tam,'BMDM.c1',tam_clusters_GSE222520)
deg_GSE227718_BMDM.c1 <-extract_DEGs(GSE227718.tam,'BMDM.c1',tam_clusters_GSE227718)
#########cluster2  GSE89567_MG.c1 GSE152277_MG.c1 GSE270109_MG.c1
deg_GSE89567_MG.c1 <- extract_DEGs(GSE89567.tam, 'MG.c1', tam_clusters_GSE89567)
deg_GSE152277_MG.c1 <- extract_DEGs(GSE152277.tam, 'MG.c1', tam_clusters_GSE152277)
deg_GSE270109_MG.c1 <- extract_DEGs(GSE270109.tam, 'MG.c1', tam_clusters_GSE270109)
######cluster3  GSE89567_MG.c2 GSE2225520_MG.c3 GSE270109_MG.c3
deg_GSE89567_MG.c2 <- extract_DEGs(GSE89567.tam, 'MG.c2', tam_clusters_GSE89567)
deg_GSE222520_MG.c3 <- extract_DEGs(GSE222520.tam, 'MG.c3', tam_clusters_GSE222520)
deg_GSE270109_MG.c3 <- extract_DEGs(GSE270109.tam, 'MG.c3', tam_clusters_GSE270109)
######cluster4  GSE152277_MG.c3 GSE202096_MG.c1
deg_GSE152277_MG.c3 <- extract_DEGs(GSE152277.tam, 'MG.c3', tam_clusters_GSE152277)
deg_GSE202096_MG.c1 <- extract_DEGs(GSE202096.tam, 'MG.c1', tam_clusters_GSE202096)
######cluster5  GSE152277_MG.c2 GSE222520_MG.c2
deg_GSE152277_MG.c2 <- extract_DEGs(GSE152277.tam, 'MG.c2', tam_clusters_GSE152277)
deg_GSE222520_MG.c2 <- extract_DEGs(GSE222520.tam, 'MG.c2', tam_clusters_GSE222520)

########绘制所有差异基因火山图
###差异基因火山图
deglist <- list(GSE152277_BMDM.c1 = deg_GSE152277_BMDM.c1,GSE202096_BMDM.c1=deg_GSE202096_BMDM.c1,
                GSE222520_BMDM.c1=deg_GSE222520_BMDM.c1,GSE227718_BMDM.c1=deg_GSE227718_BMDM.c1,
                GSE89567_MG.c1=deg_GSE89567_MG.c1,GSE152277_MG.c1=deg_GSE152277_MG.c1,GSE270109_MG.c1=deg_GSE270109_MG.c1,
                GSE89567_MG.c2=deg_GSE89567_MG.c2,GSE222520_MG.c3=deg_GSE222520_MG.c3,GSE270109_MG.c3=deg_GSE270109_MG.c3,
                GSE152277_MG.c3=deg_GSE152277_MG.c3,GSE202096_MG.c1=deg_GSE202096_MG.c1,
                GSE152277_MG.c2=deg_GSE152277_MG.c2,GSE222520_MG.c2=deg_GSE222520_MG.c2)
saveRDS(deglist,'Test/DEG/data/deglist.rds')
deglist <- readRDS('Test/DEG/data/deglist.rds')
source("VolcanoPlot.R")
for (i in 1:14){
  
  print(i)
  print(names(deglist)[i])
  
  deg <- deglist[[i]]
  deg.dif = data.frame(symbol = rownames(deg),
                       log2FoldChange=deg$avg_log2FC,
                       padj=deg$p_val_adj
  )
  png(filename = paste('Test/DEG/VolcanoPlot/',names(deglist)[i],'.png', sep = ''),height = 600,width = 800)
  VolcanoPlot(deg.dif, padj=0.05, title="cluster vs others in the same dataset", label.max = 50)+
    xlim(-4,4)
  dev.off()
}

#######筛选显著上调的基因avg_log2FC > 0.584 & p_val_adj<0.05,差异基因取交集和并集各保存一组DEGS
#cluster1
reup_deg_GSE152277_BMDM.c1 <- subset(deg_GSE152277_BMDM.c1,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE202096_BMDM.c1 <- subset(deg_GSE202096_BMDM.c1,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE222520_BMDM.c1 <- subset(deg_GSE222520_BMDM.c1,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE227718_BMDM.c1 <- subset(deg_GSE227718_BMDM.c1,avg_log2FC > 0 & p_val_adj<0.05)
#cluster2
reup_deg_GSE89567_MG.c1 <- subset(deg_GSE89567_MG.c1,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE152277_MG.c1 <- subset(deg_GSE152277_MG.c1,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE270109_MG.c1 <- subset(deg_GSE270109_MG.c1,avg_log2FC > 0 & p_val_adj<0.05)
#cluster3
reup_deg_GSE89567_MG.c2 <- subset(deg_GSE89567_MG.c2,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE222520_MG.c3 <- subset(deg_GSE222520_MG.c3,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE270109_MG.c3 <- subset(deg_GSE270109_MG.c3,avg_log2FC > 0 & p_val_adj<0.05)
#cluster4
reup_deg_GSE152277_MG.c3 <- subset(deg_GSE152277_MG.c3,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE202096_MG.c1 <- subset(deg_GSE202096_MG.c1,avg_log2FC > 0 & p_val_adj<0.05)
#cluster5
reup_deg_GSE152277_MG.c2 <- subset(deg_GSE152277_MG.c2,avg_log2FC > 0 & p_val_adj<0.05)
reup_deg_GSE222520_MG.c2 <- subset(deg_GSE222520_MG.c2,avg_log2FC > 0 & p_val_adj<0.05)

####把显著的存到列表中
reup_deglist <- list(GSE152277_BMDM.c1 = reup_deg_GSE152277_BMDM.c1,GSE202096_BMDM.c1=reup_deg_GSE202096_BMDM.c1,
                     GSE222520_BMDM.c1=reup_deg_GSE222520_BMDM.c1,GSE227718_BMDM.c1=reup_deg_GSE227718_BMDM.c1,
                     GSE89567_MG.c1=reup_deg_GSE89567_MG.c1,GSE152277_MG.c1=reup_deg_GSE152277_MG.c1,GSE270109_MG.c1=reup_deg_GSE270109_MG.c1,
                     GSE89567_MG.c2=reup_deg_GSE89567_MG.c2,GSE222520_MG.c3=reup_deg_GSE222520_MG.c3,GSE270109_MG.c3=reup_deg_GSE270109_MG.c3,
                     GSE152277_MG.c3=reup_deg_GSE152277_MG.c3,GSE202096_MG.c1=reup_deg_GSE202096_MG.c1,
                     GSE152277_MG.c2=reup_deg_GSE152277_MG.c2,GSE222520_MG.c2=reup_deg_GSE222520_MG.c2)
saveRDS(reup_deglist,'Test/DEG/data/reup_deglist.rds')
reup_deglist = readRDS('Test/DEG/data/reup_deglist.rds')


#cluster1
degs_1_1 <- rownames(reup_deg_GSE152277_BMDM.c1)
degs_1_2 <- rownames(reup_deg_GSE202096_BMDM.c1)
degs_1_3 <- rownames(reup_deg_GSE222520_BMDM.c1)
degs_1_4 <- rownames(reup_deg_GSE227718_BMDM.c1)
#cluster2
degs_2_1 <- rownames(reup_deg_GSE89567_MG.c1)
degs_2_2 <- rownames(reup_deg_GSE152277_MG.c1)
degs_2_3 <- rownames(reup_deg_GSE270109_MG.c1)
#cluster3
degs_3_1 <- rownames(reup_deg_GSE89567_MG.c2)
degs_3_2 <- rownames(reup_deg_GSE222520_MG.c3)
degs_3_3 <- rownames(reup_deg_GSE270109_MG.c3)
#cluster4
degs_4_1 <- rownames(reup_deg_GSE152277_MG.c3)
degs_4_2 <- rownames(reup_deg_GSE202096_MG.c1)
#cluster5
degs_5_1 <- rownames(reup_deg_GSE152277_MG.c2)
degs_5_2 <- rownames(reup_deg_GSE222520_MG.c2)


#交集marker
com_degs_cluster1 <- intersect(intersect(intersect(degs_1_1,degs_1_2),degs_1_3),degs_1_4) #各cluster差异基因取交集
com_degs_cluster2 <- intersect(intersect(degs_2_1,degs_2_2),degs_2_3) #各cluster差异基因取交集
com_degs_cluster3 <- intersect(intersect(degs_3_1,degs_3_2),degs_3_3) #各cluster差异基因取交集
com_degs_cluster4 <- intersect(degs_4_1,degs_4_2) #各cluster差异基因取交集
com_degs_cluster5 <- intersect(degs_5_1,degs_5_2) #各cluster差异基因取交集

########差异基因取交集
com_degs_up0 <- list(cluster1 = com_degs_cluster1, cluster2 = com_degs_cluster2, cluster3 = com_degs_cluster3,
                     cluster4 = com_degs_cluster4,cluster5 = com_degs_cluster5)
saveRDS(com_degs_up0, file = 'Test/DEG/data/com_degs_up0.rds')
com_degs <- com_degs_up0
# com_degs <- readRDS('Test/DEG/data/com_degs.rds')

# #并集marker
# union_degs_cluster1 <- union(union(union(degs_1_1,degs_1_2),degs_1_3),degs_1_4) #各cluster差异基因取并集
# union_degs_cluster2 <- union(union(degs_2_1,degs_2_2),degs_2_3) #各cluster差异基因取并集
# union_degs_cluster3 <- union(union(degs_3_1,degs_3_2),degs_3_3) #各cluster差异基因取并集
# union_degs_cluster4 <- union(degs_4_1,degs_4_2) #各cluster差异基因取并集
# union_degs_cluster5 <- union(degs_5_1,degs_5_2) #各cluster差异基因取并集
# ########差异基因取并集
# union_degs <- list(cluster1 = union_degs_cluster1, cluster2 = union_degs_cluster2, cluster3 = union_degs_cluster3,
#                    cluster4 = union_degs_cluster4, cluster5 = union_degs_cluster5)
# saveRDS(union_degs, file = 'DEG/data/union_degs.rds')
# # union_degs <- readRDS('DEG/data/union_degs.rds')


#####用差异基因做GO分析
#degs <- readRDS('Test/DEG/data/union_degs.rds')
degs <- readRDS('Test/DEG/data/com_degs_up0.rds')
indices=c(1,2,4,5)
degs=degs[indices]
comdeg_go_results <- list()
for(cluster in names(degs)){
  print(cluster)
  degs <- degs[[cluster]]
  go <- enrichGO(gene = degs,
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 readable = TRUE)
  
  comdeg_go_results[[cluster]] <- go@result
}
saveRDS(comdeg_go_results,file="Test/DEG/common/GO/goresult.rds")
comdeg_go_results <- readRDS("Test/DEG/common/GO/goresult.rds")
#显著的
sign_go_result <- list()
for (cluster in names(degs)){
  print(cluster)
  sign_go_result[[cluster]] <- comdeg_go_results[[cluster]][comdeg_go_results[[cluster]]$pvalue < 0.05,]
}
saveRDS(sign_go_result,file="Test/DEG/common/GO/sign_go_result.rds")
#存cytoscape
for(i in 1:length(sign_go_result)){
  print(i)
  go <- sign_go_result[[i]]
  go_term<-go[,c('ID','Description','p.adjust','qvalue','geneID','Count')]
  for(j in 1:nrow(go_term)){
    genes <-go_term[j,5]
    genes<-gsub('/',',',genes)
    go_term[j,5]<-genes
  }
  Count <- go_term[,6] # 示例数据，根据实际情况替换
  geneID <- go_term[,5]
  go_term[,5] <- Count
  go_term[,6] <-geneID
  colnames(go_term)<-c('Term','Description','P-value','Adjusted P-value','Gene Count','Gene list')
  write.table(go_term,file=paste('Test/DEG/common/GO/cytoscape/go_single_cluster',i,'_enrichmentmap.txt',sep=''), sep='\t',row.names = FALSE,quote=FALSE)
}

