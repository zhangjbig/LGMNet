rm(list = ls())
setwd('/media/ubuntu/XYZ/mine/data/LGG')

library(Seurat)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(harmony)
library(ggplot2)
library(clustree)
library(copykat)


seuratlist<-list()
h5file <- c("GSE270109/GSM8334604_zhongliu_2_lukai_processed/filtered_feature_bc_matrix.h5",
         "GSE270109/GSM8334605_zhongliu_1_xianghongmei_processed/filtered_feature_bc_matrix.h5",
         "GSE270109/GSM8334606_zhongliu_2_xianghongmei_processed/filtered_feature_bc_matrix.h5")
samples_name <- c('GSM8334604','GSM8334605','GSM8334606')

for(i in 1:length(samples_name)){
  seurat_data <- Read10X_h5(file = h5file[i])
  seuratlist[[i]] <- CreateSeuratObject(counts = seurat_data,
                                   project = samples_name[i],
                                   min.features = 200,
                                   min.cells = 3)
  seuratlist[[i]]<-RenameCells(seuratlist[[i]],add.cell.id=samples_name[i])
  seuratlist[[i]][['percent.mt']]<-PercentageFeatureSet(seuratlist[[i]],pattern="^MT-")
}
names(seuratlist)<-samples_name


#####去除双细胞
library(DoubletFinder)
for(sample_name in samples_name){
  print(sample_name)
  data<-seuratlist[[sample_name]]
  data<-NormalizeData(data)
  data<-FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  data <- RunPCA(data, features = VariableFeatures(object = data))
  data <- FindNeighbors(data, dims = 1:15)
  data <- FindClusters(data, resolution = 0.7)
  data <- RunUMAP(data, dims = 1:15)
  sweep.res.list_myloid_seurat <- paramSweep(data, PCs = 1:15, sct = FALSE)
  #head(sweep.res.list_myloid_seurat)
  sweep.stats_myloid_seurat <- summarizeSweep(sweep.res.list_myloid_seurat, GT = FALSE)
  bcmvn_myloid_seurat <- find.pK(sweep.stats_myloid_seurat) #可以看到最佳参数的点
  ## 所以最佳的参数是：
  mpK<-as.numeric(as.vector(bcmvn_myloid_seurat$pK[which.max(bcmvn_myloid_seurat$BCmetric)]))
  
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(data)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
  DoubletRate = 0.042104
  #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
  nExp_poi <- round(DoubletRate*length(data$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
  # 计算双细胞比例
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  data <- doubletFinder(data, PCs = 1:15, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  data <- doubletFinder(data, PCs = 1:15, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  data@meta.data["DF_hi.lo"] <- data@meta.data[8]
  data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet" & data@meta.data[10] == "Singlet")] <- "Doublet-Low Confidience"
  data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
  meta.data<-data@meta.data %>% select(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,seurat_clusters,DF_hi.lo)
  data@meta.data<-meta.data
  seuratlist[[sample_name]]<-data
}

for(i in 1:length(seuratlist)){
  seuratlist[[i]]$sample<-samples_name[i]
}
scRNA<-merge(seuratlist[[1]],seuratlist[2:length(seuratlist)])
scRNA$orig.ident<-'GSE270109'
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GSE270109.qc<-subset(scRNA,subset=nFeature_RNA>300&nFeature_RNA<7000&percent.mt<15&DF_hi.lo=='Singlet')
save(GSE270109.qc,file = 'GSE270109/data/GSE270109.df.Rdata')
load('GSE270109/data/GSE270109.df.Rdata')

####seurat标准流程
GSE270109.qc <- NormalizeData(GSE270109.qc)
GSE270109.qc <- FindVariableFeatures(GSE270109.qc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE270109.qc), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE270109.qc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE270109.qc)
GSE270109.qc <- ScaleData(GSE270109.qc, features = all.genes)
GSE270109.qc <- RunPCA(GSE270109.qc, features = VariableFeatures(object = GSE270109.qc))
## harmony
GSE270109.qc <- GSE270109.qc %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE270109.qc, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE270109.qc, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE270109.qc, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE270109.qc,reduction = 'harmony')
GSE270109.qc <- RunUMAP(GSE270109.qc, dims = 1:15,reduction = 'harmony')
# GSE270109.qc <- RunTSNE(GSE270109.qc, dims = 1:20,reduction = 'harmony')
GSE270109.qc <- FindNeighbors(GSE270109.qc, dims = 1:15,reduction = 'harmony')
GSE270109.qc <- FindClusters(GSE270109.qc, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE270109.qc) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE270109.qc$clusters<-GSE270109.qc$RNA_snn_res.0.05
GSE270109.qc$celltype<- GSE270109.qc$RNA_snn_res.0.05
png(filename = 'GSE270109/marker/umap_res.0.05.png',width = '600',height = '600')
DimPlot(GSE270109.qc,group.by = 'RNA_snn_res.0.05',label = T)
dev.off()

save(GSE270109.qc, file = 'GSE270109/data/GSE270109.qc.Rdata')
#load(file = 'GSE270109/data/GSE270109.qc.Rdata')

#copykat
setwd('/media/ubuntu/XYZ/mine/data/LGG/GSE270109/copykat')
seurat_list<-list()
samples<-names(table(GSE270109.qc$sample))
for(i in 1:length(samples)){
  seurat_list[[i]]<-subset(GSE270109.qc,sample==samples[i])
}
copykat.test_list<-list()
for(i in 1:length(samples)){
  seurat<-seurat_list[[i]]
  print(seurat$orig.ident[1])
  print(i)
  counts <-  GetAssayData(seurat, slot = 'counts')
  copykat.test_list[[i]] <- copykat(rawmat=counts, 
                                    id.type="S", 
                                    cell.line="no", 
                                    ngene.chr=5, 
                                    win.size=25, 
                                    KS.cut=0.15, 
                                    sam.name=samples[i], 
                                    distance="euclidean", 
                                    n.cores=48)
  pred.test <- data.frame(copykat.test_list[[i]]$prediction)
  seurat$copykat.pred<-pred.test$copykat.pred
  seurat_list[[i]]<-seurat
}
sc<-merge(seurat_list[[1]],seurat_list[2:length(samples)])
GSE270109.qc$copykat.pred<-sc$copykat.pred
setwd('/media/ubuntu/XYZ/mine/data/LGG')
save(copykat.test_list,file=paste0('GSE270109/copykat/copykat.all.result.Rdata'))
save(GSE270109.qc,file = 'GSE270109/data/GSE270109.qc.Rdata')


#####singleR
#准备工作
library(celldex)
library(SingleR)
library(scater)
library(viridis)
library(ggsci)
DimPlot(GSE270109.qc, reduction = "umap", group.by = "orig.ident")  
DimPlot(GSE270109.qc, reduction = "umap",label.size = 5,label = T,pt.size = 0.5)
table(GSE270109.qc$clusters)
a <- table(GSE270109.qc$orig.ident,GSE270109.qc@active.ident)
gplots::balloonplot(a,main="harmony")
saveRDS(GSE270109.qc,'GSE270109/data/GSE270109.qc_adj.rds')
#      SingleR classType
hpca.se <- celldex::HumanPrimaryCellAtlasData()
GSE270109.rds <-readRDS("GSE270109/data/GSE270109.qc_adj.rds")
clusters=GSE270109.rds@meta.data$clusters
pred.hesc <- SingleR(GSE270109.rds@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
GSE270109.rds@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
png(filename = 'GSE270109/marker/SingleR_png.png',width = '600',height = '600')
DimPlot(GSE270109.rds, group.by = 'singleR',pt.size=0.1,label = T,reduction = 'umap')
dev.off()

#####scibet
library(scibet)
exp_tpm <- calculateTPM(GSE270109.rds@assays$RNA@counts)
exp_tpm <- t(exp_tpm)
model <- readr::read_csv("major_human_cell_types.csv")
# model <-model[,-1]
model <- pro.core(model)
prd <- LoadModel(model)
scibet.type <- prd(exp_tpm)
GSE270109.rds$scibet <- scibet.type
png("GSE270109/marker/Scibet_png.png",width = 600, height = 600)
DimPlot(GSE270109.rds, group.by = "scibet", pt.size=0.1, label = T, reduction = 'umap')
dev.off()

celltype1 <- data.frame(
  Barcodes=colnames(GSE270109.rds@assays$RNA),
  ClusterID=GSE270109.rds$clusters,
  celltype=scibet.type, 
  stringsAsFactors = F
)
celltype_scibet<-celltype1[!duplicated(celltype1[,2:3]),2:3]
# saveRDS(GSE270109/data/GSE270109.rds,"GSE270109.qc_adj.rds")


##   FindAllMarkers 
DEGS<-FindAllMarkers(GSE270109.rds,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE270109/marker/Heatmap.png",width = 900, height = 900)
DoHeatmap(GSE270109.rds, features = top5$gene)
dev.off()
DotPlot(GSE270109.rds,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


GSE270109.qc <- GSE270109.rds
#####细胞类型确认
png(filename = 'GSE270109/marker/immune.png',width =350,height = 350 )
FeaturePlot(GSE270109.qc,features =c('PTPRC'),pt.size=0.5,ncol = 1,max.cutoff = 1)
dev.off()
png(filename = 'GSE270109/marker/TAM.png',width =1200,height = 350 )
FeaturePlot(GSE270109.qc,features =c('CD68','CD14','AIF1','C1QB'),pt.size=0.5,ncol = 4,max.cutoff = 1)
dev.off()
png(filename = 'GSE270109/marker/Tumor.png',width =900,height = 350 )
FeaturePlot(GSE270109.qc,features =c('FABP7','SOX2','EGFR'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE270109/marker/BMDM.png',width =1200,height = 350 )
FeaturePlot(GSE270109.qc,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 2)
dev.off()
png(filename = 'GSE270109/marker/MG.png',width =900,height = 350 )
FeaturePlot(GSE270109.qc,features =c('TMEM119','CH25H','P2RY12'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE270109/marker/T cell.png',width =900,height = 350 )
FeaturePlot(GSE270109.qc,features =c('CD3D','CD3E','CD3G'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE270109/marker/B cell.png',width =900,height = 350 )
FeaturePlot(GSE270109.qc,features =c('IGKC','CD79A','JCHAIN'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE270109/marker/NK cell.png',width =900,height = 350 )
FeaturePlot(GSE270109.qc,features =c('NKG7','GNLY','GZMB'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
new.cluster.ids <- c('not immume','not immume','not immume','not immume','TAM','not immume',
                     'not immume','not immume','not immume','T cell','not immume','not immume','TAM')
levels(GSE270109.qc@meta.data$celltype) <- new.cluster.ids#将celltype确定
DimPlot(GSE270109.qc,group.by = 'celltype',label = T)
save(GSE270109.qc,file = 'GSE270109/data/GSE270109.qc.Rdata')
# load(file = 'GSE270109/data/GSE270109.qc.Rdata')


#####GSVA验证
library(GSVA)
load('curatedMarkers.Rdata')
Idents(GSE270109.qc)<-'clusters'
expr<-AverageExpression(GSE270109.qc,assays = 'RNA',slot='data')[[1]]    
#expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr[rowSums(expr)>0,])
gsva.scdata<-gsva(expr,curatedMarkers,verbose=T, parallel.sz=10,min.sz=3)



####TAM单独处理
load(file = 'GSE270109/data/GSE270109.qc.Rdata')
GSE270109.tam <- subset(GSE270109.qc, celltype == 'TAM')
GSE270109.tam <- subset(GSE270109.tam, singleR == 'Macrophage')
GSE270109.tam <- subset(GSE270109.tam, scibet == 'Macrophage'|scibet=='Microglia')
GSE270109.tam <- NormalizeData(GSE270109.tam)
GSE270109.tam <- FindVariableFeatures(GSE270109.tam, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE270109.tam), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE270109.tam)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE270109.tam)
GSE270109.tam <- ScaleData(GSE270109.tam, features = all.genes)
GSE270109.tam <- RunPCA(GSE270109.tam, features = VariableFeatures(object = GSE270109.tam))
## harmony
GSE270109.tam <- GSE270109.tam %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE270109.tam, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE270109.tam, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE270109.tam, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE270109.tam,reduction = 'harmony')
GSE270109.tam <- RunUMAP(GSE270109.tam, dims = 1:15,reduction = 'harmony')
# GSE270109.tam <- RunTSNE(GSE270109.tam, dims = 1:20,reduction = 'harmony')
GSE270109.tam <- FindNeighbors(GSE270109.tam, dims = 1:15,reduction = 'harmony')
GSE270109.tam <- FindClusters(GSE270109.tam, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE270109.tam) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE270109.tam$clusters<-GSE270109.tam$RNA_snn_res.0.05
GSE270109.tam$celltype<-GSE270109.tam$RNA_snn_res.0.05
png(filename = 'GSE270109/marker/BMDM-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE270109.tam,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 2)
dev.off()
png(filename = 'GSE270109/marker/MG-TAM.png',width =900,height = 350 )
FeaturePlot(GSE270109.tam,features =c('TMEM119','CH25H','P2RY12'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
DimPlot(GSE270109.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)

##   FindAllMarkers 
Idents(GSE270109.tam) <- 'clusters'
DEGS<-FindAllMarkers(GSE270109.tam,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE270109/marker/TAMHeatmap.png",width = 900, height = 900)
DoHeatmap(GSE270109.tam, features = top5$gene)
dev.off()
DotPlot(GSE270109.tam,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

new.cluster.ids <- c('MG.c1','MG.c2','MG.c3','MG.c4')
levels(GSE270109.tam$clusters) <- new.cluster.ids
GSE270109.tam$celltype <- 'MG-TAM'
png("GSE270109/marker/TAM-celltype-clusters.png",width = 1200, height = 600)
DimPlot(GSE270109.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
dev.off()
save(GSE270109.tam, file = 'GSE270109/data/GSE270109.tam.Rdata')


GSE270109.MG2 <- subset(GSE270109.tam, clusters == 'MG.c2')
GSE270109.MG2 <- NormalizeData(GSE270109.MG2)
GSE270109.MG2 <- FindVariableFeatures(GSE270109.MG2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE270109.MG2), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE270109.MG2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE270109.MG2)
GSE270109.MG2 <- ScaleData(GSE270109.MG2, features = all.genes)
GSE270109.MG2 <- RunPCA(GSE270109.MG2, features = VariableFeatures(object = GSE270109.MG2))
## harmony
GSE270109.MG2 <- GSE270109.MG2 %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE270109.MG2, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE270109.MG2, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE270109.MG2, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE270109.MG2,reduction = 'harmony')
GSE270109.MG2 <- RunUMAP(GSE270109.MG2, dims = 1:20,reduction = 'harmony')
# GSE270109.MG2 <- RunTSNE(GSE270109.MG2, dims = 1:20,reduction = 'harmony')
GSE270109.MG2 <- FindNeighbors(GSE270109.MG2, dims = 1:20,reduction = 'harmony')
GSE270109.MG2 <- FindClusters(GSE270109.MG2, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE270109.MG2) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE270109.MG2$clusters <- GSE270109.MG2$RNA_snn_res.0.7
GSE270109.MG2$celltype <- GSE270109.MG2$RNA_snn_res.0.7

DimPlot(GSE270109.MG2,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
Idents(GSE270109.MG2) <- 'clusters'
DEGS<-FindAllMarkers(GSE270109.MG2,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
DotPlot(GSE270109.MG2,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

png(filename = 'GSE270109/marker/BMDM-MG2.png',width =900,height = 350 )
FeaturePlot(GSE270109.MG2,features =c('LSP1','VIM','LGALS3'),pt.size=0.5,ncol = 3)
dev.off()
png(filename = 'GSE270109/marker/MG-MG2.png',width =900,height = 350 )
FeaturePlot(GSE270109.MG2,features =c('HSPA1B','CXCL8','SPP1'),pt.size=0.5,ncol = 3)
dev.off()

new.cluster.ids <- c('MG.c2','MG.c3','BMDM.c1')
levels(GSE270109.MG2$clusters) <- new.cluster.ids
new.celltype.ids <- c('MG-TAM','MG-TAM','BMDM-TAM')
levels(GSE270109.MG2$celltype) <- new.celltype.ids

MG2.cell.ids <- Cells(GSE270109.MG2)
MG2.clusters <- GSE270109.MG2$clusters
MG2.celltypes <- GSE270109.MG2$celltype

GSE270109.tam.test <- GSE270109.tam
GSE270109.tam.test$cellname <- colnames(GSE270109.tam.test)
GSE270109.tam.test$clusters <- as.character(GSE270109.tam.test$clusters) 
GSE270109.tam.test$celltype <- as.character(GSE270109.tam.test$celltype)
GSE270109.tam.test$clusters[GSE270109.tam.test$cellname %in% MG2.cell.ids]<-MG2.clusters
GSE270109.tam.test$celltype[GSE270109.tam.test$cellname %in% MG2.cell.ids]<-MG2.celltypes

GSE270109.tam.test$clusters <- as.character(GSE270109.tam.test$clusters)
GSE270109.tam.test$clusters[GSE270109.tam.test$clusters == '1'] <- 'MG.c2'
GSE270109.tam.test$clusters[GSE270109.tam.test$clusters == '2'] <- 'MG.c3'
GSE270109.tam.test$clusters[GSE270109.tam.test$clusters == '3'] <- 'BMDM.c1'

GSE270109.tam.test$celltype <- as.character(GSE270109.tam.test$celltype)
GSE270109.tam.test$celltype[GSE270109.tam.test$celltype == '1'] <- 'MG-TAM'
GSE270109.tam.test$celltype[GSE270109.tam.test$celltype == '2'] <- 'BMDM-TAM'

meta.data <- GSE270109.tam.test@meta.data %>% dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,
                                                     seurat_clusters,DF_hi.lo,sample,RNA_snn_res.0.05,
                                                     RNA_snn_res.0.1,RNA_snn_res.0.2,RNA_snn_res.0.3,RNA_snn_res.0.4,
                                                     RNA_snn_res.0.5,RNA_snn_res.0.7,RNA_snn_res.0.9,clusters,
                                                     celltype,singleR,scibet)
GSE270109.tam.test@meta.data <- meta.data

save(GSE270109.tam,file = 'GSE270109/data/GSE270109.tam.test.Rdata')
