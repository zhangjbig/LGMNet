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

#GSE202096  sample x 1
assays<-dir('GSE202096/')
dir<-paste0('GSE202096/',assays)
seuratlist<-list()
samples_name<-assays
i <-1
seurat_data<- read.table(gzfile("GSE102130_RAW.txt.gz"), row.names = 1, header = TRUE, sep = "\t")

# 使用CreateSeuratObject()函数创建Seurat对象，并在此处指定项目名称
seuratlist[[i]] <- CreateSeuratObject(counts = seurat_data,
                                 min.features = 200,
                                 min.cells = 3,
                                 project = "102130")
for(i in 1:length(samples_name)){
  seurat_data <- Read10X(data.dir = dir[i])
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
scRNA$orig.ident<-'GSE202096'
Idents(scRNA) <- 'orig.ident'
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GSE202096.qc<-subset(scRNA,subset=nFeature_RNA>300&nFeature_RNA<7000&percent.mt<20&DF_hi.lo=='Singlet')
save(GSE202096.qc,file = 'data/GSE202096.df.Rdata')
load('GSE230389/data/GSE230389.df.Rdata')

####seurat标准流程
GSE237779.qc <- NormalizeData(GSE237779.qc)
GSE237779.qc <- FindVariableFeatures(GSE237779.qc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE237779.qc), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE237779.qc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE237779.qc)
GSE237779.qc <- ScaleData(GSE237779.qc, features = all.genes)
GSE237779.qc <- RunPCA(GSE237779.qc, features = VariableFeatures(object = GSE237779.qc))
## harmony
GSE237779.qc <- GSE237779.qc %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE237779.qc, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE237779.qc,pt.size = 1)
DimHeatmap(GSE237779.qc, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE237779.qc)
GSE237779.qc <- RunUMAP(GSE237779.qc, dims = 1:15)
# GSE237779.qc <- RunTSNE(GSE237779.qc, dims = 1:20,reduction = 'harmony')
GSE237779.qc <- FindNeighbors(GSE237779.qc, dims = 1:15)
GSE237779.qc <- FindClusters(GSE237779.qc, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE237779.qc) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE237779.qc$clusters<-GSE237779.qc$RNA_snn_res.0.1
GSE237779.qc$celltype<- GSE237779.qc$RNA_snn_res.0.1 
png(filename = 'GSE237779/marker/umap_res.0.2.png',width = '600',height = '600')
DimPlot(GSE237779.qc,group.by = 'SCT_snn_res.0.1',label = T)
dev.off()

save(GSE230389.qc, file = 'GSE230389/data/GSE230389.qc.Rdata')
#load(file = 'GSE230389/data/GSE230389.qc.Rdata')


#####singleR
#准备工作
library(celldex)
library(SingleR)
library(scater)
library(viridis)
library(ggsci)
DimPlot(GSE237779.qc, reduction = "umap", group.by = "orig.ident")  
DimPlot(GSE237779.qc, reduction = "umap",label.size = 5,label = T,pt.size = 0.5)
table(GSE237779.qc$clusters)
a <- table(GSE237779.qc$orig.ident,GSE237779.qc@active.ident)
gplots::balloonplot(a,main="harmony")
saveRDS(GSE237779.qc,'GSE237779/data/GSE237779.qc_adj.rds')
#      SingleR classType
hpca.se <- celldex::HumanPrimaryCellAtlasData()
GSE237779.rds <-readRDS("GSE237779/data/GSE237779.qc_adj.rds")
GSE237779.rds <- GSE237779.qc
clusters=GSE237779.rds@meta.data$clusters
pred.hesc <- SingleR(GSE237779.rds@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
GSE237779.rds@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
png(filename = 'GSE237779/marker/SingleR_png.png',width = '600',height = '600')
DimPlot(GSE237779.rds, group.by = 'singleR',pt.size=0.1,label = T,reduction = 'umap')
dev.off()

#####scibet
library(scibet)
exp_tpm <- calculateTPM(GSE230389.rds@assays$RNA@counts)
exp_tpm <- t(exp_tpm)
model <- readr::read_csv("major_human_cell_types.csv")
# model <-model[,-1]
model <- pro.core(model)
prd <- LoadModel(model)
scibet.type <- prd(exp_tpm)
GSE230389.rds$scibet <- scibet.type
png("GSE230389/marker/Scibet_png.png",width = 600, height = 600)
DimPlot(GSE230389.rds, group.by = "scibet", pt.size=0.1, label = T, reduction = 'umap')
dev.off()

celltype1 <- data.frame(
  Barcodes=colnames(GSE178116.rds@assays$RNA),
  ClusterID=GSE178116.rds$clusters,
  celltype=scibet.type, 
  stringsAsFactors = F
)
celltype_scibet<-celltype1[!duplicated(celltype1[,2:3]),2:3]
# saveRDS(GSE178116/data/GSE178116.rds,"GSE178116.qc_adj.rds")


##   FindAllMarkers 
DEGS<-FindAllMarkers(GSE178116.rds,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE178116/marker/Heatmap.png",width = 900, height = 900)
DoHeatmap(GSE178116.rds, features = top5$gene)
dev.off()
DotPlot(GSE178116.rds,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


GSE230389.qc <- GSE230389.rds
#####细胞类型确认
png(filename = 'GSE230389/marker/immune.png',width =350,height = 350 )
FeaturePlot(GSE230389.qc,features =c('PTPRC'),pt.size=0.5,ncol = 1,max.cutoff = 1)
dev.off()
png(filename = 'GSE230389/marker/TAM.png',width =1200,height = 350 )
FeaturePlot(GSE230389.qc,features =c('CD68','CD14','AIF1','C1QB'),pt.size=0.5,ncol = 4,max.cutoff = 1)
dev.off()
png(filename = 'GSE230389/marker/Tumor.png',width =900,height = 350 )
FeaturePlot(GSE230389.qc,features =c('FABP7','SOX2','EGFR'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE230389/marker/BMDM.png',width =1200,height = 350 )
FeaturePlot(GSE230389.qc,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 2)
dev.off()
png(filename = 'GSE230389/marker/MG.png',width =900,height = 350 )
FeaturePlot(GSE230389.qc,features =c('TMEM119','CH25H','P2RY12'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE230389/marker/T cell.png',width =900,height = 350 )
FeaturePlot(GSE230389.qc,features =c('CD3D','CD3E','CD3G'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE230389/marker/B cell.png',width =900,height = 350 )
FeaturePlot(GSE230389.qc,features =c('IGKC','CD79A','JCHAIN'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE230389/marker/NK cell.png',width =900,height = 350 )
FeaturePlot(GSE230389.qc,features =c('NKG7','GNLY','GZMB'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
new.cluster.ids <- c('not immune','not immune','not immune','not immune','not immune','not immune',
                     'TAM','TAM','not immune')
levels(GSE230389.qc@meta.data$celltype) <- new.cluster.ids#将celltype确定
DimPlot(GSE230389.qc,group.by = 'celltype',label = T)
save(GSE230389.qc,file = 'GSE230389/data/GSE230389.qc.Rdata')
# load(file = 'GSE178116/data/GSE178116.qc.Rdata')


#####GSVA验证
library(GSVA)
load('curatedMarkers.Rdata')
Idents(GSE178116.qc)<-'clusters'
expr<-AverageExpression(GSE178116.qc,assays = 'RNA',slot='data')[[1]]    
#expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr[rowSums(expr)>0,])
gsva.scdata<-gsva(expr,curatedMarkers,verbose=T, parallel.sz=10,min.sz=3)



####TAM单独处理
load(file = 'GSE202096/data/GSE202096.qc.Rdata')
GSE202096.tam <- subset(GSE202096.qc, celltype=='TAM')
GSE202096.tam <- subset(GSE202096.tam, singleR=='Macrophage'|singleR=='Monocyte')
GSE202096.tam <- subset(GSE202096.tam, scibet=='Microglia'|scibet=='Macrophage'|scibet=='Monocyte')
GSE202096.tam <- NormalizeData(GSE202096.tam)
GSE202096.tam <- FindVariableFeatures(GSE202096.tam, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE202096.tam), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE202096.tam)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE202096.tam)
GSE202096.tam <- ScaleData(GSE202096.tam, features = all.genes)
GSE202096.tam <- RunPCA(GSE202096.tam, features = VariableFeatures(object = GSE202096.tam))
##
DimPlot(GSE202096.tam, pt.size = 1)
DimHeatmap(GSE202096.tam, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE202096.tam)
GSE202096.tam <- RunUMAP(GSE202096.tam, dims = 1:11)
# GSE202096.tam <- RunTSNE(GSE202096.tam, dims = 1:20,reduction = 'harmony')
GSE202096.tam <- FindNeighbors(GSE202096.tam, dims = 1:11)
GSE202096.tam <- FindClusters(GSE202096.tam, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE202096.tam) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE202096.tam$clusters<-GSE202096.tam$RNA_snn_res.0.1
GSE202096.tam$celltype <- GSE202096.tam$RNA_snn_res.0.1
DimPlot(GSE202096.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
png(filename = 'GSE202096/marker/BMDM-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE202096.tam,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 3)
dev.off()
png(filename = 'GSE202096/marker/MG-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE202096.tam,features =c('TMEM119','CH25H','P2RY12','CRYBB1'),pt.size=0.5,ncol = 4,max.cutoff = 3)
dev.off()
new.cluster.ids <- c('MG.c1','MG.c2','BMDM.c1')
levels(GSE202096.tam$clusters) <- new.cluster.ids
new.celltype.ids <- c('MG-TAM','MG-TAM','BMDM-TAM')
levels(GSE202096.tam$celltype) <- new.celltype.ids
png("GSE202096/marker/TAM-celltype-clusters.png",width = 1200, height = 600)
DimPlot(GSE202096.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
dev.off()
save(GSE202096.tam, file = 'GSE202096/data/GSE202096.tam.Rdata')
