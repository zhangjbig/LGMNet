rm(list = ls())
setwd('E:/mine/newlgg')

library(Seurat)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(harmony)
library(ggplot2)
library(clustree)
library(copykat)

#GSE89567  IMPx4
seurat_data<- read.table(gzfile("GSE89567/GSE89567_RAW.txt.gz"), row.names = 1, header = TRUE, sep = "\t")
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 200,
                                   min.cells = 3,
                                   project = "GSE89567")

seurat_obj<-RenameCells(seurat_obj,add.cell.id='GSE89567')
seurat_obj[['percent.mt']]<-PercentageFeatureSet(seurat_obj,pattern="^MT-")
cellname <- colnames(seurat_obj)
split_elements <- function(element) { # 使用正则表达式分隔符（_ 或 .） 
  strsplit(element, "[_.]")[[1]] } # 对第六列应用分隔函数 
split_list <- lapply(cellname, split_elements) # 提取四个新向量 
seurat_obj$sample <- sapply(split_list, function(x) x[2]) 
seurat_obj$sample<-gsub('mgh103','MGH103', seurat_obj$sample)
seurat_obj$sample<-gsub('MGH107neg','MGH107', seurat_obj$sample)
seurat_obj$sample<-gsub('MGH107pos','MGH107', seurat_obj$sample)
seurat_obj$sample<-gsub('X57','MGH57', seurat_obj$sample)


#####去除双细胞
library(DoubletFinder)
data<-seurat_obj
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
data@meta.data["DF_hi.lo"] <- data@meta.data[9]
data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet" & data@meta.data[11] == "Singlet")] <- "Doublet-Low Confidience"
data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
meta.data<-data@meta.data %>% select(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,seurat_clusters,DF_hi.lo,sample)
data@meta.data<-meta.data
seurat_obj<-data


scRNA<-seurat_obj
scRNA$orig.ident<-'GSE89567'
Idents(scRNA) <- 'orig.ident'
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GSE89567.qc<-subset(scRNA,subset=nFeature_RNA>300&nFeature_RNA<7000&percent.mt<10&DF_hi.lo=='Singlet')
save(GSE89567.qc,file = 'GSE89567/data/GSE89567.df.Rdata')
load('GSE89567/data/GSE89567.df.Rdata')


####seurat标准流程
GSE89567.qc <- NormalizeData(GSE89567.qc)
GSE89567.qc <- FindVariableFeatures(GSE89567.qc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE89567.qc), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE89567.qc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE89567.qc)
GSE89567.qc <- ScaleData(GSE89567.qc, features = all.genes)
GSE89567.qc <- RunPCA(GSE89567.qc, features = VariableFeatures(object = GSE89567.qc))
## harmony
GSE89567.qc <- GSE89567.qc %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE89567.qc, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE89567.qc, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE89567.qc, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE89567.qc,reduction = 'harmony')
GSE89567.qc <- RunUMAP(GSE89567.qc, dims = 1:10,reduction = 'harmony')
# GSE89567.qc <- RunTSNE(GSE89567.qc, dims = 1:20,reduction = 'harmony')
GSE89567.qc <- FindNeighbors(GSE89567.qc, dims = 1:10,reduction = 'harmony')
GSE89567.qc <- FindClusters(GSE89567.qc, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE89567.qc) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE89567.qc$clusters<-GSE89567.qc$RNA_snn_res.0.2
GSE89567.qc$celltype<- GSE89567.qc$RNA_snn_res.0.2
png(filename = 'GSE89567/marker/umap_res.0.2_1.png',width = '600',height = '600')
DimPlot(GSE89567.qc,group.by = 'RNA_snn_res.0.2',label = F)
dev.off()
save(GSE89567.qc, file = 'GSE89567/data/GSE89567.qc.Rdata')
#load(file = 'GSE89567/data/GSE89567.qc.Rdata')

#copykat
setwd('/media/ubuntu/XYZ/mine/data/LGG/GSE89567/copykat')
seurat_list<-list()
samples<-names(table(GSE89567.qc$sample))
for(i in 1:length(samples)){
  seurat_list[[i]]<-subset(GSE89567.qc,sample==samples[i])
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
GSE89567.qc$copykat.pred<-sc$copykat.pred
setwd('/media/ubuntu/XYZ/mine/data/LGG')
save(copykat.test_list,file=paste0('GSE89567/copykat/copykat.all.result.Rdata'))
save(GSE89567.qc,file = 'GSE89567/data/GSE89567.qc.Rdata')


#####singleR
#准备工作
library(celldex)
library(SingleR)
library(scater)
library(viridis)
library(ggsci)
DimPlot(GSE89567.qc, reduction = "umap", group.by = "orig.ident")  
DimPlot(GSE89567.qc, reduction = "umap",label.size = 5,label = T,pt.size = 0.5)
table(GSE89567.qc$clusters)
a <- table(GSE89567.qc$orig.ident,GSE89567.qc@active.ident)
gplots::balloonplot(a,main="harmony")
saveRDS(GSE89567.qc,'GSE89567/data/GSE89567.qc_adj.rds')
#      SingleR classType
hpca.se <- celldex::HumanPrimaryCellAtlasData()
GSE89567.rds <-readRDS("GSE89567/data/GSE89567.qc_adj.rds")
GSE89567.rds <- GSE89567.qc
clusters=GSE89567.rds@meta.data$clusters
pred.hesc <- SingleR(GSE89567.rds@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
GSE89567.rds@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
png(filename = 'GSE89567/marker/SingleR_png.png',width = '600',height = '600')
DimPlot(GSE89567.rds, group.by = 'singleR',pt.size=0.1,label = T,reduction = 'umap')
dev.off()

#####scibet
library(scibet)
exp_tpm <- calculateTPM(GSE89567.rds@assays$RNA@counts)
exp_tpm <- t(exp_tpm)
model <- readr::read_csv("major_human_cell_types.csv")
# model <-model[,-1]
model <- pro.core(model)
prd <- LoadModel(model)
scibet.type <- prd(exp_tpm)
GSE89567.rds$scibet <- scibet.type
png("GSE89567/marker/Scibet_png.png",width = 600, height = 600)
DimPlot(GSE89567.rds, group.by = "scibet", pt.size=0.1, label = T, reduction = 'umap')
dev.off()

celltype1 <- data.frame(
  Barcodes=colnames(GSE89567.rds@assays$RNA),
  ClusterID=GSE89567.rds$clusters,
  celltype=scibet.type, 
  stringsAsFactors = F
)
celltype_scibet<-celltype1[!duplicated(celltype1[,2:3]),2:3]
# saveRDS(GSE89567/data/GSE89567.rds,"GSE89567.qc_adj.rds")


##   FindAllMarkers 
DEGS<-FindAllMarkers(GSE89567.rds,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE89567/marker/Heatmap.png",width = 900, height = 900)
DoHeatmap(GSE89567.rds, features = top5$gene)
dev.off()
DotPlot(GSE89567.rds,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


GSE89567.qc <- GSE89567.rds
#####细胞类型确认
png(filename = 'GSE89567/marker/immune.png',width =350,height = 350 )
FeaturePlot(GSE89567.qc,features =c('PTPRC'),pt.size=0.5,ncol = 1,max.cutoff = 1)
dev.off()
png(filename = 'GSE89567/marker/TAM.png',width =1200,height = 350 )
FeaturePlot(GSE89567.qc,features =c('CD68','CD14','AIF1','C1QB'),pt.size=0.5,ncol = 4,max.cutoff = 1)
dev.off()
png(filename = 'GSE89567/marker/Tumor.png',width =900,height = 350 )
FeaturePlot(GSE89567.qc,features =c('FABP7','SOX2','EGFR'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE89567/marker/BMDM.png',width =1200,height = 350 )
FeaturePlot(GSE89567.qc,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 2)
dev.off()
png(filename = 'GSE89567/marker/MG.png',width =900,height = 350 )
FeaturePlot(GSE89567.qc,features =c('TMEM119','CH25H','P2RY12'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE89567/marker/T cell.png',width =900,height = 350 )
FeaturePlot(GSE89567.qc,features =c('CD3D','CD3E','CD3G'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE89567/marker/B cell.png',width =900,height = 350 )
FeaturePlot(GSE89567.qc,features =c('IGKC','CD79A','JCHAIN'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE89567/marker/NK cell.png',width =900,height = 350 )
FeaturePlot(GSE89567.qc,features =c('NKG7','GNLY','GZMB'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
new.cluster.ids <- c('not immune','TAM','not immune','not immune','not immune','not immune')
levels(GSE89567.qc@meta.data$celltype) <- new.cluster.ids#将celltype确定
DimPlot(GSE89567.qc,group.by = 'celltype',label = T)
save(GSE89567.qc,file = 'GSE89567/data/GSE89567.qc.Rdata')
# load(file = 'GSE89567/data/GSE89567.qc.Rdata')


#####GSVA验证
library(GSVA)
load('curatedMarkers.Rdata')
Idents(GSE89567.qc)<-'clusters'
expr<-AverageExpression(GSE89567.qc,assays = 'RNA',slot='data')[[1]]    
#expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr[rowSums(expr)>0,])
gsva.scdata<-gsva(expr,curatedMarkers,verbose=T, parallel.sz=10,min.sz=3)



####TAM单独处理
load(file = 'GSE89567/data/GSE89567.qc.Rdata')
GSE89567.tam <- subset(GSE89567.qc, celltype=='TAM' & singleR=='Monocyte')
GSE89567.tam <- subset(GSE89567.tam, scibet=='Macrophage'|scibet=='Microglia'|scibet=='Monocyte')
GSE89567.tam <- NormalizeData(GSE89567.tam)
GSE89567.tam <- FindVariableFeatures(GSE89567.tam, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE89567.tam), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE89567.tam)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE89567.tam)
GSE89567.tam <- ScaleData(GSE89567.tam, features = all.genes)
GSE89567.tam <- RunPCA(GSE89567.tam, features = VariableFeatures(object = GSE89567.tam))
##
## harmony
GSE89567.tam <- GSE89567.tam %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE89567.tam, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE89567.tam, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE89567.tam, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE89567.tam,reduction = 'harmony')
GSE89567.tam <- RunUMAP(GSE89567.tam, dims = 1:15,reduction = 'harmony')
# GSE89567.tam <- RunTSNE(GSE89567.tam, dims = 1:20,reduction = 'harmony')
GSE89567.tam <- FindNeighbors(GSE89567.tam, dims = 1:15,reduction = 'harmony')
GSE89567.tam <- FindClusters(GSE89567.tam, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
png('GSE89567/marker/TAM.clustree.png', width = 1000, height = 800)
clus.tree.out <- clustree(GSE89567.tam) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()
clus.tree.out

GSE89567.tam$clusters<-GSE89567.tam$RNA_snn_res.0.4
GSE89567.tam$celltype<-GSE89567.tam$RNA_snn_res.0.4
DimPlot(GSE89567.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)


##   FindAllMarkers
Idents(GSE89567.tam) <- 'clusters'
DEGS<-FindAllMarkers(GSE89567.tam,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE89567/marker/TAMHeatmap.png",width = 900, height = 900)
DoHeatmap(GSE89567.tam, features = top5$gene)
dev.off()
png("GSE89567/marker/TAMDotmap.png", width = 900, height = 900)
DotPlot(GSE89567.tam,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

png(filename = 'GSE89567/marker/BMDM-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE89567.tam,features =c('TGFBI','S100A8','LYZ','CD163'),pt.size=1,ncol = 4)
dev.off()
png(filename = 'GSE89567/marker/MG-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE89567.tam,features =c('TMEM119','CH25H','P2RY12','CRYBB1'),pt.size=1,ncol = 4)
dev.off()

new.cluster.ids <- c('MG.c1','MG.c2','BMDM.c1')
levels(GSE89567.tam$clusters) <- new.cluster.ids
GSE89567.tam$celltype<-GSE89567.tam$RNA_snn_res.0.4
new.cluster.ids <- c('MG-TAM','MG-TAM','BMDM-TAM')
levels(GSE89567.tam$celltype) <- new.cluster.ids
png("GSE89567/marker/TAM-celltype-clusters.png",width = 1200, height = 600)
DimPlot(GSE89567.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
dev.off()
save(GSE89567.tam, file = 'GSE89567/data/GSE89567.tam.Rdata')


