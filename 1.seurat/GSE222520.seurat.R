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

#GSE222520  IMPx4
assays<-dir('GSE222520/')
dir<-paste0('GSE222520/',assays)
seuratlist<-list()
samples_name<-assays
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
scRNA$orig.ident<-'GSE222520'
Idents(scRNA) <- 'orig.ident'
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
GSE222520.qc<-subset(scRNA,subset=nFeature_RNA>300&nFeature_RNA<7000&percent.mt<20&DF_hi.lo=='Singlet')
save(GSE222520.qc,file = 'data/GSE222520.df.Rdata')
load('GSE222520/data/GSE222520.df.Rdata')

####seurat标准流程
GSE222520.qc <- NormalizeData(GSE222520.qc)
GSE222520.qc <- FindVariableFeatures(GSE222520.qc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE222520.qc), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE222520.qc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE222520.qc)
GSE222520.qc <- ScaleData(GSE222520.qc, features = all.genes)
GSE222520.qc <- RunPCA(GSE222520.qc, features = VariableFeatures(object = GSE222520.qc))
## harmony
GSE222520.qc <- GSE222520.qc %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE222520.qc, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE222520.qc, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE222520.qc, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE222520.qc,reduction = 'harmony')
GSE222520.qc <- RunUMAP(GSE222520.qc, dims = 1:16,reduction = 'harmony')
# GSE222520.qc <- RunTSNE(GSE222520.qc, dims = 1:20,reduction = 'harmony')
GSE222520.qc <- FindNeighbors(GSE222520.qc, dims = 1:16,reduction = 'harmony')
GSE222520.qc <- FindClusters(GSE222520.qc, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE222520.qc) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE222520.qc$clusters<-GSE222520.qc$RNA_snn_res.0.2
GSE222520.qc$celltype<- GSE222520.qc$RNA_snn_res.0.2 
png(filename = 'GSE222520/marker/umap_res.0.2.png',width = '600',height = '600')
DimPlot(GSE222520.qc,group.by = 'RNA_snn_res.0.2',label = T)
dev.off()

save(GSE222520.qc, file = 'GSE222520/data/GSE222520.qc.Rdata')
#load(file = 'GSE222520/data/GSE222520.qc.Rdata')

#copykat
setwd('/media/ubuntu/XYZ/mine/data/LGG/GSE222520/copykat')
seurat_list<-list()
samples<-names(table(GSE222520.qc$sample))
for(i in 1:length(samples)){
  seurat_list[[i]]<-subset(GSE222520.qc,sample==samples[i])
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
GSE222520.qc$copykat.pred<-sc$copykat.pred
setwd('/media/ubuntu/XYZ/mine/data/LGG')
save(copykat.test_list,file=paste0('GSE222520/copykat/copykat.all.result.Rdata'))
save(GSE222520.qc,file = 'GSE222520/data/GSE222520.qc.Rdata')


#####singleR
#准备工作
library(celldex)
library(SingleR)
library(scater)
library(viridis)
library(ggsci)
DimPlot(GSE222520.qc, reduction = "umap", group.by = "orig.ident")  
DimPlot(GSE222520.qc, reduction = "umap",label.size = 5,label = T,pt.size = 0.5)
table(GSE222520.qc$clusters)
a <- table(GSE222520.qc$orig.ident,GSE222520.qc@active.ident)
gplots::balloonplot(a,main="harmony")
saveRDS(GSE222520.qc,'GSE222520/data/GSE222520.qc_adj.rds')
#      SingleR classType
hpca.se <- celldex::HumanPrimaryCellAtlasData()
GSE222520.rds <-readRDS("GSE222520/data/GSE222520.qc_adj.rds")
clusters=GSE222520.rds@meta.data$clusters
pred.hesc <- SingleR(GSE222520.rds@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
GSE222520.rds@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
png(filename = 'GSE222520/marker/SingleR_png.png',width = '600',height = '600')
DimPlot(GSE222520.rds, group.by = 'singleR',pt.size=0.1,label = T,reduction = 'umap')
dev.off()

#####scibet
library(scibet)
exp_tpm <- calculateTPM(GSE222520.rds@assays$RNA@counts)
exp_tpm <- t(exp_tpm)
model <- readr::read_csv("major_human_cell_types.csv")
# model <-model[,-1]
model <- pro.core(model)
prd <- LoadModel(model)
scibet.type <- prd(exp_tpm)
GSE222520.rds$scibet <- scibet.type
png("GSE222520/marker/Scibet_png.png",width = 600, height = 600)
DimPlot(GSE222520.rds, group.by = "scibet", pt.size=0.1, label = T, reduction = 'umap')
dev.off()

celltype1 <- data.frame(
  Barcodes=colnames(GSE222520.rds@assays$RNA),
  ClusterID=GSE222520.rds$clusters,
  celltype=scibet.type, 
  stringsAsFactors = F
)
celltype_scibet<-celltype1[!duplicated(celltype1[,2:3]),2:3]
# saveRDS(GSE222520/data/GSE222520.rds,"GSE222520.qc_adj.rds")


##   FindAllMarkers 
DEGS<-FindAllMarkers(GSE222520.rds,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE222520/marker/Heatmap.png",width = 900, height = 900)
DoHeatmap(GSE222520.rds, features = top5$gene)
dev.off()
DotPlot(GSE222520.rds,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


GSE222520.qc <- GSE222520.rds
#####细胞类型确认
png(filename = 'GSE222520/marker/immune.png',width =350,height = 350 )
FeaturePlot(GSE222520.qc,features =c('PTPRC'),pt.size=0.5,ncol = 1,max.cutoff = 1)
dev.off()
png(filename = 'GSE222520/marker/TAM.png',width =1200,height = 350 )
FeaturePlot(GSE222520.qc,features =c('CD68','CD14','AIF1','C1QB'),pt.size=0.5,ncol = 4,max.cutoff = 1)
dev.off()
png(filename = 'GSE222520/marker/Tumor.png',width =900,height = 350 )
FeaturePlot(GSE222520.qc,features =c('FABP7','SOX2','EGFR'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE222520/marker/BMDM.png',width =1200,height = 350 )
FeaturePlot(GSE222520.qc,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 2)
dev.off()
png(filename = 'GSE222520/marker/MG.png',width =900,height = 350 )
FeaturePlot(GSE222520.qc,features =c('TMEM119','CH25H','P2RY12'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE222520/marker/T cell.png',width =900,height = 350 )
FeaturePlot(GSE222520.qc,features =c('CD3D','CD3E','CD3G'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE222520/marker/B cell.png',width =900,height = 350 )
FeaturePlot(GSE222520.qc,features =c('IGKC','CD79A','JCHAIN'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
png(filename = 'GSE222520/marker/NK cell.png',width =900,height = 350 )
FeaturePlot(GSE222520.qc,features =c('NKG7','GNLY','GZMB'),pt.size=0.5,ncol = 3,max.cutoff = 2)
dev.off()
new.celltype.ids <- c('TAM','TAM','others','TAM','TAM','others','TAM')
levels(GSE222520.qc@meta.data$celltype) <- new.celltype.ids#将celltype确定
DimPlot(GSE222520.qc,group.by = 'celltype',label = T)
save(GSE222520.qc,file = 'GSE222520/data/GSE222520.qc.Rdata')
#load(file = 'GSE222520/data/GSE222520.qc.Rdata')


#####GSVA验证
library(GSVA)
load('curatedMarkers.Rdata')
Idents(GSE222520.qc)<-'clusters'
expr<-AverageExpression(GSE222520.qc,assays = 'RNA',slot='data')[[1]]    
#expr<-expr[rowSums(expr)>0,]
expr<-as.matrix(expr[rowSums(expr)>0,])
gsva.scdata<-gsva(expr,curatedMarkers,verbose=T, parallel.sz=10,min.sz=3)



####TAM单独处理
GSE222520.tam <- subset(GSE222520.qc, celltype == 'TAM')
GSE222520.tam <- subset(GSE222520.tam, singleR=='Macrophage'|singleR=='Monocyte')
GSE222520.tam <- subset(GSE222520.tam, scibet=='Macrophage'|scibet=='Monocyte'|scibet=='Microglia')
GSE222520.tam <- NormalizeData(GSE222520.tam)
GSE222520.tam <- FindVariableFeatures(GSE222520.tam, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE222520.tam), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(GSE222520.tam)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(GSE222520.tam)
GSE222520.tam <- ScaleData(GSE222520.tam, features = all.genes)
GSE222520.tam <- RunPCA(GSE222520.tam, features = VariableFeatures(object = GSE222520.tam))
## harmony
GSE222520.tam <- GSE222520.tam %>% RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(GSE222520.tam, 'harmony')
harmony_embeddings[1:5, 1:5]
##
DimPlot(GSE222520.tam, reduction = "harmony",pt.size = 1)
DimHeatmap(GSE222520.tam, reduction="harmony",dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(GSE222520.tam,reduction = 'harmony')
GSE222520.tam <- RunUMAP(GSE222520.tam, dims = 1:18,reduction = 'harmony')
# GSE222520.tam <- RunTSNE(GSE222520.tam, dims = 1:20,reduction = 'harmony')
GSE222520.tam <- FindNeighbors(GSE222520.tam, dims = 1:18,reduction = 'harmony')
GSE222520.tam <- FindClusters(GSE222520.tam, resolution = c(0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9))
clus.tree.out <- clustree(GSE222520.tam) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
GSE222520.tam$clusters<-GSE222520.tam$RNA_snn_res.0.1
GSE222520.tam$celltype<-GSE222520.tam$RNA_snn_res.0.1
DimPlot(GSE222520.tam,group.by = 'clusters',pt.size = 0.01, reduction = 'umap', label = T)

##   FindAllMarkers
Idents(GSE222520.tam) <- 'clusters'
DEGS<-FindAllMarkers(GSE222520.tam,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
DEGS<-data.frame(gene=rownames(DEGS),DEGS)
# saveRDS(DEGS,"FindAllMarkers_result.rds")
# DEGS <- readRDS("./FindAllMarkers_result.rds")
top10 = DEGS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = DEGS %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top5=top5[!duplicated(top5$gene),]
select_gene <- split(top5$gene,top5$cluster)
png("GSE222520/marker/TAMHeatmap.png",width = 900, height = 900)
DoHeatmap(GSE222520.tam, features = top5$gene)
dev.off()
png("GSE222520/marker/TAMDotmap.png",width = 900, height = 900)
DotPlot(GSE222520.tam,features = top5$gene, assay = "RNA") + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

#singleR
Idents(GSE222520.tam)<-'clusters'
table(GSE222520.tam$clusters)
a <- table(GSE222520.tam$orig.ident,GSE222520.tam@active.ident)
gplots::balloonplot(a,main="harmony")
clusters=GSE222520.tam@meta.data$clusters
pred.hesc <- SingleR(GSE222520.tam@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main,clusters = clusters)
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F)
GSE222520.tam@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
png(filename = 'GSE222520/marker/SingleR_TAM.png',width = '600',height = '600')
DimPlot(GSE222520.tam, group.by = 'singleR',pt.size=0.1,label = T,reduction = 'umap')
dev.off()
#####scibet
exp_tpm <- calculateTPM(GSE222520.tam@assays$RNA@counts)
exp_tpm <- t(exp_tpm)
model <- readr::read_csv("major_human_cell_types.csv")
# model <-model[,-1]
model <- pro.core(model)
prd <- LoadModel(model)
scibet.type <- prd(exp_tpm)
GSE222520.tam$scibet <- scibet.type
png("GSE222520/marker/Scibet_TAM.png",width = 600, height = 600)
DimPlot(GSE222520.tam, group.by = "scibet", pt.size=0.1, label = T, reduction = 'umap')
dev.off()


png(filename = 'GSE222520/marker/BMDM-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE222520.tam,features =c('TGFBI','LYZ','S100A8','CD163'),pt.size=0.5,ncol = 4,max.cutoff = 'q99')
dev.off()
png(filename = 'GSE222520/marker/MG-TAM.png',width =1200,height = 350 )
FeaturePlot(GSE222520.tam,features =c('TMEM119','CH25H','P2RY12','CRYBB1'),pt.size=0.5,ncol = 4,max.cutoff = 'q99')
dev.off()
DimPlot(GSE222520.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
new.cluster.ids <- c('MG.c1','MG.c2','BMDM.c1','MG.c3')
levels(GSE222520.tam$clusters) <- new.cluster.ids
new.celltype.ids <- c('MG-TAM','MG-TAM','BMDM-TAM','MG-TAM')
levels(GSE222520.tam$celltype) <- new.celltype.ids
png("GSE222520/marker/TAM-celltype-clusters.png",width = 1200, height = 600)
DimPlot(GSE222520.tam,group.by = c('celltype','clusters'),pt.size=0.01,reduction = 'umap',label = T)
dev.off()
save(GSE222520.tam, file = 'GSE222520/data/GSE222520.tam.Rdata')
