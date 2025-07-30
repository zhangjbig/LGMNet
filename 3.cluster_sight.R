rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggpubr)

load(file = 'GSE89567/data/GSE89567.tam.Rdata')
load(file = 'GSE152277/data/GSE152277.tam.Rdata')
load(file = 'GSE202096/data/GSE202096.tam.Rdata')
load(file = 'GSE222520/data/GSE222520.tam.Rdata')
load(file = 'GSE270109/data/GSE270109.tam.Rdata')
load(file = 'GSE227718/data/GSE227718.tam.Rdata')

#计算不同簇的比例
TAM_prop_GSE89567 <- GSE89567.tam@meta.data %>%
  group_by(clusters, celltype) %>%
  summarize(count = n())
TAM_prop_GSE89567 <- TAM_prop_GSE89567 %>%
  mutate(proportion = count/ sum(TAM_prop_GSE89567$count)) %>%
  mutate(data_set = "GSE89567")

TAM_prop_GSE152277 <- GSE152277.tam@meta.data %>%
  group_by(clusters, celltype) %>%
  summarize(count = n())
TAM_prop_GSE152277 <- TAM_prop_GSE152277 %>%
  mutate(proportion = count/ sum(TAM_prop_GSE152277$count)) %>%
  mutate(data_set = "GSE152277")

TAM_prop_GSE202096 <- GSE202096.tam@meta.data %>%
  group_by(clusters, celltype) %>%
  summarize(count = n())
TAM_prop_GSE202096 <- TAM_prop_GSE202096 %>%
  mutate(proportion = count/sum(TAM_prop_GSE202096$count)) %>%
  mutate(data_set = "GSE202096")

TAM_prop_GSE222520 <- GSE222520.tam@meta.data %>%
  group_by(clusters, celltype) %>%
  summarize(count = n())
TAM_prop_GSE222520 <- TAM_prop_GSE222520 %>%
  mutate(proportion = count/ sum(TAM_prop_GSE222520$count)) %>%
  mutate(data_set = "GSE222520")

TAM_prop_GSE227718 <- GSE227718.tam@meta.data %>%
  group_by(clusters, celltype) %>%
  summarize(count = n())
TAM_prop_GSE227718 <- TAM_prop_GSE227718 %>%
  mutate(proportion = count/ sum(TAM_prop_GSE227718$count)) %>%
  mutate(data_set = "GSE227718")

TAM_prop_GSE270109 <- GSE270109.tam@meta.data %>%
  group_by(clusters, celltype) %>%
  summarize(count = n())
TAM_prop_GSE270109 <- TAM_prop_GSE270109 %>%
  mutate(proportion = count/ sum(TAM_prop_GSE270109$count)) %>%
  mutate(data_set = "GSE270109")

TAM_prop <- bind_rows(TAM_prop_GSE89567, TAM_prop_GSE152277, TAM_prop_GSE202096, 
                      TAM_prop_GSE222520, TAM_prop_GSE227718, TAM_prop_GSE270109)


### 画每个cluster在sample里分布的饼图
pie_charts_list <- list()
for (i in 1:nrow(TAM_prop)) {
  row_data <- TAM_prop[i, ]
  data_set <- row_data$data_set
  cluster <- row_data$clusters
  celltype <- row_data$celltype
  proportion <- row_data$proportion
  other_proportion <- 1-proportion
  # 包含当前行的临时数据框
  temp_data <- data.frame(celltype = c(celltype, "others"), 
                          proportion = c(proportion, other_proportion))
  temp_data$percentage <- paste0(round(temp_data$proportion * 100, 1), "%")
  mycols <- c("#0073C2FF", "#EFC000FF")
  # 绘制饼图
  p <- ggplot(temp_data, aes(x = "", y = proportion, fill = celltype)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y") + #将条形图转换为饼图
    theme_void() +
    labs(title = paste(data_set, "-", cluster), fill = "celltype") +
    scale_fill_manual(values = mycols) +
    geom_text(aes(label=percentage), position = position_stack(vjust = 0.5))
  ggsave(paste('compare/pie/',data_set,'_',cluster,'.png',sep = ''), width = 3, height = 2, bg = 'white')
  # 将生成的饼图添加到列表中
  pie_charts_list[[paste(data_set, cluster, sep = "_")]] <- p
}
saveRDS(pie_charts_list, file = "compare/pie/pie_list.rds")


##提取样本/簇/细胞类型
sample_GSE89567 <- GSE89567.tam@meta.data %>%
  group_by(clusters, sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(data_set = "GSE89567")
sample_GSE152277 <- GSE152277.tam@meta.data %>%
  group_by(clusters, sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(data_set = "GSE152277")
sample_GSE202096 <- GSE202096.tam@meta.data %>%
  group_by(clusters, sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(data_set = "GSE202096")
sample_GSE222520 <- GSE222520.tam@meta.data %>%
  group_by(clusters, sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(data_set = "GSE222520")
sample_GSE227718 <- GSE227718.tam@meta.data %>%
  group_by(clusters, sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(data_set = "GSE227718")
sample_GSE270109 <- GSE270109.tam@meta.data %>%
  group_by(clusters, sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(data_set = "GSE270109")
TAM_sample <- bind_rows(sample_GSE89567, sample_GSE152277, sample_GSE202096, 
                        sample_GSE222520, sample_GSE227718, sample_GSE270109)
saveRDS(TAM_sample, "compare/TAM_sample.rds")



sample_GSE152277_BMDM.c1 <- sample_GSE152277 %>% filter(clusters == 'BMDM.c1')
sample_GSE202096_BMDM.c1 <- sample_GSE202096 %>% filter(clusters == 'BMDM.c1')
sample_GSE222520_BMDM.c1 <- sample_GSE222520 %>% filter(clusters == 'BMDM.c1')
sample_GSE227718_BMDM.c1 <- sample_GSE227718 %>% filter(clusters == 'BMDM.c1')

TAM_sample_cluster1 <- bind_rows(sample_GSE152277_BMDM.c1,sample_GSE202096_BMDM.c1,
                                 sample_GSE222520_BMDM.c1,sample_GSE227718_BMDM.c1)

sample_GSE89567_MG.c1 <- sample_GSE89567 %>% filter(clusters == 'MG.c1')
sample_GSE152277_MG.c1 <- sample_GSE152277 %>% filter(clusters == 'MG.c1')
sample_GSE270109_MG.c1 <- sample_GSE270109 %>% filter(clusters == 'MG.c1')

TAM_sample_cluster2 <- bind_rows(sample_GSE89567_MG.c1,sample_GSE152277_MG.c1,sample_GSE270109_MG.c1)

sample_GSE89567_MG.c2 <- sample_GSE89567 %>% filter(clusters == 'MG.c2')
sample_GSE222520_MG.c3 <- sample_GSE222520 %>% filter(clusters == 'MG.c3')
sample_GSE270109_MG.c3 <- sample_GSE270109 %>% filter(clusters == 'MG.c3')

TAM_sample_cluster3 <- bind_rows(sample_GSE89567_MG.c2,sample_GSE222520_MG.c3,sample_GSE270109_MG.c3)

sample_GSE152277_MG.c3 <- sample_GSE152277 %>% filter(clusters == 'MG.c3')
sample_GSE202096_MG.c1 <- sample_GSE202096 %>% filter(clusters == 'MG.c1')

TAM_sample_cluster4 <- bind_rows(sample_GSE152277_MG.c3,sample_GSE202096_MG.c1)

sample_GSE152277_MG.c2 <- sample_GSE152277 %>% filter(clusters == 'MG.c2')
sample_GSE222520_MG.c2 <- sample_GSE222520 %>% filter(clusters == 'MG.c2')

TAM_sample_cluster5 <- bind_rows(sample_GSE152277_MG.c2, 
                                  sample_GSE222520_MG.c2)



######
p1 <- ggboxplot(TAM_sample_cluster1, x="data_set", y="count", color="data_set",
                palette = "jco", add = "jitter", xlab = "Data set", ylab = "Count", legend = "right")
my_comparisons = list(c("GSE152277", "GSE202096"), c("GSE152277", "GSE222520"),c("GSE152277", "GSE227718"), 
                      c("GSE202096", "GSE222520"), c("GSE202096", "GSE227718"),c("GSE222520", "GSE227718"))
png('compare/boxplot/cluster1.png', width =600,height = 400)
p1+stat_compare_means(comparisons = my_comparisons)
dev.off()

p2 = ggboxplot(TAM_sample_cluster2, x = 'data_set', y = 'count', color = 'data_set',
               palette = 'jco', add = 'jitter', xlab = 'Data set', ylab = 'Count', legend = 'right')
my_comparisons = list(c('GSE89567','GSE152277'),c('GSE89567','GSE270109'),c('GSE152277','GSE270109'))
png('compare/boxplot/cluster2.png', width =600,height = 400)
p2+stat_compare_means(comparisons = my_comparisons)
dev.off()

p3 = ggboxplot(TAM_sample_cluster3, x = 'data_set', y = 'count', color = 'data_set',
               palette = 'jco', add = 'jitter', xlab = 'Data set', ylab = 'Count', legend = 'right')
my_comparisons = list(c('GSE89567','GSE222520'),c('GSE89567','GSE270109'),c('GSE222520','GSE270109'))
png('compare/boxplot/cluster3.png', width =600,height = 400)
p3+stat_compare_means(comparisons = my_comparisons)
dev.off()

p4 = ggboxplot(TAM_sample_cluster4, x = 'data_set', y = 'count', color = 'data_set',
               palette = 'jco', add = 'jitter', xlab = 'Data set', ylab = 'Count', legend = 'right')
my_comparisons = list(c('GSE152277','GSE202096'))
png('compare/boxplot/cluster4.png', width =600,height = 400)
p4 + stat_compare_means(comparisons = my_comparisons)
dev.off()

p5 = ggboxplot(TAM_sample_cluster5, x = 'data_set', y = 'count', color = 'data_set',
               palette = 'jco', add = 'jitter', xlab = 'Data set', ylab = 'Count', legend = 'right')
my_comparisons = list(c('GSE152277','GSE222520'))
png('compare/boxplot/cluster5.png', width = 600, height = 400)
p5 + stat_compare_means(comparisons = my_comparisons)
dev.off()


##绘制每个cluster箱线图
datasets <- list('GSE89567' = sample_GSE89567, 'GSE152277' = sample_GSE152277, 
                 'GSE202096' = sample_GSE202096, 'GSE222520' = sample_GSE222520, 
                 'GSE227718' = sample_GSE227718, 'GSE270109' = sample_GSE270109)
boxplots_list <- list()
data_sets <- unique(TAM_sample$data_set)

# 遍历每个 data_set
for (data_set in data_sets) {
  # 过滤出当前 data_set 的数据
  df_data_set <- TAM_sample[TAM_sample$data_set == data_set, ]
  
  # 获取当前 data_set 中所有唯一的 seurat_clusters
  clusters <- unique(df_data_set$clusters)
  
  # 对每个 cluster 生成箱线图
  for (cluster in clusters) {
    # 过滤数据
    df_subset <- df_data_set[df_data_set$clusters == cluster, ]
    
    # 检查数据子集是否为空
    if (nrow(df_subset) > 0) {
      # 生成箱线图
      p <- ggplot(df_subset, aes(x = factor(clusters), y = count, color = sample)) +
        geom_boxplot(fill = "lightblue", color = "black", outlier.shape = NA) + # 绘制箱线图，隐藏默认的离群点
        geom_jitter(position = position_jitter(0.4), alpha = 0.8, size = 1) + # 添加样本点并根据样本名着色
        labs(x = "cluster", # 设置x轴标签为 "cluster"
             y = "Count", # 设置y轴标签
             color = "sample") + # 设置图例标题为 "sample"
        scale_y_continuous(breaks = seq(0, max(df_subset$count), by=100)) +
        scale_y_continuous(limits = c(0,NA))+
        theme_minimal() + # 设置主题 
        theme(panel.background = element_blank(),
              plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color="black", fill=NA, linewidth=1),
              axis.title.y = element_text(angle = 90, vjust = 0.5),
              legend.title = element_text(size = 10), # 修改图例标题的字体大小 
              legend.text = element_text(size = 8), # 修改图例标签的字体大小 
              legend.key.size = unit(0.5, "cm"), # 调整图例项的大小 
              legend.box.spacing = unit(0.5, "cm"))
      ggsave(paste('compare/boxplot/', data_set, '_', cluster, '.png', sep = ''), 
             width = 6, height = 4, bg = 'white')
      # 将图形添加到列表中
      boxplots_list[[paste(data_set, cluster, sep = "_")]] <- p
    }
  }
}

saveRDS(boxplots_list, file = 'compare/boxplot/boxplots_list.rds')
