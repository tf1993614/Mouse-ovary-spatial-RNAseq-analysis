#load libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(readbitmap)
library(hdf5r)
library(grid)
library(Matrix)
library(rjson)
library(ggsci)
library(ggrepel)
library(cowplot)

# load integrated dataset

slide.integrated = readRDS("slide.integrated.04.01.23.rds")

# Proceed with downstream analysis on the ”integrated" dataset
slide.integrated = RunPCA(object = slide.integrated, npcs = 50)
slide.integrated = FindNeighbors(object = slide.integrated, reduction = "pca", dims = 1:50,
                                 k.param = 10)
slide.integrated = FindClusters(slide.integrated, resolution = 0.9)
slide.integrated = RunUMAP(slide.integrated, reduction = "pca", dims = 1:50)

# prepare colors for visualization
colors = c(pal_lancet("lanonc")(9), pal_npg("nrc")(10)[c(1:7,9,10)], pal_nejm("default")(8))

# set "seurat_clusters" as default Ident
Idents(slide.integrated) = slide.integrated$seurat_clusters

# UMAP to show clusters
p1 = UMAPPlot(slide.integrated, cols = colors,
              group.by = "seurat_clusters",label = T, label.size = 6, repel = T) +
              theme_UMAP() +
              theme(legend.position = c(0.6,0.9),
              legend.direction = "horizontal")

# UMAP to show the distribution of young and old smaples
p2 = UMAPPlot(slide.integrated, cols = colors, group.by = "group") +
  scale_color_manual(values = colors, labels = c("Aged", "Young")) +
  theme_UMAP() +
  theme(legend.position = c(0.6,0.9),
        plot.title = element_blank(),
        legend.direction = "horizontal")

# integrate p1 and p2
p1 = p1 + theme(axis.text.x.bottom  = element_blank(),
                      axis.text.y.left  = element_blank(),
                      axis.title.x.bottom  = element_blank(),
                      axis.title.y.left = element_blank(),
                legend.position = "none",
                text = element_text(family = "sans"))

p2 = p2 + 
  scale_color_manual(values = c(aged = colors[2], young = colors[1])) +
  theme(axis.text.x.bottom  = element_blank(),
                axis.text.y.left  = element_blank(),
                axis.title.x.bottom  = element_blank(),
                axis.title.y.left = element_blank(),
                text = element_text(family = "sans"))
p3 = p1 + p2
ggsave(p3, filename = "figur1.tiff", height = 2000, width = 4000, units = "px")


# summarize the distribution of each cluster in different samples (A1-4, Y1-4)
cluster_sum = slide.integrated@meta.data %>%
              group_by(mouse, seurat_clusters) %>%
              summarise(n = n()) %>%
              group_by(mouse, .add = T) %>%
              summarise(prop = n/sum(n), cluster = seurat_clusters)

# stack bar plot to show the percentage of each cluster per sample
p4 = cluster_sum %>% ggplot(aes(mouse, prop, fill = cluster)) + geom_bar(stat = "identity") + theme_classic() +
theme(text = element_text(family = "serif"),
      axis.title.y.left = element_text(face = "bold", size = 18),
      axis.line.x.bottom = element_line(size = 0.8),
      axis.line.y.left = element_line(size = 0.8),
      axis.title.x.bottom = element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(colour = "black", size =15),
      axis.text.y = element_text(colour = "black")) +
  scale_fill_manual(values = colors) + labs(y = "Proportion") +
  scale_y_continuous(labels = ~ scales::percent(.x))

# cluster abundance analysis
fisher = fisher_test(slide.integrated@meta.data,
                              group.by.1 = "seurat_clusters",
                              group.by.2 = "group")


# To replace Inf/-Inf with +2/-2 if Inf/-Inf is observed in the table
fisher = fisher %>% mutate(across(where(is.numeric), ~ case_when(.x == Inf ~ 2,
                                                                 .x == -Inf ~ -2,
                                                                 T ~ as.double(.x))))

# Differences of cluster abundance plotted as log10 odds ratio against
# adjusted p value –log10 FDR
p5 = plotOR2FDR(fisher, cols = colors, label = T)

# export figure 4
figure4 = wrap_plots(list(p1,p2,p4,p5), ncol = 2)
cowplot::save_plot(filename = "figure1.tiff", figure4, base_height = 12,
                   base_width = 22 )

# find cluster markers (Note that the default assay is "Integrated" with only 3000 features.
#However when doing DE analysis, we need to go back original "Spatial" assay)
Idents(slide.integrated) = "seurat_clusters"
allMarkersSpatial = FindAllMarkers(slide.integrated, assay = "Spatial",
                                   only.pos = T, latent.vars = "sample",
                                   test.use = "LR")


saveRDS(allMarkersSpatial, "global_marker_genes_310523.rds")

# select top 10 marker genes in each cluster
top10_marker = allMarkersSpatial %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
DoHeatmap(slide.integrated, features = top10_marker$gene)

# find cluster markers (FindConservedMarkers() funtion ) in an alternative way show in
# https://satijalab.org/seurat/articles/integration_introduction.html
Idents(slide.integrated) = slide.integrated$seurat_clusters
cluster_marker_list  = map(0:23, ~ FindConservedMarkers(slide.integrated,
                                                        ident.1 = .x,
                                                        grouping.var = "group",
                                                        assay = "Spatial"))

cluster_marker_list = map(cluster_marker_list, ~ .x %>% rownames_to_column("Gene_name"))

# name the list
cluster_names = str_c("cluster", 0:23)
names(cluster_marker_list) = cluster_names

saveRDS(cluster_marker_list, "global_consevred_marker_genes_regardless_of_age_310523.rds")

# integrate cluster_marker_list into one dataframe
cluster_marker_list = readRDS("./global_consevred_marker_genes_regardless_of_age_310523.rds")
cluster_marker_df = bind_rows(cluster_marker_list, .id = "cluster")

# select top10 marker genes (ranked by minimump_p_val) for each cluster
cluster_marker_top10 = cluster_marker_df %>% group_by(cluster) %>%
 top_n(n=10, wt = -minimump_p_val)


# visualize top10 DE genes in different clusters by heatmap
Idents(slide.integrated) = "seurat_clusters"
# since the data stored in the "integrated" slot only contains 3000 most
# variable genes that maybe doesn't include cluster maker genes we found here.
# so we use the ScaleData function to scale all gene expression store in
# "Spatial" slot and then used for drawing heatmap.
all.genes = rownames(slide.integrated@assays$Spatial@data)
slide.integrated = ScaleData(slide.integrated, features = all.genes,
                             assay = "Spatial")
slide.integrated = NormalizeData(slide.integrated, assay = "Spatial")

p6 = DoHeatmap(slide.integrated, features = cluster_marker_top10$Gene_name, group.colors = colors, angle = 90 ,size = 4, assay = "integrated",
               group.by = "seurat_clusters") +
  guides(col = "none") +
  theme(axis.text.y.left = element_blank())

saveRDS(p6, "p6.rds")
