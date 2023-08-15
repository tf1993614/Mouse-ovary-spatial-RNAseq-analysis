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
library(pheatmap)
library(clusterProfiler)
library(ggprism)

# read Seurat object in
slide.integrated = readRDS("./slide.integrated.04.01.23.rds")

# read heatmap showing top20 marker genes in each cluster_rows
p6 = readRDS("./p6.rds")

# assign new cluster IDs
new.cluster.ids = c("CL-R", "Stroma", "Stroma", "CL-R", "Stroma", "Stroma",
                    "CL-R", "Follicle", "CL-R", "Follicle", "CL-P",
                    "Epithelium 1", "CL-P", "Stroma", "Stroma",
                    "Follicle",
                    "Follicle", "CL-P", "CL-R", "Follicle", "Follicle",
                    "CL-P", "Epithelium 2", "Epithelium 3")

# rename the Idents of Seurat object
names(new.cluster.ids) = levels(slide.integrated)
slide.integrated = RenameIdents(slide.integrated, new.cluster.ids)

# add new cluster IDs into metedata table of Seurat object
slide.integrated[["cell_type"]] = Idents(slide.integrated)


# find marker genes in per newly combined cluster
Idents(slide.integrated) = slide.integrated$cell_type

cluster_marker_after_combing_cluster = map(unique(slide.integrated$cell_type),
                                           ~ FindConservedMarkers(slide.integrated, ident.1 = .x,
                                          grouping.var = "group", assay = "integrated"))

# transform row index into a new column
cluster_marker_after_combing_cluster = map(cluster_marker_after_combing_cluster,
                                           ~ .x %>% rownames_to_column("Gene_name"))

# name the list
names(cluster_marker_after_combing_cluster) = unique(slide.integrated$cell_type)

saveRDS(cluster_marker_after_combing_cluster, "consevred_marker_genes_in_newly_combined_cluster.rds")

# select top 20 marker genes (ranked by minimump_p_val) for each cluster
cluster_marker_after_combing_cluster = readRDS("./consevred_marker_genes_in_newly_combined_cluster_310523.rds")

top20_cluster_gene = cluster_marker_after_combing_cluster %>%
  bind_rows(.id = "cluster") %>% group_by(cluster) %>%
  filter(Gene_name %in% rownames(slide.integrated)) %>%
  slice_min(n = 20, order_by = minimump_p_val)

# show the top 20 genes in the newly combined cluster by heatmap
p7 = DoHeatmap(slide.integrated, features = top20_cluster_gene$Gene_name,
               group.colors = colors, angle = 90 ,size = 4) +
  guides(col = "none") +
  theme(text = element_text(family = "sans", face = "bold"),
        axis.text.y.left = element_text(colour = "black"))

# export figure 5
figure5 = wrap_plots(list(p6, p7))
cowplot::save_plot(filename = "figure2.tiff", figure5, base_height = 15,
                   base_width = 18 )

# UMAP plot to show the distriution of newly-defined cell type
Idents(slide.integrated) = slide.integrated$cell_type
p8 = DimPlot(slide.integrated, cols = colors) + theme_UMAP() +
  theme(legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

#summarize the aged and young spots in each cell type
cellType_sum = slide.integrated@meta.data %>%
               group_by(mouse, cell_type) %>%
               summarise(count = n()) %>%
               group_by(mouse) %>%
               mutate(percent = count/sum(count))

age_sum = slide.integrated@meta.data %>%
          group_by(cell_type, group) %>%
          summarise(count = n()) %>%
          mutate(cell_type = factor(cell_type,
                                    levels = c("CL-P",
                                               "CL-R",
                                               "Stroma",
                                               "Follicle",
                                               str_c("Epithelium ", 1:3)
                                              )
                                   )
                )

# bar plot to show the proportion of cell type
p9 = age_sum %>% ggplot(aes(cell_type,count, fill = group)) +
                 geom_bar(stat = "identity", position = "fill") +
                 theme_bar() +
                 scale_fill_manual(values = colors,
                                   labels = c("Aged", "Young")) +
                 labs(y = "Proportion") +
                 scale_y_continuous(labels = ~ scales::percent(.x)) +
                 theme(legend.text = element_text(size = 12))


# summarize the proportion of spots in each newly-labeled cluster
spots_sum = slide.integrated@meta.data %>% group_by(cell_type) %>%
  tally()
spots_sum = spots_sum %>% mutate(percentage = n/sum(n))

# barplot to show the proportion of young/old spots in each cell type
p10 = spots_sum %>%
      ggplot(aes(cell_type,percentage, fill = cell_type)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bar() +
      scale_fill_manual(values = colors) +
      labs(y = "Proportion")
  #     theme(legend.position = "none",
  #           axis.text.x.bottom = element_blank(),
  #           axis.text.y.left = element_blank(),
  #           axis.title.x.bottom = element_blank(),
  #           axis.title.y.left = element_blank()) +
  #     scale_y_continuous(labels = ~ scales::percent(.x)) #+
  # #theme(axis.text.x.bottom = element_text(angle = 45))


# check cell type abundance
fisher_cellType = fisher_test(slide.integrated@meta.data,
                              group.by.1 = "cell_type",
                              group.by.2 = "group")


# visualize the cell type abundance result
p11 = plotOR2FDR(fisher_cellType, cols = colors, fill.by = "cell_type",
                 point.size = 3)

ggsave(p11, filename = "combined cluster abundance.tiff",
       height = 2000, width = 2000, units = "px")


# export figure 3
supple_figure1 = wrap_plots(list(p8, p10))
ggsave(supple_figure1, filename = "supple_figure1.tiff", height = 2000,
       width = 4000, units = "px"
      )
