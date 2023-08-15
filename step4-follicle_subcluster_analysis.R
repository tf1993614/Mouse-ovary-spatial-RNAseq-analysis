# load necessary packages
library(Seurat)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
library(ggprism)
library(ggrepel)
library(patchwork)

# prepare colors for visualization
colors = c(pal_lancet("lanonc")(9), 
           pal_npg("nrc")(10)[c(1:7,9,10)], 
           pal_nejm("default")(8))

# read Seurat object in
slide.integrated = readRDS("./slide.integrated.after.cell.type.annotation.rds")

# subset data in to follicle subclass
Idents(slide.integrated) = "cell_type"
follicle = subset(slide.integrated, idents = "Follicle")

saveRDS(follice, "./follicle.rds")

# process with downstream analysis on the integrated dataset
DefaultAssay(follicle) == "integrated"
follicle = RunPCA(follicle, npcs = 30, assay = "integrated")
follicle = FindNeighbors(follicle, reduction = "pca", dims = 1:30, k.param = 20)
follicle = FindClusters(follicle, resolution = 0.5)
follicle = RunUMAP(follicle, reduction = "pca", dims = 1:30)

p14 = UMAPPlot(follicle, cols = colors) + theme_UMAP() +
  theme(legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

p15 = UMAPPlot(follicle, cols = colors, group.by = "group") +
  theme_UMAP() +
  scale_color_manual(values = c(aged = colors[2], young = colors[1])) +
  theme(legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

figure3 = wrap_plots(p14, p15)

ggsave(figure3, filename = "figure3.tiff", height = 2000,
       width = 4000, units = "px")

# check cluster abudance in the newly-defined cluster in follicle identity
fisher_follicle = fisher_test(follicle@meta.data,
                              group.by.1 = "seurat_clusters",
                              group.by.2 = "group")

# visualize the cluster abundance result
p16 = plotOR2FDR(fisher_follicle, cols = colors, point.size = 3,
                 label.size = 10)
p16 = p16 + theme(legend.position = "none",
                  axis.text.x.bottom = element_blank(),
                  axis.text.y.left = element_blank(),
                  axis.title.x.bottom = element_blank(),
                  axis.title.y.left = element_blank())

ggsave(p16, filename = "cluster abundance in follicle identity.tiff",
       height = 2000, width = 2000, units = "px")

test = plotOR2FDR(fisher,  cols = colors, point.size = 3,
                  label.size = 10) + theme(legend.position = "none",
                                           axis.text.x.bottom = element_blank(),
                                           axis.text.y.left = element_blank(),
                                           axis.title.x.bottom = element_blank(),
                                           axis.title.y.left = element_blank())

ggsave(test, filename = "global cluster abundance.tiff",
       height = 2000, width = 2000, units = "px")

# map cell-type spots onto hex-staining image
map_plot_2 = map(c("D655", "D656", "D657", "D658"),
                  ~ SpatialDimPlot(follicle, images = .x,
                                   crop = F,
                                   pt.size.factor = 1.5)
                )

# get the original tissue image
map_plot_2_bg = map(c("D655", "D656", "D657", "D658"),
                 ~ SpatialDimPlot(follicle, images = .x,
                                  crop = F,
                                  alpha = 0,
                                  pt.size.factor = 1.5)
)

# remove legend
map_plot_2_bg = map(map_plot_2_bg, ~ .x + theme(legend.position =  "none") +
                      scale_fill_manual(values = colors))

# change the theme settings of those plots
map_plot_2[[1]] = map_plot_2[[1]] +
                 theme(legend.position = "none") +
                 scale_fill_manual(values = colors)

map_plot_2[[2]] = map_plot_2[[2]] +
                  theme(legend.position = "none") +
                  scale_fill_manual(values = colors)

map_plot_2[[3]] = map_plot_2[[3]] +
                  theme(legend.position = "none") +
                  scale_fill_manual(values = colors)

map_plot_2[[4]] = map_plot_2[[4]] +
                  theme(legend.position = "none") +
                  scale_fill_manual(values = colors)

walk(seq(4), ~ ggsave(map_plot_2[[.x]], filename = str_c("supple figure 2A", .x, ".tiff"),
                      height = 2000, width = 2000, units = "px"))

walk(seq(4), ~ ggsave(map_plot_2_bg[[.x]], filename = str_c("supple figure 2B", .x, ".tiff"),
                      height = 2000, width = 2000, units = "px"))

# find marker genes in each cluster
Idents(follicle) = "seurat_clusters"
follicle_markers = FindAllMarkers(follicle, only.pos = T,
                                  min.pct = 0.3,
                                  logfc.threshold = 0)

follicle_markers = follicle_markers %>%
                   group_by(cluster) %>%
                   group_split(.keep = T)

follicle_markers = map(follicle_markers, ~ .x %>% arrange(desc(avg_log2FC)))

names(follicle_markers) = str_c("cluster_", seq(0,5))

# select top60 marker genes for each cluster
top60_follicle_markers = follicle_markers %>%
                         group_by(cluster) %>%
                         slice_head(n = 60) %>%
                         group_split(.keep = T)

# name the list
names(top60_follicle_markers) = str_c("cluster_", seq(0,5))

# export whole marker genes list
writexl::write_xlsx(follicle_markers,
                    path = "marker genes for each subcluster in the follicle identity.xlsx")

# define the marker genes for each cell type
granulosa_markers = c("Cyp19a1", "Fshr", "Nr5a2", "Mgarp",
                      "Gldc", "Chst8", "Csn2", "Gpx3", "Slc35g1",
                      "Ca8", "Clgn", "Fam78a", "Slc16a3")

theca_markers = c("Mgp", "Dcn", "Aspn", "Aldh1a1", "Col1a2", "Fn1",
                  "Col3a1", "Ogn", "Apod", "Col5a2", "Igf2", "Nid1",
                  "Lhfp", "Acta2", "Dusp12", "Actg2", "Sparcl1",
                  "Filip1l", "Egflam", "Adamdec1", "Hpcgd", "Col12a1",
                  "Fbln5", "Ramp2", "Col15a1", "Plk2", "Col6a3",
                  "Loxl1", "Rarres1", "Fli1", "Lama2", "Insl3")

oocyte_markers = c("Gdf9", "Slc39a10", "Ddx4", "Nalp9", "Nalp5", "Zp3",
                   "Bmp15", "Zp2", "Sycp3", "Sox30", "Zar1", "Dazl",
                   "Ybx2", "Lhx8")

germ_markers = c("Oct4", "Blimp1", "Dazl", "Figla", "Fragilis",
                 "Mvh", "Par6", "Pou5f1", "Stat3", "Afp",
                 "Mater")

combination_markers = c(granulosa_markers, theca_markers,
                        oocyte_markers, germ_markers)

# get the average gene expression matrix for each cluster
Idents(follicle) = follicle$seurat_clusters
follicle_average = as.data.frame(AverageExpression(follicle,
                                                   assays = "Spatial",
                                                   slot = "counts"))

# do log1p transformation
follicle_average = log1p(follicle_average)

# filter out zero-expression genes in all clusters
# and only select genes we manually define above
follicle_combination = follicle_average %>%
  subset(rownames(.) %in% combination_markers) %>%
  dplyr::filter(!if_all(everything(), ~ .x == 0)) %>%
  rename_with(str_replace, pattern = "Spatial\\.",
              replacement = "cluster_", .cols = everything())

# heatmap to show those genes expression
tiff(filename = "figure3C.tiff", width = 2000, height = 2000, units = "px")
p17 = pheatmap::pheatmap(follicle_combination,
                         cluster_rows = T,
                         angle_col = 0,
                         cluster_cols = T,
                         scale = "row",
                         legend_breaks = c(-1.5, 0, 1.5),
                         colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(-1.5, 1.5, by= 0.1))),
                         breaks = seq(-1.5, 1.5, by= 0.1))
p17
dev.off()

# do normalization on "spatial" assay
follicle = NormalizeData(follicle, assay = "Spatial")

# retrieve spots for each subcluster
seurat_list = map(seq(0, 5), ~ subset(follicle, seurat_clusters == .x))

# do DE analysis between young and old samples in 
# each subcluster(0-5)
subcluster_DE_list = map(seq(0, 5),
                         function(x){
                           data = subset(follicle, seurat_clusters == x)
                           DE_res = FindMarkers(
                                                data, 
                                                ident.1 = "aged",
                                                assay = "Spatial",
                                                ident.2 = "young",
                                                group.by = "group",
                                                min.pct = 0.3, 
                                                logfc.threshold = 0)
                           DE_res = DE_res %>% 
                                    mutate(FDR = p.adjust(p_val, "fdr")) %>%
                                    rownames_to_column("GeneName")
                         })

subcluster_DE_list[[5]] %>% 
  mutate(DE = case_when(avg_log2FC > 0.5 & FDR < 0.05 ~ "Up",
                        avg_log2FC < -0.5 & FDR < 0.05 ~ "Down",
                        T ~ "No difference")) %>%
  count(DE)

# export DE analysis results
names(subcluster_DE_list) = str_c("cluster_", seq(0,5))
writexl::write_xlsx(subcluster_DE_list, 
                    path = "DE genes comparing old versus young samples in each subcluster of follicle identity.xlsx")
  
wrap_plots(p20_list[[2]][2:5])

# visualize the DE result by volcano plot
p18_list = map(subcluster_DE_list, ~ plot_volcano(.x, save = F, label = F))

p18_list = map(p18_list, ~ .x + 
                 theme(legend.position = "none",
                       axis.text.x.bottom = element_blank(),
                       axis.text.y.left = element_blank(),
                       axis.title.x.bottom = element_blank(),
                       axis.title.y.left = element_blank()))

walk(seq(6), ~ ggsave(p18_list[[.x]], 
                      filename = str_c("cluster_", .x-1, ".png"),
                      height = 2000,
                      width = 2000,
                      units = "px"))

# create vectors containing antioxidnat genes
interest_genes.1 = c("ACTB","RPL19", "RPL32", "FDXR", "FDX1", "CYP11A1",
                      "POR", "CYP19A1")

interest_genes.2 = c("ACTB", "RPL19", "RPL32", "GPX1", "GPX2", "GPX3", "GPX4",
                     "GPX5", "GPX6", "GPX7", "GPX8", "LRP8", "LRP2")

interest_genes.3 = c("ACTB", "RPL19", "RPL32", "SOD1", "SOD2", "CAT", "GSR")

interest_genes.4 = c("ACTB", "RPL19", "RPL32", "TXN", "TXN2",
                    "TXNRD1", "TXNRD2", "TXNRD3")

interest_genes.5 = c("ACTB", "RPL19", "RPL32", "PRDX1", "PRDX2",
                    "PRDX3", "PRDX4","PRDX5", "PRDX6")

interest_genes.6 = c("ACTB", "RPL19", "RPL32", "NOX1", "NOX2",
                     "NOX3", "NOX4", "NOX5")

interest_all = union(interest_genes.1,
                    c(interest_genes.2, interest_genes.3,
                      interest_genes.4, interest_genes.5,
                      interest_genes.6)) %>% str_to_sentence()

# get the counts matrix
ROS_matrix = map2(seurat_list, subcluster_DE_list, 
                 ~ .x@assays$Spatial@data %>% as.data.frame() %>%
                 dplyr::filter(rownames(.) %in% interest_all) %>%
                 rownames_to_column("GeneName") %>%
                 dplyr::filter(GeneName %in% .y$GeneName)
                 )

# prepare cell labels
labels = map(seurat_list,
             ~ tibble(x = .x$group) %>%
             mutate(y = row_number()) %>%
             ungroup() %>%
             mutate(label = str_c(x, "_", y)) %>%
             .$label
            )

# substitute cell id with cell label
ROS_matrix = map2(ROS_matrix, labels, function(x, y){
                data = x
                names(data) = c("GeneName", y)
                return(data)
                })

# calculate the mean and SEM of antioxidant genes expression in
# either age or young group
ROS_matrix = map(ROS_matrix,
                 ~ .x %>%
                    mutate(mean_age = pmap_dbl(.[str_which(names(.x), "aged")], ~ mean(c(...))),
                           SEM_age  = pmap_dbl(.[str_which(names(.x), "aged")], ~ sd(c(...))/sqrt(length(c(...)))),
                           mean_young = pmap_dbl(.[str_which(names(.x), "young")], ~ mean(c(...))),
                           SEM_young  = pmap_dbl(.[str_which(names(.x), "young")], ~ sd(c(...))/sqrt(length(c(...)))
                            )) %>%
                   dplyr::select(GeneName, ends_with("age"), ends_with("young")) %>%
                   mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>%
                   arrange(GeneName)
                 )

# tidy up the matrix for plotting
ROS_matrix = map(ROS_matrix,
                 ~ .x %>%
                 mutate(GeneName = factor(GeneName,
                        levels = str_sort(.$GeneName))) %>%
                 pivot_longer(-GeneName,
                              names_to = c(".value", "type"),
                              names_sep = "_")
                )

# prepare significant p value table
ROS_statiscal = map(subcluster_DE_list,
                    ~ .x %>%
                    subset(GeneName %in% interest_all) %>%
                    arrange(GeneName)
                   )

# define the p value position on the plot
y_position = map2(ROS_matrix, ROS_statiscal,
                 ~ .x %>% dplyr::filter(GeneName %in% .y$GeneName)%>%
                 group_by(GeneName) %>%
                 summarise(max_y = max(mean), SEM = SEM) %>%
                 arrange(GeneName) %>%
                 mutate(max_y = max_y + SEM + 0.1) %>%
                 group_by(GeneName) %>%
                 arrange(desc(max_y), .by_group = T) %>%
                 distinct(GeneName, .keep_all = T)
                )

# prepare significance table
df_p_val = map2(ROS_statiscal, 
                y_position,
                 ~ tibble(
                   group1 = "age",
                   group2 = "young",
                   p.adj = .x$FDR,
                   avg_log2FC = .x$avg_log2FC,
                   y.position = .y$max_y,
                   GeneName = .y$GeneName) %>%
                   mutate(p.signif = case_when(p.adj < 0.0001 ~ "****",
                                               p.adj < 0.001  ~ "***",
                                               p.adj < 0.01  ~ "**",
                                               p.adj < 0.05  ~ "*"))
               )


# visualize interested genes expression in each subcluster
p20_list = map2(ROS_matrix, 
                df_p_val, 
                function(x, y){
                  data = x
                  p_val_table = y
                  map(
                      list(
                           interest_all[c(1:8)],
                           interest_all[c(1:3, 9:16)],
                           interest_all[c(1:3, 19:22)],
                           interest_all[c(1:3, 23:27)],
                           interest_all[c(1:3, 28:38)]),
                      ~ plot_bar(data, .x, colors = c("#3C5488FF", "#E64B35FF"))
                     )
                }
)

p20_list = map(p20_list, function(data){
  map(data, ~ .x + 
        theme(legend.position = "none",
              axis.text.x.bottom = element_blank(),
              axis.text.y.left = element_blank(),
              axis.title.x.bottom = element_blank(),
              axis.title.y.left = element_blank()))
})

walk(seq(6), function(index){
    walk(seq(5), ~ ggsave(p20_list[[index]][[.x]], 
               filename = str_c("cluster_", index-1, "_P", .x , ".png"), 
               width = 2000, height = 2000,
               units = "px"))
  }
)

ROS_statiscal[[3]] %>% dplyr::filter(GeneName %in% interest_all[25:35]) %>%
  dplyr::select(-p_val_adj) %>% View()


# viusalize marker genes for newly-combined cluster in the follicle
# identity
cluster_marker_after_combing_cluster = readRDS("./consevred_marker_genes_in_newly_combined_cluster.rds")

top20_cluster_gene = cluster_marker_after_combing_cluster %>%
  bind_rows(.id = "cluster") %>% group_by(cluster) %>%
  filter(Gene_name %in% rownames(slide.integrated)) %>%
  slice_min(n = 20, order_by = minimump_p_val)

p4 = DoHeatmap(follicle, features = top20_cluster_gene$Gene_name)
ggsave(p4, filename = "newly_combined marker genes expression in follicle identity.png",
       height = 4000, width = 3000, units = "px")
