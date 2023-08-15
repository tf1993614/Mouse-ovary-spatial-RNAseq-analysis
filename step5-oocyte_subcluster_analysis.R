# load necessary packages
library(ggprism)
library(Seurat)
library(tidyverse)
library(edgeR)
library(AnnotationHub)
library(biomaRt)
library(scploid)

# Subset out spots expressing Gdf9 and Zp3 in follicle subcluster
oocyte = subset(x = follicle, subset = Gdf9 > 0 & Zp3 > 0 & seurat_clusters == 3)


# show young oocytes and old oocytes on
# the UMAP of follicle identity
young_cells = oocyte@meta.data %>% subset(group == "young") %>%
  rownames()

old_cells = oocyte@meta.data %>% subset(group == "aged") %>%
  rownames()

UMAP_matrix = follicle@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% 
  rownames_to_column("cell") %>%
  mutate(
    type = case_when(
      cell %in% young_cells ~ "Young oocyte",
      cell %in% old_cells ~ "Old oocyte",
      T ~ "No oocyte"
    )
  )

p0 = scatter_plot(UMAP_matrix, x = UMAP_1, y = UMAP_2, group = type,
             colors = c("No oocyte" = "lightgray", 
                        "Young oocyte" ="#00468BFF",
                        "Old oocyte" = "#ED0000FF"
                        )) +
  theme(legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank())

ggsave(p0, filename = "figure6A.tiff", height = 2000, 
       width = 2000, units = "px")

# Proceed with downstream analysis on the integrated dataset
DefaultAssay(oocyte) = "integrated"
oocyte = RunPCA(object = oocyte, verbose = FALSE, npcs = 10, assay = "integrated")
oocyte = FindNeighbors(oocyte, reduction = "pca", dims = 1:10, k.param = 5)
oocyte = FindClusters(oocyte, verbose = FALSE, resolution = 0.3)
oocyte = RunUMAP(oocyte, reduction = "pca", dims = 1:10)

# set "seurat_clusters" as default idents
Idents(oocyte) = "seurat_clusters"

# do normalization on "Spatial" assay
oocyte = NormalizeData(oocyte, assay = "Spatial")
oocyte = FindVariableFeatures(oocyte, assay = "Spatial")


# find DE genes between young and old oocyte spots
DE_oocyte_wilicox =  FindMarkers(oocyte,
                                 ident.1 = "aged",
                                 ident.2 = "young",
                                 assay = "Spatial",
                                 group.by = "group",
                                 min.pct = 0.3,
                                 logfc.threshold = 0)

# add FDR-corrected p value
DE_oocyte_wilicox = DE_oocyte_wilicox %>%
                    rownames_to_column("GeneName") %>%
                    mutate(FDR = p.adjust(p_val, "fdr"))

# volcano plot to show DE genes
plot_volcano(DE_oocyte_wilicox, save = T, title = "oocyteDE")


# create a list containing antioxidnat genes
interest_genes.1 = c("ACTB","RPL19", "RPL32", "FDXR", "FDX1", "CYP11A1", "POR", "CYP19A1")
interest_genes.2 = c("ACTB","RPL19", "RPL32", "GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8", "LRP8", "LRP2")
interest_genes.3 = c("ACTB","RPL19", "RPL32", "SOD1", "SOD2", "CAT", "GSR")
interest_genes.4 = c("ACTB","RPL19", "RPL32", "TXN", "TXN2", "TXNRD1", "TXNRD2", "TXNRD3")
interest_genes.5 = c("ACTB","RPL19", "RPL32", "PRDX1", "PRDX2", "PRDX3", "PRDX4","PRDX5", "PRDX6")
interest_genes.6 = c("ACTB","RPL19", "RPL32", "NOX1", "NOX2", "NOX3","NOX4", "NOX5")

interest_all = union(interest_genes.1, c(interest_genes.2, interest_genes.3,
                                         interest_genes.4, interest_genes.5, interest_genes.6))[-c(1,2,3)] %>% str_to_sentence()

# get the counts matrix for candidate genes
ROS_matrix = oocyte@assays$Spatial@data %>% as.data.frame() %>%
  dplyr::filter(rownames(.) %in% interest_all) %>%
  rownames_to_column("GeneName") %>%
  dplyr::filter(GeneName %in% DE_oocyte_wilicox$GeneName)

# substitute cell id with cell label
lables = tibble(x = oocyte$group) %>% group_by(x) %>% mutate(y = row_number()) %>%
  ungroup() %>% mutate(label = str_c(x, "_", y)) %>% .$label

names(ROS_matrix) = c("GeneName", lables)

# calculate the mean and SEM of antioxidant genes expression in
# either age or young group
ROS_matrix = ROS_matrix %>%
             mutate(mean_age = pmap_dbl(.[str_which(names(ROS_matrix), "aged")], ~ mean(c(...))),
                   SEM_age = pmap_dbl(.[str_which(names(ROS_matrix), "aged")], ~ sd(c(...))/sqrt(length(c(...)))),
                   mean_young = pmap_dbl(.[str_which(names(ROS_matrix), "young")], ~ mean(c(...))),
                   SEM_young = pmap_dbl(.[str_which(names(ROS_matrix), "young")], ~ sd(c(...))/sqrt(length(c(...))))) %>%
             dplyr::select(GeneName,ends_with("age"), ends_with("young")) %>%
             mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>%
             arrange(GeneName)

# tidy up the matrix for plotting
ROS_matrix = ROS_matrix  %>%
  mutate(GeneName = factor(GeneName, levels = str_sort(.$GeneName))) %>%
  pivot_longer(-GeneName, names_to = c(".value", "type"), names_sep = "_")


# prepare significant p value table
ROS_statiscal = DE_oocyte_wilicox %>%
  subset(GeneName %in% interest_all) %>% arrange(GeneName)

# define the p value position on the plot
y_position = ROS_matrix %>%
             dplyr::filter(GeneName %in% ROS_statiscal$GeneName) %>%
             group_by(GeneName) %>%
             summarise(max_y = max(mean), SEM = SEM) %>%
             arrange(GeneName) %>%
             mutate(max_y = max_y + SEM + 0.1) %>%
             group_by(GeneName) %>%
             arrange(desc(max_y), .by_group = T) %>%
             distinct(GeneName, .keep_all = T)

# prepare p value table
df_p_val = tibble(group1 = "age",
                  group2 = "young",
                  p.adj = ROS_statiscal$FDR,
                  avg_log2FC = ROS_statiscal$avg_log2FC,
                  y.position = y_position$max_y,
                  GeneName = y_position$GeneName) %>%
           mutate(p.signif = case_when(p.adj < 0.0001 ~ "****",
                              p.adj < 0.001  ~ "***",
                              p.adj < 0.01 ~ "**",
                              p.adj < 0.05 ~ "*")) %>%
          dplyr::select(-p.adj)


# viusalize genes expression in oocyte identity
p1 = plot_bar(ROS_matrix, features = interest_all[6:15], p_val_df = df_p_val)
p2 = plot_bar(ROS_matrix, features = interest_all[16:24], p_val_df = df_p_val)
p3 = plot_bar(ROS_matrix, features = interest_all[25:35], p_val_df = df_p_val)

# export plots
walk2(list(p1, p2, p3), c("oocyte_p1", "oocyte_p2", "oocyte_p3"),
      ~ ggsave(.x,
               filename = str_c(.y, ".png"),
               width = 3000,
               height = 3000,
               units = "px"))


## following script is used to infer aneuploidy 
## rate in young and old oocytes using the
## algorithm published in 2017.
## More details can be found in this link
## https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4253-x

# Split oocyte identity by age group
d_young = subset(oocyte, group == "young")
d_old = subset(oocyte, group == "aged")


# get the counts matrix for young oocytes
oocyte_young_counts = GetAssayData(d_young, assay = "Spatial", slot = "counts") %>%
                      as.data.frame() %>%
                      rownames_to_column("GeneName")

# rename the column names
colnames(oocyte_young_counts) = c("GeneName",
                                 str_c("young_cell_",
                                 seq(ncol(oocyte_young_counts)-1)))


# get the counts matrix for old oocytes
oocyte_old_counts = GetAssayData(d_old, assay = "Spatial", slot = "counts") %>%
                      as.data.frame() %>%
                      rownames_to_column("GeneName")

# rename the column names
colnames(oocyte_old_counts) = c("GeneName",
                                str_c("old_cell_",
                                seq(ncol(oocyte_old_counts)-1)))

# prepare file for converting gene symbol to ensembl_id
ah = AnnotationHub()
ah = ah %>% subset(species == "Mus musculus" & rdataclass == "EnsDb")
ensdb = ah[["AH100674"]]
geneGrange = genes(ensdb)
ensembl_id = mcols(geneGrange)[,c("gene_name", "gene_id")] %>% as.data.frame()

# converst gene symbol to ensembl_id
oocyte_young_counts = oocyte_young_counts %>%
  left_join(ensembl_id, by = c("GeneName" = "gene_name")) %>%
  dplyr::select(GeneName, gene_id, everything()) %>%
  dplyr::filter(! is.na(gene_id)) %>%
  dplyr::select(-GeneName) %>%
  column_to_rownames("gene_id")

oocyte_old_counts = oocyte_old_counts %>%
  left_join(ensembl_id, by = c("GeneName" = "gene_name")) %>%
  dplyr::select(GeneName, gene_id, everything()) %>%
  dplyr::filter(! is.na(gene_id)) %>%
  dplyr::select(-GeneName) %>%
  column_to_rownames("gene_id")


# identify the chromosome that each gene lies on
mouse_ensembl = useMart("ensembl")
mouse_ensembl = useDataset("mmusculus_gene_ensembl",
                            mart = mouse_ensembl)

gene_table = getBM(
    attributes = c("ensembl_gene_id", "chromosome_name"),
    mart = mouse_ensembl,
    values = as.character(rownames(oocyte_old_counts)),
    filters = "ensembl_gene_id"
)


# get genes on autosomal chromosomes
gene_table = gene_table[gene_table$chromosome_name %in% 1:19, ]

# only keep genes on autosomal chromosomes
oocyte_young_counts = oocyte_young_counts[gene_table$ensembl_gene_id, ]
oocyte_old_counts = oocyte_old_counts[gene_table$ensembl_gene_id, ]


# select most variable genes (the median TMM expression in all
# samples should be above 50) for PCA
pca_young = oocyte_young_counts %>%
            DGEList() %>%
            cpm() %>%
            as.data.frame() %>%
            dplyr::filter(pmap_lgl(.[1:ncol(oocyte_young_counts)],
                          ~ median(c(...)) >= 50))

# get a matrix with most variable genes
pca_young = oocyte_young_counts %>%
  dplyr::filter(rownames(.) %in% rownames(pca_young)) %>% DGEList()

# do PCA plot
PCA_analysis(pca_young, plot = T)

# return PCA matrix
pca_young = PCA_analysis(pca_young)

# Discard variable cells
pca_young_keep = pca_young %>% 
  subset(PC1 > -5 & PC2 < 5 & PC2 > -10)

pca_young = pca_young %>% 
  mutate(
    keep = if_else(
    sampleName %in% pca_young_keep$sampleName,
    "Keep",
    "Discard"
  )
)

# show remaining cells for down
# analysis on PCA plot
p1 = scatter_plot(pca_young, PC1, PC2, group = keep, 
             colors = c("Discard" = "lightgray", 
                        "Keep" ="#ED0000FF")) +
  labs(x = "PC1 (9%)",
       y = "PC2 (5%)") +
  theme(legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

ggsave(p1, filename = "figure6B.tiff", height = 2000, 
       width = 2000, units = "px")

# discard overdispersed cells based on PCA analysis (young oocytes)
young_keep = pca_young %>% dplyr::filter(PC1 > -5 & PC2 < 5 & PC2 > -10) %>% .$sampleName


# select most variable genes (the median TMM expression in all
# samples should be above 50) for PCA
pca_old = oocyte_old_counts %>%
          DGEList() %>%
          cpm() %>%
          as.data.frame() %>%
          dplyr::filter(pmap_lgl(.[1:ncol(oocyte_old_counts)],
                        ~ median(c(...)) >= 50))

# get a matrix with most variable genes
pca_old = oocyte_old_counts %>%
          dplyr::filter(rownames(.) %in% rownames(pca_old)) %>% DGEList()

# do PCA plot
PCA_analysis(pca_old, plot = T)

# return PCA matrix
pca_old = PCA_analysis(pca_old)

# Discard variable cells
pca_old_keep = pca_old %>% 
  subset(PC1 < 10)

pca_old = pca_old %>% 
  mutate(
    keep = if_else(
      sampleName %in% pca_old_keep$sampleName,
      "Keep",
      "Discard"
    )
  )

p2 = scatter_plot(pca_old, PC1, PC2, group = keep) +
  theme(legend.position = "none",
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

ggsave(p2, filename = "figure6C.tiff", height = 2000, 
       width = 2000, units = "px")

# discard overdispersed cells based on PCA analysis (old oocytes)
old_keep = pca_old %>% dplyr::filter(PC1 < 10) %>% .$sampleName


# keep clean cells for downstream analysis
oocyte_young_counts = oocyte_young_counts %>%
                      dplyr::select(any_of(young_keep))

oocyte_old_counts = oocyte_old_counts %>%
                    dplyr::select(any_of(old_keep))

# make ploidytest object for young and old cells respectively
ploidytest_list = map2(list(young = oocyte_young_counts,
                            old = oocyte_old_counts),
                      c("young", "old"),
                       ~ makeAneu(counts = as.matrix(.x),
                                  genes = gene_table$ensembl_gene_id,
                                  chrs = gene_table$chromosome_name,
                                  cellNames = colnames(.x),
                                  cellGroups = rep(.y, ncol(.x))))


# check cell staus by PCA
plotPCA(ploidytest_list[[1]])
plotPCA(ploidytest_list[[2]])

# get the cutoff for removing differential genes
cutoff1 = getMaxA(ploidytest_list[[1]]) %>% head() %>% .[6] %>% floor()
cutoff2 = getMaxA(ploidytest_list[[2]]) %>% head() %>% .[6] %>% floor()


# remove Differential genes within group
ploidytest_list = map2(ploidytest_list,
                       c(cutoff1, cutoff2),
                       ~ setParam(.x, param_name = "extreme.gene.thresh",
                                  param_value = .y, print = T)
)

# set p value threshold
ploidytest_list = map(ploidytest_list,  ~ setParam(.x,
                                                   param_name = "p.thresh",
                                                   param_value = 0.05))

# run the aneuploidy calculations
ploidytest_list = map(ploidytest_list, ~ doAneu(.x))


# get aneuploidy score results
score_young = getScores(ploidytest_list[[1]])
score_old = getScores(ploidytest_list[[2]])

# merge results for young and old samples
score_all = purrr::reduce(list(score_young, score_old), rbind)

# add FDR-corrected p value
score_all$p.adj = p.adjust(score_all$p, "fdr")

# tidy up the dataframe for plotting
score_all = score_all %>%
  mutate(aneuploidy = case_when(score > 1.2 & p.adj < 0.05 ~ T,
                                score < 0.8 & p.adj < 0.05 ~ T,
                                T ~ F)) %>% rownames_to_column("type") %>%
  mutate(type = str_extract(type, "young|old")) %>%
  mutate(chr = as.character(chr)) %>%
  mutate(chr = factor(chr, levels = as.character(seq(19))))


# viusalize the aneuploidy result for old oocytes
old_fig = score_all %>% subset(type == "old") %>%
  ggplot(aes(x= cell, y = chr, fill = aneuploidy)) +
  geom_tile(color = "black", size = 0.5) +
  theme_classic() +
  coord_fixed(0.75) +
  scale_fill_manual(values = c("#ADB6B699", "#ED000099"))+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y.left = element_text(colour = "black"),
        axis.text.x.bottom = element_text(angle = 45, vjust =1,
                                          hjust = 1,
                                          colour = "black"),
axis.ticks = element_line(colour = "black", linewidth = 0.5),
axis.ticks.length = unit(0.07, "inch")) +
  theme(legend.position = "none",
        axis.ticks.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())

ggsave(old_fig, filename = "old_oocyte_aneuploidy.png", units = "px",
       width = 4000, height = 3000)

# viusalize the aneuploidy result for young oocytes
young_fig = score_all %>% subset(type == "young") %>%
 ggplot(aes(x= cell, y = chr, fill = aneuploidy)) +
 geom_tile(color = "black", size = 0.5) +
 theme_classic() +
 coord_fixed(0.75) +
 scale_fill_manual(values = c("#ADB6B699", "#ED000099"))+
 theme(axis.line = element_blank(),
       axis.title = element_blank(),
       axis.text.y.left = element_text(colour = "black"),
       axis.text.x.bottom = element_text(angle = 45, vjust =1,
                                         hjust = 1,
                                         colour = "black"),
       axis.ticks = element_line(colour = "black", linewidth = 0.5),
       axis.ticks.length = unit(0.07, "inch")) +
  theme(legend.position = "none",
        axis.ticks.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank())
  
  

ggsave(young_fig, filename = "young_granulosa_aneuploidy.png", units = "px",
      width = 4000, height = 3000)
