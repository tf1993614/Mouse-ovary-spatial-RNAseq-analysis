#load libraries
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(readbitmap)
library(hdf5r)
library(grid)
library(Matrix)
library(rjson)


# load Visium data into R as showed in 10x genomics website
# https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/rkit

sample_names = str_c("Sample", as.character(seq_len(4)))

# list paths to corresponding files
image_path = list.files(pattern = ".png$", recursive = T, full.names = T)

scalefactor_path = list.files(pattern = ".json$", recursive = T, full.names = T)

tissuse_path = list.files(pattern = "positions_list.csv$", recursive = T, full.names = T)

cluster_path = list.files(pattern = "clusters.csv$", recursive = T, full.names = T)

matrix_path = list.files(pattern = ".h5$", recursive = T, full.names = T)

# read in down sampled images
images_cl = map(image_path, ~ read.bitmap(.x))

height = purrr::reduce(map(images_cl, function(x) data.frame(height = nrow(x))), rbind)

width = purrr::reduce(map(images_cl, function(x) data.frame(width = ncol(x))), rbind)

# convert images to grobs
grobs = map(images_cl, ~ rasterGrob(.x, width = unit(1, "npc"), height = unit(1, "npc")))

images_tibble = tibble(sample = factor(sample_names), grob = grobs) %>% 
  mutate(height = height$height, width = width$width)

scales = map(seq_len(4), ~ rjson::fromJSON(file = scalefactor_path[.x]))

# combine tissue information 
bcs = map(seq_len(length(sample_names)), function(x){
  tissue_position = read.csv(tissuse_path[x], col.names = c("barcode", "tissue", 
                                                            "row", "col", "imagerow", 
                                                            "imagecol"))
  
  tissue_position$imagerow = tissue_position$imagerow * scales[[x]]$tissue_lowres_scalef
  
  tissue_position$imagecol= tissue_position$imagecol * scales[[x]]$tissue_lowres_scalef
  
  tissue_position$tissue = as.factor(tissue_position$tissue)
  
  tissue_position$height = height$height[x]
  
  tissue_position$width = width$width[x]
  
  return(tissue_position)
  
})

names(bcs) = sample_names

# read in the matrix, barcodes and genes
matrix = map(seq_len(length(sample_names)), ~ as.data.frame(t(Read10X_h5(matrix_path[.x]))))

# make summary of the matrix

# total umi per spot
umi_sum = map(seq_len(length(sample_names)), 
              ~ data.frame(barcode = rownames(matrix[[.x]]),
                sum_umi = Matrix::rowSums(matrix[[.x]])))

names(umi_sum) = sample_names

umi_sum = bind_rows(umi_sum, .id = "sample")

# total genes per spot
gene_sum = map(seq_len(length(sample_names)), 
               ~ data.frame(barcode = rownames(matrix[[.x]]),
                            sum_gene = Matrix::rowSums(matrix[[.x]] != 0)))

names(gene_sum) = sample_names

gene_sum = bind_rows(gene_sum, .id = "sample")

# merge all necessary data
bcs_merge = bind_rows(bcs, .id = "sample")
bcs_merge = bcs_merge %>% inner_join(umi_sum, by = c("sample" = "sample", 
                                                     "barcode" = "barcode"))
bcs_merge = bcs_merge %>% inner_join(gene_sum, by = c("sample" = "sample", 
                                                     "barcode" = "barcode"))

saveRDS(bcs_merge, file = "bcs_merge.rds")

# divide each image by young and  aged sections
# for sample1(ends with 655): young is on the top and aged is on the bottom
# for sample2-4 (656,657,658): youngd is on the left amd aged is on the right
bcs_merge = bcs_merge %>% mutate(group = case_when(sample == "Sample1" & row < 37 ~ "aged",
                                                   sample == "Sample2" & col > 64 ~ "aged",
                                                   sample == "Sample3" & col > 74 ~ "aged",
                                                   sample == "Sample4" & row %in% seq(0, 30, by = 2) & col > 56 ~ "aged",
                                                   sample == "Sample4" & row %in% seq(1, 29, by = 2) & col > 55 ~ "aged",
                                                   sample == "Sample4" & row %in% seq(32, 76, by = 2) & col > 54 ~ "aged",
                                                   sample == "Sample4" & row %in% seq(31, 77, by = 2) & col > 55 ~ "aged",
                                                   T ~ "young"
                                                   ))

# add in metadata information about mouse
bcs_merge = bcs_merge %>% mutate(mouse = case_when(sample == "Sample1" & group == "aged" ~ "A1",
                                                   sample == "Sample2" & group == "aged" ~ "A2",
                                                   sample == "Sample3" & group == "aged" ~  "A3",
                                                   sample == "Sample4" & group == "aged" ~ "A4",
                                                   sample == "Sample1" & group == "young" ~ "Y1",
                                                   sample == "Sample2" & group == "young" ~ "Y2",
                                                   sample == "Sample3" & group == "young" ~ "Y3",
                                                   sample == "Sample4" & group == "young" ~ "Y4",))

# delineate each metadata set for later use
D21_655_meta = bcs_merge %>% filter(sample == "Sample1")
D21_656_meta = bcs_merge %>% filter(sample == "Sample2")
D21_657_meta = bcs_merge %>% filter(sample == "Sample3")
D21_658_meta = bcs_merge %>% filter(sample == "Sample4")
meta = list(D21_655_meta, D21_656_meta, D21_657_meta, D21_658_meta) 

# load each Visium slide
slide.list = map2(str_c("sample", 1:4), 
                  c("D655", "D656", "D657", "D658"), function(x, y){
                    data = Load10X_Spatial(str_c("GSE188257_RAW/", x),
                                   filename = "filtered_feature_bc_matrix.h5",
                                   slice = y)
                    data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^mt-")
                    
                    return(data)
                  }
                 )

names(slide.list) = c("D655", "D656", "D657", "D658")

# add in metadata to Seurat object
new.slide.list = map(1:4, function(x){
  seurat = slide.list[[x]]@meta.data
  seurat = seurat %>% rownames_to_column("barcode")
  metadata = meta[[x]] %>% left_join(seurat, by = "barcode") %>% 
    column_to_rownames("barcode")
  
  slide.list[[x]] = AddMetaData(slide.list[[x]], metadata)
  
})

names(new.slide.list) = c("D655", "D656", "D657", "D658")

# subset out specific non-tissue spots (< 25% tissue coverage)
clean.slide.list = map(names(new.slide.list), function(x){
  if(x == "D655"){
    data = new.slide.list[[x]] %>% subset(((row != 46 | col != 52) &
                                             (row != 43 | col != 15) &
                                             (row != 47 | col != 53) &
                                             (row != 49 | col != 55) &
                                             (row != 48 | col != 58) &
                                             (row != 48 | col != 60) &
                                             (row != 48 | col != 54) &
                                             (row != 47 | col != 59) &
                                             (row != 46 | col != 54) &
                                             (row != 48 | col != 52)))
  }
  else if(x == "D656"){
    data = new.slide.list[[x]] %>% subset(((row != 22 | col != 22) &
                                             (row != 60 | col != 18) &
                                             (row != 60 | col != 20) &
                                             (row != 72 | col != 20) &
                                             (row != 64 | col != 30) &
                                             (row != 40 | col != 30) &
                                             (row != 41 | col != 31) &
                                             (row != 49 | col != 41) &
                                             (row != 52 | col != 42) &
                                             (row != 52 | col != 16) &
                                             (row != 52 | col != 18) &
                                             (row != 53 | col != 23) &
                                             (row != 53 | col != 25) &
                                             (row != 57 | col != 35) &
                                             (row != 40 | col != 96) &
                                             (row != 69 | col != 99) &
                                             (row != 44 | col != 102) &
                                             (row != 40 | col != 94) &
                                             (row != 38 | col != 98) &
                                             (row != 39 | col != 99) &
                                             (row != 70 | col != 100) &
                                             (row != 39 | col != 93) &
                                             (row != 9 | col != 101) &
                                             (row != 39 | col != 95) &
                                             (row != 38 | col != 92) &
                                             (row != 38 | col != 96) &
                                             (row != 40 | col != 98) &
                                             (row != 40 | col != 108) &
                                             (row != 38 | col != 94) &
                                             (row != 39 | col != 97) &
                                             (row != 69 | col != 101) &
                                             (row != 0 | col != 4)))
  }
  
  else if(x == "D657"){
    data = new.slide.list[[x]] %>% subset(((row != 41 | col != 107) &
                                      (row != 10 | col != 126) &
                                      (row != 9 | col != 125) &
                                      (row != 7 | col != 127) &
                                      (row != 9 | col != 127) &
                                      (row != 8 | col != 126) &
                                      (row != 42 | col != 106) &
                                      (row != 41 | col != 105) &
                                      (row != 52 | col != 48) &
                                      (row != 47 | col != 63) &
                                      (row != 77 | col != 63) &
                                      (row != 51 | col != 47) &
                                      (row != 32 | col != 40) &
                                      (row != 77 | col != 67) &
                                      (row != 9 | col != 21)))
  }
  else{
    data = new.slide.list[[x]] %>% subset(((row != 62 | col != 98) &
                                      (row != 5 | col != 103) &
                                      (row != 45 | col != 73) &
                                      (row != 63 | col != 99) &
                                      (row != 42 | col != 68) &
                                      (row != 44 | col != 70) &
                                      (row != 16 | col != 28) &
                                      (row != 16 | col != 32) &
                                      (row != 51 | col != 47) &
                                      (row != 16 | col != 30) &
                                      (row != 65 | col != 13) &
                                      (row != 29 | col != 23) &
                                      (row != 15 | col != 31) &
                                      (row != 10 | col != 14) &
                                      (row != 20 | col != 0)))
  }
  
})

names(clean.slide.list) = c("D655", "D656", "D657", "D658")

# intergrate datasets together as show in 
# https://satijalab.org/seurat/articles/integration_introduction.html

# setup the seurat object list, and run SCTransform on each object individually
options(future.globals.maxSize = 4000 * 1024 ^2)

final.slide.list = map(clean.slide.list, 
                       ~ SCTransform(.x, assay = "Spatial", verbose = F,
                                     vars.to.regress = "percent.mt"))

# select features for downstream integration, and prepare to integrate datasets
slide.features = SelectIntegrationFeatures(object.list = final.slide.list, nfeatures = 3000)
final.slide.list = PrepSCTIntegration(object.list = final.slide.list, anchor.features = slide.features,
                                      verbose = T)

# Identify anchors and integrate the datasets
slide.anchors = FindIntegrationAnchors(object.list = final.slide.list,
                                       normalization.method = "SCT",
                                       anchor.features = slide.features)

slide.integrated = IntegrateData(anchorset = slide.anchors, normalization.method = "SCT")

# save slide.intergrated object
saveRDS(slide.integrated, "slide.integrated.rds")
