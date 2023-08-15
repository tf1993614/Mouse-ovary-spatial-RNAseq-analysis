library(ggprism)
library(tidyverse)

# create a theme type for UMAP plotted
theme_UMAP =  function(){
  theme(
    axis.title.x.bottom = element_text(face = "bold", size = 20),
    text = element_text(family = "serif"),
    axis.title.y.left = element_text(face = "bold", size = 18),
    axis.line.x.bottom = element_line(size = 0.8),
    axis.line.y.left = element_line(size = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    axis.ticks.length = unit(0.1, "inch"),
    axis.text.x.bottom = element_text(colour = "black", face = "bold"),
    axis.text.y.left = element_text(colour = "black", face = "bold"),
   plot.title = element_blank())
  }

# create a theme type for barplot
theme_bar = function(){
    theme_classic() +
      theme(
        text = element_text(family = "sans"),
        axis.title.y.left = element_text(face = "bold", size = 18),
        axis.line.x.bottom = element_line(size = 0.8),
        axis.line.y.left = element_line(size = 0.8),
        axis.title.x.bottom = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(0.1, "inch"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 15),
        axis.text.y = element_text(colour = "black", size = 12))
  }

scale_bar = function(){
    scale_fill_manual(values = colors) + labs(y = "Proportion") +
    scale_y_continuous(labels = ~ scales::percent(.x))
  }

# create a fisher exact function to check cluster abundance
fisher_test = function(data, group.by.1, group.by.2){

  group.by.1 = sym(group.by.1)
  group.by.2 = sym(group.by.2)

  data = data %>% group_by(!!group.by.1, !!group.by.2, .drop = F) %>%
    summarise(n = n()) %>% group_by(!!group.by.2) %>% mutate(prop = n/sum(n))

  aging_cluster = data %>% subset(group == "aged") %>% mutate(remaining = sum(n) - n)
  young_cluster = data %>% subset(group == "young") %>% mutate(remaining = sum(n) - n)

  clusters = aging_cluster %>% left_join(young_cluster, by = as.character(group.by.1))

  clusters = clusters %>% mutate(!! group.by.1 := as.character(!!group.by.1))


  fisher = pmap_dfr(clusters, function(...){

    data = c(...)
    name = names(data)

    tbl = matrix(as.numeric(c(data[3], data[5], data[7], data[9])), byrow = T, ncol = 2)


    fisher_res = fisher.test(tbl, alternative = "two.sided")

    data = tibble(!! sym(name[1]) := data[1],
                  group.x = data[2],
                  n.x = as.numeric(data[3]),
                  prop.x = as.numeric(data[4]),
                  remaining.x = as.numeric(data[5]),
                  group.y = data[6],
                  n.y = as.numeric(data[7]),
                  prop.y = as.numeric(data[8]),
                  remaining.y = as.numeric(data[9]),
                  oddRatio = fisher_res$estimate,
                  pValue = fisher_res$p.value)

  })

  fisher = fisher%>%
    mutate(FDR = p.adjust(pValue, "fdr"),
           logOddRatio = log10(oddRatio),
           logFDR = -log10(FDR))

}



# create a function to visualize cluster abundce result
plotOR2FDR = function(data, cols = NULL, point.size = 2, fill.by = "seurat_clusters", legend = F,
                      label = T, label.size = 10, label.font = "bold"){
  fill.by = sym(fill.by)
  p = data %>% ggplot(aes(logOddRatio, logFDR, col = !!fill.by)) +
    geom_point(position = "jitter", size = point.size) + theme_classic() +         theme_UMAP() +
    theme(axis.text.x.bottom = element_text(size = 13),
          axis.text.y.left = element_text(size = 13)) +
    labs(x = "Log10(odds ratio)", y = "-Log10(FDR)") +
    geom_vline(xintercept = c(-0.176, 0.176), linetype = 2) +
    geom_hline(yintercept = 1.3, linetype = 2)


  if(!is.null(cols)){
    p = p +scale_color_manual(values = cols)
  }
  if(!legend){
    p = p + theme(legend.position = "none")
  }
  if(label){
    p = p + geom_label_repel(aes(label = !!fill.by, size = label.size, fontface = label.font
                              ), max.overlaps = 80)
  }
  return(p)
}

# create a bar plot function
plot_bar = function(data, features,
                    p_val_df = NULL,
                    colors = NULL,
                    show.legend = F){

  data = data %>% dplyr::filter(GeneName %in% features) %>%
    mutate(type = factor(type, levels = c("young", "age"))) %>%
    mutate(GeneName = factor(GeneName, levels = features))


  p = data %>%
    ggplot(aes(GeneName, mean)) +
    geom_bar(stat = "identity", size = 1,
             width = 0.7,
             colour = "black",
             position = "dodge",
             aes(fill = type)) +
    geom_errorbar(aes(ymin = mean, ymax = mean + SEM, fill = type),
                  size = 1,
                  width =0.4,
                  position = position_dodge(width=0.7)) +
    #scale_y_continuous(limits = c(0, NA)) +
    theme_classic() +
    labs(y = "Normalized expression") +
    theme(axis.title.x.bottom = element_blank(),
          legend.title = element_blank(),
          axis.text.y.left = element_text(colour = "black"),
          axis.text.x.bottom = element_text(face = "bold.italic", colour = "black", ),
          axis.title.y.left = element_text(face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
          axis.ticks = element_line(colour = "black", linewidth = 1),
          axis.ticks.length = unit(0.1, "inch"),
          legend.position = "none",
          axis.line.y.left = element_line(size = 0.8),
          axis.line.x.bottom = element_line(size = 0.8))

  if(! is.null(colors)){
    p = p + scale_fill_manual(values = colors)
  }

  if(!is.null(p_val_df)){
    p = p + add_pvalue(p_val_df)
  }

  if(! show.legend){
    p = p + theme(legend.position = "none")
  }

  return(p)
}

# create volcano plot function
plot_volcano = function(data,
                        label = F,
                        colors = c(Down = "#00468B99", `No difference` = "#ADB6B699", Up = "#ED000099"),
                        title = NULL,
                        save = F,
                        xintercept = c(-0.5, 0.5), yintercept = -log10(0.05),
                        format = ".png",
                        fig.wdith = 2000,
                        fig.height = 2000
){

  if(is.null(rownames(data))){
    data = data %>% rownames_to_column("GeneName")
  }

  data = data %>% mutate(FDR = p.adjust(p_val,"fdr"), DE = case_when(avg_log2FC > 0.5 & FDR < 0.05  ~ "Up",
                                                                     avg_log2FC < -0.5 & FDR < 0.05  ~ "Down",
                                                                     TRUE ~ "No difference"))



  p = data %>%
    ggplot(aes(avg_log2FC, -log10(FDR), col = DE)) +
    geom_point(position = "jitter") +
    geom_hline(yintercept = yintercept, linetype = "dashed", size =0.8) +
    geom_vline(xintercept = xintercept, linetype = "dashed", size = 0.8) +
    theme_classic() +
    scale_color_manual(values = colors) +
    theme(text = element_text(family = "serif"),
          axis.title.x.bottom = element_text(size =13, face = "bold"),
          axis.text.x.bottom = element_text(face = "bold"),
          axis.title.y.left = element_text(size =13, face = "bold"),
          axis.text.y.left = element_text(face = "bold"),
          axis.line.x.bottom = element_line(size = 0.8),
          axis.line.y.left = element_line(size = 0.8),
          legend.title = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 1),
          axis.ticks.length = unit(0.1, "inch"),
          legend.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 3,
                                    face = "bold", size = 15)) +
    ggtitle(title) +
    labs(y = "-Log10(FDR)", x = "Log2(Fold Change)")

  if(label){
    p  = p + geom_label_repel(data = data %>% subset(DE != "No difference"), aes(label = GeneName),
            max.overlaps = 20)
  }
  if(save){
    ggsave(p, width = fig.wdith, height = fig.height, filename = str_c(title, format),
           units = "px")
  }
  return(p)
}

# scatter plot function
scatter_plot = function(data, x, y,
                        colors = c("Discard" = "lightgray",
                                    "Keep" ="#ED0000FF"),
                        point.size = 3,
                        group = NULL){
  data %>%
    ggplot(aes({{ x }}, {{ y }}, color = {{ group }})) +
    geom_point(size = point.size) +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(
      axis.title.x.bottom = element_text(size =13, face = "bold"),
      axis.text.x.bottom = element_text(face = "bold", colour = "black"),
      axis.title.y.left = element_text(size =13, face = "bold"),
      axis.text.y.left = element_text(face = "bold", colour = "black"),
      axis.line.x.bottom = element_line(size = 0.8),
      axis.line.y.left = element_line(size = 0.8),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, vjust = 3,
                                    face = "bold", size = 15),
      panel.grid.major = element_line(colour = "grey", linewidth = 0.5,
                                          linetype = 2),
      panel.grid.minor = element_line(colour = "grey", linewidth = 0.5,
                                          linetype = 2),
      axis.ticks = element_line(colour = "black", linewidth = 1),
      axis.ticks.length = unit(0.1, "inch")
    )
}


# create a function to show the PCA result
PCA_analysis = function(DGElistObject,
                        transpose = TRUE,
                        features = NULL,
                        scale = FALSE,
                        plot = FALSE,
                        group = NULL,
                        label = FALSE,
                        point.size = 4,
                        text.size = 6,
                        reName_samples = NULL,
                        colors = c("#E64B35FF","#4DBBD5FF", "#00A087FF"),
                        save = FALSE,
                        format = ".png",
                        width = 2000,
                        height = 2000,
                        fileName = "pca",
                        show.legend = T){

  if(is.null(features)){
    features = rownames(DGElistObject)
  }

  if(transpose){
    pca_data = DGElistObject  %>% cpm(log = T) %>%
      subset(rownames(.) %in% features) %>%
      t() %>% prcomp(scale. = scale, center = T)
  }
  else{
    pca_data = DGElistObject %>% cpm(log = T) %>%
      subset(rownames(.) %in% features) %>%
      prcomp(scale. = scale, center = T)
  }

  pca_summary = summary(pca_data)$importance %>% as.data.frame()


  samples_attr = DGElistObject$samples %>% rownames_to_column("sampleName") %>%
    dplyr::select(sampleName, group)


  if(transpose){
    pca_data = pca_data$x %>% as.data.frame() %>%
      rownames_to_column("sampleName") %>%
      left_join(samples_attr, by = "sampleName")

    if(! is.null(reName_samples)){
      pca_data[["sampleName"]] = reName_samples
    }

    pca_f =  pca_data %>%
      ggplot(aes(PC1,PC2, color = group)) + geom_point(size = point.size) +
      labs(x = str_c("PC1 (", as.character(scales::percent(pca_summary[2,1])), ")"),
           y = str_c("PC2 (", as.character(scales::percent(pca_summary[2,2])), ")")) +
      scale_color_manual(values = colors) +
      theme_classic() +
      theme(axis.title.x.bottom = element_text(size =13, face = "bold"),
            axis.text.x.bottom = element_text(face = "bold", colour = "black"),
            axis.title.y.left = element_text(size =13, face = "bold"),
            axis.text.y.left = element_text(face = "bold", colour = "black"),
            axis.line.x.bottom = element_line(size = 0.8),
            axis.line.y.left = element_line(size = 0.8),
            legend.title = element_blank(),
            legend.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, vjust = 3,
                                      face = "bold", size = 15),
            panel.grid.major = element_line(colour = "grey", linewidth = 0.5,
                                            linetype = 2),
            panel.grid.minor = element_line(colour = "grey", linewidth = 0.5,
                                            linetype = 2))
    if(! show.legend){
      pca_f = pca_f + theme(legend.position = "none")
    }
  }
  else{
    pca_data = pca_data$x %>% as.data.frame() %>% rownames_to_column("sampleName")

    pca_f =  pca_data %>%
      ggplot(aes(PC1,PC2)) + geom_point() +
      labs(x = str_c("PC1 (", as.character(scales::percent(pca_summary[2,1])), ")"),
           y = str_c("PC2 (", as.character(scales::percent(pca_summary[2,2])), ")"))
  }

  if(label){
    pca_f = pca_f + geom_label_repel(aes(label = sampleName), max.overlaps = 20,
                                     size = text.size)
  }

  if(save){
    ggsave(pca_f, width = width, height = height, units = "px",
           filename = str_c(fileName, format))
  }

  if(plot){
    return(pca_f)
  }
  else{
    return(pca_data)
  }

}
