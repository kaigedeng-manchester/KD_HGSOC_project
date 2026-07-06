# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("02.RData") 

# Load basic packages
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)
library(ggplot2)

# My basic style
theme_kd <- function(base_size = 7, base_family = "Arial", legend_pos = "bottom") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 1
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(), 
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      # 2
      axis.title = element_text(color = "black", face = "plain", size = 7, margin = margin(t = 4, r = 4)),
      axis.text = element_text(color = "black", size = 6),
      # 3
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8, margin = margin(b = 6)),
      # 4
      legend.position = legend_pos,
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(color = "black", face = "bold", size = 7),
      legend.text = element_text(color = "black", size = 6),
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(t = 0, b = 0),
      # 5
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.5),
      strip.text = element_text(color = "black", face = "bold", size = 7, margin = margin(t = 3, b = 3))
    )
}

# Load font Arial
windowsFonts(Arial = windowsFont("Arial"))

# General GSEA-----------------------------------------------------------------

# Get ranked gene list
degs_carbo <- readRDS("Data/degs_carbo.rds")
degs_taxol <- readRDS("Data/degs_taxol.rds")
function_gsea_list <- function(deg_df) {
  gene_list <- deg_df$t # using t value instead of LogFC
  names(gene_list) <- deg_df$GeneSymbol
  gene_list <- na.omit(gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}
ranked_carbo <- function_gsea_list(degs_carbo)
ranked_taxol <- function_gsea_list(degs_taxol)

# Get genesets
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
geneset_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# Run GSEA-hallmark
set.seed(202604) 
gsea_carbo <- GSEA(ranked_carbo, 
                   TERM2GENE = geneset_hallmark, 
                   pvalueCutoff  = 0.10,
                   pAdjustMethod = "BH",
                   minGSSize     = 3,      # for small genesets
                   maxGSSize     = 1000,   # for big genesets
                   verbose       = FALSE)
set.seed(202604)
gsea_taxol <- GSEA(ranked_taxol, 
                   TERM2GENE = geneset_hallmark, 
                   pvalueCutoff  = 0.10,
                   pAdjustMethod = "BH",
                   minGSSize     = 3,      # for small genesets
                   maxGSSize     = 1000,   # for big genesets
                   verbose       = FALSE)

# Clean names
library(stringr)
function_clean_gsea_names <- function(gsea_obj) {
  if (nrow(as.data.frame(gsea_obj)) > 0) {
    gsea_obj@result$Description <- gsub("HALLMARK_", "", gsea_obj@result$Description)
    gsea_obj@result$Description <- gsub("_", " ", gsea_obj@result$Description)
  }
  return(gsea_obj)
}
gsea_carbo_clean <- function_clean_gsea_names(gsea_carbo)
gsea_taxol_clean <- function_clean_gsea_names(gsea_taxol)

# Ridge plot
function_gsea_ridge <- function(gsea_obj, title_text, pathways_to_show = 15) {
  if (nrow(as.data.frame(gsea_obj)) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste("No enrichment:", title_text)))
  }
  p <- ridgeplot(gsea_obj, 
                 showCategory = pathways_to_show, 
                 fill = "NES")
  p$layers[[1]]$aes_params$alpha <- 0.85  
  p <- p + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5, linewidth = 0.5) +
    scale_fill_gradient2(
      low = "navy",        
      mid = "white",       
      high = "firebrick",  
      midpoint = 0,        
      name = "NES",
      limits = c(-3, 3)
    ) +
    labs(title = title_text, x = "Log2 Fold Change", y = NULL) + 
    theme_kd() +
    theme(axis.text.y = element_text(size = 7, face = "bold", color = "black", lineheight = 0.8)) + # for pathway names, use larger size
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  return(p)
}
p_ridge_carbo <- function_gsea_ridge(gsea_carbo_clean, "Hallmark GSEA (Carboplatin Resistance)")
p_ridge_taxol <- function_gsea_ridge(gsea_taxol_clean, "Hallmark GSEA (Taxol Resistance)")
library(patchwork)
combined_ridge <- p_ridge_carbo + p_ridge_taxol +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
cairo_pdf("Figures_raw/02_01_GSEA_Carbo_Taxol_hallmark.pdf", 
          width = 240 / 25.4, 
          height = 100 / 25.4) 
print(combined_ridge)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("02.RData")