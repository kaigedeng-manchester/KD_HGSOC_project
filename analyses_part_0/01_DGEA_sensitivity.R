# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("01.RData") 

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

# My color
colors_kd <- c(
  "up" = "firebrick", 
  "down" = "navy", 
  "not_significant" = "gray80"   
)

# Load font Arial
windowsFonts(Arial = windowsFont("Arial"))

# Preparing RNAseq data--------------------------------------------------------

# Import
reads <- read.csv(file="Data/OCM_reads.csv", row.names=1, check.names = FALSE)
anno <- reads[, 1:8]
reads <- reads[, -c(1, 3:8)]

# Duplicate gene names
reads <- reads %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), sum)) %>%  # sum counts per gene
  column_to_rownames("GeneSymbol") # raw counts of genes, no duplicates (of different versions of genes) # 59366 genes

# Keep HGSOC samples
HGS <- readRDS("Data/Clinical_OCM_HGS.rds")
HGS_screening <- readRDS("Data/Clinical_merged_OCM_HGS.rds")
all(rownames(HGS) %in% colnames(reads)) # QC
all(rownames(HGS_screening) %in% colnames(reads)) # QC
reads_HGS <- reads[, rownames(HGS)] # n = 101

# Keep protein coding genes for future analysis
protein_coding_genes <- anno %>% filter(Class == "protein_coding") %>% pull(GeneSymbol) # 19965 genes
reads_HGS <- reads_HGS[rownames(reads_HGS) %in% protein_coding_genes, ] # 19931 genes kept

# Intersect with screening data
reads_HGS_screening <- reads[, rownames(HGS_screening)] # n = 67

# DGEA-------------------------------------------------------------------------

library(limma)
library(edgeR)
library(ggrepel)

# Limma-voom pipeline
function_dgea_tertile <- function(clin_data = HGS_screening, expr_data = reads_HGS_screening, auc_var, drug_name) {
  cutoffs <- quantile(clin_data[[auc_var]], probs = c(1/3, 2/3), na.rm = TRUE)
  clin_filtered <- clin_data %>%
    mutate(
      Group = case_when(
        .data[[auc_var]] <= cutoffs[1] ~ "Sensitive",
        .data[[auc_var]] >= cutoffs[2] ~ "Resistant",
        TRUE ~ "Intermediate"
      )
    ) %>%
    filter(Group != "Intermediate") %>%
    mutate(Group = factor(Group, levels = c("Sensitive", "Resistant")))
  cat(paste("  - Sensitive N =", sum(clin_filtered$Group == "Sensitive"), "\n"))
  cat("    [Names]:", paste(rownames(clin_filtered)[clin_filtered$Group == "Sensitive"], collapse = ", "), "\n\n")
  cat(paste("  - Resistant N =", sum(clin_filtered$Group == "Resistant"), "\n"))
  cat("    [Names]:", paste(rownames(clin_filtered)[clin_filtered$Group == "Resistant"], collapse = ", "), "\n\n")
  expr_filtered <- expr_data[, rownames(clin_filtered)]
  design <- model.matrix(~ 0 + Group, data = clin_filtered)
  colnames(design) <- c("Sensitive", "Resistant")
  dge <- DGEList(counts = expr_filtered)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design = design, plot = FALSE)
  fit <- lmFit(v, design)
  contr <- makeContrasts(Resistant_vs_Sensitive = Resistant - Sensitive, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)
  de_res <- topTable(fit2, coef = "Resistant_vs_Sensitive", number = Inf, adjust = "fdr")
  de_res <- tibble::rownames_to_column(de_res, var = "GeneSymbol")
  de_res$colour <- case_when(
    de_res$P.Value < 0.05 & de_res$logFC > 1 ~ "up", # if use adj.P, no gene significant
    de_res$P.Value < 0.05 & de_res$logFC < -1 ~ "down", # so this satge using original P instead
    TRUE ~ "not_significant"
  )
  cat(paste("  - Up-regulated N =", sum(de_res$colour == "up"), "\n"))
  cat(paste("  - Down-regulated N =", sum(de_res$colour == "down"), "\n"))
  top20_genes <- bind_rows(
    de_res %>% filter(colour == "up") %>% arrange(P.Value) %>% slice_head(n = 10),
    de_res %>% filter(colour == "down") %>% arrange(P.Value) %>% slice_head(n = 10)
  )
  p_volcano <- ggplot(de_res, aes(x = logFC, y = -log10(P.Value), color = colour)) +
    geom_point(alpha = 0.8, size = 1) + 
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") + 
    scale_color_manual(values = colors_kd) +
    geom_text_repel(data = top20_genes, aes(label = GeneSymbol), 
                    fontface = "bold", color = "black", size = 2.5, # 7.11pt (significant, easy to read)
                    max.overlaps = 20, box.padding = 0.5, point.padding = 0.5, min.segment.length = 0) +
    labs(title = paste("DEGs:", drug_name, "Resistant vs Sensitive"),
         x = expression("Log"[2]*" Fold Change"),
         y = expression("-log"[10]*"(P-value)"),
         color = "Significance") +
    theme_kd() +
    theme(aspect.ratio = 1)
  print(p_volcano)
  return(list(data = de_res, plot = p_volcano))
}
res_carbo <- function_dgea_tertile(auc_var = "carboplatin_AUC", drug_name = "Carboplatin")
res_taxol <- function_dgea_tertile(auc_var = "taxol_AUC", drug_name = "Taxol")
degs_carbo <- res_carbo$data
degs_taxol <- res_taxol$data

# Volcano plot
library(patchwork)
combined_volcano <- res_carbo$plot + res_taxol$plot + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", aspect.ratio = 1)
cairo_pdf("Figures_raw/01_01_DGEA_Carbo_Taxol.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4) 
print(combined_volcano)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("01.RData")
saveRDS(degs_carbo, file = "Data/degs_carbo.rds")
saveRDS(degs_taxol, file = "Data/degs_taxol.rds")
saveRDS(reads_HGS, file = "Data/reads_OCM_HGS.rds")