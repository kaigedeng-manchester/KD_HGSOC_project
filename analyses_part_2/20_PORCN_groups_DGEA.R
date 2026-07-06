# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("20.RData") 

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

# Get PORCN+HRD groups---------------------------------------------------------

clinical_OCM_HGS_groups <- readRDS("Data/clinical_OCM_HGS_groups.rds")
clinical_TCGA_MA_HGS_groups <- readRDS("Data/clinical_TCGA_MA_HGS_groups.rds") # Not the best for DGEA?
clinical_TCGA_RNAseq_HGS <- readRDS("Data/clinical_TCGA_RNAseq_HGS.rds")
clinical_TCGA_RNAseq_HGS_groups <- clinical_TCGA_RNAseq_HGS %>%
  filter(!is.na(HRD_Status)) %>%
  filter(PORCN_Group %in% c("Low", "High")) %>%
  mutate(
    Subgroup = paste0(ifelse(HRD_Status == "HRD-", "HRP", "HRD"), "_", paste0("PORCN_", PORCN_Group)),
    Subgroup = factor(Subgroup, levels = c(
      "HRP_PORCN_Low",  
      "HRP_PORCN_High", 
      "HRD_PORCN_Low", 
      "HRD_PORCN_High"
    ))
  )

# Get original reads-----------------------------------------------------------

reads_OCM_HGS <- readRDS("Data/reads_OCM_HGS.rds") # 19931 protein coding genes
reads_TCGA_RNAseq_HGS <- readRDS("Data/reads_TCGA_RNAseq_HGS.rds") # 19938 genes

# DGEA-------------------------------------------------------------------------

library(scales)
library(limma)
library(edgeR)
library(ggrepel)
function_dgea_subgroups <- function(clin_data, expr_data, cohort_name, p_cutoff = 0.05, lfc_cutoff = 1) {
  common_samples <- intersect(clin_data$Sample.ID, colnames(expr_data))
  clin_filtered <- clin_data[match(common_samples, clin_data$Sample.ID), ]
  expr_filtered <- expr_data[, common_samples]
  rownames(clin_filtered) <- clin_filtered$Sample.ID
  cat(paste("✅", cohort_name, "Aligned! Total samples used:", length(common_samples), "\n"))
  print(summary(clin_filtered$Subgroup))
  # Design matrix and DGEList
  design <- model.matrix(~ 0 + Subgroup, data = clin_filtered)
  colnames(design) <- gsub("Subgroup", "", colnames(design))
  dge <- DGEList(counts = expr_filtered)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  # Voom and fit
  v <- voom(dge, design = design, plot = FALSE)
  fit <- lmFit(v, design)
  # Contrast (Updated to include HRD PORCN High vs HRP PORCN High)
  contr <- makeContrasts(
    Overall_PORCN_Effect = (HRD_PORCN_High + HRP_PORCN_High)/2 - (HRD_PORCN_Low + HRP_PORCN_Low)/2,
    HRD_PORCN_Effect     = HRD_PORCN_High - HRD_PORCN_Low,
    HRP_PORCN_Effect     = HRP_PORCN_High - HRP_PORCN_Low,
    Interaction          = (HRD_PORCN_High - HRD_PORCN_Low) - (HRP_PORCN_High - HRP_PORCN_Low),
    levels = design
  )
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)
  # Output helper function
  get_res_and_plot <- function(coef_name, plot_title) {
    de_res <- topTable(fit2, coef = coef_name, number = Inf, adjust = "fdr")
    de_res <- tibble::rownames_to_column(de_res, var = "GeneSymbol")
    de_res$colour <- case_when(
      de_res$P.Value < p_cutoff & de_res$logFC > lfc_cutoff ~ "up",
      de_res$P.Value < p_cutoff & de_res$logFC < -lfc_cutoff ~ "down",
      TRUE ~ "not_significant"
    )
    n_up <- sum(de_res$colour == "up")
    n_down <- sum(de_res$colour == "down")
    de_res$logP <- -log10(de_res$P.Value)
    top_genes <- bind_rows(
      de_res %>% filter(colour == "up") %>% arrange(P.Value) %>% slice_head(n = 10),
      de_res %>% filter(colour == "down") %>% arrange(P.Value) %>% slice_head(n = 10)
    )
    p_volcano <- ggplot(de_res, aes(x = logFC, y = logP, color = colour)) +
      geom_point(alpha = 0.6, size = 1) + 
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey50") + 
      geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey50") + 
      scale_color_manual(values = c("up" = "firebrick", "down" = "navy", "not_significant" = "grey80")) +
      geom_text_repel(data = top_genes, aes(label = GeneSymbol), 
                      fontface = "bold", color = "black", size = 2.5,
                      max.overlaps = 20, box.padding = 0.5, point.padding = 0.5, min.segment.length = 0) +
      annotate("text", x = Inf, y = -Inf, label = paste("UP:", n_up), 
               hjust = 1.2, vjust = -1, color = "firebrick", fontface = "bold", size = 3) +
      annotate("text", x = -Inf, y = -Inf, label = paste("DOWN:", n_down), 
               hjust = -0.2, vjust = -1, color = "navy", fontface = "bold", size = 3) +
      scale_x_continuous(breaks = seq(-6, 6, by = 3), limits = c(-6, 6), oob = scales::squish) +
      scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10), oob = scales::squish) +
      labs(title = plot_title,
           x = expression("Log"[2]*" Fold Change"),
           y = expression("-log"[10]*"(P-value)"),
           color = "Significance") +
      theme_kd() +
      theme(aspect.ratio = 1)
    return(list(data = de_res, plot = p_volcano))
  }
  # Generate results
  res_overall     <- get_res_and_plot("Overall_PORCN_Effect", paste0(cohort_name, ": Overall PORCN High vs. Low"))
  res_hrd         <- get_res_and_plot("HRD_PORCN_Effect", paste0(cohort_name, ": HRD PORCN High vs. Low"))
  res_hrp         <- get_res_and_plot("HRP_PORCN_Effect", paste0(cohort_name, ": HRP PORCN High vs. Low"))
  res_int         <- get_res_and_plot("Interaction", paste0(cohort_name, ": Interaction Effect (HRD Specific)"))
  # Return final list
  return(list(
    Overall_Effect  = res_overall,
    HRD_Effect      = res_hrd,
    HRP_Effect      = res_hrp,
    Interaction     = res_int
  ))
}
res_OCM <- function_dgea_subgroups(
  clin_data = clinical_OCM_HGS_groups, 
  expr_data = reads_OCM_HGS, 
  cohort_name = "OCM"
)
res_TCGA <- function_dgea_subgroups(
  clin_data = clinical_TCGA_RNAseq_HGS_groups, 
  expr_data = reads_TCGA_RNAseq_HGS, 
  cohort_name = "TCGA"
) # just for future analyses
library(patchwork)
p_ocm_grid <- wrap_plots(
  res_OCM$Overall_Effect$plot,
  res_OCM$HRD_Effect$plot,
  res_OCM$HRP_Effect$plot,
  res_OCM$Interaction$plot,
  ncol = 2
) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", aspect.ratio = 1)
cairo_pdf("Figures_raw/20_01_DGEA_OCM_4Panels.pdf", 
          width = 170 / 25.4, 
          height = 285 / 25.4)
print(p_ocm_grid)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("20.RData")
saveRDS(clinical_TCGA_RNAseq_HGS_groups, file = "Data/clinical_TCGA_RNAseq_HGS_groups.rds")
saveRDS(res_OCM, file = "Data/DGEA_res_OCM_group.rds")
saveRDS(res_TCGA, file = "Data/DGEA_res_TCGA_group.rds") # for future analyses
