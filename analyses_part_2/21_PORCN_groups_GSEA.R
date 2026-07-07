# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("21.RData") 

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

# Get DEG results--------------------------------------------------------------

res_OCM <- readRDS("Data/DGEA_res_OCM_group.rds")

# Check overlap----------------------------------------------------------------

library(ggVennDiagram)
function_plot_venn <- function(res_obj, cohort_name, direction = c("up", "down")) {
  direction <- match.arg(direction)
  genes_HRD <- res_obj$HRD_Effect$data %>% filter(colour == direction) %>% pull(GeneSymbol)
  genes_HRP <- res_obj$HRP_Effect$data %>% filter(colour == direction) %>% pull(GeneSymbol)
  genes_Int <- res_obj$Interaction$data %>% filter(colour == direction) %>% pull(GeneSymbol)
  venn_list <- list(
    "HRD" = genes_HRD,
    "HRP" = genes_HRP,
    "HRD Specific"  = genes_Int
  )
  if (direction == "up") {
    low_col <- "white"
    high_col <- "firebrick"
    title_text <- "Up-regulated Genes"
  } else {
    low_col <- "white"
    high_col <- "navy"
    title_text <- "Down-regulated Genes"
  }
  p_venn <- ggVennDiagram(venn_list, label_alpha = 0, edge_lty = "solid", edge_size = 0.5, set_size = 4, label_size = 3.5) + # special big font for Venn
    scale_fill_gradient(low = low_col, high = high_col) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
    labs(title = title_text)
  print(p_venn)
  return(p_venn)
}
library(patchwork)
p_ocm_up <- function_plot_venn(res_obj = res_OCM, cohort_name = "OCM", direction = "up")
p_ocm_down <- function_plot_venn(res_obj = res_OCM, cohort_name = "OCM", direction = "down")
p_ocm_combined <- p_ocm_up | p_ocm_down
cairo_pdf("Figures_raw/21_01_Venn_Combined_OCM.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4)
print(p_ocm_combined)
dev.off()

# GSEA scan--------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(forcats)
library(stringr)
function_run_gsea_hallmark <- function(res_obj, cohort_name, out_dir = "Results") {
  prep_ranked_list <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df_clean <- df %>% 
      filter(!is.na(GeneSymbol) & !is.na(t)) %>% 
      arrange(desc(t))
    ranks <- df_clean$t
    names(ranks) <- df_clean$GeneSymbol
    return(ranks)
  }
  cat("Preparing ranked DEG list...\n")
  ranked_lists <- list(
    "Overall"      = prep_ranked_list(res_obj$Overall_Effect$data),
    "HRD"          = prep_ranked_list(res_obj$HRD_Effect$data),
    "HRP"          = prep_ranked_list(res_obj$HRP_Effect$data),
    "HRD Specific" = prep_ranked_list(res_obj$Interaction$data)
  )
  ranked_lists <- ranked_lists[!sapply(ranked_lists, is.null)]
  cat("Loading Hallmark dictionary...\n")
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol) %>%
    mutate(gs_name = gsub("HALLMARK_", "", gs_name))
  cat("Running multi-group CompareCluster GSEA (no cut-off applied)...\n")
  comp_res <- compareCluster(
    geneCluster = ranked_lists, 
    fun = "GSEA", 
    TERM2GENE = m_t2g, 
    pvalueCutoff = 1.0,
    pAdjustMethod = "BH", 
    eps = 1e-10
  )
  if (is.null(comp_res) || nrow(as.data.frame(comp_res)) == 0) {
    message("No enrichment!")
    return(NULL)
  }
  df_combined <- as.data.frame(comp_res)
  if(!dir.exists(out_dir)) dir.create(out_dir)
  write.csv(df_combined, sprintf("%s/21_02_GSEA_Hallmark_Results_%s.csv", out_dir, cohort_name), row.names = FALSE)
  saveRDS(comp_res, sprintf("%s/21_02_GSEA_Hallmark_Object_%s.rds", out_dir, cohort_name))
  return(df_combined)
}
df_gsea_raw <- function_run_gsea_hallmark(res_OCM, "OCM")

# Scan and plot the results
plot_hallmark_global_landscape <- function(gsea_df, cohort_name) {
  df_padded <- gsea_df %>%
    mutate(Is_Sig = (p.adjust < 0.10 & abs(NES) > 1.0)) %>%
    complete(
      Description = unique(gsea_df$Description), 
      fill = list(NES = 0, p.adjust = 1, Is_Sig = FALSE)
    )
  pathway_order <- df_padded %>%
    group_by(Description) %>%
    summarise(
      overall_NES = sum(NES[Cluster == "Overall"], na.rm = TRUE),
      max_abs_NES = max(abs(NES), na.rm = TRUE)
    ) %>%
    arrange(overall_NES, max_abs_NES) %>%
    pull(Description)
  df_plot <- df_padded %>%
    mutate(
      Description = factor(Description, levels = pathway_order),
      logP = -log10(p.adjust)
    )
  p_landscape <- ggplot(df_plot, aes(x = Cluster, y = Description)) +
    geom_hline(yintercept = df_plot$Description, color = "grey96", linewidth = 0.3) +
    geom_point(
      aes(fill = NES, size = logP, alpha = Is_Sig, stroke = ifelse(Is_Sig, 0.6, 0.15)), 
      shape = 21, color = "black"
    ) +
    scale_fill_gradient2(
      low = "navy", mid = "white", high = "firebrick", 
      midpoint = 0, name = "NES",
      limits = c(-3, 3), oob = scales::squish
    ) +
    scale_size_continuous(name = expression("-log"[10]*"(FDR)") , range = c(1.5, 5)) +
    scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.18), guide = "none") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(gsub("_", " ", x), width = 40)) +
    theme_kd() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold", size = 7, color = "black"),
      axis.text.y = element_text(size = 5.5, face = "bold", color = "black", lineheight = 0.8),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.size = unit(0.3, "cm")
    ) +
    labs(
      title = paste0("PORCN-Related Hallmark GSEA Landscape (", cohort_name, ")"),
      x = "", y = ""
    )
  return(p_landscape)
}
p_hallmark_final <- plot_hallmark_global_landscape(df_gsea_raw, "OCM")
cairo_pdf("Figures_raw/21_02_GSEA_Hallmark_Landscape_OCM.pdf", 
          width = 170 / 25.4,
          height = 240 / 25.4)
print(p_hallmark_final)
dev.off()

# TF activity------------------------------------------------------------------

library(decoupleR)
library(OmnipathR)
TF_net <- get_collectri(organism = 'human', split_complexes = FALSE)
function_extract_tf_input <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df_clean <- df %>%
    filter(!is.na(GeneSymbol) & !is.na(t)) %>%
    group_by(GeneSymbol) %>%
    summarise(score = mean(t, na.rm = TRUE)) %>%
    distinct(GeneSymbol, .keep_all = TRUE)
  tf_mat <- matrix(df_clean$score, 
                   ncol = 1, 
                   dimnames = list(df_clean$GeneSymbol, "t_value"))
  return(tf_mat)
}
mat_lists <- list(
  "Overall"      = function_extract_tf_input(res_OCM$Overall_Effect$data),
  "HRD"          = function_extract_tf_input(res_OCM$HRD_Effect$data),
  "HRP"          = function_extract_tf_input(res_OCM$HRP_Effect$data),
  "HRD Specific" = function_extract_tf_input(res_OCM$Interaction$data)
)
mat_lists <- mat_lists[!sapply(mat_lists, is.null)]
res_tf_list <- lapply(names(mat_lists), function(grp) {
  res <- run_ulm(
    mat = mat_lists[[grp]],
    network = TF_net,
    .source = "source", 
    .target = "target",
    .mor = "mor",
    minsize = 5
  ) %>%
    mutate(Cluster = grp) 
  return(res)
})
res_tf_final <- bind_rows(res_tf_list) %>%
  mutate(Cluster = factor(Cluster, levels = c("Overall", "HRD", "HRP", "HRD Specific"))) %>%
  group_by(Cluster) %>% 
  mutate(FDR = p.adjust(p_value, method = "BH")) %>% 
  ungroup() %>%
  mutate(
    neg_log10_FDR = -log10(FDR + 1e-10), 
    Significance = case_when(
      FDR < 0.10 & p_value < 0.05 & score > 2 ~ "Activated",
      FDR < 0.10 & p_value < 0.05 & score < -2 ~ "Repressed",
      TRUE ~ "Not Sig"
    )
  )
res_tf_top <- res_tf_final %>%
  filter(Significance != "Not Sig") %>%
  group_by(Cluster, Significance) %>%
  slice_max(order_by = abs(score), n = 10, with_ties = FALSE) %>% 
  ungroup()
library(ggrepel)
library(patchwork)
plot_single_tf_volcano <- function(df_final, df_top, target_cluster) {
  df_sub <- df_final %>% filter(Cluster == target_cluster)
  top_sub <- df_top %>% filter(Cluster == target_cluster)
  p <- ggplot(df_sub, aes(x = score, y = neg_log10_FDR)) +
    geom_point(data = filter(df_sub, Significance == "Not Sig"), 
               color = "grey80", size = 1.0, alpha = 0.5) +
    geom_point(data = filter(df_sub, Significance != "Not Sig"), 
               aes(fill = Significance), shape = 21, color = "black", stroke = 0.3, size = 1.5, alpha = 0.85) +
    geom_text_repel(data = top_sub, aes(label = source, color = Significance), 
                    show.legend = FALSE, fontface = "bold", size = 2.5,
                    max.overlaps = 30, box.padding = 0.6, point.padding = 0.4, 
                    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3) +
    scale_fill_manual(values = c("Activated" = "firebrick", "Repressed" = "navy"), 
                      drop = FALSE, name = "TF Status") +
    scale_color_manual(values = c("Activated" = "firebrick", "Repressed" = "navy"), 
                       drop = FALSE) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(0.10), linetype = "dashed", color = "grey30") +
    labs(
      title = target_cluster, 
      x = "TF Activity Score (ULM)", 
      y = expression("-log"[10]*"(FDR)")
    ) +
    theme_kd() +
    theme(
      aspect.ratio = 1
    )
  return(p)
}
p_tf_overall      <- plot_single_tf_volcano(res_tf_final, res_tf_top, "Overall")
p_tf_hrd          <- plot_single_tf_volcano(res_tf_final, res_tf_top, "HRD")
p_tf_hrp          <- plot_single_tf_volcano(res_tf_final, res_tf_top, "HRP")
p_tf_hrd_specific <- plot_single_tf_volcano(res_tf_final, res_tf_top, "HRD Specific")
combined_tf_volcano <- (p_tf_overall | p_tf_hrd) / (p_tf_hrp | p_tf_hrd_specific) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.box.margin = margin(t = -10))
cairo_pdf("Figures_raw/21_03_TF_Landscape_OCM.pdf", 
          width = 170 / 25.4, 
          height = 190 / 25.4)
print(combined_tf_volcano)
dev.off()

# Targeted WNT GSEA-----------------------------------------------------------

m_t2g_go <- msigdbr(species = "Homo sapiens", category = "C5", subcollection = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)
target_wnt_pathways <- c(
  "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GOBP_WNT_SIGNALING_PATHWAY_PLANAR_CELL_POLARITY_PATHWAY",
  "GOBP_WNT_SIGNALING_PATHWAY_CALCIUM_MODULATING_PATHWAY"
)
m_t2g_wnt <- m_t2g_go %>% filter(gs_name %in% target_wnt_pathways)
prep_ranked_list_global <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df_clean <- df %>% 
    filter(!is.na(GeneSymbol) & !is.na(t)) %>% 
    arrange(desc(t))
  ranks <- df_clean$t
  names(ranks) <- df_clean$GeneSymbol
  return(ranks)
}
ranked_lists <- list(
  "Overall"      = prep_ranked_list_global(res_OCM$Overall_Effect$data),
  "HRD"          = prep_ranked_list_global(res_OCM$HRD_Effect$data),
  "HRP"          = prep_ranked_list_global(res_OCM$HRP_Effect$data),
  "HRD Specific" = prep_ranked_list_global(res_OCM$Interaction$data)
)
ranked_lists <- ranked_lists[!sapply(ranked_lists, is.null)]
comp_wnt <- compareCluster(
  geneCluster = ranked_lists, 
  fun = "GSEA", 
  TERM2GENE = m_t2g_wnt, 
  minGSSize     = 3,      # for small genesets-calcium
  maxGSSize     = 1000,   # for big genesets
  pvalueCutoff = 1.0,
  pAdjustMethod = "BH", 
  eps = 1e-10
)
df_wnt_plot <- as.data.frame(comp_wnt) %>%
  mutate(
    Description = gsub("GOBP_", "", ID),
    Description = gsub("_", " ", Description),
    Description = factor(Description, levels = c(
      "WNT SIGNALING PATHWAY CALCIUM MODULATING PATHWAY",
      "WNT SIGNALING PATHWAY PLANAR CELL POLARITY PATHWAY",
      "NON CANONICAL WNT SIGNALING PATHWAY",
      "CANONICAL WNT SIGNALING PATHWAY"
    )),
    Is_Sig = (pvalue < 0.05 & abs(NES) > 1.0),
    logP = -log10(pvalue + 1e-10) 
  )
p_wnt_bubble <- ggplot(df_wnt_plot, aes(x = Cluster, y = Description)) +
  geom_hline(yintercept = df_wnt_plot$Description, color = "grey96", linewidth = 0.3) +
  geom_point(
    aes(fill = NES, size = logP, alpha = Is_Sig, stroke = ifelse(Is_Sig, 0.6, 0.15)), 
    shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "firebrick", 
    midpoint = 0, name = "NES",
    limits = c(-3, 3), oob = scales::squish
  ) +
  scale_size_continuous(name = expression("-log"[10]*"(P-value)"), range = c(2, 6)) +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.18), guide = "none") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE, labels = function(x) stringr::str_wrap(x, width = 35)) +
  theme_kd() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7, color = "black"),
    axis.text.y = element_text(size = 6.5, face = "bold", color = "black", lineheight = 0.8)
  ) +
  labs(
    title = "Targeted GO-BP GSEA: WNT Sub-pathways",
    x = "", y = ""
  )
cairo_pdf("Figures_raw/21_04_WNT_Subpathways_Bubble.pdf", 
          width = 170 / 25.4,   
          height = 80 / 25.4)   
print(p_wnt_bubble)
dev.off()

# HRD PORCN-high core signature------------------------------------------------

library(pheatmap)
library(stringr)
library(pROC)

# Expression matrix OCM
expr_OCM_HGS <- readRDS("Data/expr_OCM_HGS.rds")
valid_expr_genes <- rownames(expr_OCM_HGS)

# Get core pathways GSEA
target_cluster <- "HRD"
top_gsea_paths <- df_gsea_raw %>%
  filter(Cluster == target_cluster, p.adjust < 0.10, abs(NES) > 1) %>% #same as screening
  arrange(desc(abs(NES))) %>%
  ungroup()

# Get leading edges GSEA
core_genes_raw <- paste(top_gsea_paths$core_enrichment, collapse = "/")
core_genes_list <- unique(unlist(strsplit(core_genes_raw, "/")))

# Filter with t value DGEA
hrd_deg_df <- res_OCM$HRD_Effect$data 
core_genes_candidate <- hrd_deg_df %>%
  filter(GeneSymbol %in% core_genes_list) %>%  
  filter(GeneSymbol %in% valid_expr_genes) %>%
  filter(P.Value < 0.05) %>% #looser than DGEA                
  arrange(desc(abs(t)))                        

# Group label
clinical_OCM_HGS_groups <- readRDS("Data/clinical_OCM_HGS_groups.rds")
hrd_clinical <- clinical_OCM_HGS_groups %>% filter(grepl("HRD", Subgroup))
hrd_ids <- hrd_clinical$Sample.ID
setdiff(hrd_ids, colnames(expr_OCM_HGS)) #check

# Expression matrix-candidate
expr_candidate_mat <- expr_OCM_HGS[core_genes_candidate$GeneSymbol, hrd_ids]

# Labels
true_labels <- ifelse(hrd_clinical$Subgroup == "HRD_PORCN_High", 1, 0)
names(true_labels) <- hrd_clinical$Sample.ID

# Set size
max_test_genes <- nrow(core_genes_candidate)

# Test different set sizes 
library(decoupleR)
mor_list <- lapply(10:max_test_genes, function(n) { #min 10 genes
  sub_genes <- core_genes_candidate[1:n, ]
  data.frame(
    source = paste0("Sig_", sprintf("%02d", n)), 
    target = sub_genes$GeneSymbol,
    weight = sub_genes$t                         
  )
})
mor_omnibus <- bind_rows(mor_list)
ulm_internal <- run_ulm(
  mat = expr_candidate_mat, 
  network = mor_omnibus,
  .source = "source",
  .target = "target",
  .mor = "weight",
  minsize = 3
)
cv_metrics_df <- ulm_internal %>%
  mutate(Sample.ID = condition) %>%
  inner_join(data.frame(Sample.ID = names(true_labels), True_Label = true_labels), by = "Sample.ID") %>%
  group_by(source) %>%
  summarise(
    Gene_Count = as.numeric(gsub("Sig_", "", source[1])),
    AUROC = as.numeric(pROC::roc(True_Label, score, quiet = TRUE)$auc),
    .groups = "drop"
  ) %>%
  arrange(Gene_Count)

# Plot original candidate heatmap
final_golden_genes <- core_genes_candidate$GeneSymbol
expr_hrd_golden <- expr_OCM_HGS[rownames(expr_OCM_HGS) %in% final_golden_genes, hrd_ids]
annotation_col <- data.frame(
  Subgroup = hrd_clinical$Subgroup,
  row.names = hrd_clinical$Sample.ID
)
ann_colors <- list(
  Subgroup = c("HRD_PORCN_High" = "#B22222", "HRD_PORCN_Low" = "#4682B4")
)
expr_hrd_z <- t(scale(t(expr_hrd_golden)))
expr_hrd_z[expr_hrd_z > 2] <- 2
expr_hrd_z[expr_hrd_z < -2] <- -2
pheatmap_obj <- pheatmap(expr_hrd_z,
                         scale = "none",             
                         cluster_cols = TRUE,        
                         cluster_rows = TRUE, 
                         show_colnames = TRUE,
                         fontsize_col = 4,
                         angle_col = 45,
                         show_rownames = TRUE,       
                         fontsize_row = 3,
                         fontsize = 5,
                         treeheight_col = 25,        
                         treeheight_row = 25,
                         cutree_cols = 2,            
                         cutree_rows = 2,
                         border_color = NA,
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(c("navy", "white", "firebrick"))(100))
cairo_pdf("Figures_raw/21_05_Heatmap_HRD.pdf", width = 120 / 25.4, height = 170 / 25.4)
print(pheatmap_obj)
dev.off()

# Plot internal validation
p_cv_curve <- ggplot(cv_metrics_df, aes(x = Gene_Count, y = AUROC)) +
  geom_line(color = "grey40", linewidth = 0.3) +
  geom_point(shape = 21, size = 0.5, color = "black", stroke = 0.15, fill = "firebrick") +
  labs(title = "Signature Size Performance in Internal Validation",
       x = "Number of Top Genes in Signature", y = "AUROC") +
  theme_kd() +
  theme(legend.position = "none", aspect.ratio = 1/1.5)
cairo_pdf("Figures_raw/21_06_Core_Signature_Internal.pdf", width = 85 / 25.4, height = 70 / 25.4)
print(p_cv_curve)
dev.off()

# Test core feature------------------------------------------------------------

batch_ulm_scoring <- function(expr_mat, mor_net, cohort_name) {
  score_res <- run_ulm(
    mat = expr_mat,
    network = mor_net,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = 3 
  ) %>%
    mutate(
      Gene_Count = as.numeric(gsub("Sig_", "", source)),
      Cohort = cohort_name
    ) %>%
    rename(Sample.ID = condition, Signature = source, ULM_Score = score) %>%
    select(Cohort, Sample.ID, Signature, Gene_Count, ULM_Score)
  return(score_res)
}
scores_OCM <- batch_ulm_scoring(expr_OCM_HGS, mor_omnibus, "OCM")
expr_CCLE_HGS <- readRDS("Data/expr_CCLE_HGS.rds")
scores_CCLE <- batch_ulm_scoring(expr_CCLE_HGS, mor_omnibus, "CCLE")
expr_TCGA_RNAseq_HGS <- readRDS("Data/expr_TCGA_RNAseq_HGS.rds")
scores_TCGA <- batch_ulm_scoring(expr_TCGA_RNAseq_HGS, mor_omnibus, "TCGA")

# Clinical labels
clinical_OCM <- readRDS("Data/clinical_OCM_HGS.rds")
clinical_TCGA <- readRDS("Data/clinical_TCGA_RNAseq_HGS.rds")
groups_ocm <- readRDS("Data/groups_ocm.rds")
groups_tcga <- readRDS("Data/groups_tcga_rna.rds")
groups_ccle <- readRDS("Data/groups_ccle.rds")
function_annotate_scores <- function(scores_df, clinical_df, groups_list, has_hrd = TRUE) {
  df <- scores_df %>%
    mutate(
      PORCN_Tertile = case_when(
        Sample.ID %in% groups_list$High ~ "High",
        Sample.ID %in% groups_list$Low ~ "Low",
        TRUE ~ "Medium"
      ),
      PORCN_Tertile = factor(PORCN_Tertile, levels = c("Low", "Medium", "High"))
    )
  if (has_hrd && !is.null(clinical_df)) {
    df <- df %>%
      left_join(clinical_df %>% select(Sample.ID, HRD_Status), by = "Sample.ID") %>%
      mutate(
        HRD_Status = as.character(HRD_Status),
        HRD_Status = case_when(
          HRD_Status == "HRD+" ~ "HRD",
          HRD_Status == "HRD-" ~ "HRP",
          is.na(HRD_Status) ~ "Unknown",
          TRUE ~ HRD_Status
        ),
        HRD_Status = factor(HRD_Status, levels = c("HRP", "HRD", "Unknown"))
      )
  } else {
    df <- df %>% mutate(HRD_Status = factor("Unknown", levels = c("HRP", "HRD", "Unknown")))
  }
  return(df)
}
anno_OCM  <- function_annotate_scores(scores_OCM, clinical_OCM, groups_ocm, TRUE)
anno_TCGA <- function_annotate_scores(scores_TCGA, clinical_TCGA, groups_tcga, TRUE)
anno_CCLE <- function_annotate_scores(scores_CCLE, NULL, groups_ccle, FALSE)

# External validation
calculate_validation_auroc <- function(df, cohort_label) {
  df_clean <- df %>%
    filter(PORCN_Tertile %in% c("High", "Low")) %>%
    mutate(True_Label = ifelse(PORCN_Tertile == "High", 1, 0))
  if(nrow(df_clean) == 0 || length(unique(df_clean$True_Label)) < 2) {
    return(data.frame(Gene_Count = 5:max_test_genes, AUROC = NA, Cohort = cohort_label))
  }
  df_clean %>%
    group_by(Gene_Count) %>%
    summarise(
      AUROC = as.numeric(pROC::roc(True_Label, ULM_Score, levels = c(0, 1), direction = "<", quiet = TRUE)$auc),
      .groups = "drop"
    ) %>%
    mutate(Cohort = cohort_label)
}
df_val_ocm <- anno_OCM %>% filter(HRD_Status == "Unknown")
res_auroc_ocm <- calculate_validation_auroc(df_val_ocm, "OCM (Subset HR Unknown)")
df_val_tcga <- anno_TCGA %>% filter(HRD_Status == "HRD")
res_auroc_tcga <- calculate_validation_auroc(df_val_tcga, "TCGA (HRD)")
res_auroc_ccle <- calculate_validation_auroc(anno_CCLE, "CCLE (All HR Unknown)")
external_auroc_df <- bind_rows(res_auroc_ocm, res_auroc_tcga, res_auroc_ccle)
mean_auroc_df <- external_auroc_df %>%
  pivot_wider(names_from = Cohort, values_from = AUROC) %>%
  rowwise() %>%
  mutate(Mean_AUROC = mean(c_across(starts_with(c("OCM", "TCGA", "CCLE"))), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Gene_Count)

# Plot external validation
cohort_colors <- c("OCM (Subset HR Unknown)" = "#E69F00", "TCGA (HRD)" = "firebrick", "CCLE (All HR Unknown)" = "#56B4E9")
p_external_landscape <- ggplot(external_auroc_df, aes(x = Gene_Count, y = AUROC, color = Cohort)) +
  geom_vline(xintercept = 51, color = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.3) +
  geom_point(aes(fill = Cohort), shape = 21, size = 0.5, stroke = 0.15, color = "black") +
  geom_smooth(data = mean_auroc_df, aes(x = Gene_Count, y = Mean_AUROC), 
              inherit.aes = FALSE, method = "loess", span = 0.3, 
              color = "black", linewidth = 0.3, se = FALSE, linetype = "dashed") +
  scale_color_manual(values = cohort_colors) +
  scale_fill_manual(values = cohort_colors) +
  theme_kd(legend_pos = "bottom") +
  labs(
    title = "Signature Generalization Across Independent Cohorts",
    x = "Number of Top Genes in Signature",
    y = "Validation AUROC"
  ) +
  theme(
    aspect.ratio = 1/1.5,
    legend.title = element_blank()
  )
cairo_pdf("Figures_raw/21_07_Core_Signature_External.pdf", width = 85 / 25.4, height = 70 / 25.4)
print(p_external_landscape)
dev.off()

# Comparisons
library(patchwork)
library(ggpubr)
target_gene_count <- 51
df_combined <- bind_rows(
  anno_OCM %>% mutate(Cohort = "OCM"),
  anno_TCGA %>% mutate(Cohort = "TCGA"),
  anno_CCLE %>% mutate(Cohort = "CCLE")
) %>% 
  filter(Gene_Count == target_gene_count) %>%
  mutate(
    Cohort = factor(Cohort, levels = c("OCM", "TCGA", "CCLE")),
    HRD_Status = factor(HRD_Status, levels = c("HRD", "HRP", "Unknown")),
    PORCN_Tertile = factor(PORCN_Tertile, levels = c("Low", "Medium", "High"))
  )
df_zscored <- df_combined %>%
  group_by(Cohort) %>%
  mutate(ULM_Z = as.numeric(scale(ULM_Score))) %>%
  ungroup()

# level 1
comp_lvl1 <- list(c("OCM", "TCGA"), c("TCGA", "CCLE"), c("OCM", "CCLE"))
p_level1 <- ggplot(df_combined, aes(x = Cohort, y = ULM_Score)) +
  geom_violin(aes(fill = Cohort), trim = FALSE, alpha = 0.4, color = NA, scale = "width", width = 0.3) +
  geom_jitter(aes(color = Cohort), width = 0.1, size = 1.0, alpha = 0.5) +
  geom_boxplot(width = 0.15, fill = alpha("white", 0.3), color = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("OCM" = "#E69F00", "TCGA" = "#009E73", "CCLE" = "#56B4E9")) +
  scale_color_manual(values = c("OCM" = "#E69F00", "TCGA" = "#009E73", "CCLE" = "#56B4E9")) +
  stat_compare_means(comparisons = comp_lvl1, method = "wilcox.test", label = "p.format", tip.length = 0.02, size = 2.5) +
  stat_compare_means(label.y.npc = "top", fontface = "bold", size = 2.5) +
  theme_kd(legend_pos = "none") +
  theme(axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black")) +
  labs(title = "Inter-Cohort Baseline Comparison", x = "Cohort", y = "Raw Signature ULM Score") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# level 2
df_level2 <- df_zscored %>% filter(Cohort %in% c("OCM", "TCGA"))
comp_lvl2 <- list(c("HRP", "HRD"), c("HRD", "Unknown"), c("HRP", "Unknown"))
p_level2 <- ggplot(df_level2, aes(x = HRD_Status, y = ULM_Z)) +
  geom_violin(aes(fill = HRD_Status), trim = FALSE, alpha = 0.4, color = NA, scale = "width", width = 0.3) +
  geom_jitter(aes(color = HRD_Status), width = 0.1, size = 1.0, alpha = 0.5) +
  geom_boxplot(width = 0.15, fill = alpha("white", 0.3), color = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("HRD" = "firebrick", "HRP" = "steelblue", "Unknown" = "grey60")) +
  scale_color_manual(values = c("HRD" = "firebrick", "HRP" = "steelblue", "Unknown" = "grey60")) +
  stat_compare_means(comparisons = comp_lvl2, method = "wilcox.test", label = "p.format", tip.length = 0.02, size = 2.5) +
  stat_compare_means(label.y.npc = "top", fontface = "bold", size = 2.5) +
  facet_wrap(~ Cohort, ncol = 1) + 
  theme_kd(legend_pos = "none") +
  theme(axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
        strip.text = element_text(size = 8, face = "bold", margin = margin(t = 4, b = 4))) +
  labs(title = "Intra-Cohort HRD Impact", x = "HRD Status", y = "Signature ULM Score (Z-Standardised)") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# level 3
df_level3 <- df_zscored %>%
  filter(!(Cohort == "CCLE" & HRD_Status %in% c("HRD", "HRP"))) %>%
  mutate(
    Facet_Name = factor(paste(Cohort, HRD_Status, sep = " | "),
                        levels = c(
                          "OCM | HRD", "OCM | HRP", "OCM | Unknown",
                          "TCGA | HRD", "TCGA | HRP", "TCGA | Unknown",
                          "CCLE | Unknown"
                        ))
  )
comp_lvl3 <- list(c("Low", "Medium"), c("Medium", "High"), c("Low", "High"))
color_tertiles <- c("Low" = "steelblue", "Medium" = "grey70", "High" = "firebrick")
p_level3 <- ggplot(df_level3, aes(x = PORCN_Tertile, y = ULM_Z)) +
  geom_violin(aes(fill = PORCN_Tertile), trim = FALSE, alpha = 0.4, color = NA, scale = "width", width = 0.3) +
  geom_jitter(aes(color = PORCN_Tertile), width = 0.1, size = 1.0, alpha = 0.5) +
  geom_boxplot(width = 0.15, fill = alpha("white", 0.3), color = "black", outlier.shape = NA) +
  scale_fill_manual(values = color_tertiles) +
  scale_color_manual(values = color_tertiles) +
  stat_compare_means(comparisons = comp_lvl3, method = "wilcox.test", label = "p.format", tip.length = 0.02, size = 2.5) +
  stat_compare_means(label.y.npc = "top", fontface = "bold", size = 2.5) +
  facet_wrap(~ Facet_Name, ncol = 3) + 
  theme_kd(legend_pos = "none") +
  theme(axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
        strip.text = element_text(size = 8, face = "bold", margin = margin(t = 4, b = 4)),
        panel.spacing.x = unit(4, "mm"), panel.spacing.y = unit(4, "mm")) +
  labs(title = "Intra-Cohort PORCN Correlation Across HRD Context", x = "PORCN Expression Tertile", y = "Signature ULM Score (Z-Standardised)") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Combine
layout_design <- "
  ##A
  BCC
"
layout_final <- p_level1 + p_level2 + p_level3 + 
  plot_layout(design = layout_design, heights = c(1, 2))
cairo_pdf("Figures_raw/21_08_Signature_PORCN.pdf", 
          width = 240 / 25.4, 
          height = 200 / 25.4)
print(layout_final)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("21.RData")
saveRDS(df_zscored, "Data/df_51genes_zscored.rds")