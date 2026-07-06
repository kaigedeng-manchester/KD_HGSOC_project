# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("04.RData") 

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

# Risk gene list--------------------------------------------------------------

trans_degs <- readRDS("Data/trans_degs.rds")
risk_genes <- (trans_degs %>% filter(Translational_Status == "Resistant & Risk"))$GeneSymbol

# TCGA matrix and clinical data----------------------------------------------

# TCGA_RNAseq
reads_TCGA_RNAseq <- read.csv(file="Data/matrix_TCGA_RNAseq_sum.csv", row.names=1, check.names = FALSE) # N=434
clinical_GDC <- read.csv("Data/Clinical_GDC_total.csv")
clinical_TCGA_RNAseq <- read.csv("Data/Clinical_TCGA_RNAseq.csv")
clinical_TCGA_RNAseq_HGS <- clinical_TCGA_RNAseq %>% filter(
  Data.Type == "Gene Expression Quantification" &
    Tissue.Type == "Tumor" &
    Inferred.identity %in% c(
      "HGSOC (good confidence)",
      "HGSOC (high confidence)",
      "inferred HGSOC (good confidence)",
      "inferred HGSOC (high confidence)"
    )
)
reads_TCGA_RNAseq_HGS <- reads_TCGA_RNAseq[, colnames(reads_TCGA_RNAseq) %in% clinical_TCGA_RNAseq_HGS$Sample.ID] # N=418 HGSOC samples

# Filter low expression (keep genes expressed in >=20% samples with CPM > 1)
library(limma)
library(edgeR)
dge <- DGEList(counts = reads_TCGA_RNAseq_HGS)
keep <- rowSums(cpm(dge) > 1) >= (0.2 * ncol(reads_TCGA_RNAseq))
dge <- dge[keep, , keep.lib.sizes = FALSE] # keep 14884 genes from 19938 coding genes

# Calculate Normalization Factors (TMM)
dge <- calcNormFactors(dge, method = "TMM")

# Convert to log2-CPM (Continuous values needed for correction)
cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE)
expr <- log2(cpm_mat + 1) #log2(cpm+1) transformation

# Output 
expr_TCGA_RNAseq_HGS <- expr
rm(expr, cpm_mat, dge, keep)

# TCGA_MA
reads_TCGA_MA <- read.csv(file="Data/matrix_TCGA_MA_mean.csv", row.names=1, check.names = FALSE) # N=540
clinical_TCGA_MA <- read.csv("Data/Clinical_TCGA_MA.csv")
clinical_TCGA_MA_HGS <- clinical_TCGA_MA %>% filter(
  Data.Type == "Raw Intensities" &
    Tissue.Type == "Tumor" &
    Inferred.identity %in% c(
      "HGSOC (good confidence)",
      "HGSOC (high confidence)",
      "inferred HGSOC (good confidence)",
      "inferred HGSOC (high confidence)"
    )
)
reads_TCGA_MA_HGS <- reads_TCGA_MA[, colnames(reads_TCGA_MA) %in% clinical_TCGA_MA_HGS$Sample.ID] # N=516 HGSOC samples
expr_TCGA_MA_HGS <- as.matrix(reads_TCGA_MA_HGS) # 13039 genes all kept (MA background noise already been handled in previous processing)

# Risk genes quantification in TCGA RNAseq vs. MA------------------------------

valid_genes <- intersect(risk_genes, rownames(expr_TCGA_RNAseq_HGS)) # N = 12/13
valid_genes <- intersect(valid_genes, rownames(expr_TCGA_MA_HGS)) # N = 6/13
common_samples <- intersect(colnames(expr_TCGA_RNAseq_HGS), colnames(expr_TCGA_MA_HGS)) # N = 407 (418*516)
mat_rna <- expr_TCGA_RNAseq_HGS[valid_genes, common_samples]
mat_ma <- expr_TCGA_MA_HGS[valid_genes, common_samples]
risk_gene_cor <- lapply(valid_genes, function(gene) {
  expr_r <- as.numeric(mat_rna[gene, ])
  expr_m <- as.numeric(mat_ma[gene, ])
  if(sd(expr_r, na.rm=TRUE) == 0 | sd(expr_m, na.rm=TRUE) == 0) {
    return(data.frame(GeneSymbol = gene, R = 0, P.Value = 1))
  }
  res <- cor.test(expr_r, expr_m, method = "pearson")
  data.frame(
    GeneSymbol = gene,
    R = res$estimate,
    P.Value = res$p.value
  )
})
risk_gene_cor_df <- bind_rows(risk_gene_cor) %>%
  arrange(R) %>%
  mutate(
    GeneSymbol = factor(GeneSymbol, levels = GeneSymbol),
    Consistency = case_when(
      R >= 0.6 ~ "High (R ≥ 0.6)",
      R >= 0.3 ~ "Moderate (0.3 ≤ R < 0.6)",
      TRUE ~ "Low (R < 0.3)"
    )
  )
robust_genes <- risk_gene_cor_df %>% filter(R >= 0.3) %>% pull(GeneSymbol) %>% as.character() # N = 5/13
p_lollipop <- ggplot(risk_gene_cor_df, aes(x = GeneSymbol, y = R, color = Consistency)) +
  geom_segment(aes(x = GeneSymbol, xend = GeneSymbol, y = 0, yend = R), 
               linewidth = 1, alpha = 0.7) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "#2166AC", alpha = 0.5) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "#B2182B", alpha = 0.5) +
  scale_color_manual(values = c(
    "High (R ≥ 0.6)" = "#2166AC",      
    "Moderate (0.3 ≤ R < 0.6)" = "#F4A582",
    "Low (R < 0.3)" = "#B2182B"         
  )) +
  coord_flip() +
  labs(
    title = "TCGA RNAseq vs. Microarray Quantification of Risk Genes",
    x = "Genes",
    y = "Pearson Correlation Coefficient (R)"
  ) +
  theme_kd() +
  theme(axis.text.y = element_text(size = 7, face = "bold", color = "black", lineheight = 0.8)) # for gene names, use larger size
cairo_pdf("Figures_raw/04_S1_DEGS_cross_platform_quantification.pdf", 
          width = 100 / 25.4, 
          height = 85 / 25.4) 
print(p_lollipop)
dev.off()

# Risk genes validation in TCGA-MA---------------------------------------------

library(survival)
library(broom)

# Match patient survival
survival_TCGA <- clinical_GDC %>%
  filter(in.TCGA.OV == "Y") %>%
  dplyr::select(Patient_ID = Master.Unique.Case.ID, 
                OS_Time_TCGA = OS..days..GDC.clinical.label., 
                OS_Status_TCGA = Death..GDC.clinical.label., 
                DFS_Time_TCGA = DFS..months..cBioPortal.clinical.label.,
                DFS_Status_TCGA = Disease.status..cBioPortal.clinical.label.) %>%
  distinct(Patient_ID, .keep_all = TRUE) %>%
  mutate(
    OS_Time_TCGA = as.numeric(OS_Time_TCGA),
    OS_Status_TCGA = as.integer(substr(OS_Status_TCGA, 1, 1)),
    DFS_Time_TCGA = as.integer(as.numeric(DFS_Time_TCGA) * 30), # months in cBioPortal data
    DFS_Status_TCGA = as.integer(substr(DFS_Status_TCGA, 1, 1))
  )
clinical_TCGA_RNAseq_HGS <- clinical_TCGA_RNAseq_HGS %>%
  mutate(Patient_ID = substr(Sample.ID, 1, 12)) %>%
  left_join(survival_TCGA, by = "Patient_ID")
clinical_TCGA_MA_HGS <- clinical_TCGA_MA_HGS %>%
  mutate(Patient_ID = substr(Sample.ID, 1, 12)) %>%
  left_join(survival_TCGA, by = "Patient_ID")
function_cox_batch_tcga <- function(
    gene_list = robust_genes,               
    expr_mat = expr_TCGA_MA_HGS,                
    clin_df = clinical_TCGA_MA_HGS,                 
    time_col,    
    status_col        
) {
  expr_samples <- colnames(expr_mat)
  clin_samples <- clin_df$Sample.ID
  common_samples <- intersect(expr_samples, clin_samples)
  target_genes <- intersect(gene_list, rownames(expr_mat))
  expr_sub <- expr_mat[target_genes, common_samples, drop = FALSE]
  clin_sub <- clin_df[match(common_samples, clin_df$Sample.ID), , drop = FALSE] # make sure sample order matched
  cat("Calculating Cox regression for", length(target_genes), "target degs...\n")
  res_list <- lapply(target_genes, function(gene) {
    tmp_data <- data.frame(
      Time = as.numeric(as.character(clin_sub[[time_col]])),
      Status = as.numeric(as.character(clin_sub[[status_col]])),
      Expression = as.numeric(expr_sub[gene, ])
    )
    tmp_data <- na.omit(tmp_data)
    if (sd(tmp_data$Expression, na.rm = TRUE) == 0) return(NULL)
    fit <- try(coxph(Surv(Time, Status) ~ Expression, data = tmp_data), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL) 
    tidy_fit <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      select(
        HR = estimate, 
        LCI = conf.low, 
        UCI = conf.high, 
        P.Value = p.value
      ) %>%
      mutate(GeneSymbol = gene)
    return(tidy_fit)
  })
  final_res <- bind_rows(res_list) %>%
    mutate(
      adj.P.Val = p.adjust(P.Value, method = "fdr"),
      log2HR = log2(HR),
      Survival_Type = case_when(
        P.Value < 0.05 & HR > 1 ~ "Risk in TCGA",
        P.Value < 0.05 & HR < 1 ~ "Protective in TCGA",
        TRUE ~ "Not Significant"
      )
    ) %>%
    relocate(GeneSymbol, HR, LCI, UCI, log2HR, P.Value, adj.P.Val, Survival_Type) %>%
    arrange(P.Value) 
  cat("Done! Processed", nrow(final_res), "genes.\n")
  return(final_res)
}
cox_TCGA_OS <- function_cox_batch_tcga(time_col = "OS_Time_TCGA", status_col = "OS_Status_TCGA")
cox_TCGA_DFS <- function_cox_batch_tcga(time_col = "DFS_Time_TCGA", status_col = "DFS_Status_TCGA")

# Visualisation
library(grid)
library(forestploter)
tcga_plot_data <- full_join(
  cox_TCGA_OS %>% select(GeneSymbol, HR, LCI, UCI, P.Value),
  cox_TCGA_DFS %>% select(GeneSymbol, HR, LCI, UCI, P.Value),
  by = "GeneSymbol",
  suffix = c("_OS", "_DFS")
) %>%
  arrange(desc(HR_OS)) 
format_hr <- function(hr, lci, uci) {
  ifelse(is.na(hr), "", sprintf("%.2f (%.2f-%.2f)", hr, lci, uci))
}
format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001 ***",
    p < 0.01 ~ sprintf("%.3f **", p),
    p < 0.05 ~ sprintf("%.3f *", p),
    TRUE ~ sprintf("%.3f", p)
  )
}
tcga_plot_data <- tcga_plot_data %>%
  mutate(
    HR_text_OS  = format_hr(HR_OS, LCI_OS, UCI_OS),
    P_text_OS   = format_p(P.Value_OS),
    HR_text_DFS = format_hr(HR_DFS, LCI_DFS, UCI_DFS),
    P_text_DFS  = format_p(P.Value_DFS)
  )
space_placeholder <- paste(rep(" ", 20), collapse = " ")
plot_df <- data.frame(
  Variable = tcga_plot_data$GeneSymbol,
  N = as.character(length(common_samples)), 
  `TCGA OS` = space_placeholder,            
  `HR (95% CI) ` = tcga_plot_data$HR_text_OS, 
  `P Value ` = tcga_plot_data$P_text_OS,
  `TCGA DFS` = space_placeholder,                 
  `HR (95% CI)` = tcga_plot_data$HR_text_DFS,
  `P Value` = tcga_plot_data$P_text_DFS,
  check.names = FALSE
)
est_list   <- list(tcga_plot_data$HR_OS, tcga_plot_data$HR_DFS)
lower_list <- list(tcga_plot_data$LCI_OS, tcga_plot_data$LCI_DFS)
lower_list <- lapply(lower_list, function(x) ifelse(x <= 0, 0.05, x)) 
upper_list <- list(tcga_plot_data$UCI_OS, tcga_plot_data$UCI_DFS)
upper_list <- lapply(upper_list, function(x) ifelse(x > 10, 10, x)) 
tm <- forest_theme(
  base_size = 7,
  base_family = "Arial",
  ci_pch = 15,           
  ci_col = "black",  
  ci_fill = "navy",      
  ci_alpha = 1,        
  ci_lty = 1,            
  ci_lwd = 0.75,         
  ci_Theight = 0.2,      
  refline_lwd = 0.5,     
  refline_lty = "dashed",  
  refline_col = "grey50",
  core = list(padding = unit(c(2, 2), "mm")) 
)
p_forest_tcga <- forest(
  plot_df,
  est = est_list,
  lower = lower_list,
  upper = upper_list,
  ci_column = c(3, 6),                             
  ref_line = 1, 
  sizes = 0.4,
  x_trans = "log",
  xlim = list(c(0.5, 3), c(0.5, 3)),
  ticks_at = c(0.5, 1, 1.5, 2, 3), 
  arrow_lab = c("Protective", "Risk"),
  theme = tm
) %>% 
  add_border(part = "header", where = "top") %>%
  add_border(part = "header", where = "bottom") %>%
  add_border(row = nrow(plot_df), where = "bottom") %>% 
  edit_plot(
    part = "header",       
    col = c(3, 6), 
    hjust = 0.5,           
    x = unit(0.5, "npc")                
  )
cairo_pdf(file = "Figures_raw/04_01_Cox_regression_TCGA.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4) 
print(p_forest_tcga)
dev.off()

# Risk stratification by PORCN-----------------------------------------------

library(survminer)
clinical_TCGA_MA_HGS <- clinical_TCGA_MA_HGS %>%
  mutate(
    PORCN_Expr = as.numeric(expr_TCGA_MA_HGS["PORCN", as.character(Sample.ID)]),
    PORCN_Group = case_when(
      is.na(PORCN_Expr) ~ NA_character_,
      PORCN_Expr >= quantile(PORCN_Expr, probs = 2/3, na.rm = TRUE) ~ "High",
      PORCN_Expr <= quantile(PORCN_Expr, probs = 1/3, na.rm = TRUE) ~ "Low",
      TRUE ~ "Intermediate"),
    PORCN_Group = factor(PORCN_Group, levels = c("Low", "Intermediate", "High"))
  )
clinical_TCGA_RNAseq_HGS <- clinical_TCGA_RNAseq_HGS %>%
  mutate(
    PORCN_Expr = as.numeric(expr_TCGA_RNAseq_HGS["PORCN", as.character(Sample.ID)]),
    PORCN_Group = case_when(
      is.na(PORCN_Expr) ~ NA_character_,
      PORCN_Expr >= quantile(PORCN_Expr, probs = 2/3, na.rm = TRUE) ~ "High",
      PORCN_Expr <= quantile(PORCN_Expr, probs = 1/3, na.rm = TRUE) ~ "Low",
      TRUE ~ "Intermediate"),
    PORCN_Group = factor(PORCN_Group, levels = c("Low", "Intermediate", "High"))
  ) # just for further analysis
expr_OCM_HGS <- readRDS("Data/expr_OCM_HGS.rds")
clinical_OCM_HGS <- readRDS("Data/Clinical_OCM_HGS.rds") %>%
  rownames_to_column(var = "Sample.ID") %>%
  mutate(
    PORCN_Expr = as.numeric(expr_OCM_HGS["PORCN", as.character(Sample.ID)]),
    PORCN_Group = case_when(
      is.na(PORCN_Expr) ~ NA_character_,
      PORCN_Expr >= quantile(PORCN_Expr, probs = 2/3, na.rm = TRUE) ~ "High",
      PORCN_Expr <= quantile(PORCN_Expr, probs = 1/3, na.rm = TRUE) ~ "Low",
      TRUE ~ "Intermediate"),
    PORCN_Group = factor(PORCN_Group, levels = c("Low", "Intermediate", "High"))
  ) # just for further analysis
function_km_porcn_extremes <- function(df, time_col, status_col, title_text) {
  df_clean <- df %>%
    filter(!is.na(.data[[time_col]]) & !is.na(.data[[status_col]])) %>%
    filter(PORCN_Group %in% c("Low", "High")) %>%
    mutate(PORCN_Group = droplevels(PORCN_Group))
  form <- as.formula(paste("Surv(", time_col, ",", status_col, ") ~ PORCN_Group"))
  fit <- surv_fit(form, data = df_clean)
  p <- ggsurvplot(
    fit,
    data = df_clean,
    pval = TRUE,
    pval.size = 2.5, # 7.11pt as basic
    conf.int = FALSE,
    palette = c("navy", "firebrick"),
    title = title_text,
    xlab = "Time (Days)",
    ylab = "Survival Probability",
    legend.title = "PORCN Tertiles",
    legend.labs = c("Low", "High"),
    risk.table = TRUE,
    risk.table.fontsize = 2.5, # 7.11pt as basic
    risk.table.height = 0.25,
    ggtheme = theme_kd(),
    tables.theme = theme_kd()
  )
  return(p)
}
p_km_os <- function_km_porcn_extremes(
  df = clinical_TCGA_MA_HGS, 
  time_col = "OS_Time_TCGA", 
  status_col = "OS_Status_TCGA", 
  title_text = "TCGA OS vs. PORCN"
)
p_km_dfs <- function_km_porcn_extremes(
  df = clinical_TCGA_MA_HGS, 
  time_col = "DFS_Time_TCGA", 
  status_col = "DFS_Status_TCGA", 
  title_text = "TCGA DFS vs. PORCN"
)
plot_list <- list(p_km_os, p_km_dfs)
p_km__combined <- arrange_ggsurvplots(
  plot_list,
  print = TRUE,
  ncol = 2,
  nrow = 1
)
ggsave(
  filename = "Figures_raw/04_02_KM_PORCN_TCGA.pdf", 
  plot = p_km__combined, 
  width = 170 / 25.4, 
  height = 100 / 25.4, 
  device = cairo_pdf
)

# Keep RDS--------------------------------------------------------------------

save.image("04.RData")
saveRDS(reads_TCGA_RNAseq_HGS, file = "Data/reads_TCGA_RNAseq_HGS.rds") # for limma pipeline
saveRDS(expr_TCGA_RNAseq_HGS, file = "Data/expr_TCGA_RNAseq_HGS.rds")
saveRDS(expr_TCGA_MA_HGS, file = "Data/expr_TCGA_MA_HGS.rds")
saveRDS(clinical_GDC, file = "Data/clinical_GDC.rds")
saveRDS(clinical_OCM_HGS, file = "Data/clinical_OCM_HGS_PORCN.rds")
saveRDS(clinical_TCGA_RNAseq_HGS, file = "Data/clinical_TCGA_RNAseq_HGS_PORCN.rds")
saveRDS(clinical_TCGA_MA_HGS, file = "Data/clinical_TCGA_MA_HGS_PORCN.rds")