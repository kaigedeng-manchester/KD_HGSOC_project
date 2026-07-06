# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("10.RData") 

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

# Define groups----------------------------------------------------------------

# Get expression data
expr_OCM_HGS <- readRDS("Data/expr_OCM_HGS.rds")
expr_TCGA_RNAseq_HGS <- readRDS("Data/expr_TCGA_RNAseq_HGS.rds")
expr_TCGA_MA_HGS <- readRDS("Data/expr_TCGA_MA_HGS.rds")

# Add new_CCLE data
library(data.table)
reads_CCLE <- fread("Data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", check.names = FALSE) %>% as.data.frame() # N=1754, fread() better for very big data
clinical_CCLE <- read.csv("Data/Clinical_CCLE.csv", row.names=1, na.strings = "")
clinical_CCLE_HGS <- clinical_CCLE %>% filter(NMF.subtype == "HGSOC" | Inferred.identity..Genomic.QC.KD.18.03.26. == "HGSOC") # N=28
reads_CCLE_HGS <- reads_CCLE[reads_CCLE$ModelID %in% rownames(clinical_CCLE_HGS), ] # N=26
reads_CCLE_HGS <- reads_CCLE_HGS %>%
  remove_rownames() %>%
  mutate(Cell_Line_Name = clinical_CCLE_HGS$cellLineDisplayName[match(ModelID, rownames(clinical_CCLE_HGS))]) %>%
  column_to_rownames(var = "Cell_Line_Name") %>%
  select(-c(1:6)) %>% # annotations from CCLE data
  rename_with(~ gsub(" \\(.*\\)", "", .x))

# Duplicates
colnames(reads_CCLE_HGS) %>% length() #19215 genes in total (coding genes)
colnames(reads_CCLE_HGS) %>% unique() %>% length() #19215 No duplicate

# Get matrix
expr_CCLE_HGS <- reads_CCLE_HGS %>% as.matrix() %>% t()

# Filter low expression (keep genes expressed in >=20% samples with TPM > 1)
library(limma)
library(edgeR)
keep <- rowSums(expr_CCLE_HGS > 1) >= (0.2 * ncol(expr_CCLE_HGS))
expr <- expr_CCLE_HGS[keep, ] # keep 13372 genes from 19215 coding genes

# Output 
expr_CCLE_HGS <- expr
rm(expr, keep)

# PORCN expression groups
function_get_extreme_groups <- function(expr_mat, gene_name, dataset_name) {
  if (!(gene_name %in% rownames(expr_mat))) {
    cat("⚠️ Warning:", gene_name, "not found in", dataset_name, "\n")
    return(NULL)
  }
  expr_vals <- as.numeric(expr_mat[gene_name, ])
  names(expr_vals) <- colnames(expr_mat)
  p33 <- quantile(expr_vals, 1/3, na.rm = TRUE)
  p66 <- quantile(expr_vals, 2/3, na.rm = TRUE)
  low_samples <- names(expr_vals)[expr_vals <= p33]
  high_samples <- names(expr_vals)[expr_vals >= p66]
  cat(sprintf("✅ %-25s | Low (<= %.2f): %3d samples | High (>= %.2f): %3d samples\n", 
              dataset_name, p33, length(low_samples), p66, length(high_samples)))
  return(list(Low = low_samples, High = high_samples))
}
cat("\n📊 PORCN Extreme Tertiles Sample Sizes:\n")
cat("--------------------------------------------------------------------------------\n")
groups_ocm <- function_get_extreme_groups(expr_OCM_HGS, "PORCN", "OCM RNAseq")
groups_tcga_rna <- function_get_extreme_groups(expr_TCGA_RNAseq_HGS, "PORCN", "TCGA RNAseq")
groups_tcga_ma <- function_get_extreme_groups(expr_TCGA_MA_HGS, "PORCN", "TCGA Microarray")
groups_ccle <- function_get_extreme_groups(expr_CCLE_HGS, "PORCN", "CCLE Cell Lines")
cat("--------------------------------------------------------------------------------\n")

# Clinical data preparation----------------------------------------------------

# Get clinical data
clinical_OCM_HGS <- readRDS("Data/clinical_OCM_HGS_PORCN.rds")
clinical_TCGA_RNAseq_HGS <- readRDS("Data/clinical_TCGA_RNAseq_HGS_PORCN.rds")
clinical_TCGA_MA_HGS <- readRDS("Data/clinical_TCGA_MA_HGS_PORCN.rds")
clinical_GDC <- readRDS("Data/clinical_GDC.rds")

# (1) Match TCGA baseline data
library(stringr)
baseline_TCGA <- clinical_GDC %>%
  filter(in.TCGA.OV == "Y") %>%
  select(
    Patient_ID = Master.Unique.Case.ID,
    Age = Age..GDC.clincal.label.,                                 
    Stage = Stage..GDC.clinical.label.                 
  ) %>%
  distinct(Patient_ID, .keep_all = TRUE) %>%
  mutate(Age = as.numeric(Age),
         Stage = str_remove_all(Stage, "(?i)Stage\\s*"),
         Stage = str_replace_all(Stage, c("IV" = "4", 
                                          "III" = "3", 
                                          "II" = "2", 
                                          "I" = "1"))
  )
clinical_TCGA_RNAseq_HGS <- clinical_TCGA_RNAseq_HGS %>%
  left_join(baseline_TCGA, by = "Patient_ID")
clinical_TCGA_MA_HGS <- clinical_TCGA_MA_HGS %>%
  left_join(baseline_TCGA, by = "Patient_ID")

# (2) Match TCGA HRD data
HRD_TCGA <- read_csv("Data/hrd_features_TCGA_OV_clean.csv", show_col_types = FALSE)
HRD_TCGA <- HRD_TCGA %>%
  mutate(Patient_ID = substr(Sample_ID, 1, 12)) %>%
  arrange(Patient_ID, desc(HRD_Sum)) %>%
  distinct(Patient_ID, .keep_all = TRUE) %>%
  select(Patient_ID, HRD_LOH, LST, TAI, scarHRD_score = HRD_Sum)

# BRCA1/2 somatic mutation
maf <- fread("Data/mc3.v0.2.8.PUBLIC.maf.gz")
maf_TCGA_combined <- maf[substr(Tumor_Sample_Barcode, 1, 16) %in% 
                           union(clinical_TCGA_RNAseq_HGS$Sample.ID, clinical_TCGA_MA_HGS$Sample.ID)]
maf_TCGA_combined$Variant_Classification %>% factor() %>% summary()
non_silence <- c(
  "Frame_Shift_Del", 
  "Frame_Shift_Ins", 
  "Nonsense_Mutation", 
  "Splice_Site", 
  "Translation_Start_Site", 
  "Nonstop_Mutation",
  "In_Frame_Del",        
  "In_Frame_Ins",        
  "Missense_Mutation" 
)
maf_brca_somatic <- maf_TCGA_combined[
  Hugo_Symbol %in% c("BRCA1", "BRCA2") & 
    Variant_Classification %in% non_silence
]
somatic_mutated_samples <- unique(substr(maf_brca_somatic$Tumor_Sample_Barcode, 1, 16)) # N=28
somatic_mutated_patients <- unique(substr(somatic_mutated_samples, 1, 12)) # N=28
rm(maf)

# BRCA1/2 germline mutation-Access control

# BRCA1 hypermethylation

# Composite Element REF
probe_list <- fread("Data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv", select = 1)
brca1_alt_probes <- c(
  "cg08902545",
  "cg04658354", 
  "cg19531713", 
  "cg16630982", 
  "cg19088651",
  "cg19515236"
)
matched_probes <- intersect(probe_list[[1]], brca1_alt_probes)
print(matched_probes)
if(length(matched_probes) > 0) {
  target_probe <- matched_probes[1]
  row_idx <- which(probe_list[[1]] == target_probe)
  header <- fread("Data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv", nrows = 0)
  meth_brca1_data <- fread("Data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv", skip = row_idx, nrows = 1, header = FALSE)
  colnames(meth_brca1_data) <- colnames(header)
  rm(probe_list)
  gc()
  cat("✅ Extraction perfect for probe [", target_probe, "]! You can now use `meth_brca1_data` to filter for hypermethylation (Beta > 0.3).\n", sep = "")
} else {
  cat("❌ Extremely rare: All alternative BRCA1 probes are missing from the dataset.\n")
}
meth_long <- meth_brca1_data %>%
  pivot_longer(cols = -1, names_to = "Sample_ID", values_to = "Beta_Value") %>%
  mutate(Beta_Value = as.numeric(Beta_Value)) %>%
  filter(!is.na(Beta_Value))
meth_patients <- meth_long %>%
  filter(Beta_Value > 0.3) %>%
  mutate(Patient_ID = substr(Sample_ID, 1, 12)) %>%
  pull(Patient_ID) %>%
  unique()
pure_meth_patients <- setdiff(meth_patients, somatic_mutated_patients) # should be mutual exclusive

# HRD score distribution
df_hrd_validation <- clinical_TCGA_RNAseq_HGS %>%
  left_join(HRD_TCGA, by = "Patient_ID") %>%
  mutate(
    BRCA_Status = case_when(
      Patient_ID %in% somatic_mutated_patients ~ "BRCA1/2_Somatic_Mutant",
      Patient_ID %in% pure_meth_patients ~ "BRCA1_Hypermethylated",
      TRUE ~ "WT_or_Unknown"
    ),
    BRCA_Status = factor(BRCA_Status, levels = c("WT_or_Unknown", "BRCA1/2_Somatic_Mutant", "BRCA1_Hypermethylated"))
  ) %>%
  filter(!is.na(scarHRD_score))
p_density <- ggplot(df_hrd_validation, aes(x = scarHRD_score, fill = BRCA_Status, color = BRCA_Status)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  geom_vline(xintercept = 42, color = "black", linetype = "dashed", linewidth = 0.5) + # for import lines, use 0.5pt
  annotate("text", x = 44, y = 0.026, label = "Clinical Cutoff: 42", angle = 90, fontface = "italic", color = "black") +
  geom_rug(alpha = 0.7, linewidth = 0.5, length = unit(0.04, "npc")) +
  scale_fill_manual(values = c("steelblue", "firebrick", "purple")) + 
  scale_color_manual(values = c("steelblue", "firebrick", "purple")) + 
  theme_kd() +
  labs(
    title = "TCGA scarHRD Landscape by BRCA Status",
    x = "Total scarHRD Score (LOH + TAI + LST)",
    y = "Density",
    fill = "BRCA Status",
    color = "BRCA Status"
  )

# HRD score distribution comparison
library(ggpubr)
comp_brca <- list(
  c("BRCA1/2_Somatic_Mutant", "WT_or_Unknown"),
  c("BRCA1/2_Somatic_Mutant", "BRCA1_Hypermethylated"),
  c("BRCA1_Hypermethylated", "WT_or_Unknown")
)
color_brca <- c("WT_or_Unknown" = "steelblue", 
                "BRCA1/2_Somatic_Mutant" = "firebrick", 
                "BRCA1_Hypermethylated" = "purple")
p_hrd_comparison <- ggplot(df_hrd_validation, aes(x = BRCA_Status, y = scarHRD_score)) +
  geom_violin(aes(fill = BRCA_Status), trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = BRCA_Status), width = 0.15, size = 1, alpha = 0.5) + # large number of dots, use 1pt
  geom_hline(yintercept = 42, color = "black", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 0.8, y = 46, label = "Clincal Cutoff: 42", fontface = "italic") +
  scale_fill_manual(values = color_brca) +
  scale_color_manual(values = color_brca) +
  scale_x_discrete(labels = c(
    "WT_or_Unknown" = "WT or\nUnknown",
    "BRCA1/2_Somatic_Mutant" = "BRCA1/2\nSomatic Mutant",
    "BRCA1_Hypermethylated" = "BRCA1\nHypermethylated"
  )) +
  stat_compare_means(comparisons = comp_brca, method = "wilcox.test", 
                     label = "p.format", tip.length = 0.02, size = 3) +
  theme_kd() + 
  theme(
    axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    title = "TCGA scarHRD Score by BRCA Status",
    x = "",
    y = "Total scarHRD Score (LOH + TAI + LST)",
    fill = "BRCA Status",
    color = "BRCA Status"
  )
plot_brca_sub <- function(data, y_col, y_label, title_text) {
  df_clean <- data %>% filter(!is.na(.data[[y_col]]), !is.na(BRCA_Status))
  p <- ggplot(df_clean, aes(x = BRCA_Status, y = .data[[y_col]])) +
    geom_violin(aes(fill = BRCA_Status), trim = FALSE, alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.5) +
    geom_jitter(aes(color = BRCA_Status), width = 0.15, size = 1, alpha = 0.5) +
    scale_fill_manual(values = color_brca) +
    scale_color_manual(values = color_brca) +
    scale_x_discrete(labels = c(
      "WT_or_Unknown" = "WT or\nUnknown",
      "BRCA1/2_Somatic_Mutant" = "BRCA1/2\nSomatic Mutant",
      "BRCA1_Hypermethylated" = "BRCA1\nHypermethylated"
    )) +
    stat_compare_means(comparisons = comp_brca, method = "wilcox.test", 
                       label = "p.format", tip.length = 0.02, size = 2.8) +
    theme_kd() +
    theme(axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"), 
          legend.position = "none") +
    labs(title = title_text, x = "", y = y_label)
  return(p)
}
p_loh_brca <- plot_brca_sub(df_hrd_validation, "HRD_LOH", "LOH Score", "LOH score by BRCA")
p_lst_brca <- plot_brca_sub(df_hrd_validation, "LST", "LST Score", "LST score by BRCA")
p_tai_brca <- plot_brca_sub(df_hrd_validation, "TAI", "TAI Score", "TAI score by BRCA")
library(patchwork)
p_hrd <- (p_hrd_comparison | p_loh_brca) / 
  (p_lst_brca | p_tai_brca) +
  plot_layout(guides = "collect") & 
  theme(
    aspect.ratio = 1, 
    legend.position = "none",      
    axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black")
  ) &
  guides(
    fill = guide_legend(
      title.position = "top",          
      title.hjust = 0.5,               
      nrow = 2,                        
      byrow = TRUE                     
    ),
    color = "none"
  )
cairo_pdf("Figures_raw/10_S1_TCGA_scarHRD.pdf", 
          width = 170 / 25.4, 
          height = 190 / 25.4) 
print(p_hrd)
dev.off()

# Clinical comparison---------------------------------------------------------

library(gtsummary)
clinical_OCM_HGS_comparison <- clinical_OCM_HGS %>%
  filter(PORCN_Group %in% c("Low", "High")) %>%
  mutate(PORCN_Group = droplevels(PORCN_Group)) %>%
  mutate(
    Age = as.numeric(Age),
    Stage = factor(Stage),
    HRD_Status = ifelse(Seqone_score >= 0.5, "HRD+", "HRD-"),
    HRD_Status = factor(HRD_Status, levels = c("HRD-", "HRD+"))
  )
clinical_TCGA_MA_HGS_comparison <- clinical_TCGA_MA_HGS %>%
  filter(PORCN_Group %in% c("Low", "High")) %>%
  mutate(PORCN_Group = droplevels(PORCN_Group)) %>%
  left_join(HRD_TCGA %>% select("Patient_ID", "scarHRD_score"), by = "Patient_ID") %>%
  mutate(
    Age = as.numeric(Age),
    Stage = factor(Stage),
    HRD_Status = ifelse(scarHRD_score >= 42, "HRD+", "HRD-"),
    HRD_Status = factor(HRD_Status, levels = c("HRD-", "HRD+"))
  )
table_ocm <- clinical_OCM_HGS_comparison %>%
  select(PORCN_Group, Age, Stage, Seqone_score, HRD_Status) %>%
  tbl_summary(
    by = PORCN_Group, 
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})", 
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age (Years)",
      Stage ~ "Clinical Stage",
      Seqone_score ~ "HRD Score (SeqOne)",
      HRD_Status ~ "HRD Status (Cutoff: 0.5)"
    ),
    missing_text = "Unknown"
  ) %>%
  add_p() %>% 
  bold_p(t = 0.05) %>% 
  modify_header(label = "**Baseline Characteristics**") %>%
  modify_caption("**Clinical characteristics of HGS OCMs stratified by PORCN expression**")
table_ocm
table_tcga <- clinical_TCGA_MA_HGS_comparison %>%
  select(PORCN_Group, Age, Stage, scarHRD_score, HRD_Status) %>%
  tbl_summary(
    by = PORCN_Group,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})", 
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age (Years)",
      Stage ~ "Clinical Stage",
      scarHRD_score ~ "HRD Score (scarHRD)",
      HRD_Status ~ "HRD Status (Cutoff: 42)"
    ),
    missing_text = "Unknown"
  ) %>%
  add_p() %>%
  bold_p(t = 0.05) %>% 
  modify_header(label = "**Baseline Characteristics**") %>%
  modify_caption("**Clinical characteristics of HGSOC TCGA samples stratified by PORCN expression**")
table_tcga

# Adding new HRD data into clinical data
clinical_OCM_HGS <- clinical_OCM_HGS %>%
  mutate(
    Age = as.numeric(Age),
    Stage = factor(Stage),
    HRD_Status = ifelse(Seqone_score >= 0.5, "HRD+", "HRD-"),
    HRD_Status = factor(HRD_Status, levels = c("HRD-", "HRD+"))
  )
clinical_TCGA_RNAseq_HGS <- clinical_TCGA_RNAseq_HGS %>%
  left_join(HRD_TCGA %>% select("Patient_ID", "scarHRD_score"), by = "Patient_ID") %>%
  mutate(
    Age = as.numeric(Age),
    Stage = factor(Stage),
    HRD_Status = ifelse(scarHRD_score >= 42, "HRD+", "HRD-"),
    HRD_Status = factor(HRD_Status, levels = c("HRD-", "HRD+"))
  )
clinical_TCGA_MA_HGS <- clinical_TCGA_MA_HGS %>%
  left_join(HRD_TCGA %>% select("Patient_ID", "scarHRD_score"), by = "Patient_ID") %>%
  mutate(
    Age = as.numeric(Age),
    Stage = factor(Stage),
    HRD_Status = ifelse(scarHRD_score >= 42, "HRD+", "HRD-"),
    HRD_Status = factor(HRD_Status, levels = c("HRD-", "HRD+"))
  )

# Keep RDS--------------------------------------------------------------------

save.image("10.RData")
saveRDS(expr_CCLE_HGS, file = "Data/expr_CCLE_HGS.rds") # new RNA dataset
saveRDS(clinical_OCM_HGS, file = "Data/clinical_OCM_HGS.rds") # updated clinical data
saveRDS(clinical_TCGA_RNAseq_HGS, file = "Data/clinical_TCGA_RNAseq_HGS.rds")
saveRDS(clinical_TCGA_MA_HGS, file = "Data/clinical_TCGA_MA_HGS.rds") 
saveRDS(groups_ocm, file = "Data/groups_ocm.rds") # PORCN group labels for further analyses
saveRDS(groups_ccle, file = "Data/groups_ccle.rds")
saveRDS(groups_tcga_ma, file = "Data/groups_tcga_ma.rds")
saveRDS(groups_tcga_rna, file = "Data/groups_tcga_rna.rds")
saveRDS(df_hrd_validation, file = "Data/df_hrd_validation.rds")
