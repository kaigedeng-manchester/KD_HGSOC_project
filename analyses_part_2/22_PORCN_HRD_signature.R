# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("22.RData") 

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

# Get signature score----------------------------------------------------------

df_scores <- readRDS("Data/df_51genes_zscored.rds")

# Get expression---------------------------------------------------------------

expr_OCM <- readRDS("Data/expr_OCM_HGS.rds")
expr_TCGA <- readRDS("Data/expr_TCGA_RNAseq_HGS.rds")
expr_CCLE <- readRDS("Data/expr_CCLE_HGS.rds")

# Get clinical----------------------------------------------------------------

clinical_OCM <- readRDS("Data/clinical_OCM_HGS.rds")
clinical_TCGA <- readRDS("Data/clinical_TCGA_RNAseq_HGS.rds")
clinical_CCLE <- read.csv("Data/Clinical_CCLE.csv", row.names=1, na.strings = "")
clinical_CCLE_HGS <- clinical_CCLE %>% filter(NMF.subtype == "HGSOC" | Inferred.identity..Genomic.QC.KD.18.03.26. == "HGSOC") # N=28

# Master df-------------------------------------------------------------------

porcn_ocm_df <- data.frame(Sample.ID = colnames(expr_OCM), 
                           PORCN_Raw = as.numeric(expr_OCM["PORCN", ]), 
                           Cohort = "OCM")
porcn_tcga_df <- data.frame(Sample.ID = colnames(expr_TCGA), 
                            PORCN_Raw = as.numeric(expr_TCGA["PORCN", ]), 
                            Cohort = "TCGA")
porcn_ccle_df <- data.frame(Sample.ID = colnames(expr_CCLE), 
                            PORCN_Raw = as.numeric(expr_CCLE["PORCN", ]), 
                            Cohort = "CCLE")
porcn_all_df <- bind_rows(porcn_ocm_df, porcn_tcga_df, porcn_ccle_df)
master_df <- df_scores %>%
  left_join(porcn_all_df, by = c("Sample.ID", "Cohort")) %>%
  group_by(Cohort) %>%
  mutate(
    sig_q33 = quantile(ULM_Score, 1/3, na.rm = TRUE),
    sig_q67 = quantile(ULM_Score, 2/3, na.rm = TRUE),
    porcn_q33 = quantile(PORCN_Raw, 1/3, na.rm = TRUE),
    porcn_q67 = quantile(PORCN_Raw, 2/3, na.rm = TRUE),
    Signature_Class = case_when(
      ULM_Score >= sig_q67 ~ "High",
      ULM_Score <= sig_q33 ~ "Low",
      TRUE ~ "Medium"
    ),
    PORCN_Class = case_when(
      PORCN_Raw >= porcn_q67 ~ "High",
      PORCN_Raw <= porcn_q33 ~ "Low",
      TRUE ~ "Medium"
    ),
    PORCN_Z = as.numeric(scale(PORCN_Raw))
  ) %>%
  ungroup() %>%
  mutate(
    Signature_Class = factor(Signature_Class, levels = c("High", "Medium", "Low")),
    PORCN_Class = factor(PORCN_Class, levels = c("High", "Medium", "Low"))
  ) %>%
  select(-sig_q33, -sig_q67, -porcn_q33, -porcn_q67)

# Add CIN data-----------------------------------------------------------------

HRD_TCGA <- read_csv("Data/hrd_features_TCGA_OV_clean.csv", show_col_types = FALSE) %>%
  mutate(Patient_ID = substr(Sample_ID, 1, 12)) %>%
  arrange(Patient_ID, desc(HRD_Sum)) %>%
  distinct(Patient_ID, .keep_all = TRUE) %>%
  select(Patient_ID, scarHRD_score = HRD_Sum, LOH_score = HRD_LOH, LST_score = LST, TAI_score = TAI)
master_df <- master_df %>%
  left_join(
    clinical_OCM %>% select(Sample.ID, Seqone_score), 
    by = "Sample.ID"
  ) %>%
  mutate(Patient_ID = ifelse(Cohort == "TCGA", substr(Sample.ID, 1, 12), NA)) %>%
  left_join(HRD_TCGA, by = "Patient_ID") %>%
  select(-Patient_ID)

# Add origin data--------------------------------------------------------------

origin_df <- read_csv("Data/signature_scores.csv", show_col_types = FALSE)
setdiff(origin_df$Sample.ID, master_df$Sample.ID) #check
master_df <- master_df %>%
  left_join(origin_df, by = "Sample.ID")

# Add subtyping data-----------------------------------------------------------

ocm_clusters <- read_csv("Data/Cluster_membership_for_lab.csv", show_col_types = FALSE)
setdiff(ocm_clusters$Sample.ID, master_df$Sample.ID) #check
ocm_subtype_df <- ocm_clusters %>%
  select(Sample.ID, Label_k3) %>%
  mutate(Value = 1) %>%
  pivot_wider(
    names_from = Label_k3,
    values_from = Value,
    names_prefix = "Subtype_OCM_",
    values_fill = list(Value = 0)
  )
protype_mor_df <- read_csv("Data/PrOTYPE.csv")
protype_mor_df <- protype_mor_df %>%
  group_by(source, target) %>%
  filter(n() == 1) %>% 
  ungroup()
library(decoupleR)
tcga_ulm_res <- run_ulm(
  mat = expr_TCGA,
  network = protype_mor_df,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)
tcga_subtype_df <- tcga_ulm_res %>%
  select(Sample.ID = condition, source, score) %>%
  pivot_wider(
    names_from = source,
    values_from = score,
    names_prefix = "Subtype_TCGA_"
  )
master_df <- master_df %>%
  left_join(ocm_subtype_df, by = "Sample.ID") %>%
  left_join(tcga_subtype_df, by = "Sample.ID")

# Add response data-----------------------------------------------------------

sensitivity_OCM_HGS <- readRDS("Data/Clinical_merged_OCM_HGS.rds")
ocm_chemo_df <- sensitivity_OCM_HGS %>%
  rownames_to_column("Sample.ID") %>%
  select(Sample.ID, carboplatin_log10_IC50 = log10_carboplatin_IC50, taxol_log10_IC50 = log10_taxol_IC50)
setdiff(ocm_chemo_df$Sample.ID, master_df$Sample.ID) #check
ocm_parp_df <- read_csv("Data/OCM_parp.csv", show_col_types = FALSE) %>%
  select(Sample.ID, PARPi_IC50 = PARPi, PARGi_IC50 = PARGi) %>%
  mutate(
    PARPi_log10_IC50 = log10(PARPi_IC50),
    PARGi_log10_IC50 = log10(PARGi_IC50)
  ) %>%
  select(-PARPi_IC50, -PARGi_IC50)
setdiff(ocm_parp_df$Sample.ID, master_df$Sample.ID) #check
library(data.table)
prism_mat <- fread("Data/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv") %>% as.data.frame()
prism_info <- fread("Data/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv") %>% as.data.frame()
prism_dict <- prism_info %>%
  select(Broad_ID = IDs, Drug_Name = Drug.Name, MOA, Target = repurposing_target) %>%
  distinct(Broad_ID, .keep_all = TRUE)
clean_drug_names <- c(
  "carboplatin"  = "carboplatin",
  "paclitaxel"   = "taxol",
  "olaparib"     = "olaparib",
  "niraparib"    = "niraparib",
  "pdd-00017273" = "PARGi_PDD"
)
target_drugs_regex <- "(?i)^(carboplatin|paclitaxel|olaparib|niraparib|PDD-00017273)$"
prism_selected_drugs <- prism_dict %>%
  filter(grepl(target_drugs_regex, Drug_Name)) %>%
  mutate(Lower_Name = tolower(Drug_Name)) %>%
  group_by(Lower_Name) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Standard_Name = clean_drug_names[Lower_Name]) %>%
  select(Broad_ID, Standard_Name)
ccle_bridge <- clinical_CCLE_HGS %>%
  rownames_to_column("ModelID") %>%
  select(ModelID, Sample.ID = cellLineDisplayName)
ccle_prism_df <- prism_mat %>%
  rename(Broad_ID = 1) %>% 
  filter(Broad_ID %in% prism_selected_drugs$Broad_ID) %>% 
  pivot_longer(cols = -Broad_ID, names_to = "ModelID", values_to = "Viability") %>%
  drop_na(Viability) %>%
  left_join(prism_selected_drugs, by = "Broad_ID") %>%
  inner_join(ccle_bridge, by = "ModelID") %>%
  mutate(Feature_Name = paste0(Standard_Name, "_PRISM", "log2_viab")) %>%
  select(Sample.ID, Feature_Name, Viability) %>%
  pivot_wider(
    names_from = Feature_Name,
    values_from = Viability,
    values_fn = mean 
  )
setdiff(ccle_prism_df$Sample.ID, master_df$Sample.ID) #check
master_df <- master_df %>%
  left_join(ocm_chemo_df, by = "Sample.ID") %>%
  left_join(ocm_parp_df, by = "Sample.ID") %>%
  left_join(ccle_prism_df, by = "Sample.ID")

# Sample-resolution annotation and comparison---------------------------------

all_features <- c(
  # Level 1: CIN / HRD Scars
  "Seqone_score",   
  "scarHRD_score",  
  "LOH_score",
  "LST_score",
  "TAI_score",
  # Level 2: Origin
  "FNE_score",
  "OSE_score",
  # Level 3: Subtyping
  "Subtype_OCM_Alpha",
  "Subtype_OCM_Beta",
  "Subtype_OCM_Gamma",
  "Subtype_TCGA_C1.MES",
  "Subtype_TCGA_C2.IMM",
  "Subtype_TCGA_C4.DIFF",
  "Subtype_TCGA_C5.PRO",
  # Level 4: Pharmacogenomics (Drug Response)
  # --- OCM Cohort (log10 IC50) ---
  "carboplatin_log10_IC50",
  "taxol_log10_IC50",
  "PARPi_log10_IC50",
  "PARGi_log10_IC50",
  # --- CCLE Cohort (PRISM log2 Viability) ---
  "carboplatin_PRISMlog2_viab",
  "taxol_PRISMlog2_viab",
  "olaparib_PRISMlog2_viab",
  "niraparib_PRISMlog2_viab",
  "PARGi_PDD_PRISMlog2_viab"
)

# Groups
plot_df <- master_df %>%
  filter(
    (Cohort == "OCM" & HRD_Status %in% c("HRD", "HRP")) |
      (Cohort == "TCGA" & HRD_Status %in% c("HRD", "HRP")) |
      (Cohort == "CCLE")
  ) %>%
  mutate(
    Split_Group = case_when(
      Cohort == "OCM" & HRD_Status == "HRD" ~ "OCM\n(HRD)",
      Cohort == "OCM" & HRD_Status == "HRP" ~ "OCM\n(HRP)",
      Cohort == "TCGA" & HRD_Status == "HRD" ~ "TCGA\n(HRD)",
      Cohort == "TCGA" & HRD_Status == "HRP" ~ "TCGA\n(HRP)",
      Cohort == "CCLE" ~ "CCLE\n(All)",
      TRUE ~ "Other"
    ),
    Split_Group = factor(Split_Group, levels = c("OCM\n(HRD)", "OCM\n(HRP)", "TCGA\n(HRD)", "TCGA\n(HRP)", "CCLE\n(All)"))
  ) %>%
  arrange(Split_Group, desc(ULM_Z))
all_features %in% colnames(plot_df) #check

# Annotations
sig_colors <- case_when(plot_df$Signature_Class == "High" ~ "firebrick", plot_df$Signature_Class == "Low" ~ "steelblue", TRUE ~ "grey80")
porcn_colors <- case_when(plot_df$PORCN_Class == "High" ~ "firebrick", plot_df$PORCN_Class == "Low" ~ "steelblue", TRUE ~ "grey80")
top_anno <- HeatmapAnnotation(
  `Signature Z` = anno_barplot(plot_df$ULM_Z, bar_width = 0.8, gp = gpar(fill = sig_colors, col = NA), height = unit(2.5, "cm")),
  `PORCN Z` = anno_barplot(plot_df$PORCN_Z, bar_width = 0.8, gp = gpar(fill = porcn_colors, col = NA), height = unit(2.5, "cm")),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8, fontface = "bold"), gap = unit(2, "mm")
)

# Main plot
ht_list <- NULL
for (i in seq_along(all_features)) {
  feat <- all_features[i]
  if (!(feat %in% colnames(plot_df))) next
  feat_vals <- plot_df[[feat]]
  is_binary_feature <- length(unique(na.omit(feat_vals))) <= 2 && all(na.omit(feat_vals) %in% c(0, 1))
  aux_val <- NULL
  if (feat == "Seqone_score") {
    aux_val <- 0.5
  } else if (feat == "scarHRD_score") {
    aux_val <- 42
  }
  y_range_vals <- na.omit(c(feat_vals, aux_val))
  if (length(y_range_vals) == 0) {
    y_min <- 0; y_max <- 1
  } else {
    y_min <- min(y_range_vals)
    y_max <- max(y_range_vals)
  }
  y_range <- y_max - y_min
  if (y_range == 0) y_range <- 1
  y_min_pad <- y_min - y_range * 0.1
  y_max_pad <- y_max + y_range * 0.15
  row_label <- feat
  mat <- matrix(feat_vals, nrow = 1, dimnames = list(row_label, plot_df$Sample.ID))
  if (is_binary_feature || length(y_range_vals) == 0) {
    y_axis_anno <- NULL 
  } else {
    y_breaks <- pretty(c(y_min, y_max), n = 3) 
    y_breaks <- y_breaks[y_breaks >= y_min_pad & y_breaks <= y_max_pad]
    y_axis_anno <- rowAnnotation(
      Scale = function(index) {
        pushViewport(viewport(xscale = c(0, 1), yscale = c(y_min_pad, y_max_pad)))
        grid.segments(x0 = unit(0, "npc"), y0 = unit(y_breaks, "native"), 
                      x1 = unit(1.0, "mm"), y1 = unit(y_breaks, "native"), 
                      gp = gpar(col = "black", lwd = 0.5))
        grid.text(y_breaks, x = unit(1.5, "mm"), y = unit(y_breaks, "native"), 
                  just = "left", gp = gpar(fontsize = 4, col = "black"))
        popViewport()
      },
      width = unit(5, "mm"),
      show_annotation_name = FALSE
    )
  }
  lollipop_layer <- function(j, i, x, y, w, h, fill) {
    vals <- pindex(mat, i, j)
    valid <- !is.na(vals)
    if (any(valid)) {
      grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8)) # 底层黑框 
      x_v <- x[valid]
      v_v <- as.numeric(vals[valid]) 
      w_v <- w[valid]
      if (is_binary_feature) {
        idx_yes <- which(v_v == 1)
        if (length(idx_yes) > 0) {
          grid.rect(
            x = x_v[idx_yes], 
            y = unit(0.5, "npc"), 
            width = w_v[idx_yes], 
            height = unit(1, "npc"), 
            gp = gpar(fill = "firebrick", col = "white", lwd = 0.5)
          )
        }
      } else {
        block_mean <- mean(v_v, na.rm = TRUE)
        y_scaled <- (v_v - y_min_pad) / (y_max_pad - y_min_pad)
        y_ref_scaled <- (block_mean - y_min_pad) / (y_max_pad - y_min_pad)
        local_min <- min(v_v, na.rm = TRUE)
        local_max <- max(v_v, na.rm = TRUE)
        local_mid <- median(v_v, na.rm = TRUE)
        if (local_max > local_min) {
          local_col_fun <- colorRamp2(c(local_min, local_mid, local_max), c("navy", "grey80", "firebrick"))
          pt_cols <- local_col_fun(v_v)
        } else {
          pt_cols <- rep("grey50", length(v_v))
        }
        if (!is.null(aux_val)) {
          aux_scaled <- (aux_val - y_min_pad) / (y_max_pad - y_min_pad)
          grid.lines(x = c(0, 1), y = c(aux_scaled, aux_scaled), gp = gpar(col = "firebrick", lty = 3, lwd = 1))
        }
        grid.lines(x = c(0, 1), y = c(y_ref_scaled, y_ref_scaled), gp = gpar(col = "grey30", lty = 2, lwd = 0.5))
        grid.segments(x0 = x_v, y0 = y_ref_scaled, x1 = x_v, y1 = y_scaled, gp = gpar(col = "grey50", lwd = 0.5))
        grid.points(x = x_v, y = y_scaled, pch = 16, size = unit(2.5, "pt"), gp = gpar(col = pt_cols))
      }
    } else {
      grid.rect(gp = gpar(fill = "grey92", col = "black", lwd = 0.8))
      grid.text("N/A", gp = gpar(fontsize = 6, col = "grey60", fontface = "italic", fontfamily = "sans"))
    }
  }
  if (i == 1) {
    ht <- Heatmap(
      mat, name = feat,
      cluster_columns = FALSE, show_column_names = FALSE, 
      column_split = plot_df$Split_Group, 
      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
      top_annotation = top_anno, 
      right_annotation = y_axis_anno, 
      row_names_side = "left", row_names_gp = gpar(fontsize = 8, fontface = "bold"),
      height = unit(1.0, "cm"), 
      rect_gp = gpar(type = "none"), 
      layer_fun = lollipop_layer,
      show_heatmap_legend = FALSE
    )
  } else {
    ht <- Heatmap(
      mat, name = feat,
      cluster_columns = FALSE, show_column_names = FALSE, 
      right_annotation = y_axis_anno, 
      row_names_side = "left", row_names_gp = gpar(fontsize = 8, fontface = "bold"),
      height = unit(0.8, "cm"), 
      rect_gp = gpar(type = "none"), 
      layer_fun = lollipop_layer,
      show_heatmap_legend = FALSE
    )
  }
  if (is.null(ht_list)) ht_list <- ht else ht_list <- ht_list %v% ht
}
cairo_pdf("Figures_raw/22_01_Dynamic_Lollipop_Landscape.pdf", width = 240 / 25.4, height = 260 / 25.4)
draw(ht_list, 
     ht_gap = unit(1.0, "mm"),
     column_title = "Pan-Cohort Signature-PORCN Landscape", 
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# Comparisons------------------------------------------------------------------

library(dplyr)
library(tidyr)
calc_stats <- function(df, feat, group_col) {
  sub_df <- df %>% 
    filter(!is.na(.data[[feat]]), .data[[group_col]] %in% c("High", "Low"))
  
  if (nrow(sub_df) < 6 || length(unique(sub_df[[group_col]])) < 2) {
    return(NULL)
  }
  val_vec <- as.numeric(sub_df %>% pull(feat))
  clean_vals <- na.omit(val_vec)
  is_binary <- length(unique(clean_vals)) <= 2 && all(clean_vals %in% c(0, 1))
  if (is_binary) {
    # === Binary Fisher Exact Test ===
    tab <- table(sub_df[[group_col]], val_vec)
    if(nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
    p_val <- fisher.test(tab)$p.value
    prop_h <- mean(val_vec[sub_df[[group_col]] == "High"], na.rm = TRUE)
    prop_l <- mean(val_vec[sub_df[[group_col]] == "Low"], na.rm = TRUE)
    # Keep trend
    direction <- ifelse(prop_h > prop_l, "Up", ifelse(prop_h < prop_l, "Down", "NS"))
    if (is.na(p_val)) direction <- "NS"
    return(data.frame(
      N_High = sum(sub_df[[group_col]] == "High", na.rm = TRUE),
      N_Low = sum(sub_df[[group_col]] == "Low", na.rm = TRUE),
      Stat_High = prop_h,  
      Stat_Low = prop_l,   
      P_Value = p_val,
      Direction = direction,
      Stat_Type = "Binary_Proportion",
      stringsAsFactors = FALSE
    ))
  } else {
    # === Continuous Wilcoxon Rank-Sum Test ===
    val_h <- val_vec[sub_df[[group_col]] == "High"]
    val_l <- val_vec[sub_df[[group_col]] == "Low"]
    p_val <- suppressWarnings(wilcox.test(val_h, val_l)$p.value)
    med_h <- median(val_h, na.rm = TRUE)
    med_l <- median(val_l, na.rm = TRUE)
    # Keep trend
    direction <- ifelse(med_h > med_l, "Up", ifelse(med_h < med_l, "Down", "NS"))
    if (is.na(p_val)) direction <- "NS"
    return(data.frame(
      N_High = length(na.omit(val_h)), 
      N_Low = length(na.omit(val_l)),
      Stat_High = med_h, 
      Stat_Low = med_l,
      P_Value = p_val, 
      Direction = direction,
      Stat_Type = "Continuous_Median",
      stringsAsFactors = FALSE
    ))
  }
}
results_list <- list()
analysis_grid <- list(
  list(cohort = "OCM",  hrd = "Overall"),
  list(cohort = "OCM",  hrd = "HRD"),
  list(cohort = "OCM",  hrd = "HRP"),
  list(cohort = "TCGA", hrd = "Overall"),
  list(cohort = "TCGA", hrd = "HRD"),
  list(cohort = "TCGA", hrd = "HRP"),
  list(cohort = "CCLE", hrd = "Overall")
)
for (feat in all_features) {
  if (!(feat %in% colnames(master_df))) next
  for (grid in analysis_grid) {
    c_cohort <- grid$cohort
    c_hrd <- grid$hrd
    target_df <- master_df %>% filter(Cohort == c_cohort)
    if (c_hrd != "Overall") {
      target_df <- target_df %>% filter(HRD_Status == c_hrd)
    }
    # Signature
    res_sig <- calc_stats(target_df, feat, "Signature_Class")
    if (!is.null(res_sig)) {
      results_list[[length(results_list) + 1]] <- cbind(
        Feature = feat, Cohort = c_cohort, HRD_Subgroup = c_hrd,
        Grouping = "Signature_Class", res_sig
      )
    }
    # PORCN
    res_porcn <- calc_stats(target_df, feat, "PORCN_Class")
    if (!is.null(res_porcn)) {
      results_list[[length(results_list) + 1]] <- cbind(
        Feature = feat, Cohort = c_cohort, HRD_Subgroup = c_hrd,
        Grouping = "PORCN_Class", res_porcn
      )
    }
  }
}
df_stats_master <- bind_rows(results_list) %>%
  mutate(
    P_Value_Format = case_when(
      is.na(P_Value)  ~ "NA",
      P_Value < 0.001 ~ "<0.001***", 
      P_Value < 0.01  ~ sprintf("%.3f**", P_Value),
      P_Value < 0.05  ~ sprintf("%.3f*", P_Value),
      TRUE            ~ sprintf("%.3f", P_Value)
    ),
    Trend = case_when(
      Direction == "Up"   ~ paste0("↑", P_Value_Format),
      Direction == "Down" ~ paste0("↓", P_Value_Format),
      TRUE                ~ P_Value_Format
    )
  )

# Add survival data and COX----------------------------------------------------

library(survival)
library(survminer)
library(broom)
surv_tcga <- clinical_TCGA %>%
  select(
    Sample.ID, 
    OS_Time = OS_Time_TCGA, OS_Status = OS_Status_TCGA,
    DFS_Time = DFS_Time_TCGA, DFS_Status = DFS_Status_TCGA
  )
surv_ocm <- clinical_OCM %>%
  select(
    Sample.ID, 
    OS_Time = OS_days_from_diagnosis, OS_Status = OS_status,
    DFS_Time = OS_days_from_OCM, DFS_Status = OS_status #actually OS from OCM
  )
surv_combined <- bind_rows(surv_tcga, surv_ocm) %>%
  distinct(Sample.ID, .keep_all = TRUE)
master_df <- master_df %>%
  left_join(surv_combined, by = "Sample.ID")
run_cox_analysis <- function(df, time_col, status_col, var_col) {
  surv_df <- df %>% 
    select(Time = all_of(time_col), Status = all_of(status_col), Var = all_of(var_col)) %>%
    filter(!is.na(Time), !is.na(Status), !is.na(Var), Var %in% c("High", "Low"))
  if (nrow(surv_df) < 10) return(NULL)
  # 1. 强制转为数值，避免 factor 干扰
  surv_df$Status <- as.numeric(as.character(surv_df$Status))
  # 2. 🌟 修复核心 Bug：不再看 unique，而是看“是否有足够的死亡/进展事件”
  # 只要事件数 >= 3 个，Cox 回归就有意义 (就算 100% 是 1 也没问题)
  if (sum(surv_df$Status == 1, na.rm = TRUE) < 3) {
    return(NULL)
  }
  # 3. 🌟 新增防御：确保这个亚组内同时存在 High 和 Low 两个组，否则无法做对比
  if (length(unique(surv_df$Var)) < 2) {
    return(NULL)
  }
  surv_df$Var <- factor(surv_df$Var, levels = c("Low", "High"))
  form <- as.formula("Surv(Time, Status) ~ Var")
  fit <- tryCatch(coxph(form, data = surv_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "VarHigh") %>%
    mutate(
      N_Total = nrow(surv_df),
      N_Events = sum(surv_df$Status == 1, na.rm = TRUE),
      HR_Format = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_Value = p.value
    ) %>%
    select(N_Total, N_Events, HR = estimate, HR_Format, P_Value)
  return(res)
}
cox_results_list <- list()
endpoints <- list(
  list(time = "OS_Time", status = "OS_Status", name = "OS"),
  list(time = "DFS_Time", status = "DFS_Status", name = "DFS")
)
for (grid in analysis_grid) {
  c_cohort <- grid$cohort
  c_hrd <- grid$hrd
  if (c_cohort == "CCLE") next
  target_df <- master_df %>% filter(Cohort == c_cohort)
  if (c_hrd != "Overall") {
    target_df <- target_df %>% filter(HRD_Status == c_hrd)
  }
  for (ep in endpoints) {
    #Signature_Class
    res_sig <- run_cox_analysis(target_df, ep$time, ep$status, "Signature_Class")
    if (!is.null(res_sig)) {
      cox_results_list[[length(cox_results_list) + 1]] <- cbind(
        Cohort = c_cohort, HRD_Subgroup = c_hrd, Endpoint = ep$name, 
        Grouping = "Signature_Class", res_sig
      )
    }
    #PORCN_Class
    res_porcn <- run_cox_analysis(target_df, ep$time, ep$status, "PORCN_Class")
    if (!is.null(res_porcn)) {
      cox_results_list[[length(cox_results_list) + 1]] <- cbind(
        Cohort = c_cohort, HRD_Subgroup = c_hrd, Endpoint = ep$name, 
        Grouping = "PORCN_Class", res_porcn
      )
    }
  }
}
df_cox_master <- bind_rows(cox_results_list)
df_cox_master <- df_cox_master %>%
  mutate(
    # Format the P-value with significance stars
    P_Value_Format = case_when(
      P_Value < 0.001 ~ "<0.001***",
      P_Value < 0.01  ~ sprintf("%.3f**", P_Value),
      P_Value < 0.05  ~ sprintf("%.3f*", P_Value),
      TRUE            ~ sprintf("%.3f", P_Value)
    ),
    # Determine direction based on Hazard Ratio
    Direction_Arrow = case_when(
      HR > 1 ~ "↑",
      HR < 1 ~ "↓",
      TRUE   ~ ""
    ),
    # Combine everything into the final Trend string
    Trend = paste0(Direction_Arrow, HR_Format, ", ", P_Value_Format)
  )

# Keep RDS---------------------------------------------------------------------

save.image("22.RData")