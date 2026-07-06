# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("12.RData") 

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

clinical_OCM_HGS <- readRDS("Data/clinical_OCM_HGS.rds")
clinical_TCGA_MA_HGS <- readRDS("Data/clinical_TCGA_MA_HGS.rds")
sensitivity_OCM_HGS <- readRDS("Data/Clinical_merged_OCM_HGS.rds")
clinical_OCM_HGS_groups <- clinical_OCM_HGS %>%
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
clinical_TCGA_MA_HGS_groups <- clinical_TCGA_MA_HGS %>%
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

# Screening results------------------------------------------------------------

library(ggpubr)
library(tibble)
df_drug <- clinical_OCM_HGS_groups %>%
  left_join(
    sensitivity_OCM_HGS %>% rownames_to_column("Sample.ID"), 
    by = "Sample.ID"
  )
color_groups <- c(
  "HRP_PORCN_Low"  = "#8BA3C1",  
  "HRP_PORCN_High" = "#D39695",  
  "HRD_PORCN_Low"  = "#2166AC",  
  "HRD_PORCN_High" = "#B2182B"
)
my_comparisons <- list(
  c("HRD_PORCN_High", "HRD_PORCN_Low"), 
  c("HRP_PORCN_High", "HRP_PORCN_Low"),
  c("HRD_PORCN_High", "HRP_PORCN_High"), # the HRD-high group, comparable results with HRP
  c("HRD_PORCN_High", "HRP_PORCN_Low")
)
function_drug_auc <- function(df, drug_col, drug_name) {
  df_clean <- df %>% filter(!is.na(.data[[drug_col]]))
  p <- ggplot(df_clean, aes(x = Subgroup, y = .data[[drug_col]])) +
    geom_violin(aes(fill = Subgroup), trim = FALSE, alpha = 0.5, color = NA) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(aes(color = Subgroup), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = color_groups) +
    scale_color_manual(values = color_groups) +
    scale_x_discrete(labels = c(
      "HRP_PORCN_Low"  = "HRP\nPORCN Low",
      "HRP_PORCN_High" = "HRP\nPORCN High",
      "HRD_PORCN_Low"  = "HRD\nPORCN Low",
      "HRD_PORCN_High" = "HRD\nPORCN High"
    )) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                       label = "p.format", tip.length = 0.02, size = 2.5) +
    stat_compare_means(label.y.npc = "top", fontface = "bold", size = 2.5) +
    theme_kd(legend_pos = "none") +
    theme(
      axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black")
    ) +
    labs(
      title = paste(drug_name, " Sensitivity"),
      x = "Subgroups",
      y = paste(drug_name, " AUC")
    )
  return(p)
}
p_carbo <- function_drug_auc(df_drug, "carboplatin_AUC", "Carboplatin")
p_taxol <- function_drug_auc(df_drug, "taxol_AUC", "Taxol")
library(patchwork)
p_sensitivity <- p_carbo + p_taxol +
  plot_layout(guides = "collect") & 
  theme(
    aspect.ratio = 1
  )
cairo_pdf("Figures_raw/12_01_group_chemo.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4) 
print(p_sensitivity)
dev.off()

# Survival---------------------------------------------------------------------

library(survival)
library(survminer)
function_four_panel_surv <- function(df) {
  df_hrp <- df %>% filter(HRD_Status == "HRD-") 
  df_hrd <- df %>% filter(HRD_Status == "HRD+")
  make_plot <- function(data, time_col, status_col, sub_title) {
    form <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ PORCN_Group"))
    fit <- surv_fit(form, data = data)
    pal <- c("navy", "firebrick")
    ggsurvplot(
      fit, data = data, pval = TRUE, pval.size = 2.5, # 7.11pt as basic
      conf.int = FALSE,
      risk.table = TRUE, risk.table.fontsize = 2.5, # 7.11pt as basic
      palette = pal, legend.title = "PORCN Tertiles", legend.labs = c("Low", "High"),
      xlab = "Time (days)", ylab = "Survival Probability",
      ggtheme = theme_kd(), tables.theme = theme_kd(), title = sub_title
    )
  }
  p_os_hrp <- make_plot(df_hrp, "OS_Time_TCGA", "OS_Status_TCGA", "TCGA OS vs. PORCN (HRP)")
  p_os_hrd <- make_plot(df_hrd, "OS_Time_TCGA", "OS_Status_TCGA", "TCGA OS vs. PORCN (HRD)")
  p_dfs_hrp <- make_plot(df_hrp, "DFS_Time_TCGA", "DFS_Status_TCGA", "TCGA DFS vs. PORCN (HRP)")
  p_dfs_hrd <- make_plot(df_hrd, "DFS_Time_TCGA", "DFS_Status_TCGA", "TCGA DFS vs. PORCN (HRD)")
  arrange_ggsurvplots(
    list(p_os_hrp, p_os_hrd, p_dfs_hrp, p_dfs_hrd), 
    print = TRUE, 
    ncol = 2, 
    nrow = 2
  )
}
cairo_pdf("Figures_raw/12_02_group_survival_TCGA.pdf", 
          width = 170 / 25.4, 
          height = 200 / 25.4)
function_four_panel_surv(clinical_TCGA_MA_HGS_groups)
dev.off()

# Survival-pairwise------------------------------------------------------------

cairo_pdf("Figures_raw/12_03_group4_survival_TCGA.pdf", 
          width = 120 / 25.4, 
          height = 110 / 25.4)
ggsurvplot(surv_fit(Surv(OS_Time_TCGA, OS_Status_TCGA) ~ Subgroup, 
                    data = clinical_TCGA_MA_HGS_groups),
           data = clinical_TCGA_MA_HGS_groups,
           pval = TRUE, 
           pval.size = 2.5, # 7.11pt as basic
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.fontsize = 2.5, # 7.11pt as basic
           palette = c(
             "HRP_PORCN_Low"  = "#8BA3C1",  
             "HRP_PORCN_High" = "#D39695",  
             "HRD_PORCN_Low"  = "#2166AC",  
             "HRD_PORCN_High" = "#B2182B"  
           ),
           legend.title = "Subgroups",
           legend.labs = levels(clinical_TCGA_MA_HGS_groups$Subgroup),
           xlab = "Time (days)",
           ylab = "Survival Probability",
           ggtheme = theme_kd(),
           title = "TCGA OS Across Subgroups"
) %>% print()
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("12.RData")
saveRDS(clinical_OCM_HGS_groups, file = "Data/clinical_OCM_HGS_groups.rds")
saveRDS(clinical_TCGA_MA_HGS_groups, file = "Data/clinical_TCGA_MA_HGS_groups.rds")