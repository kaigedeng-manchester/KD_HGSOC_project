# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("11.RData") 

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

# Get PORCN groups------------------------------------------------------------

clinical_OCM_HGS <- readRDS("Data/clinical_OCM_HGS.rds")
clinical_OCM_HGS_treatment <- clinical_OCM_HGS %>%
  filter(!is.na(PORCN_Expr)) %>%
  mutate(
    Treatment_interval = suppressWarnings(as.numeric(ifelse(Treatment_interval == "naive", NA, Treatment_interval))),
    Treatment_exposure = factor(
      Treatment_exposure, 
      levels = c("naive", ">21d", "<=21d")
    )
  ) %>%
  filter(!is.na(Treatment_exposure))

# Treatment exposure vs. PORCN-------------------------------------------------

my_comparisons_exposure <- list(
  c("naive", ">21d"),
  c("naive", "<=21d"),
  c(">21d", "<=21d")
)
library(ggpubr)
p_exposure <- ggplot(clinical_OCM_HGS_treatment, aes(x = Treatment_exposure, y = PORCN_Expr)) +
  geom_violin(aes(fill = Treatment_exposure), trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = Treatment_exposure), width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("gray60", "#2166AC", "#B2182B")) +
  scale_color_manual(values = c("gray60", "#2166AC", "#B2182B")) +
  stat_compare_means(comparisons = my_comparisons_exposure, method = "wilcox.test", 
                     label = "p.format", tip.length = 0.02, size = 2.5) +
  stat_compare_means(method = "kruskal.test", 
                     label.y = max(clinical_OCM_HGS_treatment$PORCN_Expr, na.rm = TRUE) * 1.15, 
                     color = "black", fontface = "bold", size = 2.5) +
  scale_x_discrete(labels = c(
    "naive" = "Treatment-naïve",
    ">21d"  = "Post-chemo\n(> 21 days)",
    "<=21d" = "Recent chemo\n(\u2264 21 days)"
  )) +
  theme_kd() +
  theme(axis.text.x = element_text(face = "bold", color = "black")) +
  labs(
    title = "PORCN Expression across Treatment Exposures",
    x = "Treatment Exposure Status",
    y = "PORCN Expression Level"
  )
df_ocm_treated_only <- clinical_OCM_HGS_treatment %>% 
  filter(!is.na(Treatment_interval))
p_interval <- ggplot(df_ocm_treated_only, aes(x = Treatment_interval, y = PORCN_Expr, color = HRD_Status)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(aes(fill = HRD_Status), method = "lm", linetype = "dashed", linewidth = 0.5, alpha = 0.15) +
  stat_cor(method = "spearman", 
           label.x = max(df_ocm_treated_only$Treatment_interval, na.rm = TRUE) * 0.3, 
           label.y.npc = "top", 
           fontface = "bold", show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = NA, linetype = "dotted", linewidth = 0.5, inherit.aes = FALSE, aes(x = Treatment_interval, y = PORCN_Expr)) +
  stat_cor(inherit.aes = FALSE, aes(x = Treatment_interval, y = PORCN_Expr),
           method = "spearman", color = "black", 
           label.x = max(df_ocm_treated_only$Treatment_interval, na.rm = TRUE) * 0.3, 
           label.y.npc = 0.6, 
           fontface = "bold", size = 2.5) +
  scale_color_manual(values = c("HRD+" = "firebrick", "HRD-" = "navy"), na.value = "gray50", labels = c("HRD+" = "HRD", "HRD-" = "HRP")) +
  scale_fill_manual(values = c("HRD+" = "firebrick", "HRD-" = "navy"), na.value = "gray50", labels = c("HRD+" = "HRD", "HRD-" = "HRP")) +
  theme_kd() +
  labs(
    title = "Chemo-Free Interval vs. PORCN by HRD Status",
    x = "Treatment Interval (Days since last chemotherapy)",
    y = "PORCN Expression Level",
    color = "HRD Status",
    fill = "HRD Status"
  )
library(patchwork)
p_treatment <- p_exposure + p_interval +
  theme(aspect.ratio = 1)
cairo_pdf("Figures_raw/11_01_treatment_effect.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4) 
print(p_treatment)
dev.off()

# Longitudinal samples vs. PORCN-----------------------------------------------

df_longitudinal <- clinical_OCM_HGS %>%
  group_by(Patient_ID) %>%
  mutate(Sample_Count = n()) %>%
  ungroup() %>%
  filter(Sample_Count > 1) %>%
  mutate(
    Treatment_exposure = as.character(Treatment_exposure),
    Treatment_exposure = ifelse(is.na(Treatment_exposure), "Unknown", Treatment_exposure),
    Treatment_exposure = factor(Treatment_exposure, levels = c("naive", "<=21d", ">21d", "Unknown"))
  ) %>%
  arrange(Patient_ID, Days_to_OCM)
consistency_stat <- df_longitudinal %>%
  group_by(Patient_ID) %>%
  summarise(
    Subtypes_Seen = list(unique(as.character(PORCN_Group))),
    Unique_Subtype_Count = length(unique(as.character(PORCN_Group))),
    Is_Stable = Unique_Subtype_Count == 1,
    .groups = "drop"
  )
stable_rate <- sum(consistency_stat$Is_Stable) / nrow(consistency_stat) * 100 # 56.25

# Trajectory
df_longitudinal_ordinal <- df_longitudinal %>%
  group_by(Patient_ID) %>%
  arrange(Days_to_OCM, .by_group = TRUE) %>%
  mutate(
    Timepoint = row_number(),
    Timepoint_Label = factor(paste0("T", Timepoint), levels = paste0("T", sort(unique(Timepoint))))
  ) %>%
  ungroup()
p_facet_lanes <- ggplot(df_longitudinal_ordinal, aes(x = Timepoint_Label, y = PORCN_Expr, group = Patient_ID)) +
  geom_line(color = "gray50", linewidth = 1) + # for trajectory, use 1
  geom_point(aes(fill = PORCN_Group, shape = Treatment_exposure), 
             size = 4, color = "black", stroke = 0.8) +
  facet_wrap(~ Patient_ID, ncol = 4) + 
  scale_fill_manual(values = c("Low" = "navy", "Intermediate" = "gray80", "High" = "firebrick")) +
  scale_shape_manual(values = c("naive" = 21, "<=21d" = 24, ">21d" = 22, "Unknown" = 23)) +
  theme_kd() +
  labs(
    title = "PORCN Expression Trajectories by Patient",
    x = "Sampling Sequence",
    y = "PORCN Expression Level",
    fill = "PORCN Subtype",
    shape = "Treatment Exposure"
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)))
cairo_pdf("Figures_raw/11_02_patient_trajectory.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4) 
print(p_facet_lanes)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("11.RData")