# set working directory to parent folder of the current script file
# myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) #this.path prints the path to the current script while dirname prints the folder in which it is in
# setwd(paste0(myDir[-length(myDir)], collapse = '/')) #this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() #just checking
setwd("C:/Temp processing") #for long path
#load("CCLE.RData") #load data
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggvenn)
library(fgsea)
library(limma)
library(edgeR)

# Sample check-----------------------------------------------------------------

models_OV <- read.csv("CCLE_OV_samples.csv", check.names = FALSE)
hgsoc_known_df <- models_OV %>%
  filter(`NMF subtype` == "HGSOC") %>%
  select(depmapId, cellLineDisplayName)
target_df <- models_OV %>%
  filter(`NMF subtype` == "not found") %>%
  select(depmapId, cellLineDisplayName)
hgsoc_known_ids <- hgsoc_known_df$depmapId
target_ids <- target_df$depmapId

# check availability of QC source data
library(data.table)
cnv_data <- fread("DepMap Public 25Q3 Files/OmicsCNSegmentsWGS.csv") %>% as.data.frame()
setnames(cnv_data, 1, "Number")
cnv_ids <- unique(cnv_data$ModelID) # N=1095
mut_data <- fread("DepMap Public 25Q3 Files/OmicsSomaticMutations.csv") %>% as.data.frame()
setnames(mut_data, 1, "Number")
mut_ids <- unique(mut_data$ModelID) # N=1955
sig_data <- fread("DepMap Public 25Q3 Files/OmicsGlobalSignatures.csv") %>% as.data.frame()
setnames(sig_data, 1, "Number")
sig_ids <- unique(sig_data$ModelID) # N=1955
check_availability <- function(group_name, id_list) {
  total_in_group <- length(id_list)
  # 
  has_cnv <- sum(id_list %in% cnv_ids)
  has_mut <- sum(id_list %in% mut_ids)
  has_sig <- sum(id_list %in% sig_ids)
  # 
  has_all_three <- sum((id_list %in% cnv_ids) & 
                         (id_list %in% mut_ids) & 
                         (id_list %in% sig_ids))
  # 
  cat(sprintf("Checking for [%s] (%d models in total) \n", group_name, total_in_group))
  cat(sprintf("(1) CNV_WGS available: %d (%.1f%%)\n", has_cnv, (has_cnv/total_in_group)*100))
  cat(sprintf("(2) Mut_WGS/WES available: %d (%.1f%%)\n", has_mut, (has_mut/total_in_group)*100))
  cat(sprintf("(3) Global_WGS/WES available: %d (%.1f%%)\n", has_sig, (has_sig/total_in_group)*100))
  cat(sprintf("(4) All data available: %d (%.1f%%)\n\n", has_all_three, (has_all_three/total_in_group)*100))
}
check_availability("hgsoc_known_ids", hgsoc_known_ids)
# Checking for [hgsoc_known_ids] (16 models in total) 
# (1) CNV_WGS available: 13 (81.2%)
# (2) Mut_WGS/WES available: 16 (100.0%)
# (3) Global_WGS/WES available: 16 (100.0%)
# (4) All data available: 13 (81.2%)
check_availability("target_ids", target_ids)
# Checking for [target_ids] (31 models in total) 
# (1) CNV_WGS available: 9 (29.0%)
# (2) Mut_WGS/WES available: 29 (93.5%)
# (3) Global_WGS/WES available: 29 (93.5%)
# (4) All data available: 9 (29.0%)

# HGSOC CIN + TP53 QC ---------------------------------------------------------

all_ids <- unique(c(hgsoc_known_ids, target_ids)) # 31+16
QC_sig <- sig_data %>% 
  filter(ModelID %in% all_ids) %>%            
  filter(!is.na(CIN)) %>%                    
  arrange(ModelID, desc(IsDefaultEntryForModel)) %>% 
  distinct(ModelID, .keep_all = TRUE)        # N=40
QC_tp53 <- mut_data %>% filter(ModelID %in% all_ids) %>%
  filter(HugoSymbol == "TP53") %>%
  group_by(ModelID) %>%
  summarise(TP53_Details = paste(unique(VariantInfo), collapse = "; ")) # N=33
QC_df <- QC_sig %>%
  left_join(QC_tp53, by = "ModelID") %>%
  # label TP53 WT
  mutate(TP53_Status = ifelse(is.na(TP53_Details), "WT", "Mutant")) %>%
  # mark reference list
  mutate(Group = case_when(
    ModelID %in% hgsoc_known_ids ~ "Known_HGSOC",
    ModelID %in% target_ids ~ "Target_Unknown",
    TRUE ~ "Other"
  )) %>%
  # Model name labels
  left_join(models_OV %>% select(depmapId, cellLineDisplayName), 
            by = c("ModelID" = "depmapId")) # N=26/29(31)+14/16
write_csv(QC_df, "HGSOC_TP53_Global_Signatures_QC.csv")
library(ggplot2)
library(ggrepel)
library(patchwork)
thresholds <- QC_df %>%
  filter(Group == "Known_HGSOC") %>%
  summarise(
    LoH_min = quantile(LoHFraction, 0.05, na.rm = TRUE),
    CIN_min = quantile(CIN, 0.05, na.rm = TRUE),
    Aneuploidy_min = quantile(Aneuploidy, 0.05, na.rm = TRUE)
  )
color_map <- c("Known_HGSOC" = "#CD534CFF", "Target_Unknown" = "#0073C2FF")
shape_map <- c("Mutant" = 16, "WT" = 17)
my_theme <- function(base_size = 7, base_family = "Arial", legend_pos = "bottom") {
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
plot_CIN <- ggplot(QC_df, aes(x = LoHFraction, y = CIN)) +
  geom_vline(xintercept = thresholds$LoH_min, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = thresholds$CIN_min, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Group, shape = TP53_Status, size = MSIScore), alpha = 0.8) +
  geom_text_repel(
    data = QC_df,
    aes(label = cellLineDisplayName), 
    size = 2.5, max.overlaps = 20, box.padding = 0.5
  ) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  scale_size_continuous(range = c(2, 8), breaks = c(1, 3, 5, 10)) + 
  labs(
    title = "",
    x = "LoH Fraction", y = "CIN Score"
  ) +
  my_theme()
plot_Aneuploidy <- ggplot(QC_df, aes(x = LoHFraction, y = Aneuploidy)) +
  geom_vline(xintercept = thresholds$LoH_min, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = thresholds$Aneuploidy_min, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Group, shape = TP53_Status, size = MSIScore), alpha = 0.8) +
  geom_text_repel(
    data = QC_df,
    aes(label = cellLineDisplayName), 
    size = 2.5, max.overlaps = 20, box.padding = 0.5
  ) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  scale_size_continuous(range = c(2, 8), breaks = c(1, 3, 5, 10)) +
  labs(
    title = "",
    x = "LoH Fraction", y = "Aneuploidy Score"
  ) +
  my_theme()
plot_QC <- plot_CIN + plot_Aneuploidy + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", aspect.ratio = 1)
cairo_pdf("QC_CCLE.pdf", 
          width = 170 / 25.4, 
          height = 100 / 25.4) 
print(plot_QC)
dev.off()
