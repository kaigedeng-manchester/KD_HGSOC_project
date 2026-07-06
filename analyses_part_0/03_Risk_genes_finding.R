# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("03.RData") 

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

# DEG list--------------------------------------------------------------------

degs_carbo <- readRDS("Data/degs_carbo.rds")
degs_taxol <- readRDS("Data/degs_taxol.rds")
sig_genes_union <- union(degs_carbo %>% filter(colour != "not_significant") %>% pull(GeneSymbol),
                         degs_taxol %>% filter(colour != "not_significant") %>% pull(GeneSymbol))
merged_degs <- full_join(
  degs_carbo %>% filter(GeneSymbol %in% sig_genes_union) %>% select(GeneSymbol, logFC_Carbo = logFC, colour_Carbo = colour),
  degs_taxol %>% filter(GeneSymbol %in% sig_genes_union) %>% select(GeneSymbol, logFC_Taxol = logFC, colour_Taxol = colour),
  by = "GeneSymbol"
) %>%
  mutate(
    logFC_Carbo = ifelse(is.na(logFC_Carbo), 0, logFC_Carbo),
    logFC_Taxol = ifelse(is.na(logFC_Taxol), 0, logFC_Taxol),
    colour_Carbo = ifelse(is.na(colour_Carbo), "not_significant", colour_Carbo),
    colour_Taxol = ifelse(is.na(colour_Taxol), "not_significant", colour_Taxol)
  )
merged_degs <- merged_degs %>%
  mutate(
    Status = case_when(
      (colour_Carbo == "up" & colour_Taxol == "down") | (colour_Carbo == "down" & colour_Taxol == "up") ~ "Contradictory",
      (colour_Carbo == "up" & colour_Taxol == "up") ~ "Dual_Resistant",
      (colour_Carbo == "down" & colour_Taxol == "down") ~ "Dual_Sensitive",
      (colour_Carbo == "up" | colour_Taxol == "up") ~ "Single_Resistant",
      (colour_Carbo == "down" | colour_Taxol == "down") ~ "Single_Sensitive"
    )
  )

# Quadrant plot
porcn_data <- merged_degs %>% filter(GeneSymbol == "PORCN")
merged_degs <- merged_degs %>%
  mutate(Status = factor(Status, levels = c(
    "Single_Sensitive",
    "Single_Resistant",
    "Dual_Sensitive",
    "Dual_Resistant",
    "Contradictory"
  )))
library(ggrepel)
p_quadrant <- ggplot(merged_degs, aes(x = logFC_Carbo, y = logFC_Taxol, color = Status)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  scale_color_manual(values = c(
    "Single_Resistant" = "#F4A582", 
    "Single_Sensitive" = "#92C5DE",
    "Dual_Resistant" = "#B2182B",
    "Dual_Sensitive" = "#2166AC",  
    "Contradictory" = "#E66101"                    
  )) +
  annotate("text", x = -2, y = 2, label = "Contradictory", 
           color = "black", size = 3, fontface = "italic", lineheight = 0.9) +
  annotate("text", x = 2, y = -2, label = "Contradictory", 
           color = "black", size = 3, fontface = "italic", lineheight = 0.9) +
  geom_text_repel(data = porcn_data, aes(label = GeneSymbol), 
                  fontface = "bold", color = "black", size = 2.5,
                  max.overlaps = 20, box.padding = 0.5, point.padding = 0.5, min.segment.length = 0) +
  coord_fixed(ratio = 1) + 
  labs(
    title = "Carboplatin vs. Taxol DEGs",
    x = expression("Carboplatin log"[2]*"(Fold Change)"),
    y = expression("Taxol log"[2]*"(Fold Change)"),
    color = "Gene Shift Status"
  ) +
  theme_kd()

# Prognostic value of DEG------------------------------------------------------

# Preparing expression matrix
reads_OCM_HGS <- readRDS("Data/reads_OCM_HGS.rds")

# Filter low expression (keep genes expressed in >=20% samples with CPM > 1)
library(limma)
library(edgeR)
dge <- DGEList(counts = reads_OCM_HGS)
keep <- rowSums(cpm(dge) > 1) >= (0.2 * ncol(reads_OCM_HGS))
dge <- dge[keep, , keep.lib.sizes = FALSE] # keep 16610 genes out of 19931 total genes

# Calculate Normalization Factors (TMM)
dge <- calcNormFactors(dge, method = "TMM")

# Convert to log2-CPM (Continuous values needed for correction)
cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE)
expr <- log2(cpm_mat + 1) #log2(cpm+1) transformation

# Output
expr_OCM_HGS <- expr
rm(expr, cpm_mat, dge, keep)

# Cox regression
library(survival)
library(broom)
clinical_OCM_HGS <- readRDS("Data/Clinical_OCM_HGS.rds")
function_cox_batch <- function(
    expr_mat = expr_OCM_HGS,
    merged_df = merged_degs,
    clin_df = clinical_OCM_HGS, 
    time_col = "OS_days_from_diagnosis", 
    status_col = "OS_status"
) {
  target_genes <- intersect(merged_df$GeneSymbol, rownames(expr_mat))
  samples_to_use <- intersect(colnames(expr_mat), rownames(clin_df))
  expr_sub <- expr_mat[target_genes, samples_to_use, drop = FALSE]
  clin_sub <- clin_df[samples_to_use, , drop = FALSE]
  cat("Calculating Cox regression for", length(target_genes), "target degs...\n")
  res_list <- lapply(target_genes, function(gene) {
    tmp_data <- data.frame(
      Time = clin_sub[[time_col]],
      Status = clin_sub[[status_col]],
      Expression = as.numeric(expr_sub[gene, ])
    )
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
        adj.P.Val < 0.10 & HR > 1 ~ "Risk",
        adj.P.Val < 0.10 & HR < 1 ~ "Protective",
        TRUE ~ "NS"
      )
    ) %>%
    relocate(GeneSymbol, HR, log2HR, P.Value, adj.P.Val) %>%
    arrange(adj.P.Val)
  cat("Done! Processed", nrow(final_res), "genes.\n")
  return(final_res)
}
cox_degs <- function_cox_batch()

# Translational analysis
trans_data <- inner_join(merged_degs, cox_degs, by = "GeneSymbol") %>%
  mutate(
    Mean_logFC = (logFC_Carbo + logFC_Taxol) / 2,
    Translational_Status = case_when(
      Mean_logFC > 0 & log2HR > 0 & adj.P.Val < 0.10 ~ "Resistant & Risk",
      Mean_logFC < 0 & log2HR < 0 & adj.P.Val < 0.10 ~ "Sensitive & Protective",
      Mean_logFC > 0 & log2HR < 0 & adj.P.Val < 0.10 |
      Mean_logFC < 0 & log2HR > 0 & adj.P.Val < 0.10 ~ "Contradictory",
      TRUE ~ "Not Significant"
    )
  )
trans_data <- trans_data %>%
  mutate(Translational_Status = factor(Translational_Status, levels = c(
    "Resistant & Risk",
    "Sensitive & Protective",
    "Contradictory",
    "Not Significant"
  )))

# Translational quadrant plot
porcn_trans <- trans_data %>% filter(GeneSymbol == "PORCN")
limit_x <- max(abs(trans_data$log2HR), na.rm = TRUE) * 1.1
limit_y <- max(abs(trans_data$Mean_logFC), na.rm = TRUE) * 1.1
p_trans_quadrant <- ggplot(trans_data, aes(x = log2HR, y = Mean_logFC, color = Translational_Status)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_continuous(limits = c(-limit_x, limit_x)) +
  scale_y_continuous(limits = c(-limit_y, limit_y)) +
  annotate("text", x = limit_x * 0.6, y = limit_y * 0.8, 
           label = "Risk", color = "#B2182B", fontface = "bold", size = 3) +
  annotate("text", x = -limit_x * 0.6, y = -limit_y * 0.8, 
           label = "Protective", color = "#2166AC", fontface = "bold", size = 3) +
  annotate("text", x = -limit_x * 0.6, y = limit_y * 0.8, 
           label = "Contradictory", color = "black", fontface = "italic", size = 3) +
  annotate("text", x = limit_x * 0.6, y = -limit_y * 0.8, 
           label = "Contradictory", color = "black", fontface = "italic", size = 3) +
  scale_color_manual(values = c(
    "Not Significant" = "grey85",
    "Contradictory" = "#E66101",                
    "Resistant & Risk" = "#B2182B",         
    "Sensitive & Protective" = "#2166AC"
  )) +
  geom_text_repel(data = porcn_trans, aes(label = GeneSymbol), 
                  fontface = "bold", color = "black", size = 2.5,
                  max.overlaps = 20, box.padding = 0.5, point.padding = 0.5, min.segment.length = 0) +
  labs(
    title = "DEGs: Sensitivity vs. Survival Relavance",
    x = expression("Clinical Risk: Cox "*log[2]*"(Hazard Ratio)"),
    y = expression("In Vitro Response: Mean "*log[2]*"(Fold Change)"),
    color = "Translational Profile"
  ) +
  theme_kd()

# Combined quadrant plots
library(patchwork)
combined_quadrant <- p_quadrant + p_trans_quadrant + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",      
    legend.box = "horizontal",       
    legend.spacing.x = unit(1, "cm"),
    aspect.ratio = 1
  ) &
  guides(color = guide_legend(
    title.position = "top",          
    title.hjust = 0.5,               
    nrow = 2,                        
    byrow = TRUE                     
  ))
cairo_pdf("Figures_raw/03_01_DEGS_with_survival.pdf", 
          width = 170 / 25.4, 
          height = 120 / 25.4) 
print(combined_quadrant)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("03.RData")
saveRDS(expr_OCM_HGS, file = "Data/expr_OCM_HGS.rds")
saveRDS(trans_data, file = "Data/trans_degs.rds")