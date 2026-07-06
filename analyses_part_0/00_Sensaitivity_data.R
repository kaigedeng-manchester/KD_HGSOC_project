# Set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) # this.path prints the path to the current script while dirname prints the folder in which it is in
setwd(paste0(myDir[-length(myDir)], collapse = '/')) # this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() # checking

# Recover work if necessary
# load("00.RData") 

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
  "Carboplatin" = "firebrick",  
  "Taxol"       = "navy"      
)

# Load font Arial
windowsFonts(Arial = windowsFont("Arial"))

# Baseline comparison----------------------------------------------------------

clinical_data <- read.csv("Data/OCM_HGS_clinical.csv", row.names = 1, na.strings = "unknown") # n = 101 
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
df_ocm <- clinical_data %>%
  rownames_to_column("Sample_ID") %>%
  select(Patient_ID, Sample_ID, Age, Stage, OS_days = OS_days_from_diagnosis, OS_status) %>%
  mutate(
    Cohort = "OCM",
    OS_status = as.numeric(OS_status)
  )
process_tcga_clinical <- function(df_subset, cohort_name) {
  df_subset %>%
    left_join(clinical_GDC, by = c("Case.ID" = "Master.Unique.Case.ID")) %>%
    mutate(
      Cohort = cohort_name,
      Patient_ID = Case.ID,        
      Sample_ID = Sample.ID,       
      Age = as.numeric(Age..GDC.clincal.label.),
      Stage = Stage..GDC.clinical.label., 
      OS_days = as.numeric(OS..days..GDC.clinical.label.),
      OS_status = ifelse(Death..GDC.clinical.label. == "1:Dead", 1, 0)
    ) %>%
    select(Patient_ID, Sample_ID, Age, Stage, OS_days, OS_status, Cohort) %>%
    distinct(Patient_ID, .keep_all = TRUE)
}
df_tcga_rna <- process_tcga_clinical(clinical_TCGA_RNAseq_HGS, "TCGA_RNAseq")
df_tcga_ma  <- process_tcga_clinical(clinical_TCGA_MA_HGS, "TCGA_Microarray")
df_all <- bind_rows(df_ocm, df_tcga_rna, df_tcga_ma) %>%
  mutate(
    Cohort = factor(Cohort, levels = c("OCM", "TCGA_RNAseq", "TCGA_Microarray")),
    Stage_Raw = toupper(as.character(Stage)),
    Stage_Clean = case_when(
      grepl("IV|4", Stage_Raw) ~ "Stage IV",
      grepl("III|3", Stage_Raw) ~ "Stage III",
      grepl("II|2", Stage_Raw) ~ "Stage II",
      grepl("I|1", Stage_Raw) ~ "Stage I",
      TRUE ~ "Unknown"
    ),
    Stage_Clean = factor(Stage_Clean, levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Unknown"))
  )
summary(df_all)
library(gtsummary)
baseline_table <- df_all %>%
  select(Cohort, Age, Stage_Clean, OS_status) %>%
  tbl_summary(
    by = Cohort,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age at Diagnosis (Years)",
      Stage_Clean ~ "FIGO Tumor Stage", 
      OS_status ~ "Mortality Status (Dead)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(all_categorical() ~ "fisher.test"),
    test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE),
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) %>%
  add_overall() %>%
  bold_p() %>%
  bold_labels() %>%
  modify_header(label = "**Clinical Characteristic**")
print(baseline_table)
library(survival)
library(survminer)
surv_fit <- survfit(Surv(OS_days, OS_status) ~ Cohort, data = df_all)
p_surv <- ggsurvplot(
  surv_fit,
  data = df_all,
  pval = TRUE,
  pval.size = 2.5,
  pval.coord = c(0, 0.1),
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.y.text = FALSE,
  risk.table.fontsize = 2.5,
  risk.table.height = 0.3,
  palette = c("firebrick", "steelblue", "navy"),
  legend.title = "Study Cohorts",
  legend.labs = c("OCM", "TCGA RNAseq", "TCGA Microarray"),
  xlab = "Overall Survival (Days)",
  ylab = "Survival Probability",
  ggtheme = theme_kd(legend_pos = "top"),
  tables.theme = theme_kd(legend_pos = "none")
)
cairo_pdf("Figures_raw/00_01_Cohort_OS_comparison.pdf", 
          width = 85 / 25.4, height = 100 / 25.4)
print(p_surv)
dev.off()

# Distribution of sensitivity data---------------------------------------------

# Import data
screening_data <- read.csv("Data/OCM_mix_screening.csv", row.names = 1) # n = 83 
clinical_data <- read.csv("Data/OCM_HGS_clinical.csv", row.names = 1, na.strings = "unknown") # n = 101 

# Only keep HGSOC sensitivity data
setdiff(rownames(clinical_data), rownames(screening_data)) # check is there any missing values caused by different labels
setdiff(rownames(screening_data), rownames(clinical_data))
merged_data <- merge(clinical_data, screening_data, by = "row.names", all = FALSE) # n = 67
merged_data <- column_to_rownames(merged_data, var = "Row.names")
merged_data$log10_carboplatin_IC50 <- log10(merged_data$carboplatin_IC50)
merged_data$log10_taxol_IC50 <- log10(merged_data$taxol_IC50)

# Clean data-factors and numbers
merged_data$Stage <- as.factor(merged_data$Stage)
merged_data$HRD <- as.factor(merged_data$HRD)
merged_data$Chemonaive <- as.factor(merged_data$Chemonaive)
merged_data$Treatment_exposure <- as.factor(merged_data$Treatment_exposure)
merged_data$Treatment_interval <- as.numeric(merged_data$Treatment_interval)

# Distribution of log10(IC50) & AUC
df_ic50 <- merged_data %>%
  select(log10_carboplatin_IC50, log10_taxol_IC50) %>%
  pivot_longer(cols = everything(), 
               names_to = "Drug", 
               values_to = "logIC50") %>%
  mutate(Drug = ifelse(grepl("carboplatin", Drug), "Carboplatin", "Taxol"))
df_auc <- merged_data %>%
  select(carboplatin_AUC, taxol_AUC) %>%
  pivot_longer(cols = everything(), 
               names_to = "Drug", 
               values_to = "AUC") %>%
  mutate(Drug = ifelse(grepl("carboplatin", Drug), "Carboplatin", "Taxol"))

# Plot the distribution
library(ggpubr)
library(patchwork)
p_scatter_ic50 <- ggplot(merged_data, aes(x = log10_carboplatin_IC50, y = log10_taxol_IC50)) +
  geom_point(shape = 21, size = 2.5, fill = "#F39B7F", color = "black", alpha = 0.8) +
  theme_kd() +
  theme(aspect.ratio = 1) +
  labs(
    title = "",
    x = expression(log[10]("Carboplatin IC"[50]~"("*mu*"M)")),
    y = expression(log[10]("Taxol IC"[50]~"("*mu*"M)"))
  )
p_scatter_auc <- ggplot(merged_data, aes(x = carboplatin_AUC, y = taxol_AUC)) +
  geom_point(shape = 21, size = 2.5, fill = "#4DBBD5", color = "black", alpha = 0.8) +
  theme_kd() +
  theme(aspect.ratio = 1) +
  labs(
    title = "",
    x = "Carboplatin AUC",
    y = "Taxol AUC"
  )
library(ggExtra)
p_ic50_marginal <- ggMarginal(p_scatter_ic50, type = "density", fill = "grey80", color = "black")
p_auc_marginal <- ggMarginal(p_scatter_auc, type = "density", fill = "grey80", color = "black")
p_combined_scatter <- wrap_elements(p_ic50_marginal) | wrap_elements(p_auc_marginal)
cairo_pdf("Figures_raw/00_02_Drug_Sensitivity_Scatter.pdf", 
          width = 170 / 25.4, height = 85 / 25.4)
print(p_combined_scatter)
dev.off()

# Clinical and survival relevance of sensitivity data--------------------------

library(survival)
library(broom)

# Factors
covariates <- c("Age", "Stage", 
                "Seqone_score", "HRD",
                "log10_carboplatin_IC50", "carboplatin_AUC",
                "log10_taxol_IC50", "taxol_AUC")

# Endpoints
endpoints <- c("OS_days_from_diagnosis", "OS_days_from_OCM")

# COX regression
function_cox_batch <- function(data, covars, time_var) {
  res_list <- list()
  for (var in covars) {
    form <- as.formula(paste("Surv(", time_var, ", OS_status) ~", var))
    fit <- coxph(form, data = data)
    tidy_fit <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    n_total <- fit$n
    if (is.factor(data[[var]])) {
      lvls <- levels(data[[var]])
      counts <- table(model.frame(fit)[[var]]) 
      header <- tibble(Variable = var, term = NA, N = as.character(n_total), 
                       HR = NA, LCI = NA, UCI = NA, P = NA)
      ref_n <- as.vector(counts[lvls[1]])
      ref <- tibble(Variable = paste0("  ", lvls[1], " (Ref)"), term = NA, 
                    N = ifelse(is.na(ref_n), "0", as.character(ref_n)), 
                    HR = 1, LCI = 1, UCI = 1, P = NA)
      comp_n <- as.vector(counts[lvls[-1]])
      res <- tibble(Variable = paste0("  ", lvls[-1]), term = tidy_fit$term,
                    N = as.character(comp_n),
                    HR = tidy_fit$estimate, LCI = tidy_fit$conf.low, UCI = tidy_fit$conf.high, P = tidy_fit$p.value)
      
      res_list[[var]] <- bind_rows(header, ref, res)
      
    } else {
      res <- tibble(Variable = var, term = tidy_fit$term, 
                    N = as.character(n_total),
                    HR = tidy_fit$estimate, LCI = tidy_fit$conf.low, UCI = tidy_fit$conf.high, P = tidy_fit$p.value)
      res_list[[var]] <- res
    }
  }
  final_df <- bind_rows(res_list) %>%
    mutate(
      HR_text = case_when(
        is.na(P) & is.na(term) & HR == 1 ~ "1.00 (Reference)",
        is.na(HR) ~ "",
        TRUE ~ sprintf("%.2f (%.2f-%.2f)", HR, LCI, UCI)
      ),
      P_text = case_when(
        is.na(P) ~ "",
        P < 0.001 ~ "<0.001 ***",
        P < 0.01 ~ sprintf("%.3f **", P),
        P < 0.05 ~ sprintf("%.3f *", P),
        TRUE ~ sprintf("%.3f", P)
      )
    )
  return(final_df)
}
cox_results_diag <- function_cox_batch(merged_data, covariates, "OS_days_from_diagnosis")
cox_results_ocm <- function_cox_batch(merged_data, covariates, "OS_days_from_OCM")

# Forest plot
library(grid)
library(forestploter)
space_placeholder <- paste(rep(" ", 20), collapse = " ")
plot_df <- data.frame(
  Variable = cox_results_diag$Variable,
  N = cox_results_diag$N,
  `OS (Diagnosis)` = space_placeholder,           
  `HR (95% CI) ` = cox_results_diag$HR_text, 
  `P Value ` = cox_results_diag$P_text,
  `OS (OCM)` = space_placeholder,                 
  `HR (95% CI)` = cox_results_ocm$HR_text,
  `P Value` = cox_results_ocm$P_text,
  check.names = FALSE
)
est_list   <- list(cox_results_diag$HR, cox_results_ocm$HR)
est_list   <- lapply(est_list, function(x) ifelse(x <= 0, 0.05, x))
lower_list <- list(cox_results_diag$LCI, cox_results_ocm$LCI)
lower_list <- lapply(lower_list, function(x) ifelse(x <= 0, 0.05, x))
upper_list <- list(cox_results_diag$UCI, cox_results_ocm$UCI)
upper_list <- lapply(upper_list, function(x) ifelse(x > 100, 100, x))
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
forest_plot <- forest(
  plot_df,
  est = est_list,
  lower = lower_list,
  upper = upper_list,
  ci_column = c(3, 6),                            
  ref_line = 1, 
  sizes = 0.4,
  x_trans = "log",
  xlim = list(c(0.1, 10), c(0.1, 10)),
  ticks_at = c(0.1, 0.5, 1, 2, 5, 10), 
  arrow_lab = c("Protective", "Risk"),
  theme = tm
) %>% 
  add_border(part = "header", where = "top", gp = gpar(lwd = 1)) %>%
  add_border(part = "header", where = "bottom", gp = gpar(lwd = 1)) %>%
  add_border(row = nrow(plot_df), where = "bottom", gp = gpar(lwd = 1)) %>% 
  edit_plot(
    part = "header",       
    col = c(3, 6), 
    hjust = 0.5,           
    x = unit(0.5, "npc")                
  )
cairo_pdf(file = "Figures_raw/00_02_Cox_regression.pdf", 
          width = 170 / 25.4, 
          height = 85 / 25.4) 
print(forest_plot)
dev.off()

# Keep RDS--------------------------------------------------------------------

save.image("00.RData")
saveRDS(merged_data, file = "Data/Clinical_merged_OCM_HGS.rds")
saveRDS(clinical_data, file = "Data/Clinical_OCM_HGS.rds")
