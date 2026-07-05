# set working directory to parent folder of the current script file
# myDir <- unlist(strsplit(dirname(this.path::this.path()), '/')) #this.path prints the path to the current script while dirname prints the folder in which it is in
# setwd(paste0(myDir[-length(myDir)], collapse = '/')) #this removes the last part of the path (i.e. the /code part of the path and rebuilds the path to be the parent folder)
getwd() #just checking
setwd("C:/Temp processing") #for long path
#load("TCGA.RData") #load data
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

# Sample overlap check---------------------------------------------------------

sample_TCGA_OV_RNAseq <- read_tsv("GDC_TCGA-OV_RNA-seq_sample_sheet.2026-02-12.tsv", show_col_types = FALSE)
sample_TCGA_OV_MA <- read_tsv("GDC_TCGA-OV_MA_sample_sheet.2026-02-12.tsv", show_col_types = FALSE)

# filter
sample_TCGA_OV_RNAseq %>% 
  filter(`Data Category` == "Transcriptome Profiling") %>%
  filter(`Tissue Type` == "Tumor") -> sample_TCGA_OV_RNAseq
sample_TCGA_OV_MA %>% 
  filter(`Data Category` == "Transcriptome Profiling") %>%
  filter(`Tissue Type` == "Tumor") -> sample_TCGA_OV_MA

# check cases and samples
venn_list_samples <- list(
  "TCGA RNAseq" = sample_TCGA_OV_RNAseq$`Sample ID` %>% unique(),
  "TCGA Array"  = sample_TCGA_OV_MA$`Sample ID` %>% unique()
)
venn_list_cases <- list(
  "TCGA RNAseq" = sample_TCGA_OV_RNAseq$`Case ID` %>% unique(),
  "TCGA Array"  = sample_TCGA_OV_MA$`Case ID` %>% unique()
)

# Overlap
ggvenn(
  venn_list_samples, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, 
  set_name_size = 4,
  show_percentage = TRUE) +
  ggtitle("Overlap of Sample IDs") +
  theme(plot.title = element_text(hjust = 0.5))
ggvenn(
  venn_list_cases, 
  fill_color = c("#0073C2FF", "#CD534CFF"), 
  stroke_size = 0.5, 
  set_name_size = 4,
  show_percentage = TRUE) +
  ggtitle("Overlap of Case IDs") +
  theme(plot.title = element_text(hjust = 0.5))

# Dataset 1: TCGA RNAseq-------------------------------------------------------

# Batch retrieval
data_dir <- paste0(getwd(), "/GDC_TCGA-OV_RNA-seq_data.2026-02-12") 
read_RNAseq_counts <- function(file_id, file_name, sample_id) {
  # Construct the full path.
  # GDC nests files inside a folder named after the File ID.
  # If your files are all in one big folder, remove 'file_id' from this path.
  file_path <- file.path(data_dir, file_id, file_name)
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  # Read the TSV file
  # TCGA STAR files usually have a header row and then data
  # 'skip = 1' might be needed depending on the exact version, but usually
  # 'comment = #' handles header metadata if present.
  counts_data <- read_tsv(file_path, show_col_types = FALSE, skip = 1)
  # Clean and Select Data
  # TCGA STAR files have columns: gene_id, gene_name, gene_type, unstranded, etc.
  # We usually want 'unstranded' (raw counts) or 'tpm_unstranded'.
  # We also filter out the "N_" summary stats rows (N_unmapped, etc.) at the top.
  clean_counts <- counts_data %>%
    filter(!grepl("^N_", gene_id)) %>%         #remove summary stats rows
    select(gene_id, gene_name, gene_type, unstranded) %>% #keep ID, Name, and Raw Counts
    rename(!!sample_id := unstranded)          #rename counts column to the Sample ID
  return(clean_counts)
}
matrix_TCGA_RNAseq <- pmap(list(sample_TCGA_OV_RNAseq$`File ID`,
                                sample_TCGA_OV_RNAseq$`File Name`,
                                sample_TCGA_OV_RNAseq$`Sample ID`),
                           read_RNAseq_counts) %>% 
  reduce(inner_join, by = c("gene_id", "gene_name", "gene_type"))

# Duplicates: sum
matrix_TCGA_RNAseq_sum <- matrix_TCGA_RNAseq %>%
  filter(gene_type == "protein_coding") %>% #Keep coding genes
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum))
colnames(matrix_TCGA_RNAseq_sum)[1] <- "GeneSymbol"

# Export
write_csv(matrix_TCGA_RNAseq, "output/matrix_TCGA_RNAseq.csv")
write_csv(matrix_TCGA_RNAseq_sum, "output/matrix_TCGA_RNAseq_sum.csv") #no duplicate

# Dataset 2: TCGA MA-----------------------------------------------------------

library(affy)

# Batch retrieval
data_dir <- paste0(getwd(), "/GDC_TCGA-OV_MA_data.2026-02-12")
MA_paths <- file.path(data_dir, 
                      sample_TCGA_OV_MA$`File ID`, 
                      sample_TCGA_OV_MA$`File Name`) #get exact file paths
exists_paths <- file.exists(MA_paths)
if (all(exists_paths)) {
  cat("Perfect! All", length(MA_paths), ".CEL files found\n")
} else {
  cat("Caution! There are", sum(!exists_paths), ".CEL files not found\n")
}
print(head(MA_paths[!exists_paths])) #check

# get raw CEL data
raw_MA_data <- ReadAffy(filenames = MA_paths)

# check chip platform
cdfName(raw_MA_data) # "HT_HG-U133A"

# Normalisation: Robust Multi-Array Average expression measure
  # -Background correcting
  # -Normalizing
  # -Calculating Expression
matrix_TCGA_MA_Probe <- raw_MA_data %>% rma() %>% exprs()

# Probe to Symbol annotation
library(hgu133a.db)
probes <- rownames(matrix_TCGA_MA_Probe)
symbols <- mapIds(hgu133a.db,
                  keys = probes,
                  column = "SYMBOL",
                  keytype = "PROBEID",
                  multiVals = "first")
matrix_TCGA_MA_symbol <- matrix_TCGA_MA_Probe %>% as.data.frame()
matrix_TCGA_MA_symbol$GeneSymbol <- symbols
matrix_TCGA_MA_symbol <- matrix_TCGA_MA_symbol %>% filter(!is.na(GeneSymbol)) #remove unmatched probes

# Duplicates: mean
matrix_TCGA_MA_mean <- matrix_TCGA_MA_symbol %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), mean))

# match sample ID
id_map <- setNames(sample_TCGA_OV_MA$`Sample ID`, sample_TCGA_OV_MA$`File Name`)
colnames(matrix_TCGA_MA_mean)[-1] <- id_map[colnames(matrix_TCGA_MA_mean)[-1]]

# Export
write_csv(matrix_TCGA_MA_mean, "output/matrix_TCGA_MA_mean.csv")

# HGSOC CNV QC-----------------------------------------------------------------

# (1) Read files---------------------------------------------------------------

sample_TCGA_OV_CNV <- read_tsv("GDC_TCGA-OV_masked_copy_number_segment_sample_sheet.2026-03-04.tsv", show_col_types = FALSE)

# Batch retrieval
data_dir <- paste0(getwd(), "/[DATA] GDC_TCGA-OV_masked_copy_number_segment_data.2026-03-04") 
read_cnv <- function(file_id, file_name, sample_id) {
  # Construct the full path.
  # GDC nests files inside a folder named after the File ID.
  # If your files are all in one big folder, remove 'file_id' from this path.
  file_path <- file.path(data_dir, file_id, file_name)
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  # Read the TXT file
  cnv_data <- read_tsv(file_path, show_col_types = FALSE)
  # Clean and Select Data
  # masked copy number segment files have columns: Chromosome, Start, End, Num_Probes, Segment_Mean
  clean_cnv <- cnv_data %>%
    select(Chromosome, Start, End, Num_Probes, Segment_Mean) %>% 
    # add sample label for later group_by analyses
    mutate(Sample_ID = sample_id,
           File_Name = file_name)
  return(clean_cnv)
}
master_cnv_segments <- pmap_dfr(list(sample_TCGA_OV_CNV$`File ID`,
                                     sample_TCGA_OV_CNV$`File Name`,
                                     sample_TCGA_OV_CNV$`Sample ID`),
                                read_cnv)
nrow(master_cnv_segments) #check
head(master_cnv_segments) #check

# Duplicate files: best file
best_files_index <- master_cnv_segments %>%
  group_by(Sample_ID, File_Name) %>%
  # 
  summarise(Total_Probes = sum(Num_Probes, na.rm = TRUE), .groups = "drop") %>%
  # 
  arrange(Sample_ID, desc(Total_Probes)) %>%
  # 
  distinct(Sample_ID, .keep_all = TRUE)
master_cnv_clean <- master_cnv_segments %>%
  semi_join(best_files_index, by = c("Sample_ID", "File_Name"))

# (2) Check FGA----------------------------------------------------------------

threshold <- 0.2 

# FGA calculation
cnv_qc_summary <- master_cnv_clean %>%
  mutate(Segment_Length = End - Start,
         Is_Altered = abs(Segment_Mean) > threshold) %>%
  group_by(Sample_ID) %>%
  summarise(
    Total_Genome_Length = sum(Segment_Length, na.rm = TRUE),
    Altered_Length = sum(Segment_Length[Is_Altered == TRUE], na.rm = TRUE),
    FGA = Altered_Length / Total_Genome_Length
  ) %>%
  arrange(FGA) # sort FGA of samples

# Export
write_csv(cnv_qc_summary, "output/cnv_qc_summary.csv")

# (3) Visualise genome CNV-----------------------------------------------------

library(scales)
ordered_samples <- cnv_qc_summary$Sample_ID
plot_cnv_segments <- master_cnv_segments %>%
  filter(Chromosome %in% c(as.character(1:22),"X")) %>%
  mutate(Chromosome = factor(Chromosome, levels = c(as.character(1:22),"X"))) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = ordered_samples)) %>%
  mutate(Y_index = as.numeric(Sample_ID))

# Annotation left
library(patchwork)
excluded_not_hgsoc <- c(
  "TCGA-20-0996-01A", "TCGA-59-2349-01A", "TCGA-61-1721-01A", 
  "TCGA-61-2017-01A", "TCGA-61-2095-01A"
)
excluded_low_conf <- c(
  "TCGA-10-0934-01A", "TCGA-13-1477-01A", "TCGA-24-0966-01A", 
  "TCGA-24-1565-01A", "TCGA-24-2036-01A", "TCGA-24-2038-01A", 
  "TCGA-29-1690-01A", "TCGA-29-1785-01A", "TCGA-61-1738-01A", 
  "TCGA-61-2092-01A", "TCGA-09-1664-01A", "TCGA-25-2408-01A", 
  "TCGA-29-1771-01A", "TCGA-61-2095-02A"
)
included <- readLines("New Text Document.txt")
anno_data <- data.frame(Sample_ID = ordered_samples) %>%
  mutate(
    # 
    QC_Status = case_when(
      Sample_ID %in% excluded_not_hgsoc ~ "Excluded (Not HGSOC)",
      Sample_ID %in% excluded_low_conf ~ "Excluded (Low Confidence)",
      Sample_ID %in% included ~ "Included (HGSOC Transcriptomics)",
      TRUE ~ "Others"
    ),
    # 
    Y_index = as.numeric(factor(Sample_ID, levels = ordered_samples)),
    # 
    X_pos = 1 
  ) %>%
  # 
  mutate(QC_Status = factor(QC_Status, levels = c("Included (HGSOC Transcriptomics)", "Excluded (Low Confidence)", "Excluded (Not HGSOC)", "Others")))

# Anno plot
p_anno <- ggplot(anno_data, aes(x = X_pos, y = Y_index, fill = QC_Status)) +
  geom_tile() +
  # 
  scale_fill_manual(values = c(
    "Included (HGSOC Transcriptomics)" = "forestgreen", 
    "Excluded (Low Confidence)" = "orange", 
    "Excluded (Not HGSOC)" = "firebrick",
    "Others" = "gray90"
  )) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
  theme_void() +
  theme(
    # 
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length.y = unit(0.2, "cm"),
    axis.title.y = element_text(size = 12, face = "bold", angle = 90, margin = margin(r = 10)),
    # 
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold", vjust = 0.8), 
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"), 
    legend.key.width = unit(0.5, "cm")   
  ) +
  labs(y = "Index of tumours ordered by FGA")

# Heatmap plot
p_heatmap <- ggplot(plot_cnv_segments) +
  geom_rect(aes(xmin = Start, xmax = End, 
                ymin = Y_index - 0.5, ymax = Y_index + 0.5, 
                fill = Segment_Mean)) +
  scale_fill_gradient2(low = "navyblue", mid = "white", high = "firebrick", 
                       midpoint = 0, limits = c(-1, 1), oob = squish, 
                       name = "Log2 Ratio") +
  facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),    
    axis.ticks = element_blank(),   
    axis.title = element_blank(),
    panel.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(size = 9, face = "bold", margin = margin(b = 5)),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.3, "cm"), 
    legend.title = element_text(size = 10, face = "bold", vjust = 0.8),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 15))
  ) +
  labs(
    title = "Genome-wide Copy Number Alterations in TCGA-OV (606 tumours/589 cases)",
    y = "Index of tumours ordered by FGA" 
  )

# Final plot 
final_plot <- p_anno + p_heatmap + plot_layout(widths = c(1, 40), guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",       
    legend.box.just = "center",      
    legend.margin = margin(t = 10)   
  )
ggsave("HGSOC_QC_Annotated_Heatmap.png", final_plot, width = 17, height = 9, dpi = 600)
ggsave("HGSOC_QC_Annotated_Heatmap.pdf", final_plot, width = 17, height = 9)

# HRD ASCN inference-----------------------------------------------------------

sample_TCGA_OV_ASCN <- read_tsv("GDC_TCGA-OV_ASCN_segment_sample_sheet.2024-04-01.tsv", show_col_types = FALSE)
ascn_dir <- paste0(getwd(), "/[DATA] GDC_TCGA-OV_ASCN_segment.2024-04-01")
# install.packages("devtools")
# library(devtools)
# remotes::install_github("ShixiangWang/copynumber")
# install.packages("https://cran.r-project.org/src/contrib/Archive/sequenza/sequenza_3.0.0.tar.gz", repos = NULL, type = "source")
# devtools::install_github("sztup/scarHRD")
library(scarHRD)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
dir.create("Temp_scarHRD", showWarnings = FALSE)

# Single sample inference function
calc_hrd_with_scarHRD <- function(file_name, sample_id) {
  file_path <- file.path(ascn_dir, file_name)
  if (!file.exists(file_path)) {
    return(tibble(Sample_ID = sample_id, HRD_LOH = NA, LST = NA, TAI = NA, HRD_Sum = NA, Note = "File missing"))
  }
  df <- read_tsv(file_path, show_col_types = FALSE)
  if (nrow(df) == 0 || !all(c("Chromosome", "Start", "End", "Copy_Number", "Major_Copy_Number", "Minor_Copy_Number") %in% colnames(df))) {
    return(tibble(Sample_ID = sample_id, HRD_LOH = NA, LST = NA, TAI = NA, HRD_Sum = NA, Note = "Invalid data format"))
  }
  estimated_ploidy <- sum(df$Copy_Number * (df$End - df$Start), na.rm = TRUE) / sum(df$End - df$Start, na.rm = TRUE)
  df_scar <- df %>%
    mutate(Chromosome = gsub("chr", "", Chromosome)) %>%
    filter(Chromosome %in% as.character(1:22)) %>%
    select(
      Chromosome = Chromosome,
      Start_position = Start,
      End_position = End,
      total_cn = Copy_Number,
      A_cn = Major_Copy_Number,   
      B_cn = Minor_Copy_Number    
    ) %>%
    drop_na() %>% 
    mutate(
      SampleID = sample_id,
      ploidy = round(estimated_ploidy)
    ) %>%
    select(SampleID, Chromosome, Start_position, End_position, total_cn, A_cn, B_cn, ploidy)
  temp_file <- file.path("Temp_scarHRD", paste0(sample_id, "_temp.txt"))
  write_tsv(df_scar, temp_file)
  res <- tryCatch({
    capture.output(
      score <- scar_score(temp_file, reference = "grch38", seqz = FALSE)
    )
    as.data.frame(score, check.names = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  unlink(temp_file)
  if (is.null(res) || nrow(res) == 0) {
    return(tibble(Sample_ID = sample_id, HRD_LOH = NA, LST = NA, TAI = NA, HRD_Sum = NA, Note = "scarHRD Error"))
  }
  return(tibble(
    Sample_ID = sample_id,
    HRD_LOH = res$HRD[1],               #LPC
    LST = res$LST[1],                   #LGA
    TAI = res$`Telomeric AI`[1],        
    HRD_Sum = res$`HRD-sum`[1],         
    Note = "Success"
  ))
}

# Batch calculation
cat("Starting mass HRD feature extraction for", nrow(sample_TCGA_OV_ASCN), "samples...\n")
hrd_features_TCGA <- pmap_dfr(
  list(
    sample_TCGA_OV_ASCN$`File Name`, 
    sample_TCGA_OV_ASCN$`Sample ID`
  ), 
  calc_hrd_with_scarHRD
)

# Clean results
hrd_features_clean <- hrd_features_TCGA %>%
  filter(Note == "Success") %>%
  arrange(Sample_ID, desc(HRD_Sum)) %>% 
  distinct(Sample_ID, .keep_all = TRUE) %>%
  select(Sample_ID, HRD_LOH, LST, TAI, HRD_Sum)

# Export
write_csv(hrd_features_TCGA, "hrd_features_TCGA_OV_raw.csv")
write_csv(hrd_features_clean, "hrd_features_TCGA_OV_clean.csv")