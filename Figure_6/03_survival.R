#!/usr/bin/env Rscript

# Pan-cancer survival and vital status analysis using a 9-gene RF tumor score
#
# This script:
# 1. trains a pan-cancer random forest classifier on matched TCGA tumour-normal cohorts
# 2. scores tumour samples across TCGA projects using local SummarizedExperiment RDS files
# 3. harmonizes clinical annotations and computes overall survival
# 4. generates Kaplan-Meier plots comparing top vs bottom RF score quartiles
# 5. generates Kaplan-Meier plots comparing top vs bottom TTLL11 expression quartiles
# 6. fits a joint Cox model including TTLL11 expression and RF tumor score
# 7. generates a forest plot for the joint Cox model
#
# Input:
# - dds_*_matched_group.rds files
# - TCGA-<PROJECT>_RNAseq.rds files
#
# Notes:
# - Edit all input paths below before running
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - Overall survival is computed from days to death or last follow-up
# - Survival analysis includes only patients in the bottom and top RF score quartiles
# - Vital status analysis compares RF tumor score distributions in Alive vs Deceased patients
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(SummarizedExperiment)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(TCGAbiolinks)
  library(coin)
  library(effsize)
  library(broom)
})

set.seed(123)

# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

tcga_dds_dir <- "path-to-tcga-dds-files"
gdc_rds_dir <- "path-to-gdc-rds-files"

if (!dir.exists(tcga_dds_dir)) {
  stop("Please set 'tcga_dds_dir' to the folder containing dds_*_matched_group.rds files.")
}
if (!dir.exists(gdc_rds_dir)) {
  stop("Please set 'gdc_rds_dir' to the folder containing TCGA-<PROJECT>_RNAseq.rds files.")
}

# ------------------------------------------------------------------
# Fixed 9-gene panel
# ------------------------------------------------------------------

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

strip_ver <- function(x) sub("\\..*", "", x)

numify <- function(v) suppressWarnings(as.numeric(v))

coalesce_num <- function(df, cols) {
  present <- cols[cols %in% names(df)]
  if (!length(present)) return(rep(NA_real_, nrow(df)))
  out <- numify(df[[present[1]]])
  for (cc in present[-1]) {
    out <- dplyr::coalesce(out, numify(df[[cc]]))
  }
  out
}

map_deceased <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[y %in% c("dead", "deceased", "1", "true", "yes")] <- 1
  out[y %in% c("alive", "living", "0", "false", "no")] <- 0
  out
}

# ------------------------------------------------------------------
# Ensembl to HGNC mapping
# ------------------------------------------------------------------

message("Fetching Ensembl to HGNC symbol mapping...")
ensembl <- tryCatch(
  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"),
  error = function(e) useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# ------------------------------------------------------------------
# Train pan-cancer RF model on matched TCGA cohorts
# ------------------------------------------------------------------

message("Training pan-cancer RF model...")
dds_files <- list.files(
  tcga_dds_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No matched TCGA DESeq2 files found in 'tcga_dds_dir'.")
}

all_data <- list()

for (file in dds_files) {
  project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
  message("  processing: ", project)

  dds <- readRDS(file)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  ens_ids <- strip_ver(rownames(mat))
  syms <- gene_map$external_gene_name[match(ens_ids, gene_map$ensembl_gene_id)]

  keep <- !is.na(syms) & syms != ""
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- syms[keep]

  present <- intersect(feature_genes, rownames(mat))
  if (!length(present)) next

  expr <- t(mat[present, , drop = FALSE]) %>% as.data.frame()
  expr$sample <- rownames(expr)

  cond <- tryCatch(colData(dds)[expr$sample, "condition"], error = function(e) NULL)
  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  expr$project <- project
  all_data[[length(all_data) + 1]] <- expr
}

tcga_df <- bind_rows(all_data)

if (nrow(tcga_df) == 0) {
  stop("No training data could be assembled from the matched TCGA cohorts.")
}

tcga_clean <- tcga_df %>%
  dplyr::select(all_of(feature_genes), condition, sample, project) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

x_pan <- tcga_clean[, feature_genes, drop = FALSE]
y_pan <- tcga_clean$condition_bin

rf_cv_pan <- train(
  x = x_pan,
  y = y_pan,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

message(sprintf(
  "Pan-cancer RF trained. Best mtry = %s; CV AUC ≈ %.3f",
  rf_cv_pan$bestTune$mtry,
  as.numeric(pROC::auc(pROC::roc(rf_cv_pan$pred$obs, rf_cv_pan$pred$Tumor)))
))

# ------------------------------------------------------------------
# Load TCGA project RDS files and score tumour samples
# ------------------------------------------------------------------

message("Scoring TCGA tumour samples across projects...")

cancer_types <- c(
  "TCGA-UCEC", "TCGA-THCA", "TCGA-STAD", "TCGA-READ", "TCGA-PRAD",
  "TCGA-LUSC", "TCGA-LUAD", "TCGA-LIHC", "TCGA-KIRP", "TCGA-KIRC",
  "TCGA-KICH", "TCGA-HNSC", "TCGA-ESCA", "TCGA-COAD", "TCGA-CESC",
  "TCGA-BLCA", "TCGA-BRCA"
)

symbol_to_ensembl <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = feature_genes,
  mart = ensembl
) %>%
  filter(hgnc_symbol %in% feature_genes)

ensembl_lookup <- setNames(symbol_to_ensembl$hgnc_symbol, symbol_to_ensembl$ensembl_gene_id)

expression_list <- list()
clinical_list <- list()

for (cancer in cancer_types) {
  file_path <- file.path(gdc_rds_dir, paste0(cancer, "_RNAseq.rds"))
  if (!file.exists(file_path)) next

  message("  reading: ", cancer)
  se <- readRDS(file_path)

  assays_avail <- assayNames(se)
  counts_assay <- if ("unstranded" %in% assays_avail) "unstranded" else assays_avail[1]
  if (!counts_assay %in% assays_avail) next

  raw_counts <- assay(se, counts_assay)

  rd <- as.data.frame(rowData(se))
  gene_name_col <- intersect(c("gene_name", "hgnc_symbol", "external_gene_name"), names(rd))

  if (length(gene_name_col) > 0) {
    symbols_for_rows <- rd[[gene_name_col[1]]]
  } else {
    rn <- rownames(raw_counts)
    if (!is.null(rn) && mean(grepl("^ENSG", rn)) < 0.2) {
      symbols_for_rows <- rn
    } else {
      ens <- if ("gene_id" %in% names(rd)) rd$gene_id else rn
      ens <- strip_ver(ens)
      symbols_for_rows <- ensembl_lookup[ens]
    }
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(raw_counts),
    colData = as.data.frame(colData(se)),
    design = ~1
  )

  vsd <- vst(dds, blind = TRUE)
  expr_vst <- assay(vsd)

  if (length(symbols_for_rows) != nrow(expr_vst)) next

  keep_idx <- which(symbols_for_rows %in% feature_genes)
  if (!length(keep_idx)) next

  expr_sel <- expr_vst[keep_idx, , drop = FALSE]
  sym_sel <- symbols_for_rows[keep_idx]

  tmp <- as.data.frame(expr_sel) %>%
    mutate(.symbol = sym_sel)

  expr_collapsed <- tmp %>%
    group_by(.symbol) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

  mat_sym <- as.matrix(expr_collapsed[, -1])
  rownames(mat_sym) <- expr_collapsed$.symbol

  miss_rows <- setdiff(feature_genes, rownames(mat_sym))
  if (length(miss_rows) > 0) {
    mat_sym <- rbind(
      mat_sym,
      matrix(
        NA_real_,
        nrow = length(miss_rows),
        ncol = ncol(mat_sym),
        dimnames = list(miss_rows, colnames(mat_sym))
      )
    )
  }

  mat_sym <- mat_sym[feature_genes, , drop = FALSE]

  expr_t <- as.data.frame(t(mat_sym))
  expr_t$sample <- substr(rownames(expr_t), 1, 12)
  expr_t$Cancer_Type <- cancer
  expression_list[[cancer]] <- expr_t

  clinical_data <- as.data.frame(colData(se))
  clinical_data$bcr_patient_barcode <- substr(clinical_data$bcr_patient_barcode, 1, 12)
  clinical_list[[cancer]] <- clinical_data
}

if (!length(expression_list)) {
  stop("No TCGA expression datasets were loaded from 'gdc_rds_dir'.")
}

expr_df <- bind_rows(expression_list, .id = "source_id")

if (!"Cancer_Type" %in% names(expr_df)) {
  expr_df$Cancer_Type <- expr_df$source_id
}
if (!"sample" %in% names(expr_df)) {
  if (!is.null(rownames(expr_df))) {
    expr_df$sample <- substr(rownames(expr_df), 1, 12)
  } else {
    stop("Combined expression table lacks 'sample' and rownames.")
  }
}

expr_df$sample <- substr(expr_df$sample, 1, 12)

missing_feats <- setdiff(feature_genes, colnames(expr_df))
if (length(missing_feats) > 0) {
  for (m in missing_feats) expr_df[[m]] <- NA_real_
}

newx <- expr_df[, feature_genes, drop = FALSE]
valid_rows <- rowSums(is.na(newx)) < ncol(newx)

rf_tumor <- rep(NA_real_, nrow(expr_df))
if (any(valid_rows)) {
  pp <- predict(rf_cv_pan, newdata = newx[valid_rows, , drop = FALSE], type = "prob")
  pp_df <- as.data.frame(pp)
  cn <- tolower(colnames(pp_df))
  col_tum <- which(cn %in% c("tumor", "tumour", "1", "positive"))
  if (length(col_tum) == 0) {
    col_tum <- if (ncol(pp_df) == 1) 1L else stop("No Tumor probability column found.")
  }
  rf_tumor[valid_rows] <- as.numeric(pp_df[[col_tum[1]]])
}

rf_scores_df <- expr_df %>%
  dplyr::select(sample, Cancer_Type) %>%
  mutate(rf_tumor_score = rf_tumor)

# ------------------------------------------------------------------
# Harmonize clinical data and compute overall survival
# ------------------------------------------------------------------

clinical_list2 <- lapply(clinical_list, function(df) {
  df <- as.data.frame(df)

  if (!"bcr_patient_barcode" %in% names(df)) {
    guess <- intersect(c("patient_id", "case_submitter_id", "submitter_id", "barcode"), names(df))
    if (length(guess) > 0) df$bcr_patient_barcode <- df[[guess[1]]]
  }

  df[] <- lapply(df, function(x) {
    if (is.list(x)) {
      vapply(x, function(y) paste0(y, collapse = "; "), character(1))
    } else if (is.factor(x)) {
      as.character(x)
    } else {
      as.character(x)
    }
  })

  df$bcr_patient_barcode <- substr(df$bcr_patient_barcode, 1, 12)
  df
})

merged_clinical <- bind_rows(clinical_list2, .id = "Cancer_Type")
merged_clinical$sample <- merged_clinical$bcr_patient_barcode
merged_clinical$deceased <- map_deceased(merged_clinical$vital_status)

death_cols <- c(
  "days_to_death", "paper_days_to_death", "OS.time", "os_time", "os_days",
  "OS_Days", "death_days_to", "paper_OS.time", "paper_OS_days"
)

followup_cols <- c(
  "days_to_last_follow_up", "days_to_last_followup", "last_contact_days_to",
  "days_to_last_contact", "days_to_last_known_status", "days_to_last_known_alive",
  "followup_time", "follow_up_time", "followup_days",
  "paper_days_to_last_followup", "paper_followup_days", "paper_os_days"
)

merged_clinical$dt_death <- coalesce_num(merged_clinical, death_cols)
merged_clinical$dt_followup <- coalesce_num(merged_clinical, followup_cols)

needs_patch <- merged_clinical %>%
  group_by(Cancer_Type) %>%
  summarise(
    n_events = sum(deceased == 1, na.rm = TRUE),
    any_follow = any(!is.na(dt_followup)),
    .groups = "drop"
  ) %>%
  filter(n_events > 0 & !any_follow) %>%
  pull(Cancer_Type)

if (length(needs_patch)) {
  message("Patching follow-up from GDC for: ", paste(needs_patch, collapse = ", "))

  barcode_cols <- c("bcr_patient_barcode", "submitter_id", "case_submitter_id", "patient_id", "barcode")
  followup_patch_cols <- c(
    "days_to_last_follow_up", "days_to_last_followup", "last_contact_days_to",
    "days_to_last_contact", "days_to_last_known_status", "days_to_last_known_alive",
    "followup_time", "follow_up_time", "followup_days",
    "paper_days_to_last_followup", "paper_followup_days"
  )

  patch_df <- lapply(needs_patch, function(prj) {
    cli <- tryCatch(GDCquery_clinic(project = prj, type = "clinical"), error = function(e) NULL)
    if (is.null(cli) || !nrow(cli)) return(NULL)

    cli <- as.data.frame(cli)
    bc_col <- intersect(barcode_cols, names(cli))
    if (!length(bc_col)) return(NULL)

    cli$bcr_patient_barcode <- substr(cli[[bc_col[1]]], 1, 12)

    pf <- {
      present <- followup_patch_cols[followup_patch_cols %in% names(cli)]
      if (!length(present)) {
        rep(NA_real_, nrow(cli))
      } else {
        v <- suppressWarnings(as.numeric(cli[[present[1]]]))
        for (cc in present[-1]) {
          v <- dplyr::coalesce(v, suppressWarnings(as.numeric(cli[[cc]])))
        }
        v
      }
    }

    data.frame(
      Cancer_Type = prj,
      bcr_patient_barcode = cli$bcr_patient_barcode,
      patch_follow = pf,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  if (!is.null(patch_df) && nrow(patch_df)) {
    merged_clinical <- merged_clinical %>%
      left_join(patch_df, by = c("Cancer_Type", "bcr_patient_barcode")) %>%
      mutate(dt_followup = coalesce(dt_followup, patch_follow)) %>%
      dplyr::select(-patch_follow)
  }
}

merged_clinical$overall_survival_days <- ifelse(
  merged_clinical$deceased == 1,
  merged_clinical$dt_death,
  merged_clinical$dt_followup
)

# ------------------------------------------------------------------
# Merge RF scores with survival data
# ------------------------------------------------------------------

scores_df <- rf_scores_df %>%
  distinct(sample, Cancer_Type, .keep_all = TRUE)

merged_df <- merged_clinical %>%
  inner_join(scores_df, by = c("sample", "Cancer_Type")) %>%
  filter(!is.na(overall_survival_days), overall_survival_days > 0, !is.na(deceased)) %>%
  mutate(
    time = as.numeric(overall_survival_days),
    status = as.integer(deceased),
    time_years = time / 365.25
  )

message("N ready for survival: ", nrow(merged_df))

# ------------------------------------------------------------------
# Kaplan-Meier analysis: top vs bottom RF score quartiles
# ------------------------------------------------------------------

df_for_km <- merged_df %>%
  dplyr::select(time_years, status, rf_tumor_score) %>%
  filter(!is.na(time_years), !is.na(status), !is.na(rf_tumor_score), time_years > 0)

q1 <- quantile(df_for_km$rf_tumor_score, 0.25, na.rm = TRUE)
q3 <- quantile(df_for_km$rf_tumor_score, 0.75, na.rm = TRUE)

message(sprintf("Q1 = %.4f | Q3 = %.4f", q1, q3))

quartile_df <- df_for_km %>%
  mutate(
    RF_group = case_when(
      rf_tumor_score <= q1 ~ "Low RF Score (Q1)",
      rf_tumor_score >= q3 ~ "High RF Score (Q4)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(RF_group)) %>%
  mutate(
    RF_group = factor(
      RF_group,
      levels = c("High RF Score (Q4)", "Low RF Score (Q1)")
    )
  )

message("N used for KM (Q1 vs Q4): ", nrow(quartile_df))

fit_q <- survfit(Surv(time_years, status) ~ RF_group, data = quartile_df)

pal <- c("#B45757", "#4F7052")

km <- ggsurvplot(
  fit_q,
  data = quartile_df,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  censor = FALSE,
  surv.size = 1.4,
  legend.labs = levels(quartile_df$RF_group),
  legend.title = "Predicted Tumor Probability",
  title = "Pan-cancer survival: top vs bottom RF score quartiles",
  xlab = "Time (years)",
  ylab = "Survival Probability",
  xlim = c(0, 10),
  break.x.by = 1,
  palette = pal,
  ggtheme = theme_minimal(base_size = 14)
)

s5 <- summary(fit_q, times = 5)

fiveyr <- data.frame(strata = s5$strata, surv = s5$surv) %>%
  filter(!is.na(surv)) %>%
  mutate(label = sprintf("%.1f%%", 100 * surv))

km$plot <- km$plot +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1)
  ) +
  geom_segment(
    data = fiveyr,
    aes(x = 0, xend = 5, y = surv, yend = surv),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  geom_segment(
    data = fiveyr,
    aes(x = 5, xend = 5, y = 0, yend = surv),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  geom_text(
    data = fiveyr,
    aes(x = 5, y = surv, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4.5
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )

print(km)

# ------------------------------------------------------------------
# Pan-cancer RF tumor score by vital status
# ------------------------------------------------------------------

pan_df <- merged_df %>%
  filter(!is.na(status), !is.na(rf_tumor_score)) %>%
  mutate(
    Status = factor(
      ifelse(status == 1, "Deceased", "Alive"),
      levels = c("Alive", "Deceased")
    )
  )

wilx <- wilcox.test(rf_tumor_score ~ Status, data = pan_df, exact = FALSE)
p_overall <- wilx$p.value

cd <- effsize::cliff.delta(rf_tumor_score ~ Status, data = pan_df)
cd_est <- unname(cd$estimate)
cd_mag <- cd$magnitude

pan_df2 <- pan_df %>%
  filter(!is.na(Cancer_Type), is.finite(rf_tumor_score)) %>%
  mutate(Cancer_Type = factor(Cancer_Type))

strata_ok <- pan_df2 %>%
  group_by(Cancer_Type) %>%
  summarise(n_levels = n_distinct(Status), .groups = "drop") %>%
  filter(n_levels == 2) %>%
  pull(Cancer_Type)

pan_df_ok <- pan_df2 %>%
  filter(Cancer_Type %in% strata_ok)

wt_strat <- coin::wilcox_test(
  rf_tumor_score ~ Status | Cancer_Type,
  data = pan_df_ok,
  distribution = "approximate"
)

p_strat <- coin::pvalue(wt_strat)

cat(sprintf(
  "Overall Wilcoxon p = %.3g | Cliff's delta = %.3f (%s)\nStratified Wilcoxon p = %.3g\n",
  p_overall, cd_est, cd_mag, p_strat
))

write.csv(
  pan_df,
  file.path(out_dir, "tcga_pan_rf_scores_alive_vs_deceased.csv"),
  row.names = FALSE
)

y_min <- min(pan_df$rf_tumor_score, na.rm = TRUE)
y_max <- max(pan_df$rf_tumor_score, na.rm = TRUE)
y_rng <- y_max - y_min

ns <- pan_df %>%
  count(Status)

p_label <- paste0(
  "Wilcoxon p = ",
  ifelse(p_overall < 1e-3, "< 0.001", signif(p_overall, 3))
)

p_pan <- ggplot(pan_df, aes(x = Status, y = rf_tumor_score, fill = Status)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    alpha = 0.9,
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c("Alive" = "#4F7052", "Deceased" = "#B45757")
  ) +
  labs(
    title = "Pan-cancer RF tumor score by vital status",
    x = NULL,
    y = "RF Tumor Score (probability)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15),
    legend.position = "none",
    plot.margin = margin(15, 25, 15, 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  ) +
  coord_cartesian(ylim = c(y_min, y_max + 0.12 * y_rng), clip = "off") +
  annotate(
    "text",
    x = 1,
    y = y_min,
    label = paste0("n = ", ns$n[ns$Status == "Alive"]),
    vjust = 1.2,
    size = 5
  ) +
  annotate(
    "text",
    x = 2,
    y = y_min,
    label = paste0("n = ", ns$n[ns$Status == "Deceased"]),
    vjust = 1.2,
    size = 5
  ) +
  annotate(
    "text",
    x = 1.5,
    y = y_max + 0.08 * y_rng,
    label = p_label,
    vjust = 0,
    size = 5
  )

print(p_pan)

ggsave(
  file.path(out_dir, "PanTCGA_RFscore_Alive_vs_Deceased.png"),
  p_pan,
  width = 6,
  height = 5,
  dpi = 300
)

