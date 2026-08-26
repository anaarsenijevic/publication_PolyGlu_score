#!/usr/bin/env Rscript

# ==============================================================================
# Pan-cancer progression-free interval analysis using a 9-gene RF tumor score
# ==============================================================================
#
# This script:
# 1. trains a pan-cancer random forest classifier on matched TCGA tumour-normal cohorts
# 2. scores tumour samples across TCGA projects using local SummarizedExperiment RDS files
# 3. downloads curated TCGA PanCancer Atlas survival data from UCSC Xena
# 4. uses curated PFI and PFI.time as the progression-free interval endpoint
# 5. generates Kaplan-Meier plots comparing top vs bottom RF score quartiles
#
#
# Input:
# - dds_*_matched_group.rds files
# - TCGA-<PROJECT>_RNAseq.rds files
#
# Output:
# - RF score table
# - PFI analysis dataset
# - KM plot PDF/PNG
# - Cox model result tables
#
# Notes:
# - Edit all input paths below before running.
# - The fixed 9-gene panel is:
#   TTLL4, TTLL5, TTLL6, TTLL7, TTLL11, AGBL1, AGBL3, AGBL4, AGBL5
# - PFI is taken from the curated PanCancer Atlas survival table:
#   Survival_SupplementalTable_S1_20171025_xena_sp
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(SummarizedExperiment)
  library(biomaRt)
  library(randomForest)
  library(caret)
  library(pROC)
  library(survival)
  library(survminer)
  library(UCSCXenaTools)
  library(broom)
})

set.seed(123)

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

tcga_dds_dir <- "path-to-tcga-dds-files"
gdc_rds_dir  <- "path-to-gdc-rds-files"
out_dir      <- "results_pfi"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(tcga_dds_dir)) {
  stop("Please set 'tcga_dds_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(gdc_rds_dir)) {
  stop("Please set 'gdc_rds_dir' to the folder containing TCGA-<PROJECT>_RNAseq.rds files.")
}

feature_genes <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

cancer_types <- c(
  "TCGA-UCEC", "TCGA-THCA", "TCGA-STAD", "TCGA-READ", "TCGA-PRAD",
  "TCGA-LUSC", "TCGA-LUAD", "TCGA-LIHC", "TCGA-KIRP", "TCGA-KIRC",
  "TCGA-KICH", "TCGA-HNSC", "TCGA-ESCA", "TCGA-COAD", "TCGA-CESC",
  "TCGA-BLCA", "TCGA-BRCA"
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

strip_ver <- function(x) sub("\\..*", "", x)

collapse_gene_duplicates <- function(expr_mat, symbols) {
  as.data.frame(expr_mat) %>%
    mutate(.symbol = symbols) %>%
    group_by(.symbol) %>%
    summarise(
      across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )
}

format_p <- function(p) {
  ifelse(p < 1e-3, "<0.001", signif(p, 2))
}

get_ensembl_connection <- function() {
  tryCatch(
    useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"),
    error = function(e) useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  )
}

plot_km <- function(df, title, y_label, fit_cox, out_prefix = NULL) {

  # Low is always the reference group
  df <- df %>%
    mutate(
      grp = factor(
        grp,
        levels = c("Low", "High")
      )
    )

  # Cox model: HR = High vs Low
  fit_cox <- coxph(
    Surv(time_years, status) ~ grp,
    data = df
  )

  cox_summary <- summary(fit_cox)

  hr <- cox_summary$coef[1, "exp(coef)"]
  hr_l <- cox_summary$conf.int[1, "lower .95"]
  hr_u <- cox_summary$conf.int[1, "upper .95"]

  # Kaplan-Meier fit
  fit <- survfit(
    Surv(time_years, status) ~ grp,
    data = df
  )

  # Explicit log-rank test
  logrank <- survdiff(
    Surv(time_years, status) ~ grp,
    data = df
  )

  logrank_p <- pchisq(
    logrank$chisq,
    df = 1,
    lower.tail = FALSE
  )

  p_label <- if (logrank_p < 0.0001) {
    "< 0.0001"
  } else {
    format.pval(
      logrank_p,
      digits = 2
    )
  }

  hr_label <- sprintf(
    "HR = %.2f (95%% CI %.2f-%.2f)\nLog-rank p = %s",
    hr,
    hr_l,
    hr_u,
    p_label
  )

  km <- ggsurvplot(
    fit,
    data = df,
    pval = FALSE,
    conf.int = TRUE,
    risk.table = FALSE,
    censor = FALSE,
    surv.size = 1.4,

    legend.title = title,
    legend.labs = c(
      "Low (Q1)",
      "High (Q4)"
    ),

    xlab = "Time (years)",
    ylab = y_label,

    xlim = c(0, 10),
    break.x.by = 1,

    palette = c(
      "#4F7052",
      "#B45757"
    ),

    ggtheme = theme_classic(
      base_size = 14
    ) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(
          color = "black",
          linewidth = 0.8
        ),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.8
        ),
        axis.text = element_text(
          color = "black"
        ),
        axis.title = element_text(
          color = "black",
          face = "bold"
        ),
        legend.title = element_text(
          face = "bold"
        ),
        legend.position = "top"
      )
  )

  # 5-year PFI
  s5 <- summary(
    fit,
    times = 5
  )

  fiveyr <- data.frame(
    strata = s5$strata,
    surv = s5$surv
  ) %>%
    filter(
      !is.na(surv)
    ) %>%
    mutate(
      label = sprintf(
        "%.1f%%",
        100 * surv
      )
    )

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
      aes(
        x = 0,
        xend = 5,
        y = surv,
        yend = surv
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.6
    ) +
    geom_segment(
      data = fiveyr,
      aes(
        x = 5,
        xend = 5,
        y = 0,
        yend = surv
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.6
    ) +
    geom_text(
      data = fiveyr,
      aes(
        x = 5,
        y = surv,
        label = label
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      color = "black",
      vjust = -0.6,
      hjust = 0.5,
      size = 4.5
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.15,
      label = hr_label,
      hjust = 0,
      size = 4.5
    )

  if (!is.null(out_prefix)) {

    pdf(
      file.path(
        out_dir,
        paste0(out_prefix, ".pdf")
      ),
      width = 7,
      height = 8
    )

    print(km)
    dev.off()

    png(
      file.path(
        out_dir,
        paste0(out_prefix, ".png")
      ),
      width = 2100,
      height = 2400,
      res = 300
    )

    print(km)
    dev.off()
  }

  km
}
run_quartile_km <- function(df, var, title, y_label, out_prefix = NULL) {
  q1 <- quantile(df[[var]], 0.25, na.rm = TRUE)
  q3 <- quantile(df[[var]], 0.75, na.rm = TRUE)

  q_df <- df %>%
    mutate(
      grp = case_when(
        .data[[var]] <= q1 ~ "Low",
        .data[[var]] >= q3 ~ "High",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(grp)) %>%
    mutate(grp = factor(grp, levels = c("Low", "High")))

  message(sprintf("%s quartiles: Q1 = %.4f | Q3 = %.4f", title, q1, q3))

  q_summary <- q_df %>%
    group_by(grp) %>%
    summarise(
      n = n(),
      events = sum(status),
      event_rate = mean(status),
      median_score = median(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    )

  print(q_summary)

  fit_cox <- coxph(
  Surv(time_years, status) ~ grp,
  data = q_df
)
  print(summary(fit_cox))

  km <- plot_km(
    df = q_df,
    title = title,
    y_label = y_label,
    fit_cox = fit_cox,
    out_prefix = out_prefix
  )

  invisible(list(
    data = q_df,
    summary = q_summary,
    q1 = q1,
    q3 = q3,
    fit = fit_cox,
    km = km
  ))
}

# ------------------------------------------------------------------------------
# Ensembl to HGNC mapping
# ------------------------------------------------------------------------------

message("Fetching Ensembl to HGNC symbol mapping...")

ensembl <- get_ensembl_connection()

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
) %>%
  filter(external_gene_name != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

symbol_to_ensembl <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = feature_genes,
  mart = ensembl
) %>%
  filter(hgnc_symbol %in% feature_genes)

ensembl_lookup <- setNames(
  symbol_to_ensembl$hgnc_symbol,
  symbol_to_ensembl$ensembl_gene_id
)

# ------------------------------------------------------------------------------
# Train pan-cancer RF model on matched TCGA cohorts
# ------------------------------------------------------------------------------

message("Training pan-cancer RF model...")

dds_files <- list.files(
  tcga_dds_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No matched TCGA DESeq2 files found in 'tcga_dds_dir'.")
}

training_list <- list()

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

  expr <- t(mat[present, , drop = FALSE]) %>%
    as.data.frame()

  expr$sample <- rownames(expr)

  cond <- tryCatch(
    colData(dds)[expr$sample, "condition"],
    error = function(e) NULL
  )
  if (is.null(cond)) next

  expr$condition <- as.character(cond)
  expr$project <- project

  training_list[[project]] <- expr
}

training_df <- bind_rows(training_list)

if (nrow(training_df) == 0) {
  stop("No training data could be assembled from matched TCGA cohorts.")
}

training_clean <- training_df %>%
  select(all_of(feature_genes), condition, sample, project) %>%
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

rf_cv_pan <- train(
  x = training_clean[, feature_genes, drop = FALSE],
  y = training_clean$condition_bin,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

cv_auc <- as.numeric(
  pROC::auc(
    pROC::roc(rf_cv_pan$pred$obs, rf_cv_pan$pred$Tumor)
  )
)

message(sprintf(
  "Pan-cancer RF trained. Best mtry = %s; CV AUC = %.3f",
  rf_cv_pan$bestTune$mtry,
  cv_auc
))

# ------------------------------------------------------------------------------
# Load TCGA project RDS files and score tumour samples
# ------------------------------------------------------------------------------

message("Scoring TCGA tumour samples across projects...")

expression_list <- list()

for (cancer in cancer_types) {
  file_path <- file.path(gdc_rds_dir, paste0(cancer, "_RNAseq.rds"))

  if (!file.exists(file_path)) {
    warning("Missing file: ", file_path)
    next
  }

  message("  reading: ", cancer)

  se <- readRDS(file_path)

  assays_avail <- assayNames(se)
  counts_assay <- if ("unstranded" %in% assays_avail) "unstranded" else assays_avail[1]
  raw_counts <- assay(se, counts_assay)

  rd <- as.data.frame(rowData(se))
  gene_name_col <- intersect(
    c("gene_name", "hgnc_symbol", "external_gene_name"),
    names(rd)
  )

  if (length(gene_name_col) > 0) {
    symbols_for_rows <- rd[[gene_name_col[1]]]
  } else {
    rn <- rownames(raw_counts)

    if (!is.null(rn) && mean(grepl("^ENSG", rn)) < 0.2) {
      symbols_for_rows <- rn
    } else {
      ens <- if ("gene_id" %in% names(rd)) rd$gene_id else rn
      symbols_for_rows <- ensembl_lookup[strip_ver(ens)]
    }
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(raw_counts),
    colData = as.data.frame(colData(se)),
    design = ~ 1
  )

  vsd <- vst(dds, blind = TRUE)
  expr_vst <- assay(vsd)

  if (length(symbols_for_rows) != nrow(expr_vst)) {
    warning(cancer, ": symbol vector length mismatch; skipping.")
    next
  }

  keep_idx <- which(symbols_for_rows %in% feature_genes)
  if (!length(keep_idx)) {
    warning(cancer, ": none of the panel genes found; skipping.")
    next
  }

  expr_collapsed <- collapse_gene_duplicates(
    expr_mat = expr_vst[keep_idx, , drop = FALSE],
    symbols = symbols_for_rows[keep_idx]
  )

  mat_sym <- as.matrix(expr_collapsed[, -1])
  rownames(mat_sym) <- expr_collapsed$.symbol

  missing_genes <- setdiff(feature_genes, rownames(mat_sym))
  if (length(missing_genes) > 0) {
    mat_sym <- rbind(
      mat_sym,
      matrix(
        NA_real_,
        nrow = length(missing_genes),
        ncol = ncol(mat_sym),
        dimnames = list(missing_genes, colnames(mat_sym))
      )
    )
  }

  mat_sym <- mat_sym[feature_genes, , drop = FALSE]

  expr_t <- as.data.frame(t(mat_sym))
  expr_t$sample <- substr(rownames(expr_t), 1, 12)
  expr_t$Cancer_Type <- cancer

  expression_list[[cancer]] <- expr_t
}

if (!length(expression_list)) {
  stop("No TCGA expression datasets were loaded from 'gdc_rds_dir'.")
}

expr_df <- bind_rows(expression_list) %>%
  mutate(sample = substr(sample, 1, 12))

missing_features <- setdiff(feature_genes, colnames(expr_df))
if (length(missing_features) > 0) {
  for (gene in missing_features) {
    expr_df[[gene]] <- NA_real_
  }
}

newx <- expr_df[, feature_genes, drop = FALSE]
valid_rows <- rowSums(is.na(newx)) == 0

rf_tumor <- rep(NA_real_, nrow(expr_df))

if (any(valid_rows)) {
  pp <- predict(
    rf_cv_pan,
    newdata = newx[valid_rows, , drop = FALSE],
    type = "prob"
  )

  pp_df <- as.data.frame(pp)
  tumor_col <- which(tolower(colnames(pp_df)) %in% c("tumor", "tumour", "1", "positive"))

  if (!length(tumor_col)) {
    stop("No Tumor probability column found in RF prediction output.")
  }

  rf_tumor[valid_rows] <- as.numeric(pp_df[[tumor_col[1]]])
}

scores_df <- expr_df %>%
  select(sample, Cancer_Type) %>%
  mutate(rf_tumor_score = rf_tumor) %>%
  distinct(sample, Cancer_Type, .keep_all = TRUE)

write.csv(
  scores_df,
  file.path(out_dir, "tcga_pan_rf_scores.csv"),
  row.names = FALSE
)

message("Scored samples: ", sum(!is.na(scores_df$rf_tumor_score)))

# ------------------------------------------------------------------------------
# Download curated PanCancer Atlas PFI data from UCSC Xena
# ------------------------------------------------------------------------------

message("Downloading curated PanCancer Atlas PFI data from UCSC Xena...")

surv_query <- XenaGenerate() %>%
  XenaFilter(
    filterDatasets = "Survival_SupplementalTable_S1_20171025_xena_sp"
  ) %>%
  XenaQuery()

surv_download <- XenaDownload(surv_query)
surv_xena <- XenaPrepare(surv_download)

if (!all(c("sample", "PFI", "PFI.time") %in% names(surv_xena))) {
  stop("Downloaded Xena survival table does not contain sample, PFI, and PFI.time columns.")
}

pfi_df <- surv_xena %>%
  mutate(
    sample = substr(sample, 1, 12),
    PFI = as.integer(PFI),
    PFI.time = as.numeric(PFI.time)
  ) %>%
  filter(
    !is.na(PFI),
    !is.na(PFI.time),
    PFI.time > 0
  ) %>%
  select(sample, PFI, PFI.time) %>%
  distinct(sample, .keep_all = TRUE)

pfi_summary <- pfi_df %>%
  summarise(
    n = n(),
    events = sum(PFI == 1),
    event_rate = mean(PFI == 1),
    median_time_days = median(PFI.time)
  )

print(pfi_summary)

write.csv(
  pfi_summary,
  file.path(out_dir, "xena_pfi_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Merge RF scores with PFI endpoint
# ------------------------------------------------------------------------------

dat_pfi <- scores_df %>%
  inner_join(pfi_df, by = "sample") %>%
  mutate(
    time = PFI.time,
    status = PFI,
    time_years = time / 365.25
  ) %>%
  filter(
    is.finite(rf_tumor_score),
    is.finite(time),
    time > 0,
    status %in% c(0, 1)
  )

message("PFI analysis dataset: N = ", nrow(dat_pfi), " | Events = ", sum(dat_pfi$status))

write.csv(
  dat_pfi,
  file.path(out_dir, "tcga_pan_rf_score_pfi_dataset.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# RF score PFI KM analysis
# ------------------------------------------------------------------------------

res_rf_q_pfi <- run_quartile_km(
  df = dat_pfi,
  var = "rf_tumor_score",
  title = "Polyglutamylation score",
  y_label = "Progression-free interval probability",
  out_prefix = "PanCancer_PFI_by_polyglutamylation_score"
)

write.csv(
  res_rf_q_pfi$data,
  file.path(out_dir, "rf_score_pfi_top_bottom_quartiles.csv"),
  row.names = FALSE
)

write.csv(
  res_rf_q_pfi$summary,
  file.path(out_dir, "rf_score_pfi_top_bottom_quartile_summary.csv"),
  row.names = FALSE
)

