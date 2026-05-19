# TCGA random forest cross-validation ROC analysis
#
# This script:
# 1. reads matched TCGA tumour-normal DESeq2 objects
# 2. extracts VST-normalized expression for 15-gene and 9-gene panels
# 3. trains random forest classifiers using 5-fold cross-validation
# 4. computes fold-wise and mean ROC curves
# 5. plots ROC curves and confusion matrices
#
# Input:
# - dds_*_matched_group.rds files
#
#
# Notes:
# - Edit 'tcga_dir' and 'output_dir' before running
# - Ensembl IDs are mapped to gene symbols using org.Hs.eg.db
# - Duplicate symbols are collapsed by mean expression

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
})

# ------------------------------------------------------------------
# Input/output paths
# ------------------------------------------------------------------

tcga_dir <- "path-to-tcga-dds-files"
output_dir <- "results/tcga_rf_cv"

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------
# Gene panels
# ------------------------------------------------------------------

genes_15 <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7",
  "TTLL9", "TTLL11", "TTLL13",
  "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1"
)

genes_9 <- c(
  "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL11",
  "AGBL1", "AGBL3", "AGBL4", "AGBL5"
)

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

  if (mode == "sum") {
    rowsum(mat, group = symbols)
  } else {
    sp <- split(seq_len(nrow(mat)), symbols)
    out <- vapply(
      sp,
      function(ix) colMeans(mat[ix, , drop = FALSE]),
      numeric(ncol(mat))
    )
    out <- t(out)
    rownames(out) <- names(sp)
    out
  }
}

load_tcga_expression <- function(tcga_dir, genes) {
  dds_files <- list.files(
    tcga_dir,
    pattern = "dds_.*_matched_group\\.rds$",
    full.names = TRUE
  )

  if (length(dds_files) == 0) {
    stop("No DESeq2 files found in tcga_dir.")
  }

  tcga_data <- list()

  for (file in dds_files) {
    project <- sub("dds_(.*)_matched_group\\.rds", "\\1", basename(file))
    message("Processing ", project)

    dds <- readRDS(file)
    vsd <- vst(dds, blind = TRUE)
    expr_mat <- assay(vsd)

    ens_ids <- sub("\\..*$", "", rownames(expr_mat))

    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ens_ids,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )

    expr_sym <- aggregate_dup_symbols(expr_mat, symbols)

    matched_genes <- intersect(genes, rownames(expr_sym))
    if (length(matched_genes) == 0) next

    expr_df <- t(expr_sym[matched_genes, , drop = FALSE]) %>%
      as.data.frame()

    for (g in setdiff(genes, colnames(expr_df))) {
      expr_df[[g]] <- NA_real_
    }

    expr_df <- expr_df[, genes, drop = FALSE]
    expr_df$sample <- rownames(expr_df)

    cond <- tryCatch({
      as.character(colData(dds)[expr_df$sample, "condition"])
    }, error = function(e) NULL)

    if (is.null(cond)) next

    expr_df$condition <- cond
    expr_df$project <- project

    tcga_data[[length(tcga_data) + 1]] <- expr_df
  }

  bind_rows(tcga_data) %>%
    dplyr::select(all_of(genes), condition, sample, project) %>%
    filter(if_all(all_of(genes), ~ !is.na(.))) %>%
    mutate(
      condition_bin = factor(
        ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
        levels = c("Normal", "Tumor")
      )
    )
}

run_rf_cv_roc <- function(data, genes, model_name, seed = 123) {
  df <- data %>%
    dplyr::select(all_of(genes), condition_bin) %>%
    mutate(condition_bin = factor(condition_bin, levels = c("Normal", "Tumor")))

  x <- df[, genes, drop = FALSE]
  y <- df$condition_bin

  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )

  set.seed(seed)
  rf_cv <- train(
    x = x,
    y = y,
    method = "rf",
    metric = "ROC",
    trControl = ctrl,
    ntree = 500
  )

  cv_preds <- rf_cv$pred %>%
    filter(mtry == rf_cv$bestTune$mtry) %>%
    mutate(Resample = factor(Resample))

  roc_list <- cv_preds %>%
    split(.$Resample) %>%
    lapply(function(d) {
      roc(
        response = d$obs,
        predictor = d$Tumor,
        levels = c("Normal", "Tumor"),
        direction = "<",
        quiet = TRUE
      )
    })

  fold_auc <- tibble(
    Fold = names(roc_list),
    AUC = as.numeric(sapply(roc_list, auc)),
    Model = model_name
  )

  roc_df <- purrr::map2_dfr(
    roc_list,
    names(roc_list),
    function(r, fold_name) {
      tibble(
        FPR = 1 - r$specificities,
        TPR = r$sensitivities,
        Fold = fold_name,
        Model = model_name
      )
    }
  )

  mean_fpr <- seq(0, 1, length.out = 101)

  mean_tpr_mat <- sapply(roc_list, function(r) {
    fpr <- 1 - r$specificities
    tpr <- r$sensitivities

    ord <- order(fpr)

    approx(
      x = fpr[ord],
      y = tpr[ord],
      xout = mean_fpr,
      ties = max,
      rule = 2
    )$y
  })

  mean_roc_df <- tibble(
    FPR = mean_fpr,
    TPR = rowMeans(mean_tpr_mat, na.rm = TRUE),
    Model = model_name
  )

  list(
    model = rf_cv,
    roc_df = roc_df,
    mean_roc_df = mean_roc_df,
    fold_auc = fold_auc,
    mean_auc = mean(fold_auc$AUC),
    sd_auc = sd(fold_auc$AUC)
  )
}

plot_roc <- function(roc_obj, title_txt) {
  ggplot() +
    geom_line(
      data = roc_obj$roc_df,
      aes(x = FPR, y = TPR, group = Fold),
      linewidth = 0.7,
      alpha = 0.45,
      color = "gray70"
    ) +
    geom_line(
      data = roc_obj$mean_roc_df,
      aes(x = FPR, y = TPR),
      linewidth = 1.4,
      color = "coral1"
    ) +
    geom_abline(linetype = "dashed", color = "gray70") +
    labs(
      title = title_txt,
      subtitle = sprintf(
        "Mean AUC = %.3f ± %.3f",
        roc_obj$mean_auc,
        roc_obj$sd_auc
      ),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid = element_blank()
    )
}

plot_confusion_matrix <- function(roc_obj, title_txt) {
  cv_preds <- roc_obj$model$pred %>%
    filter(mtry == roc_obj$model$bestTune$mtry) %>%
    mutate(
      obs = factor(obs, levels = c("Normal", "Tumor")),
      pred = factor(pred, levels = c("Normal", "Tumor"))
    )

  cm <- caret::confusionMatrix(
    data = cv_preds$pred,
    reference = cv_preds$obs,
    positive = "Tumor"
  )

  print(cm)

  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Reference", "Freq")

  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "coral1") +
    labs(
      title = title_txt,
      subtitle = "5-fold cross-validation on pan-cancer TCGA samples",
      x = "True label",
      y = "Predicted label",
      fill = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid = element_blank()
    )
}

# ------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------

tcga_df <- load_tcga_expression(tcga_dir, genes_15)

saveRDS(
  tcga_df,
  file.path(output_dir, "tcga_expression_15gene.rds")
)

roc_15 <- run_rf_cv_roc(
  data = tcga_df,
  genes = genes_15,
  model_name = "15-gene model"
)

roc_9 <- run_rf_cv_roc(
  data = tcga_df,
  genes = genes_9,
  model_name = "9-gene model"
)

saveRDS(roc_15$model, file.path(output_dir, "rf_cv_15gene.rds"))
saveRDS(roc_9$model, file.path(output_dir, "rf_cv_9gene.rds"))

# ------------------------------------------------------------------
# plots
# ------------------------------------------------------------------

p_roc_15 <- plot_roc(roc_15, "ROC Curve: 15-Gene Model")
p_roc_9  <- plot_roc(roc_9, "ROC Curve: 9-Gene Model")

p_cm_15 <- plot_confusion_matrix(roc_15, "Confusion Matrix — 15-Gene RF Model")
p_cm_9  <- plot_confusion_matrix(roc_9, "Confusion Matrix — 9-Gene RF Model")
