# TCGA preprocessing, training, and cross-validation of 15-gene random forest classifier
#
# Input:
# - dds_*_matched_group.rds files
#
#
# Notes:
# - Edit 'tcga_dir' and 'output_dir' below before running
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

tcga_dir   <- "path-to-tcga-dds-files"

if (!dir.exists(tcga_dir)) {
  stop("Please set 'tcga_dir' to the folder containing dds_*_matched_group.rds files.")
}


# ------------------------------------------------------------------
# Define 15-gene panel
# ------------------------------------------------------------------

ttll_genes <- c(
  "TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL13"
)

ccp_genes <- c(
  "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5", "AGTPBP1"
)

feature_genes <- unique(c(ttll_genes, ccp_genes))

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

aggregate_dup_symbols <- function(mat, symbols, mode = "mean") {
  ok <- !is.na(symbols) & symbols != ""
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]

  if (mode == "sum") {
    out <- rowsum(mat, group = symbols)
  } else {
    sp <- split(seq_len(nrow(mat)), symbols)
    out <- vapply(
      sp,
      function(ix) colMeans(mat[ix, , drop = FALSE]),
      numeric(ncol(mat))
    )
    out <- t(out)
    rownames(out) <- names(sp)
  }

  out
}

plot_cm <- function(pred, truth, title, high_col = "#00B2EE") {
  pred <- factor(pred, levels = levels(truth))
  cm <- confusionMatrix(pred, truth)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Reference", "Freq")

  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = high_col) +
    labs(title = title, x = "Actual", y = "Predicted") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# ==================================================================
# 1. PREPROCESS TCGA DATA
# ==================================================================

dds_files <- list.files(
  tcga_dir,
  pattern = "dds_.*_matched_group\\.rds$",
  full.names = TRUE
)

if (length(dds_files) == 0) {
  stop("No TCGA DESeq2 files found in 'tcga_dir'.")
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

  expr_sym <- aggregate_dup_symbols(expr_mat, symbols, mode = "mean")

  matched_genes <- intersect(feature_genes, rownames(expr_sym))
  if (length(matched_genes) == 0) next

  expr_df <- t(expr_sym[matched_genes, , drop = FALSE]) %>% as.data.frame()
  expr_df$sample <- rownames(expr_df)

  cond <- tryCatch({
    as.character(colData(dds)[expr_df$sample, "condition"])
  }, error = function(e) NULL)

  if (is.null(cond)) next

  expr_df$condition <- cond
  expr_df$project <- project
  tcga_data[[length(tcga_data) + 1]] <- expr_df
}

tcga_df <- bind_rows(tcga_data) %>%
  dplyr::select(any_of(feature_genes), condition, sample, project) %>%
  filter(if_all(all_of(feature_genes), ~ !is.na(.))) %>%
  mutate(
    condition_bin = factor(
      ifelse(condition == "Primary Tumor", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    )
  )

if (nrow(tcga_df) == 0) {
  stop("No TCGA samples available after preprocessing.")
}

saveRDS(tcga_df, file.path(output_dir, "tcga_15gene_expression.rds"))

# ==================================================================
# 2. TRAIN THE 15-GENE MODEL
# ==================================================================

x_tcga <- tcga_df[, feature_genes, drop = FALSE]
y_tcga <- tcga_df$condition_bin

set.seed(123)
rf_15 <- randomForest(
  x = x_tcga,
  y = y_tcga,
  ntree = 500,
  importance = TRUE
)

importance_15 <- importance(rf_15, type = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  arrange(desc(MeanDecreaseAccuracy))

print(importance_15)

saveRDS(rf_15, file.path(output_dir, "rf_15gene_tcga_model.rds"))

# ==================================================================
# 3. CROSS-VALIDATION ROC ANALYSIS — 15-GENE MODEL ONLY
# ==================================================================

genes_15 <- feature_genes

run_rf_cv_roc <- function(data, genes, model_name, seed = 123) {
  
  genes_present <- intersect(genes, colnames(data))
  missing_genes <- setdiff(genes, genes_present)
  
  if (length(missing_genes) > 0) {
    message(model_name, " missing genes: ", paste(missing_genes, collapse = ", "))
  }
  
  df <- data %>%
    dplyr::select(all_of(genes_present), condition_bin) %>%
    filter(if_all(all_of(genes_present), ~ !is.na(.))) %>%
    mutate(condition_bin = factor(condition_bin, levels = c("Normal", "Tumor")))
  
  x <- df[, genes_present, drop = FALSE]
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
    trControl = ctrl,
    metric = "ROC",
    ntree = 500
  )
  
  best_mtry <- rf_cv$bestTune$mtry
  
  cv_preds <- rf_cv$pred %>%
    filter(mtry == best_mtry) %>%
    mutate(Resample = factor(Resample))
  
  roc_list <- cv_preds %>%
    split(.$Resample) %>%
    lapply(function(d) {
      pROC::roc(
        response = d$obs,
        predictor = d$Tumor,
        levels = c("Normal", "Tumor"),
        direction = "<",
        quiet = TRUE
      )
    })
  
  fold_auc <- tibble(
    Fold = names(roc_list),
    AUC = as.numeric(sapply(roc_list, pROC::auc)),
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
    fpr <- c(0, 1 - r$specificities, 1)
    tpr <- c(0, r$sensitivities, 1)
    
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
    sd_auc = sd(fold_auc$AUC),
    genes_used = genes_present
  )
}

plot_fold_rocs_tcga_style <- function(roc_obj, title_txt) {
  
  df_mean <- roc_obj$mean_roc_df %>%
    transmute(
      fpr = FPR,
      tpr = TPR,
      model = "TCGA CV"
    )
  
  ggplot() +
    geom_line(
      data = roc_obj$roc_df,
      aes(x = FPR, y = TPR, group = Fold),
      linewidth = 0.7,
      alpha = 0.45,
      color = "gray70"
    ) +
    geom_line(
      data = df_mean,
      aes(x = fpr, y = tpr, color = model),
      linewidth = 1.4
    ) +
    geom_abline(linetype = "dashed", color = "gray70") +
    scale_color_manual(
      values = c("TCGA CV" = "coral1"),
      breaks = c("TCGA CV")
    ) +
    labs(
      title = title_txt,
      subtitle = sprintf(
        "Mean AUC = %.3f ± %.3f",
        roc_obj$mean_auc,
        roc_obj$sd_auc
      ),
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = "Dataset"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid = element_blank()
    )
}

roc_15 <- run_rf_cv_roc(
  data = tcga_df,
  genes = genes_15,
  model_name = "15-gene model",
  seed = 123
)

saveRDS(roc_15$model, file.path(output_dir, "rf_15gene_tcga_cv.rds"))

p_roc_15 <- plot_fold_rocs_tcga_style(
  roc_15,
  "ROC Curve: 15-Gene Model"
)

print(p_roc_15)

ggsave(
  filename = file.path(output_dir, "rf_15gene_foldwise_ROC.pdf"),
  plot = p_roc_15,
  width = 5.5,
  height = 5
)

ggsave(
  filename = file.path(output_dir, "rf_15gene_foldwise_ROC.png"),
  plot = p_roc_15,
  width = 5.5,
  height = 5,
  dpi = 300
)

# ==================================================================
# 4. 15-GENE POOLED CROSS-VALIDATION CONFUSION MATRIX
# ==================================================================

cv_preds_15 <- roc_15$model$pred %>%
  dplyr::filter(mtry == roc_15$model$bestTune$mtry) %>%
  dplyr::mutate(
    obs = factor(obs, levels = c("Normal", "Tumor")),
    pred = factor(pred, levels = c("Normal", "Tumor"))
  )

cm_15 <- caret::confusionMatrix(
  data = cv_preds_15$pred,
  reference = cv_preds_15$obs,
  positive = "Tumor"
)

print(cm_15)

cm_15_df <- as.data.frame(cm_15$table)
colnames(cm_15_df) <- c("Prediction", "Reference", "Freq")

p_cm_15_coral <- ggplot(cm_15_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Freq), size = 6, fontface = "bold") +
  scale_fill_gradient(
    low = "white",
    high = "coral1"
  ) +
  labs(
    title = "Confusion Matrix — 15-Gene RF Model",
    subtitle = "Pooled out-of-fold predictions from 5-fold cross-validation",
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
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_cm_15_coral)

