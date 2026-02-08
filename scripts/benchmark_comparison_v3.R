#!/usr/bin/env Rscript
# IlluMeta Benchmark Comparison Script v3
# Uses a provided SampleSheet.csv for minfi/ChAMP fairness and group labels.

suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(data.table)
})

# Force single-threaded BLAS/OMP to avoid pthread errors.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

args <- commandArgs(trailingOnly = TRUE)
idat_dir <- if (length(args) >= 1) args[1] else "projects/GSE115508/idat"
output_dir <- if (length(args) >= 2) args[2] else "benchmarks/tool_comparison"
group_con <- if (length(args) >= 3) args[3] else ""
group_test <- if (length(args) >= 4) args[4] else ""
max_per_group <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else NA_integer_
sample_sheet <- if (length(args) >= 6 && nzchar(args[6])) args[6] else file.path(idat_dir, "SampleSheet.csv")
seed <- if (length(args) >= 7 && nzchar(args[7])) as.integer(args[7]) else 42

run_limma_dmp <- function(beta, groups) {
  groups <- as.factor(groups)
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- make.names(levels(groups))
  if (ncol(design) < 2) return(NULL)
  mvals <- log2(beta / (1 - beta))
  mvals[is.infinite(mvals)] <- NA
  mvals <- mvals[complete.cases(mvals), ]
  fit <- lmFit(mvals, design)
  contrast_str <- paste0(colnames(design)[2], "-", colnames(design)[1])
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, number = Inf, sort.by = "P")
  res$CpG <- rownames(res)
  res
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=", rep("=", 60), "\n", sep = "")
cat("IlluMeta Benchmark Comparison v3\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("IDAT directory:", idat_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Sample sheet:", sample_sheet, "\n")
cat("Group control:", group_con, "\n")
cat("Group test:", group_test, "\n\n")

# Load sample sheet
if (!file.exists(sample_sheet)) {
  stop("SampleSheet.csv not found: ", sample_sheet)
}
ss <- data.table::fread(sample_sheet)
if (!"Sample_Group" %in% colnames(ss)) {
  stop("SampleSheet missing Sample_Group column")
}

# Filter to selected groups if provided
if (nzchar(group_con) && nzchar(group_test)) {
  con_lower <- tolower(group_con)
  test_lower <- tolower(group_test)
  ss <- ss[tolower(Sample_Group) %in% c(con_lower, test_lower)]
  if (nrow(ss) == 0) stop("No samples match requested groups")
}

# Subset per group if requested
if (!is.na(max_per_group) && max_per_group > 0) {
  set.seed(seed)
  keep_rows <- ss[, .I[sample(.N, min(.N, max_per_group))], by = Sample_Group]$V1
  ss <- ss[keep_rows]
}

# Build basenames for minfi
if ("Basename" %in% colnames(ss) && all(nzchar(ss$Basename))) {
  basenames <- ss$Basename
} else if ("Sample_Name" %in% colnames(ss)) {
  basenames <- file.path(idat_dir, ss$Sample_Name)
} else {
  stop("SampleSheet must include Basename or Sample_Name")
}

# Groups for minfi/limma
sample_groups <- ss$Sample_Group
if (length(unique(sample_groups)) < 2) {
  stop("Need at least two groups for comparison")
}

results <- list()
array_type_detected <- NA_character_

# ============================================================
# 1. MINFI MANUAL WORKFLOW
# ============================================================
cat("\n", rep("-", 60), "\n", sep = "")
cat("1. MINFI MANUAL WORKFLOW\n")
cat(rep("-", 60), "\n", sep = "")

minfi_start <- Sys.time()

tryCatch({
  cat("  Reading IDAT files...\n")
  rgSet <- read.metharray(basenames, verbose = FALSE)

  ann <- tryCatch(minfi::annotation(rgSet), error = function(e) NULL)
  if (!is.null(ann) && length(ann) >= 1) {
    arr <- tolower(ann[1])
    if (grepl("450", arr)) array_type_detected <- "450K"
    if (grepl("epic", arr)) array_type_detected <- "EPIC"
  }

  cat("  Preprocessing (NOOB normalization)...\n")
  mSet <- preprocessNoob(rgSet)

  cat("  Getting beta/M values...\n")
  betas <- getBeta(mSet)

  cat("  Filtering probes...\n")
  detP <- detectionP(rgSet)
  keep <- rowSums(detP > 0.05) < ncol(detP) * 0.1
  betas_filt <- betas[keep, ]
  betas_filt <- betas_filt[complete.cases(betas_filt), ]

  cat("  Running limma DMP analysis...\n")
  mvals <- log2(betas_filt / (1 - betas_filt))
  mvals[is.infinite(mvals)] <- NA
  mvals <- mvals[complete.cases(mvals), ]

  design <- model.matrix(~ 0 + factor(sample_groups))
  colnames(design) <- make.names(levels(factor(sample_groups)))

  if (ncol(design) < 2) stop("Not enough groups after filtering")

  fit <- lmFit(mvals, design)
  # Explicit contrast direction: Test - Control (matches IlluMeta default).
  con_name <- make.names(group_con)
  test_name <- make.names(group_test)
  if (!(con_name %in% colnames(design)) || !(test_name %in% colnames(design))) {
    # Fallback to design order: second level minus first level.
    con_name <- colnames(design)[1]
    test_name <- colnames(design)[2]
  }
  contrast_matrix <- makeContrasts(contrasts = paste0(test_name, "-", con_name), levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  dmps_minfi <- topTable(fit2, number = Inf, sort.by = "P")
  dmps_minfi$CpG <- rownames(dmps_minfi)

  minfi_end <- Sys.time()

  results$minfi <- list(
    time = as.numeric(difftime(minfi_end, minfi_start, units = "secs")),
    n_probes_input = nrow(betas),
    n_probes_analyzed = nrow(mvals),
    n_dmps_sig = sum(dmps_minfi$adj.P.Val < 0.05, na.rm = TRUE),
    n_dmps_nominal = sum(dmps_minfi$P.Value < 0.05, na.rm = TRUE),
    top_dmps = head(dmps_minfi, 1000)
  )

  cat("  Completed in", round(results$minfi$time, 1), "seconds\n")
  cat("  Probes analyzed:", results$minfi$n_probes_analyzed, "\n")
  cat("  Significant DMPs (FDR<0.05):", results$minfi$n_dmps_sig, "\n")
  cat("  Nominal DMPs (p<0.05):", results$minfi$n_dmps_nominal, "\n")

  write.csv(dmps_minfi, file.path(output_dir, "minfi_dmps.csv"), row.names = FALSE)

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  results$minfi <<- list(error = e$message)
})

suppressWarnings({
  rm(rgSet, mSet, betas, betas_filt, mvals, fit, fit2, dmps_minfi)
})
gc(verbose = FALSE)

# ============================================================
# 2. ChAMP WORKFLOW
# ============================================================
cat("\n", rep("-", 60), "\n", sep = "")
cat("2. ChAMP WORKFLOW\n")
cat(rep("-", 60), "\n", sep = "")

champ_available <- requireNamespace("ChAMP", quietly = TRUE)

if (champ_available) {
  champ_start <- Sys.time()

  tryCatch({
    suppressPackageStartupMessages(library(ChAMP))

    # Ensure SampleSheet.csv lives inside idat_dir for champ.load
    champ_sheet <- file.path(idat_dir, "SampleSheet.csv")
    if (normalizePath(sample_sheet, winslash = "/", mustWork = TRUE) !=
        normalizePath(champ_sheet, winslash = "/", mustWork = FALSE)) {
      file.copy(sample_sheet, champ_sheet, overwrite = TRUE)
    }

    array_type <- if (!is.na(array_type_detected)) array_type_detected else "EPIC"

    cat("  Loading data with champ.load()...\n")
    myLoad <- suppressMessages(champ.load(
      directory = idat_dir,
      method = "minfi",
      methValue = "B",
      filterXY = TRUE,
      filterDetP = TRUE,
      filterBeads = FALSE,
      filterNoCG = TRUE,
      filterSNPs = TRUE,
      arraytype = array_type
    ))

    cat("  Normalization with champ.norm() (BMIQ)...\n")
    myNorm <- suppressMessages(champ.norm(
      beta = myLoad$beta,
      method = "BMIQ",
      plotBMIQ = FALSE,
      cores = 1
    ))

    cat("  DMP analysis with champ.DMP()...\n")
    myDMP <- tryCatch(
      suppressMessages(champ.DMP(
        beta = myNorm,
        pheno = myLoad$pd$Sample_Group,
        adjPVal = 0.05,
        adjust.method = "BH"
      )),
      error = function(e) e
    )

    dmps_champ <- NULL
    champ_note <- ""
    if (inherits(myDMP, "error")) {
      champ_note <- paste0("champ.DMP failed: ", myDMP$message, "; fallback to limma on BMIQ")
    } else {
      dmps_champ <- myDMP[[1]]
      if (is.null(dmps_champ) || nrow(dmps_champ) == 0) {
        champ_note <- "champ.DMP returned no significant CpGs; fallback to limma on BMIQ"
        dmps_champ <- NULL
      }
    }

    if (is.null(dmps_champ)) {
      dmps_champ <- run_limma_dmp(myNorm, myLoad$pd$Sample_Group)
    }

    champ_end <- Sys.time()

    results$champ <- list(
      time = as.numeric(difftime(champ_end, champ_start, units = "secs")),
      n_probes_input = nrow(myLoad$beta),
      n_probes_analyzed = nrow(myNorm),
      n_dmps_sig = if (!is.null(dmps_champ) && nrow(dmps_champ) > 0) sum(dmps_champ$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      n_dmps_nominal = if (!is.null(dmps_champ) && nrow(dmps_champ) > 0) sum(dmps_champ$P.Value < 0.05, na.rm = TRUE) else 0,
      top_dmps = if (!is.null(dmps_champ) && nrow(dmps_champ) > 0) head(dmps_champ, 1000) else NULL,
      note = champ_note
    )

    cat("  Completed in", round(results$champ$time, 1), "seconds\n")
    cat("  Probes analyzed:", results$champ$n_probes_analyzed, "\n")
    cat("  Significant DMPs (FDR<0.05):", results$champ$n_dmps_sig, "\n")
    cat("  Nominal DMPs (p<0.05):", results$champ$n_dmps_nominal, "\n")

    if (!is.null(dmps_champ) && nrow(dmps_champ) > 0) {
      write.csv(dmps_champ, file.path(output_dir, "champ_dmps.csv"), row.names = FALSE)
    }

  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    results$champ <<- list(error = e$message)
  })

} else {
  cat("  ChAMP not installed, skipping...\n")
  results$champ <- list(error = "Not installed")
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("BENCHMARK SUMMARY\n")
cat(rep("=", 60), "\n\n", sep = "")

summary_df <- data.frame(
  Tool = c("minfi (manual)", "ChAMP"),
  Time_sec = c(
    if (!is.null(results$minfi) && is.null(results$minfi$error)) round(results$minfi$time, 1) else NA,
    if (!is.null(results$champ) && is.null(results$champ$error)) round(results$champ$time, 1) else NA
  ),
  Probes_Analyzed = c(
    if (!is.null(results$minfi) && is.null(results$minfi$error)) results$minfi$n_probes_analyzed else NA,
    if (!is.null(results$champ) && is.null(results$champ$error)) results$champ$n_probes_analyzed else NA
  ),
  DMPs_FDR_0.05 = c(
    if (!is.null(results$minfi) && is.null(results$minfi$error)) results$minfi$n_dmps_sig else NA,
    if (!is.null(results$champ) && is.null(results$champ$error)) results$champ$n_dmps_sig else NA
  ),
  Auto_QC = c("No", "Yes"),
  Auto_Batch = c("No", "Yes"),
  Note = c("", if (!is.null(results$champ) && !is.null(results$champ$note)) results$champ$note else ""),
  stringsAsFactors = FALSE
)

print(summary_df)

cat("\nSaving results to:", output_dir, "\n")
write.csv(summary_df, file.path(output_dir, "benchmark_summary.csv"), row.names = FALSE)
saveRDS(results, file.path(output_dir, "benchmark_results.rds"))

# DMP concordance (if both tools ran successfully)
if (!is.null(results$minfi$top_dmps) && !is.null(results$champ$top_dmps)) {
  cat("\nDMP CONCORDANCE (minfi vs ChAMP):\n")
  minfi_top <- head(results$minfi$top_dmps$CpG, 500)
  champ_top <- head(rownames(results$champ$top_dmps), 500)
  overlap <- length(intersect(minfi_top, champ_top))
  cat("  Top 500 overlap:", overlap, "(", round(overlap/500*100, 1), "%)\n", sep = "")
}

cat("\nDone!\n")
