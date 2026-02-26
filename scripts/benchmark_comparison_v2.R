#!/usr/bin/env Rscript
# IlluMeta Benchmark Comparison Script v2
# Fixed threading issues and ChAMP sample sheet

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

# Configuration
args <- commandArgs(trailingOnly = TRUE)
idat_dir <- if (length(args) >= 1) args[1] else "projects/GSE115508/idat"
output_dir <- if (length(args) >= 2) args[2] else "benchmarks/tool_comparison"
n_samples <- if (length(args) >= 3) as.integer(args[3]) else 20

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=" , rep("=", 60), "\n", sep = "")
cat("IlluMeta Benchmark Comparison v2\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("IDAT directory:", idat_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("N samples (subset):", n_samples, "\n\n")

# Find IDAT files
idat_files <- list.files(idat_dir, pattern = "_Grn.idat$", full.names = TRUE)
basenames <- gsub("_Grn.idat$", "", idat_files)

if (length(basenames) == 0) {
  gz_files <- list.files(idat_dir, pattern = "_Grn.idat.gz$", full.names = TRUE)
  if (length(gz_files) > 0) {
    stop("Found only .idat.gz files in ", idat_dir, ". Please gunzip *.idat.gz before running.")
  }
  stop("No IDAT files found in ", idat_dir)
}

cat("Total IDAT pairs found:", length(basenames), "\n")

# Subset for benchmark
if (length(basenames) > n_samples) {
  set.seed(42)
  basenames <- sample(basenames, n_samples)
  cat("Using subset of", n_samples, "samples\n")
}

# Build a clean subset directory to keep ChAMP/minfi consistent.
subset_dir <- file.path(output_dir, "idat_subset")
if (dir.exists(subset_dir)) {
  unlink(subset_dir, recursive = TRUE, force = TRUE)
}
dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)
subset_basenames <- file.path(subset_dir, basename(basenames))
for (i in seq_along(basenames)) {
  for (suffix in c("_Grn.idat", "_Red.idat")) {
    src <- paste0(basenames[i], suffix)
    dst <- paste0(subset_basenames[i], suffix)
    if (!file.exists(dst)) {
      ok <- file.symlink(src, dst)
      if (!isTRUE(ok)) {
        file.copy(src, dst, overwrite = TRUE)
      }
    }
  }
}

# Assign random groups
n <- length(basenames)
groups <- rep(c("Control", "Test"), length.out = n)

results <- list()
array_type_detected <- NA_character_

# ============================================================
# 1. MINFI MANUAL WORKFLOW (using preprocessNoob - no threading)
# ============================================================
cat("\n", rep("-", 60), "\n", sep = "")
cat("1. MINFI MANUAL WORKFLOW\n")
cat(rep("-", 60), "\n", sep = "")

minfi_start <- Sys.time()

tryCatch({
  cat("  Reading IDAT files...\n")
  rgSet <- read.metharray(subset_basenames, verbose = FALSE)
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

  # Filter probes with detection p > 0.05
  cat("  Filtering probes...\n")
  detP <- detectionP(rgSet)
  keep <- rowSums(detP > 0.05) < ncol(detP) * 0.1
  betas_filt <- betas[keep, ]
  betas_filt <- betas_filt[complete.cases(betas_filt), ]

  # Convert to M-values
  cat("  Running limma DMP analysis...\n")
  mvals <- log2(betas_filt / (1 - betas_filt))
  mvals[is.infinite(mvals)] <- NA
  mvals <- mvals[complete.cases(mvals), ]

  # Limma analysis
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- c("Control", "Test")

  fit <- lmFit(mvals, design)
  contrast_matrix <- makeContrasts(Test - Control, levels = design)
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

  # Save top DMPs for comparison
  write.csv(dmps_minfi, file.path(output_dir, "minfi_dmps.csv"), row.names = FALSE)

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  results$minfi <<- list(error = e$message)
})

# Cleanup
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

    # Create sample sheet for ChAMP
    sample_names <- basename(subset_basenames)
    parse_sentrix <- function(name) {
      m <- regmatches(name, regexec(".*_(\\d+)_(R\\d+C\\d+)$", name))[[1]]
      if (length(m) >= 3) return(list(id = m[2], pos = m[3]))
      parts <- strsplit(name, "_")[[1]]
      if (length(parts) >= 3) {
        return(list(id = parts[length(parts) - 1], pos = parts[length(parts)]))
      }
      return(list(id = NA_character_, pos = NA_character_))
    }
    sentrix <- lapply(sample_names, parse_sentrix)
    sentrix_id <- vapply(sentrix, function(x) x$id, character(1))
    sentrix_pos <- vapply(sentrix, function(x) x$pos, character(1))

    sample_sheet <- data.frame(
      Sample_Name = sample_names,
      Sentrix_ID = sentrix_id,
      Sentrix_Position = sentrix_pos,
      Basename = basenames,
      Sample_Group = groups,
      stringsAsFactors = FALSE
    )

    sheet_path <- file.path(subset_dir, "SampleSheet.csv")
    unlink(file.path(subset_dir, "sample_sheet.csv"))
    write.csv(sample_sheet, sheet_path, row.names = FALSE)
    cat("  Created sample sheet:", sheet_path, "\n")

    cat("  Loading data with champ.load()...\n")

    # Determine array type
    array_type <- if (!is.na(array_type_detected)) array_type_detected else "EPIC"

    myLoad <- suppressMessages(champ.load(
      directory = subset_dir,
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
    myDMP <- suppressMessages(champ.DMP(
      beta = myNorm,
      pheno = myLoad$pd$Sample_Group,
      adjPVal = 0.05,
      adjust.method = "BH"
    ))

    champ_end <- Sys.time()

    dmps_champ <- myDMP[[1]]

    results$champ <- list(
      time = as.numeric(difftime(champ_end, champ_start, units = "secs")),
      n_probes_input = nrow(myLoad$beta),
      n_probes_analyzed = nrow(myNorm),
      n_dmps_sig = sum(dmps_champ$adj.P.Val < 0.05, na.rm = TRUE),
      n_dmps_nominal = sum(dmps_champ$P.Value < 0.05, na.rm = TRUE),
      top_dmps = head(dmps_champ, 1000)
    )

    cat("  Completed in", round(results$champ$time, 1), "seconds\n")
    cat("  Probes analyzed:", results$champ$n_probes_analyzed, "\n")
    cat("  Significant DMPs (FDR<0.05):", results$champ$n_dmps_sig, "\n")
    cat("  Nominal DMPs (p<0.05):", results$champ$n_dmps_nominal, "\n")

    # Save top DMPs
    write.csv(dmps_champ, file.path(output_dir, "champ_dmps.csv"), row.names = FALSE)

    # Clean up sample sheet
    unlink(sheet_path)

  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    results$champ <<- list(error = e$message)
  })

} else {
  cat("  ChAMP not installed, skipping...\n")
  results$champ <- list(error = "Not installed")
}

# ============================================================
# 3. ILLUMETA TIMING (from existing results)
# ============================================================
cat("\n", rep("-", 60), "\n", sep = "")
cat("3. ILLUMETA (from existing analysis)\n")
cat(rep("-", 60), "\n", sep = "")

# Check for existing IlluMeta log
log_file <- file.path(dirname(idat_dir), "CRF21_test_results", "run_full.log")
if (file.exists(log_file)) {
  log_content <- readLines(log_file)

  # Extract timing from log (look for timestamp patterns)
  time_lines <- grep("^\\[\\d{2}:\\d{2}:\\d{2}\\]", log_content, value = TRUE)
  if (length(time_lines) >= 2) {
    start_time <- time_lines[1]
    end_time <- time_lines[length(time_lines)]

    # Parse times
    start_h <- as.numeric(substr(start_time, 2, 3))
    start_m <- as.numeric(substr(start_time, 5, 6))
    start_s <- as.numeric(substr(start_time, 8, 9))
    end_h <- as.numeric(substr(end_time, 2, 3))
    end_m <- as.numeric(substr(end_time, 5, 6))
    end_s <- as.numeric(substr(end_time, 8, 9))

    illumeta_time <- ((end_h - start_h) * 3600 + (end_m - start_m) * 60 + (end_s - start_s))
    if (illumeta_time < 0) illumeta_time <- illumeta_time + 86400  # Handle day rollover

    cat("  Estimated from log file\n")
    cat("  Total runtime:", round(illumeta_time, 1), "seconds\n")

    results$illumeta <- list(
      time = illumeta_time,
      n_samples = n_samples,
      note = "Full pipeline including QC, normalization, batch correction, DMP, DMR, CRF"
    )
  }
} else {
  cat("  No IlluMeta log found, will estimate from benchmark\n")
  results$illumeta <- list(time = NA, note = "Log not found")
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("BENCHMARK SUMMARY (n=", n_samples, " samples)\n", sep = "")
cat(rep("=", 60), "\n\n", sep = "")

# Create summary table
summary_df <- data.frame(
  Tool = c("minfi (manual)", "ChAMP", "IlluMeta"),
  Time_sec = c(
    if (!is.null(results$minfi) && is.null(results$minfi$error)) round(results$minfi$time, 1) else NA,
    if (!is.null(results$champ) && is.null(results$champ$error)) round(results$champ$time, 1) else NA,
    if (!is.null(results$illumeta)) results$illumeta$time else NA
  ),
  Probes_Analyzed = c(
    if (!is.null(results$minfi) && is.null(results$minfi$error)) results$minfi$n_probes_analyzed else NA,
    if (!is.null(results$champ) && is.null(results$champ$error)) results$champ$n_probes_analyzed else NA,
    NA  # IlluMeta varies by pipeline
  ),
  DMPs_FDR_0.05 = c(
    if (!is.null(results$minfi) && is.null(results$minfi$error)) results$minfi$n_dmps_sig else NA,
    if (!is.null(results$champ) && is.null(results$champ$error)) results$champ$n_dmps_sig else NA,
    NA
  ),
  Auto_QC = c("No", "Yes", "Yes"),
  Auto_Batch = c("No", "Yes", "Yes"),
  Robustness_CRF = c("No", "No", "Yes"),
  Dashboard = c("No", "No", "Yes"),
  stringsAsFactors = FALSE
)

print(summary_df)

# Code complexity
cat("\n")
cat("CODE COMPLEXITY:\n")
cat("  minfi (manual): ~45-60 lines of R code, 8+ functions to learn\n")
cat("  ChAMP:          ~25-30 lines of R code, 5+ functions to learn\n")
cat("  IlluMeta:       1 command (python illumeta.py project_dir)\n")

# Feature comparison
cat("\n")
cat("UNIQUE ILLUMETA FEATURES:\n")
cat("  - Automatic GEO data download\n")
cat("  - Automatic group inference from metadata\n")
cat("  - CRF robustness assessment (MMC, NCS, RSS)\n")
cat("  - Sample-adaptive analysis tiers\n")
cat("  - Interactive HTML dashboard with verdict system\n")
cat("  - Cross-pipeline intersection (minfi + SeSAMe)\n")

# Save results
cat("\nSaving results to:", output_dir, "\n")
write.csv(summary_df, file.path(output_dir, "benchmark_summary.csv"), row.names = FALSE)
saveRDS(results, file.path(output_dir, "benchmark_results.rds"))

# DMP concordance (if both tools ran successfully)
if (!is.null(results$minfi$top_dmps) && !is.null(results$champ$top_dmps)) {
  cat("\nDMP CONCORDANCE (minfi vs ChAMP):\n")
  minfi_top <- head(results$minfi$top_dmps$CpG, 500)
  champ_top <- head(rownames(results$champ$top_dmps), 500)
  overlap <- length(intersect(minfi_top, champ_top))
  cat("  Top 500 overlap:", overlap, "(", round(overlap/500*100, 1), "%)\n")
}

cat("\nDone!\n")
