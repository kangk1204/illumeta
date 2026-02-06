#!/usr/bin/env Rscript
# IlluMeta Benchmark Comparison Script
# Compares IlluMeta vs ChAMP vs minfi manual workflow
# For Bioinformatics Application Note

suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(data.table)
})

# Configuration
args <- commandArgs(trailingOnly = TRUE)
idat_dir <- if (length(args) >= 1) args[1] else "projects/GSE121633/idat"
output_dir <- if (length(args) >= 2) args[2] else "benchmarks/tool_comparison"
n_samples <- if (length(args) >= 3) as.integer(args[3]) else 20  # Subset for quick test

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=" , rep("=", 60), "\n", sep = "")
cat("IlluMeta Benchmark Comparison\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("IDAT directory:", idat_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("N samples (subset):", n_samples, "\n\n")

# Find IDAT files
idat_files <- list.files(idat_dir, pattern = "_Grn.idat$", full.names = TRUE)
basenames <- gsub("_Grn.idat$", "", basename(idat_files))

if (length(basenames) == 0) {
  stop("No IDAT files found in ", idat_dir)
}

cat("Total IDAT pairs found:", length(basenames), "\n")

# Subset for benchmark
if (length(basenames) > n_samples) {
  set.seed(42)
  basenames <- sample(basenames, n_samples)
  cat("Using subset of", n_samples, "samples\n")
}

# Create targets dataframe
targets <- data.frame(
  Basename = file.path(idat_dir, basenames),
  Sample_Name = basenames,
  stringsAsFactors = FALSE
)

# Assign random groups for benchmark (half control, half test)
n <- nrow(targets)
targets$Group <- rep(c("Control", "Test"), length.out = n)

results <- list()

# ============================================================
# 1. MINFI MANUAL WORKFLOW
# ============================================================
cat("\n", rep("-", 60), "\n", sep = "")
cat("1. MINFI MANUAL WORKFLOW\n")
cat(rep("-", 60), "\n", sep = "")

minfi_start <- Sys.time()
minfi_mem_start <- gc(reset = TRUE)

tryCatch({
  # Read data
  cat("  Reading IDAT files...\n")
  rgSet <- read.metharray(targets$Basename, verbose = FALSE)

  # Preprocess
  cat("  Preprocessing (functional normalization)...\n")
  grSet <- preprocessFunnorm(rgSet)

  # Get beta values
  cat("  Extracting beta values...\n")
  betas <- getBeta(grSet)

  # Filter probes
  cat("  Filtering probes...\n")
  # Remove probes with detection p > 0.05
  detP <- detectionP(rgSet)
  keep <- rowSums(detP > 0.05) < ncol(detP) * 0.1
  betas_filt <- betas[keep, ]

  # Remove probes with NA
  betas_filt <- betas_filt[complete.cases(betas_filt), ]

  # DMP analysis
  cat("  Running limma DMP analysis...\n")
  design <- model.matrix(~ 0 + Group, data = targets)
  colnames(design) <- gsub("Group", "", colnames(design))

  # Convert to M-values
  mvals <- log2(betas_filt / (1 - betas_filt))
  mvals[is.infinite(mvals)] <- NA
  mvals <- mvals[complete.cases(mvals), ]

  fit <- lmFit(mvals, design)
  contrast_matrix <- makeContrasts(Test - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  dmps_minfi <- topTable(fit2, number = Inf, sort.by = "P")
  dmps_minfi$CpG <- rownames(dmps_minfi)

  minfi_end <- Sys.time()
  minfi_mem_end <- gc()

  results$minfi <- list(
    time = as.numeric(difftime(minfi_end, minfi_start, units = "secs")),
    memory_mb = sum(minfi_mem_end[, 2] - minfi_mem_start[, 2]),
    n_probes_input = nrow(betas),
    n_probes_analyzed = nrow(mvals),
    n_dmps_sig = sum(dmps_minfi$adj.P.Val < 0.05, na.rm = TRUE),
    n_dmps_nominal = sum(dmps_minfi$P.Value < 0.05, na.rm = TRUE),
    top_dmps = head(dmps_minfi, 100)
  )

  cat("  Completed in", round(results$minfi$time, 1), "seconds\n")
  cat("  Significant DMPs (FDR<0.05):", results$minfi$n_dmps_sig, "\n")

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  results$minfi <<- list(error = e$message)
})

# Cleanup
rm(rgSet, grSet, betas, betas_filt, mvals, fit, fit2)
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
  champ_mem_start <- gc(reset = TRUE)

  tryCatch({
    library(ChAMP)

    # ChAMP requires specific directory structure
    # Create sample sheet
    sample_sheet <- data.frame(
      Sample_Name = targets$Sample_Name,
      Sample_Group = targets$Group,
      Pool_ID = NA,
      Project = "Benchmark",
      Sample_Plate = "Plate1",
      Sample_Well = paste0("A", 1:nrow(targets)),
      Basename = targets$Basename,
      stringsAsFactors = FALSE
    )

    sheet_path <- file.path(output_dir, "champ_sample_sheet.csv")
    write.csv(sample_sheet, sheet_path, row.names = FALSE)

    cat("  Loading data with champ.load()...\n")
    myLoad <- champ.load(
      directory = idat_dir,
      method = "minfi",
      methValue = "B",
      filterXY = TRUE,
      filterDetP = TRUE,
      filterBeads = FALSE,
      filterNoCG = TRUE,
      filterSNPs = TRUE,
      arraytype = "450K"
    )

    cat("  Normalization with champ.norm()...\n")
    myNorm <- champ.norm(beta = myLoad$beta, method = "BMIQ", cores = 1)

    cat("  DMP analysis with champ.DMP()...\n")
    myDMP <- champ.DMP(
      beta = myNorm,
      pheno = myLoad$pd$Sample_Group,
      adjPVal = 0.05,
      adjust.method = "BH"
    )

    champ_end <- Sys.time()
    champ_mem_end <- gc()

    dmps_champ <- myDMP[[1]]

    results$champ <- list(
      time = as.numeric(difftime(champ_end, champ_start, units = "secs")),
      memory_mb = sum(champ_mem_end[, 2] - champ_mem_start[, 2]),
      n_probes_input = nrow(myLoad$beta),
      n_probes_analyzed = nrow(myNorm),
      n_dmps_sig = sum(dmps_champ$adj.P.Val < 0.05, na.rm = TRUE),
      n_dmps_nominal = sum(dmps_champ$P.Value < 0.05, na.rm = TRUE),
      top_dmps = head(dmps_champ, 100)
    )

    cat("  Completed in", round(results$champ$time, 1), "seconds\n")
    cat("  Significant DMPs (FDR<0.05):", results$champ$n_dmps_sig, "\n")

  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    results$champ <<- list(error = e$message)
  })

} else {
  cat("  ChAMP not installed, skipping...\n")
  results$champ <- list(error = "Not installed")
}

# ============================================================
# 3. CODE COMPLEXITY COMPARISON
# ============================================================
cat("\n", rep("-", 60), "\n", sep = "")
cat("3. CODE COMPLEXITY COMPARISON\n")
cat(rep("-", 60), "\n", sep = "")

code_complexity <- data.frame(
  Tool = c("IlluMeta", "minfi (manual)", "ChAMP"),
  Lines_of_Code = c(
    1,  # python illumeta.py project_dir
    45, # Manual minfi workflow above
    25  # ChAMP workflow above
  ),
  Functions_to_Learn = c(
    1,  # Just run illumeta.py
    8,  # read.metharray, preprocessFunnorm, getBeta, detectionP, lmFit, makeContrasts, contrasts.fit, eBayes
    5   # champ.load, champ.norm, champ.DMP, etc.
  ),
  Auto_QC = c("Yes", "No", "Yes"),
  Auto_Batch_Correction = c("Yes", "No", "Yes"),
  Robustness_Assessment = c("Yes (CRF)", "No", "No"),
  Interactive_Dashboard = c("Yes", "No", "No"),
  stringsAsFactors = FALSE
)

print(code_complexity)

# ============================================================
# SUMMARY
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("BENCHMARK SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

summary_df <- data.frame(
  Metric = c("Processing Time (sec)", "Memory (MB)", "Probes Analyzed",
             "DMPs (FDR<0.05)", "DMPs (p<0.05)"),
  stringsAsFactors = FALSE
)

if (!is.null(results$minfi) && is.null(results$minfi$error)) {
  summary_df$minfi <- c(
    round(results$minfi$time, 1),
    round(results$minfi$memory_mb, 1),
    results$minfi$n_probes_analyzed,
    results$minfi$n_dmps_sig,
    results$minfi$n_dmps_nominal
  )
} else {
  summary_df$minfi <- rep(NA, 5)
}

if (!is.null(results$champ) && is.null(results$champ$error)) {
  summary_df$ChAMP <- c(
    round(results$champ$time, 1),
    round(results$champ$memory_mb, 1),
    results$champ$n_probes_analyzed,
    results$champ$n_dmps_sig,
    results$champ$n_dmps_nominal
  )
} else {
  summary_df$ChAMP <- rep(NA, 5)
}

print(summary_df)

# Save results
cat("\nSaving results...\n")
write.csv(summary_df, file.path(output_dir, "benchmark_summary.csv"), row.names = FALSE)
write.csv(code_complexity, file.path(output_dir, "code_complexity.csv"), row.names = FALSE)
saveRDS(results, file.path(output_dir, "benchmark_results.rds"))

cat("\nResults saved to:", output_dir, "\n")
cat("Done!\n")
