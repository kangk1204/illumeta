#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(sesame))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(dmrff))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(reformulas))

option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL, 
              help="Path to configure.tsv", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="Output directory for results", metavar="character"),
  make_option(c("-m", "--max_plots"), type="integer", default=10000, 
              help="Max points for interactive plots (default: 10000)", metavar="integer"),
  make_option(c("-p", "--pval"), type="double", default=0.05, 
              help="Adjusted P-value threshold (default: 0.05)", metavar="double"),
  make_option(c("-l", "--lfc"), type="double", default=0.5, 
              help="LogFC threshold (default: 0.5)", metavar="double"),
  make_option(c("--group_con"), type="character", default=NULL, 
              help="Control group label", metavar="character"),
  make_option(c("--group_test"), type="character", default=NULL, 
              help="Test group label", metavar="character"),
  make_option(c("--disable_auto_covariates"), action="store_true", default=FALSE,
              help="Disable automatic covariate selection via PCs"),
  make_option(c("--disable_sva"), action="store_true", default=FALSE,
              help="Disable surrogate variable analysis"),
  make_option(c("--include_covariates"), type="character", default="",
              help="Comma-separated covariate names to always try to include (if present)"),
  make_option(c("--permutations"), type="integer", default=0,
              help="Number of label permutations for null DMP counts (0 to skip)"),
  make_option(c("--vp_top"), type="integer", default=5000,
              help="Number of top-variable CpGs for variancePartition (default: 5000)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Timestamped messages for all logging in this script
ts_message <- function(..., domain = NULL, appendLF = TRUE) {
  base::message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " ")),
                domain = domain, appendLF = appendLF)
}
message <- ts_message

if (is.null(opt$config) || is.null(opt$group_con) || is.null(opt$group_test)){
  print_help(opt_parser)
  stop("Configuration file, group_con, and group_test must be supplied", call.=FALSE)
}

config_file <- opt$config
out_dir <- opt$out
max_points <- opt$max_plots
pval_thresh <- opt$pval
lfc_thresh <- opt$lfc
group_con_in <- opt$group_con
group_test_in <- opt$group_test
perm_n <- opt$permutations
vp_top <- opt$vp_top
disable_auto_cov <- opt$disable_auto_covariates
disable_sva <- opt$disable_sva
include_cov <- ifelse(opt$include_covariates == "", character(0), trimws(strsplit(opt$include_covariates, ",")[[1]]))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Set ExperimentHub/AnnotationHub Cache ---
sesame_cache_dir <- file.path(getwd(), "cache") 

if (!dir.exists(sesame_cache_dir)) dir.create(sesame_cache_dir, recursive = TRUE)

# Force BiocFileCache/Sesame to use this directory
Sys.setenv(R_USER_CACHE_DIR = sesame_cache_dir)
Sys.setenv(XDG_CACHE_HOME = sesame_cache_dir)

options(ExperimentHub.cache = sesame_cache_dir)
options(AnnotationHub.cache = sesame_cache_dir)

message(paste("Setting ExperimentHub/AnnotationHub cache to:", sesame_cache_dir))

# Ensure essential Sesame data is cached
# Skip sesameDataCache on analysis if cache already exists
sesame_cache_dir <- file.path(getwd(), "cache")
if (!dir.exists(sesame_cache_dir)) {
    tryCatch({
        message("  - Running sesameDataCache() for essential data...")
        sesameDataCache() 
        message("Sesame data caching complete in project cache.")
    }, error = function(e) {
        message("Warning: Failed to run sesameDataCache in analyze.R. Analysis may try to download files at runtime.")
        message(paste("Error detail:", e$message))
    })
} else {
    message("  - Cache directory exists; skipping sesameDataCache() to save time.")
}

# --- Helper Functions ---

save_interactive_plot <- function(p, filename, dir) {
  pp <- ggplotly(p)
  saveWidget(pp, file = file.path(dir, filename), selfcontained = TRUE)
}

save_datatable <- function(df, filename, dir) {
  # Identify P.Value column index (0-based for JS)
  p_idx <- which(colnames(df) == "P.Value") - 1
  if (length(p_idx) == 0) p_idx <- 0 # Fallback
  
  dt <- datatable(df, extensions = 'Buttons', options = list(
    dom = 'Bfrtip',
    pageLength = 25,
    buttons = c('csv', 'excel'),
    order = list(list(p_idx, 'asc'))
  ))
  saveWidget(dt, file = file.path(dir, filename), selfcontained = TRUE)
}

logit_offset <- function(betas, offset = 1e-6) {
  log2((betas + offset) / (1 - betas + offset))
}

auto_detect_covariates <- function(targets, gsm_col, pca_scores, alpha=0.05, max_pcs=5) {
  # Use PC–covariate association to decide covariates to include
  meta_cols <- setdiff(colnames(targets), c("primary_group", "Basename", "filenames", gsm_col))
  keep <- c()
  log_df <- data.frame(Variable=character(0), MinP=numeric(0), stringsAsFactors = FALSE)
  
  if (length(meta_cols) == 0 || ncol(pca_scores) == 0) {
    return(list(selected = keep, log = log_df))
  }
  
  # Limit PCs checked
  pca_use <- pca_scores[, 1:min(max_pcs, ncol(pca_scores)), drop=FALSE]
  
  for (col in meta_cols) {
    val <- targets[[col]]
    if (length(unique(val)) < 2) next
    
    # Skip if almost perfectly collinear with group (chi-square check for factors)
    if (is.character(val) || is.factor(val)) {
      tbl <- table(val, targets$primary_group)
      p_chisq <- tryCatch(chisq.test(tbl)$p.value, error=function(e) 1)
      if (!is.na(p_chisq) && p_chisq < 1e-6) next
    }
    
    p_vec <- rep(NA_real_, ncol(pca_use))
    for (k in seq_len(ncol(pca_use))) {
      if (is.numeric(val)) {
        fit <- lm(pca_use[,k] ~ val)
        coefs <- summary(fit)$coefficients
        if (nrow(coefs) >= 2) p_vec[k] <- coefs[2,4]
      } else {
        p_vec[k] <- kruskal.test(pca_use[,k] ~ as.factor(val))$p.value
      }
    }
    
    min_p <- suppressWarnings(min(p_vec, na.rm = TRUE))
    if (!is.finite(min_p)) next
    
    log_df <- rbind(log_df, data.frame(Variable = col, MinP = min_p))
    if (min_p < alpha) keep <- c(keep, col)
  }
  
  return(list(selected = keep, log = log_df))
}

  drop_linear_dependencies <- function(mat, group_cols) {
  # Reorder: group columns first to favor keeping them
  all_cols <- colnames(mat)
  group_cols <- intersect(group_cols, all_cols)
  other_cols <- setdiff(all_cols, group_cols)
  mat2 <- mat[, c(group_cols, other_cols), drop=FALSE]
  
  qr_obj <- qr(mat2)
  keep_idx <- sort(qr_obj$pivot[seq_len(qr_obj$rank)])
  keep_cols <- colnames(mat2)[keep_idx]
  drop_cols <- setdiff(colnames(mat2), keep_cols)
  
  list(mat = mat2[, keep_cols, drop=FALSE], dropped = drop_cols)
}

filter_covariates <- function(targets, covariates, group_col) {
  keep <- c()
  drops <- data.frame(Variable=character(0), Reason=character(0), stringsAsFactors = FALSE)
  n <- nrow(targets)
  group_vals <- targets[[group_col]]
  
  for (cv in covariates) {
    val <- targets[[cv]]
    uniq <- length(unique(val))
    if (uniq < 2) {
      drops <- rbind(drops, data.frame(Variable=cv, Reason="constant_or_single_level"))
      next
    }
    if (uniq >= (n - 1)) {
      drops <- rbind(drops, data.frame(Variable=cv, Reason="high_cardinality_vs_samples"))
      next
    }
    # Check confounding with group (any level exclusive to one group or strong chi-square)
    if (is.character(val) || is.factor(val)) {
      tbl <- table(val, group_vals)
      if (any(tbl == 0)) {
        drops <- rbind(drops, data.frame(Variable=cv, Reason="confounded_with_group"))
        next
      }
      p_chisq <- tryCatch(chisq.test(tbl)$p.value, error=function(e) NA)
      if (!is.na(p_chisq) && p_chisq < 1e-6) {
        drops <- rbind(drops, data.frame(Variable=cv, Reason="confounded_with_group"))
        next
      }
    }
    keep <- c(keep, cv)
  }
  list(keep = keep, dropped = drops)
}

summarize_pvals <- function(pmat, exclude = character(0), alpha = 0.05) {
  if (is.null(pmat) || length(pmat) == 0) return(list(sig_count = 0, min_p = NA_real_))
  rows <- setdiff(rownames(pmat), exclude)
  if (length(rows) == 0) return(list(sig_count = 0, min_p = NA_real_))
  vals <- pmat[rows, , drop = FALSE]
  vals_vec <- as.numeric(vals)
  vals_vec <- vals_vec[is.finite(vals_vec)]
  if (length(vals_vec) == 0) return(list(sig_count = 0, min_p = NA_real_))
  list(sig_count = sum(vals_vec < alpha, na.rm = TRUE),
       min_p = suppressWarnings(min(vals_vec, na.rm = TRUE)))
}

# --- 1. Load and Prepare Data ---

set.seed(12345) # Ensure reproducibility

message("Loading configuration...")
targets <- read.delim(config_file, stringsAsFactors = FALSE)

# Force Batch/ID columns to be factors (categorical)
if ("Sentrix_ID" %in% colnames(targets)) targets$Sentrix_ID <- as.factor(targets$Sentrix_ID)
if ("Sentrix_Position" %in% colnames(targets)) targets$Sentrix_Position <- as.factor(targets$Sentrix_Position)

if (!"primary_group" %in% colnames(targets)) stop("configure.tsv missing 'primary_group' column.")
targets$primary_group <- trimws(targets$primary_group) # remove whitespace

project_dir <- dirname(config_file)
idat_dir <- file.path(project_dir, "idat")

files <- list.files(idat_dir, pattern = "_Grn.idat", full.names = TRUE)
basenames <- unique(sub("_Grn.idat.*", "", files))

gsm_col <- grep("GSM|geo_accession", colnames(targets), value = TRUE, ignore.case = TRUE)[1]
if (is.na(gsm_col)) stop("Could not identify GSM ID column in configure.tsv")

targets$Basename <- NA
for (i in 1:nrow(targets)) {
  gsm <- targets[[gsm_col]][i]
  match <- grep(gsm, basenames, value = TRUE)
  if (length(match) > 0) targets$Basename[i] <- match[1]
}

targets <- targets[!is.na(targets$Basename), ]
if (nrow(targets) == 0) stop("No matching IDAT files found for the defined samples.")

# 1. Filter for specified groups (Case-Insensitive)
# Standardize to check match
target_groups_lower <- tolower(targets$primary_group)
con_lower <- tolower(group_con_in)
test_lower <- tolower(group_test_in)

keep_idx <- which(target_groups_lower %in% c(con_lower, test_lower))
targets <- targets[keep_idx, ]

if (nrow(targets) == 0) stop(paste("No samples found matching groups:", group_con_in, "or", group_test_in))

# Check if both groups exist
existing_groups <- unique(tolower(targets$primary_group))
if (!(con_lower %in% existing_groups)) stop(paste("Control group", group_con_in, "not found in configuration file."))
if (!(test_lower %in% existing_groups)) stop(paste("Test group", group_test_in, "not found in configuration file."))

# Count samples per group
n_con <- sum(tolower(targets$primary_group) == con_lower)
n_test <- sum(tolower(targets$primary_group) == test_lower)

message(paste("Found", nrow(targets), "samples for analysis."))
message(paste("  - Control Group:", group_con_in, paste0("(n=", n_con, ")")))
message(paste("  - Test Group:   ", group_test_in, paste0("(n=", n_test, ")")))

message("\n--- Sample Preview (Top 5) ---")
preview_cols <- c("primary_group", gsm_col, "Basename")
# Check if columns exist before selecting
preview_cols <- preview_cols[preview_cols %in% colnames(targets)]
preview_out <- capture.output(print(head(targets[, preview_cols], 5)))
message(paste(preview_out, collapse = "\n"))
message("------------------------------\n")

# 2. Set Factor Levels for Contrast (Test - Control)
actual_con <- unique(targets$primary_group[tolower(targets$primary_group) == con_lower])[1]
actual_test <- unique(targets$primary_group[tolower(targets$primary_group) == test_lower])[1]

clean_con <- make.names(actual_con)
clean_test <- make.names(actual_test)

targets$primary_group[tolower(targets$primary_group) == con_lower] <- clean_con
targets$primary_group[tolower(targets$primary_group) == test_lower] <- clean_test

targets$primary_group <- factor(targets$primary_group, levels = c(clean_con, clean_test))


# --- 2. Minfi Analysis ---

message("--- Starting Minfi Analysis ---")
rgSet <- read.metharray.exp(targets = targets)

# A. Sample-Level QC
message("Performing Sample QC...")
detP <- detectionP(rgSet)
qc <- getQC(preprocessRaw(rgSet))
bad_samples_idx <- which(qc$mMed < 10.5 | qc$uMed < 10.5)
if (length(bad_samples_idx) > 0) {
    bad_samples <- rownames(qc)[bad_samples_idx]
    message(paste("WARNING: Removing", length(bad_samples), "samples due to low signal intensity (< 10.5 log2):"))
    message(paste(bad_samples, collapse=", "))
    targets <- targets[-bad_samples_idx, ]
    rgSet <- rgSet[, -bad_samples_idx]
    detP <- detP[, -bad_samples_idx]
    if (nrow(targets) < 2) stop("Too few samples remaining after QC.")
}

# B. Normalization
message("Preprocessing (Noob)...")
mSet <- preprocessNoob(rgSet)
gmSet <- mapToGenome(mSet)

# C. Probe-Level QC
message("Performing Probe QC...")
keep_detP <- rowSums(detP < 0.01) == ncol(gmSet)
gmSet <- gmSet[keep_detP, ]
message(paste("  - Removed", sum(!keep_detP), "probes with high detection P-value."))

gmSet <- dropLociWithSnps(gmSet, snps=c("SBE","CpG"), maf=0)
message(paste("  - SNP filtering applied. Current probes:", nrow(gmSet)))

ann <- getAnnotation(gmSet)
keep_sex <- !(ann$chr %in% c("chrX", "chrY"))
gmSet <- gmSet[keep_sex, ]
message(paste("  - Removed Sex Chromosome probes. Final probes:", nrow(gmSet)))

beta_minfi <- getBeta(gmSet)

anno_data <- getAnnotation(gmSet)
anno <- anno_data[, c("chr", "pos", "Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
anno_df <- as.data.frame(anno)
anno_df$CpG <- rownames(anno_df)

colnames(anno_df)[colnames(anno_df) == "UCSC_RefGene_Name"] <- "Gene"
colnames(anno_df)[colnames(anno_df) == "UCSC_RefGene_Group"] <- "Region"
colnames(anno_df)[colnames(anno_df) == "Relation_to_Island"] <- "Island_Context"


# --- 3. Sesame Analysis ---
message("--- Starting Sesame Analysis ---")
ssets <- lapply(targets$Basename, function(x) readIDATpair(x))
betas_sesame_list <- lapply(ssets, function(x) getBetas(dyeBiasCorrTypeINorm(noob(x))))
common_probes_sesame <- Reduce(intersect, lapply(betas_sesame_list, names))
beta_sesame <- do.call(cbind, lapply(betas_sesame_list, function(x) x[common_probes_sesame]))
colnames(beta_sesame) <- targets[[gsm_col]]
beta_sesame <- na.omit(beta_sesame)
message(paste("Sesame processed", nrow(beta_sesame), "probes."))

# --- 4. Intersection ---
common_cpgs <- intersect(rownames(beta_minfi), rownames(beta_sesame))
message(paste("Intersection:", length(common_cpgs), "CpGs common to Minfi and Sesame."))

# --- Pipeline Function ---

run_pipeline <- function(betas, prefix, annotation_df) {
  message(paste("Running pipeline for:", prefix))
  
  curr_anno <- annotation_df[rownames(betas), ]
  common_ids <- intersect(rownames(betas), rownames(curr_anno))
  betas <- betas[common_ids, ]
  curr_anno <- curr_anno[common_ids, ]
  
  # Drop probes with zero variance across samples to avoid PCA failures
  var_vals <- apply(betas, 1, var, na.rm = TRUE)
  keep_var <- var_vals > 0
  if (sum(keep_var) < 2) {
    stop(paste(prefix, "has fewer than 2 probes with non-zero variance; cannot proceed."))
  }
  betas <- betas[keep_var, , drop = FALSE]
  
  # A. PCA and Covariate Check
  vars <- apply(betas, 1, var)
  top_vars <- head(order(vars, decreasing=TRUE), 1000)
  beta_clust <- betas[top_vars, ]
  
  sample_dists <- as.matrix(dist(t(beta_clust)))
  dist_df <- as.data.frame(as.table(sample_dists))
  colnames(dist_df) <- c("Sample1", "Sample2", "Distance")
  
  p_dist <- ggplot(dist_df, aes(x=Sample1, y=Sample2, fill=Distance, text=Distance)) +
      geom_tile() +
      scale_fill_gradient(low="red", high="white") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6)) +
      ggtitle(paste(prefix, "Sample Distance Matrix (Top 1000 Var Probes)"))
      
  save_interactive_plot(p_dist, paste0(prefix, "_Sample_Clustering_Distance.html"), out_dir)

  pca_res <- prcomp(t(betas), scale. = TRUE)
  pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = targets$primary_group)
  p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) + geom_point(size=3) + theme_minimal() + ggtitle(paste(prefix, "PCA (Raw)"))
  save_interactive_plot(p_pca, paste0(prefix, "_PCA_Before.html"), out_dir)
  
  # Auto-detect covariates via PC association (top PCs) unless disabled
  covariates <- character(0)
  covar_log_path <- file.path(out_dir, paste0(prefix, "_AutoCovariates.csv"))
  drop_log <- data.frame(Variable=character(0), Reason=character(0))
  
  if (!disable_auto_cov) {
    covar_info <- auto_detect_covariates(targets, gsm_col, pca_res$x, alpha = 0.05, max_pcs = 5)
    if (nrow(covar_info$log) > 0) write.csv(covar_info$log, covar_log_path, row.names = FALSE)
    covariates <- covar_info$selected
    if (length(covariates) > 0) {
      filt <- filter_covariates(targets, covariates, group_col = "primary_group")
      covariates <- filt$keep
      if (nrow(filt$dropped) > 0) drop_log <- rbind(drop_log, filt$dropped)
    }
  } else {
    message("  Auto covariate selection disabled by flag.")
  }
  
  # Force-include covariates if present
  if (length(include_cov) > 0) {
    present_force <- intersect(include_cov, colnames(targets))
    missing_force <- setdiff(include_cov, present_force)
    covariates <- unique(c(covariates, present_force))
    if (length(missing_force) > 0) {
      drop_log <- rbind(drop_log, data.frame(Variable = missing_force, Reason = "not_in_config"))
    }
  }
  
  if (length(covariates) > 0) {
    message(paste("  Covariate candidates after filtering:", paste(covariates, collapse = ", ")))
  } else {
    message("  No covariate candidates retained before design drop.")
  }
  
  # Build Design Matrix
  formula_str <- "~ 0 + primary_group"
  if (length(covariates) > 0) formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = targets)
  colnames(design) <- gsub("primary_group", "", colnames(design))
  colnames(design) <- make.unique(make.names(colnames(design)))
  group_cols <- make.names(levels(targets$primary_group))
  
  sv_cols <- character(0)
  
  # --- Surrogate Variable Analysis (SVA) ---
  if (!disable_sva) {
    message("  Running SVA to detect hidden batch effects...")
    
    f_str <- "~ primary_group"
    if (length(covariates) > 0) f_str <- paste(f_str, "+", paste(covariates, collapse = " + "))
    mod_sva <- model.matrix(as.formula(f_str), data = targets)
    
    f0_str <- "~ 1"
    if (length(covariates) > 0) f0_str <- paste(f0_str, "+", paste(covariates, collapse = " + "))
    mod0_sva <- model.matrix(as.formula(f0_str), data = targets)
    
    run_sva <- function(mod, mod0, label) {
      if (nrow(mod) <= ncol(mod) + 2) {
        message(paste("    - Skipping SVA (", label, "): Sample size too small for stable estimation.", sep=""))
        return(NULL)
      }
      svobj <- sva(betas, mod, mod0)
      return(svobj)
    }
    
    tryCatch({
        svobj <- run_sva(mod_sva, mod0_sva, "full")
        if (is.null(svobj)) {
          # skip
        } else {
          n_sv <- svobj$n.sv
          if (n_sv == 0) message("    - No significant surrogate variables found.")
        }
    }, error = function(e) {
        message("    - SVA failed on full model: ", e$message, " | Retrying with group-only model...")
        mod_sva2 <- model.matrix(~ primary_group, data = targets)
        mod0_sva2 <- model.matrix(~ 1, data = targets)
        tryCatch({
          svobj <<- run_sva(mod_sva2, mod0_sva2, "group-only")
        }, error = function(e2) {
          message("    - SVA fallback failed: ", e2$message)
        })
    })
    
    if (exists("svobj") && !is.null(svobj)) {
        n_sv <- svobj$n.sv
        if (n_sv > 0) {
            max_sv <- min(n_sv, 5, floor(nrow(targets) * 0.2))
            if (max_sv < 1) max_sv <- 1
            
            if (max_sv < n_sv) {
                message(paste("    - Detected", n_sv, "SVs. Capping at", max_sv, "SVs."))
                sv_mat <- svobj$sv[, 1:max_sv, drop=FALSE]
                colnames(sv_mat) <- paste0("SV", 1:max_sv)
            } else {
                message(paste("    - Detected and using", n_sv, "SVs."))
                sv_mat <- svobj$sv
                colnames(sv_mat) <- paste0("SV", 1:n_sv)
            }
            
            design <- cbind(design, sv_mat)
            targets <- cbind(targets, as.data.frame(sv_mat))
            sv_cols <- colnames(sv_mat)
        } else {
            message("    - No significant surrogate variables found.")
        }
    }
  } else {
    message("  SVA disabled by flag.")
  }
  colnames(design) <- make.unique(colnames(design))

  design_ld <- drop_linear_dependencies(design, group_cols = group_cols)
  if (length(design_ld$dropped) > 0) {
    drop_log <- rbind(drop_log, data.frame(Variable = design_ld$dropped, Reason = "linear_dependency"))
    message(paste("  Dropped collinear terms:", paste(design_ld$dropped, collapse=", ")))
  }
  design <- design_ld$mat
  covariates <- setdiff(covariates, design_ld$dropped)
  used_covariates <- covariates
  if (length(used_covariates) > 0) {
    message(paste("  Covariates used in design:", paste(used_covariates, collapse = ", ")))
  } else {
    message("  Covariates used in design: none")
  }
  
  # Persist drop log (filtering + linear dependency)
  if (nrow(drop_log) > 0) {
    drop_path <- file.path(out_dir, paste0(prefix, "_DroppedCovariates.csv"))
    write.table(drop_log, drop_path, sep=",", row.names=FALSE, col.names=TRUE)
  }

  # PCA After Correction (covariates + SVs if present)
  all_covariates <- c(used_covariates, sv_cols)
  if (length(all_covariates) > 0) {
    cov_formula <- as.formula(paste("~ 0 +", paste(all_covariates, collapse = " + ")))
    cov_mat <- model.matrix(cov_formula, data = targets)
    
    design_keep <- model.matrix(~ primary_group, data=targets)
    
    clean_betas <- removeBatchEffect(betas, covariates=cov_mat, design=design_keep)
    
    pca_clean <- prcomp(t(clean_betas), scale. = TRUE)
    pca_clean_df <- data.frame(PC1 = pca_clean$x[,1], PC2 = pca_clean$x[,2], Group = targets$primary_group)
    p_pca_clean <- ggplot(pca_clean_df, aes(x=PC1, y=PC2, color=Group)) + 
      geom_point(size=3) + theme_minimal() + ggtitle(paste(prefix, "PCA (Corrected: covariates + SVs)"))
    save_interactive_plot(p_pca_clean, paste0(prefix, "_PCA_After_Correction.html"), out_dir)
  }
  
  groups <- levels(targets$primary_group)
  if (length(groups) < 2) {
      message("Skipping stats: Not enough groups.")
      return(NULL)
  }
  contrast_str <- paste0(groups[2], "-", groups[1])
  message(paste("  - Contrast:", contrast_str))
  
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  
  # B. Stats
  m_vals <- logit_offset(betas)
  fit <- lmFit(m_vals, design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  res$CpG <- rownames(res)
  
  # Add annotation
  res <- merge(res, curr_anno, by.x="CpG", by.y="CpG")
  
  clean_gene_names <- function(g_str) {
    if (is.na(g_str) || g_str == "") return("")
    genes <- unlist(strsplit(as.character(g_str), ";"))
    genes <- unique(genes)
    if (length(genes) > 2) genes <- genes[1:2]
    final_str <- paste(genes, collapse=";")
    if (nchar(final_str) > 15) final_str <- paste0(substr(final_str, 1, 12), "...")
    return(final_str)
  }
  
  res$Gene <- sapply(res$Gene, clean_gene_names)

  # C. Plots
  p_vals <- res$P.Value
  p_vals <- p_vals[!is.na(p_vals)]
  n_p <- length(p_vals)
  
  chisq_vals <- qchisq(1 - p_vals, 1)
  lambda_val <- median(chisq_vals) / qchisq(0.5, 1)
  message(paste("  Genomic Inflation Factor (Lambda):", round(lambda_val, 3)))
  
  expected <- -log10(ppoints(n_p))
  observed <- -log10(sort(p_vals))
  qq_df <- data.frame(Expected = expected, Observed = observed)
  
  if (nrow(qq_df) > max_points) {
      idx <- unique(c(1:1000, seq(1001, n_p, length.out=max_points)))
      qq_df <- qq_df[idx, ]
  }
  
  p_qq <- ggplot(qq_df, aes(x=Expected, y=Observed)) +
      geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
      geom_point(alpha=0.5, size=1) +
      theme_minimal() +
      labs(x="Expected -log10(P)", y="Observed -log10(P)") +
      ggtitle(paste(prefix, "Q-Q Plot (Lambda =", round(lambda_val, 3), ")"))
      
  save_interactive_plot(p_qq, paste0(prefix, "_QQPlot.html"), out_dir)

  plot_res <- res[order(res$P.Value), ]
  if (nrow(plot_res) > max_points) {
      message(paste("Subsampling top", max_points, "probes for interactive plots..."))
      plot_res <- plot_res[1:max_points, ]
  }
  
  plot_res$diffexpressed <- "NO"
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC > lfc_thresh] <- "UP"
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC < -lfc_thresh] <- "DOWN"
  
  subtitle_str <- paste0("Control: ", group_con_in, " (n=", n_con, ") vs Test: ", group_test_in, " (n=", n_test, ")")
  
  p_vol <- ggplot(plot_res, aes(x=logFC, y=-log10(P.Value), color=diffexpressed, 
                  text=paste("CpG:", CpG, "<br>Gene:", Gene, "<br>Region:", Region, "<br>Island:", Island_Context))) +
    geom_point(alpha=0.5) + theme_minimal() +
    scale_color_manual(values=c("DOWN"="blue", "NO"="grey", "UP"="red")) +
    labs(title = paste(prefix, "Volcano Plot"), 
         subtitle = subtitle_str,
         color = "Status")
  save_interactive_plot(p_vol, paste0(prefix, "_Volcano.html"), out_dir)
  
  chr_map <- c(as.character(1:22), "X", "Y", "M")
  chr_vals <- 1:25
  names(chr_vals) <- chr_map
  
  clean_chr <- gsub("chr", "", plot_res$chr)
  plot_res$chr_num <- chr_vals[clean_chr]
  
  plot_res_man <- plot_res[!is.na(plot_res$chr_num) & !is.na(plot_res$pos), ]
  plot_res_man <- plot_res_man[order(plot_res_man$chr_num, plot_res_man$pos), ]
  
  chr_len <- tapply(plot_res_man$pos, plot_res_man$chr_num, max)
  chr_len <- chr_len[!is.na(chr_len)]
  
  chrs_present <- sort(unique(plot_res_man$chr_num))
  
  if (length(chrs_present) > 0) {
      cum_lengths <- cumsum(as.numeric(chr_len))
      offsets <- c(0, cum_lengths[-length(cum_lengths)])
      names(offsets) <- names(chr_len)
      
      plot_res_man$pos_cum <- plot_res_man$pos + offsets[as.character(plot_res_man$chr_num)]
      
      axis_set <- plot_res_man %>%
        group_by(chr_num) %>%
        summarize(center = (min(pos_cum) + max(pos_cum)) / 2)
      
      plot_res_man$color_group <- factor(plot_res_man$chr_num %% 2)
      chr_labels <- names(chr_vals)[match(axis_set$chr_num, chr_vals)]
      
      p_man <- ggplot(plot_res_man, aes(x=pos_cum, y=-log10(P.Value), color=as.factor(chr_num), 
                        text=paste("CpG:", CpG, "<br>Gene:", Gene, "<br>Region:", Region, "<br>Island:", Island_Context))) +
        geom_point(alpha=0.7, size=1) +
        scale_x_continuous(label = chr_labels, breaks = axis_set$center) +
        theme_minimal() +
        theme(
          legend.position="none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
        ) +
        xlab("Chromosome") +
        ggtitle(paste(prefix, "Manhattan Plot (Top", nrow(plot_res_man), "CpG sites)"))
      
      save_interactive_plot(p_man, paste0(prefix, "_Manhattan.html"), out_dir)
  } else {
      message("Skipping Manhattan: No valid chromosome data found in subsample.")
  }
  
  message("  Generating Top 100 Heatmap...")
  top_100 <- res[order(res$P.Value), ][1:min(100, nrow(res)), ]
  valid_cpgs <- intersect(top_100$CpG, rownames(betas))
  if (length(valid_cpgs) < 2) {
      message("  Skipping heatmap: too few matching CpGs found.")
  } else {
      top_betas <- betas[valid_cpgs, , drop=FALSE]
      row_vars <- apply(top_betas, 1, var)
      top_betas <- top_betas[row_vars > 0, , drop=FALSE]
      
      if (nrow(top_betas) >= 2) {
          top_betas_scaled <- t(scale(t(top_betas)))
          gene_map <- top_100$Gene
          names(gene_map) <- top_100$CpG
          
          hc_r <- hclust(dist(top_betas_scaled))
          row_ord <- hc_r$order
          hc_c <- hclust(dist(t(top_betas_scaled)))
          col_ord <- hc_c$order
          
          hm_df <- expand.grid(CpG = rownames(top_betas_scaled), Sample = colnames(top_betas_scaled))
          hm_df$Value <- as.vector(top_betas_scaled)
          
          hm_df$Group <- targets$primary_group[match(hm_df$Sample, targets[[gsm_col]])]
          if(any(is.na(hm_df$Group))) {
               hm_df$Group <- targets$primary_group[match(hm_df$Sample, colnames(betas))]
          }
          
          hm_df$Gene <- gene_map[as.character(hm_df$CpG)]
          hm_df$Label <- paste(hm_df$CpG, hm_df$Gene, sep="_")
          
          # Cluster rows and columns for ordering
          hc_r <- hclust(dist(top_betas_scaled))
          hc_c <- hclust(dist(t(top_betas_scaled)))
          row_ord <- hc_r$order
          col_ord <- hc_c$order
          
          ordered_cpgs <- rownames(top_betas_scaled)[row_ord]
          ordered_labels <- paste(ordered_cpgs, gene_map[ordered_cpgs], sep="_")
          
          hm_df$Label <- factor(hm_df$Label, levels = ordered_labels)
          
          # Fix sample order by group (Control first, then Test), within-group keep cluster order
          samples_ordered <- colnames(top_betas_scaled)[col_ord]
          grp_vec <- targets$primary_group[match(samples_ordered, targets[[gsm_col]])]
          grp_vec[is.na(grp_vec)] <- targets$primary_group[match(samples_ordered, colnames(betas))]
          grp_levels <- unique(as.character(targets$primary_group))
          if (length(grp_levels) == 2) {
              samples_by_group <- c(samples_ordered[grp_vec == grp_levels[1]], samples_ordered[grp_vec == grp_levels[2]])
          } else {
              samples_by_group <- samples_ordered
          }
          hm_df$Sample <- factor(hm_df$Sample, levels = samples_by_group)
          
          p_heat_top <- ggplot(hm_df, aes(x=Sample, y=Label, fill=Value, text=paste("Group:", Group))) +
              geom_tile() +
              scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
                    axis.text.y = element_text(size=6),
                    axis.title.x = element_blank(),
                    legend.position = "none")
          
          anno_df <- unique(hm_df[, c("Sample", "Group")])
          p_anno <- ggplot(anno_df, aes(x=Sample, y=1, fill=Group, text=paste("Group:", Group))) +
              geom_tile() +
              theme_void() +
              theme(legend.position = "none",
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    plot.margin = margin(0, 0, 0, 0))
          
          p_heat_top_ly <- ggplotly(p_heat_top, tooltip = c("x", "y", "text"))
          p_anno_ly <- ggplotly(p_anno, tooltip = c("text"))
          
          final_heatmap <- subplot(p_anno_ly, p_heat_top_ly, nrows = 2, heights = c(0.05, 0.95), shareX = TRUE, titleX = FALSE) %>%
            layout(title = paste(prefix, "Top 100 CpG sites Heatmap"))
              
          saveWidget(final_heatmap, file = file.path(out_dir, paste0(prefix, "_Top100_Heatmap.html")), selfcontained = TRUE)
      }
  }

  col_order <- c("CpG", "Gene", setdiff(colnames(res), c("CpG", "Gene")))
  res <- res[, col_order]

  save_datatable(head(res, 2000), paste0(prefix, "_Top_DMPs.html"), out_dir)
  
  # --- D. DMR Analysis (dmrff) ---
  message("  Running DMR analysis (dmrff)...")
  
  # Prepare inputs for dmrff
  # We need: estimate (logFC), se, p.value, methylation, chr, pos
  # fit2 has: coefficients (estimate), stdev.unscaled * sigma (se), p.value
  
  # Ensure we match the rows of fit2 and betas
  use_idx <- match(rownames(fit2), rownames(betas))
  betas_dmr <- betas[use_idx, ]
  
  dmr_est <- fit2$coefficients[, 1]
  dmr_se <- fit2$stdev.unscaled[, 1] * fit2$sigma
  dmr_p <- fit2$p.value[, 1]
  
  # Annotation for dmrff
  dmr_anno <- curr_anno[match(names(dmr_est), rownames(curr_anno)), ]
  
  tryCatch({
    dmr_res <- dmrff(
      estimate = dmr_est,
      se = dmr_se,
      p.value = dmr_p,
      methylation = betas_dmr,
      chr = dmr_anno$chr,
      pos = dmr_anno$pos,
      maxgap = 500,
      p.cutoff = 0.05
    )
    
    if (nrow(dmr_res) > 0) {
        message(paste("    - Found", nrow(dmr_res), "DMRs."))
        
        # Save Table
        dmr_filename <- paste0(prefix, "_DMRs_Table.html")
        save_datatable(head(dmr_res, 3000), dmr_filename, out_dir)
        write.csv(dmr_res, file.path(out_dir, paste0(prefix, "_DMRs.csv")))
        
        # 1. DMR Volcano Plot
        # dmrff returns 'estimate' (avg LogFC) and 'p.adjust'
        dmr_res$diffexpressed <- "NO"
        dmr_res$diffexpressed[dmr_res$p.adjust < 0.05 & dmr_res$estimate > 0] <- "UP"
        dmr_res$diffexpressed[dmr_res$p.adjust < 0.05 & dmr_res$estimate < 0] <- "DOWN"
        
        p_vol_dmr <- ggplot(dmr_res, aes(x=estimate, y=-log10(p.adjust), color=diffexpressed,
                        text=paste("Region:", paste0(chr, ":", start, "-", end), "<br>nCpG:", n, "<br>P-adj:", p.adjust))) +
            geom_point(alpha=0.6, size=2) + theme_minimal() +
            scale_color_manual(values=c("DOWN"="blue", "NO"="grey", "UP"="red")) +
            labs(title = paste(prefix, "DMR Volcano Plot"), 
                 subtitle = "dmrff analysis",
                 x = "Mean Methylation Difference (Estimate)", y = "-log10 Adjusted P-value",
                 color = "Status")
        save_interactive_plot(p_vol_dmr, paste0(prefix, "_DMR_Volcano.html"), out_dir)
        
        # 2. DMR Manhattan Plot
        chr_map <- c(as.character(1:22), "X", "Y", "M")
        chr_vals <- 1:25
        names(chr_vals) <- chr_map
        
        clean_chr_dmr <- gsub("chr", "", dmr_res$chr)
        dmr_res$chr_num <- chr_vals[clean_chr_dmr]
        dmr_res$mid_pos <- (dmr_res$start + dmr_res$end) / 2
        
        plot_dmr_man <- dmr_res[!is.na(dmr_res$chr_num), ]
        plot_dmr_man <- plot_dmr_man[order(plot_dmr_man$chr_num, plot_dmr_man$start), ]
        
        if (nrow(plot_dmr_man) > 0) {
            chr_len_dmr <- tapply(plot_dmr_man$end, plot_dmr_man$chr_num, max) # approximate
            cum_lengths <- cumsum(as.numeric(chr_len_dmr))
            offsets <- c(0, cum_lengths[-length(cum_lengths)])
            names(offsets) <- names(chr_len_dmr)
            
            plot_dmr_man$pos_cum <- plot_dmr_man$mid_pos + offsets[as.character(plot_dmr_man$chr_num)]
            
            axis_set_dmr <- plot_dmr_man %>%
                group_by(chr_num) %>%
                summarize(center = (min(pos_cum) + max(pos_cum)) / 2)
            
            chr_labels_dmr <- names(chr_vals)[match(axis_set_dmr$chr_num, chr_vals)]
            
            p_man_dmr <- ggplot(plot_dmr_man, aes(x=pos_cum, y=-log10(p.adjust), color=as.factor(chr_num %% 2),
                                text=paste("Region:", paste0(chr, ":", start, "-", end), "<br>nCpG:", n))) +
                geom_point(alpha=0.7, size=2) +
                scale_color_manual(values=c("0"="#2c3e50", "1"="#3498db")) +
                scale_x_continuous(label = chr_labels_dmr, breaks = axis_set_dmr$center) +
                theme_minimal() +
                theme(legend.position="none",
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)) +
                xlab("Chromosome") +
                ggtitle(paste(prefix, "DMR Manhattan Plot"))
            
            save_interactive_plot(p_man_dmr, paste0(prefix, "_DMR_Manhattan.html"), out_dir)
        }
        
        # 3. DMR Heatmap (Top 50 DMRs - Average Methylation per Region)
        top_dmrs <- head(dmr_res[order(dmr_res$p.adjust), ], 50)
        
        # Calculate average methylation for each top DMR across all samples
        dmr_means <- matrix(NA, nrow=nrow(top_dmrs), ncol=ncol(betas_dmr))
        rownames(dmr_means) <- paste0(top_dmrs$chr, ":", top_dmrs$start, "-", top_dmrs$end)
        colnames(dmr_means) <- colnames(betas_dmr)
        
        for (i in 1:nrow(top_dmrs)) {
            # Find CpGs in this region
            # Using the annotation we prepared earlier (dmr_anno)
            # dmrff results provide chr, start, end
            # We need to find CpGs in dmr_anno that fall in this range
            
            # Filter annotation for this chr
            chr_cpgs <- dmr_anno[dmr_anno$chr == top_dmrs$chr[i], ]
            in_region <- chr_cpgs$pos >= top_dmrs$start[i] & chr_cpgs$pos <= top_dmrs$end[i]
            cpg_ids <- rownames(chr_cpgs)[in_region]
            cpg_ids <- intersect(cpg_ids, rownames(betas_dmr))
            
            if (length(cpg_ids) > 0) {
                if (length(cpg_ids) == 1) {
                    dmr_means[i, ] <- betas_dmr[cpg_ids, ]
                } else {
                    dmr_means[i, ] <- colMeans(betas_dmr[cpg_ids, , drop=FALSE], na.rm=TRUE)
                }
            }
        }
        
        dmr_means <- dmr_means[rowSums(is.na(dmr_means)) == 0, , drop=FALSE]
        
        if (nrow(dmr_means) >= 2) {
             dmr_means_scaled <- t(scale(t(dmr_means)))
             
             hm_dmr_df <- expand.grid(Region = rownames(dmr_means_scaled), Sample = colnames(dmr_means_scaled))
             hm_dmr_df$Value <- as.vector(dmr_means_scaled)
             
             hm_dmr_df$Group <- targets$primary_group[match(hm_dmr_df$Sample, targets[[gsm_col]])]
             if(any(is.na(hm_dmr_df$Group))) {
                  hm_dmr_df$Group <- targets$primary_group[match(hm_dmr_df$Sample, colnames(betas))]
             }
             
             # Clustering order
             hc_r <- hclust(dist(dmr_means_scaled))
             hc_c <- hclust(dist(t(dmr_means_scaled)))
             
             hm_dmr_df$Region <- factor(hm_dmr_df$Region, levels = rownames(dmr_means_scaled)[hc_r$order])
             
             # Sample order: Control then Test, but clustered within
             samples_ordered <- colnames(dmr_means_scaled)[hc_c$order]
             grp_vec <- targets$primary_group[match(samples_ordered, targets[[gsm_col]])]
             grp_vec[is.na(grp_vec)] <- targets$primary_group[match(samples_ordered, colnames(betas))]
             
             if (length(levels(targets$primary_group)) == 2) {
                  g_lev <- levels(targets$primary_group)
                  samples_by_group <- c(samples_ordered[grp_vec == g_lev[1]], samples_ordered[grp_vec == g_lev[2]])
             } else {
                  samples_by_group <- samples_ordered
             }
             hm_dmr_df$Sample <- factor(hm_dmr_df$Sample, levels = samples_by_group)
             
             p_heat_dmr <- ggplot(hm_dmr_df, aes(x=Sample, y=Region, fill=Value, text=paste("Group:", Group))) +
                  geom_tile() +
                  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
                        axis.text.y = element_text(size=6),
                        axis.title.x = element_blank(),
                        legend.position = "none") +
                  ggtitle(paste(prefix, "Top 50 DMRs Heatmap (Avg Methylation)"))
                  
             save_interactive_plot(p_heat_dmr, paste0(prefix, "_Top_DMRs_Heatmap.html"), out_dir)
        }
        
    } else {
        message("    - No DMRs found.")
    }
  }, error = function(e) {
      message("    - DMR analysis failed: ", e$message)
  })

  evaluate_batch <- function(data_betas, meta, label, file_suffix, allowed_vars) {
     pca <- prcomp(t(data_betas), scale.=TRUE)
     eig <- (pca$sdev)^2
     var_expl <- eig/sum(eig)
     
     n_pcs <- min(10, ncol(pca$x))
     pcs <- pca$x[, 1:n_pcs]
     
     allowed <- unique(c("primary_group", allowed_vars))
     allowed <- intersect(allowed, colnames(meta))
     keep_cols <- vapply(allowed, function(nm) {
         if (nm == "primary_group") return(TRUE)
         if (grepl("^SV", nm)) return(length(unique(meta[[nm]])) > 1)
         uniq <- length(unique(meta[[nm]]))
         return(uniq > 1 & uniq < nrow(meta))
     }, logical(1))
     test_meta <- meta[, allowed[keep_cols], drop=FALSE]
     
     pvals <- matrix(NA, nrow=ncol(test_meta), ncol=n_pcs)
     rownames(pvals) <- colnames(test_meta)
     colnames(pvals) <- paste0("PC", 1:n_pcs)
     
     for (m in colnames(test_meta)) {
         for (k in 1:n_pcs) {
             val <- test_meta[[m]]
             if (is.numeric(val)) {
                 fit <- lm(pcs[,k] ~ val)
                 coefs <- summary(fit)$coefficients
                 if (nrow(coefs) >= 2) {
                     pvals[m,k] <- coefs[2,4]
                 } else {
                     pvals[m,k] <- NA
                 }
             } else {
                 pvals[m,k] <- kruskal.test(pcs[,k] ~ as.factor(val))$p.value
             }
         }
     }
     
     csv_name <- paste0(prefix, "_Batch_Evaluation_", file_suffix, "_Table.csv")
     write.csv(pvals, file.path(out_dir, csv_name))
     
     pvals_melt <- as.data.frame(as.table(pvals))
     colnames(pvals_melt) <- c("Variable", "PC", "P_Value")
     pvals_melt$LogP <- -log10(pvals_melt$P_Value)
     
     p_heat <- ggplot(pvals_melt, aes(x=PC, y=Variable, fill=LogP, text=P_Value)) +
         geom_tile() +
         scale_fill_gradient(low="white", high="red") +
         theme_minimal() +
         ggtitle(paste(label, "- PC-Covariate Association (Top", n_pcs, "PCs)"))
         
     return(list(plot = p_heat, pvals = pvals))
  }
  
  used_covariates <- setdiff(covariates, unique(drop_log$Variable))
  allowed_vars <- unique(c(used_covariates, sv_cols))
  
  p_eval_before <- evaluate_batch(betas, targets, paste(prefix, "Before"), "Before", allowed_vars)
  save_interactive_plot(p_eval_before$plot, paste0(prefix, "_Batch_Evaluation_Before.html"), out_dir)
  pvals_before <- p_eval_before$pvals
  
  p_eval_after <- NULL
  pvals_after <- NULL
  if (length(all_covariates) > 0) {
      cov_formula <- as.formula(paste("~ 0 +", paste(all_covariates, collapse = " + ")))
      cov_mat <- model.matrix(cov_formula, data = targets)
      
      design_keep <- model.matrix(~ primary_group, data=targets)
      clean_betas <- removeBatchEffect(betas, covariates=cov_mat, design=design_keep)
      
      p_eval_after <- evaluate_batch(clean_betas, targets, paste(prefix, "After"), "After", allowed_vars)
      pvals_after <- p_eval_after$pvals
      save_interactive_plot(p_eval_after$plot, paste0(prefix, "_Batch_Evaluation_After.html"), out_dir)
  }
  
  # Permutation test on top-variable CpGs (optional)
  perm_summary <- NULL
  if (perm_n > 0) {
      message(paste("  Running", perm_n, "permutations for null DMP counts (top", vp_top, "CpGs)..."))
      top_perm <- head(order(apply(betas, 1, var), decreasing = TRUE), min(vp_top, nrow(betas)))
      m_perm <- logit_offset(betas[top_perm, , drop=FALSE])
      perm_results <- data.frame(run = integer(0), sig_count = integer(0), min_p = numeric(0))
      perm_group_cols <- make.names(levels(targets$primary_group))
      
      for (i in seq_len(perm_n)) {
          perm_labels <- sample(targets$primary_group)
          perm_df <- data.frame(perm_labels = perm_labels)
          perm_design <- model.matrix(~ 0 + perm_labels, data = perm_df)
          colnames(perm_design) <- gsub("perm_labels", "", colnames(perm_design))
          
          if (length(all_covariates) > 0) {
              cov_formula <- as.formula(paste("~ 0 +", paste(all_covariates, collapse = " + ")))
              cov_mat <- model.matrix(cov_formula, data = targets)
              perm_design <- cbind(perm_design, cov_mat)
          }
          
          perm_ld <- drop_linear_dependencies(perm_design, group_cols = perm_group_cols)
          perm_design <- perm_ld$mat
          if (ncol(perm_design) < 2) next
          
          cont_vec <- rep(0, ncol(perm_design))
          cont_vec[2] <- 1; cont_vec[1] <- -1
          cm_perm <- matrix(cont_vec, ncol = 1)
          rownames(cm_perm) <- colnames(perm_design)
          
          fitp <- lmFit(m_perm, perm_design)
          fitp2 <- contrasts.fit(fitp, cm_perm)
          fitp2 <- eBayes(fitp2)
          resp <- topTable(fitp2, coef = 1, number = Inf, adjust.method = "BH")
          
          sig <- sum(resp$adj.P.Val < pval_thresh & abs(resp$logFC) > lfc_thresh, na.rm=TRUE)
          minp <- suppressWarnings(min(resp$P.Value, na.rm=TRUE))
          perm_results <- rbind(perm_results, data.frame(run = i, sig_count = sig, min_p = minp))
      }
      perm_summary <- data.frame(
        stat = c("mean_sig", "median_sig", "max_sig", "min_p_median"),
        value = c(mean(perm_results$sig_count), median(perm_results$sig_count),
                  max(perm_results$sig_count), median(perm_results$min_p))
      )
      write.csv(perm_results, file.path(out_dir, paste0(prefix, "_Permutation_Results.csv")), row.names = FALSE)
      write.csv(perm_summary, file.path(out_dir, paste0(prefix, "_Permutation_Summary.csv")), row.names = FALSE)
  }
  
  # variancePartition on top-variable CpGs (optional)
  vp_primary <- NA
  if (vp_top > 0) {
      message(paste("  Running variancePartition on top", vp_top, "CpGs..."))
      top_vp <- head(order(apply(betas, 1, var), decreasing = TRUE), min(vp_top, nrow(betas)))
      m_vp <- logit_offset(betas[top_vp, , drop=FALSE])
      
      meta_vp <- targets
      # keep only variables in design (group + covariates + SVs)
      keep_meta <- unique(c("primary_group", used_covariates, sv_cols))
      meta_vp <- meta_vp[, intersect(colnames(meta_vp), keep_meta), drop=FALSE]
      # coerce characters to factors
      for (nm in colnames(meta_vp)) {
        if (is.character(meta_vp[[nm]])) meta_vp[[nm]] <- as.factor(meta_vp[[nm]])
      }
      # align rownames with expression matrix columns
      rownames(meta_vp) <- colnames(m_vp)
      
      form_terms <- keep_meta
      if (length(form_terms) == 0) {
        message("    - Skipping variancePartition: no terms to evaluate.")
      } else {
        form_str <- paste("~", paste(form_terms, collapse = " + "))
        form_vp <- as.formula(form_str)
        vp_res <- tryCatch({
          fitExtractVarPartModel(m_vp, form_vp, meta_vp)
        }, error = function(e) { message("    - variancePartition failed: ", e$message); NULL })
        
        if (!is.null(vp_res)) {
          # vp_res: genes x variables; summarize by mean across genes
          vp_summary <- data.frame(Variable = colnames(vp_res), Mean = colMeans(vp_res, na.rm=TRUE))
          write.csv(vp_summary, file.path(out_dir, paste0(prefix, "_VariancePartition.csv")), row.names = FALSE)
          if ("primary_group" %in% vp_summary$Variable) {
            vp_primary <- vp_summary$Mean[vp_summary$Variable == "primary_group"][1]
          }
        }
      }
  }
  
  # Metrics summary for reviewers
  exclude_batch <- c("primary_group", grep("^SV", rownames(pvals_before), value=TRUE))
  batch_before <- summarize_pvals(pvals_before, exclude = exclude_batch)
  batch_after <- summarize_pvals(pvals_after, exclude = exclude_batch)
  grp_min_before <- tryCatch(suppressWarnings(min(pvals_before["primary_group",], na.rm=TRUE)), error=function(e) NA)
  grp_min_after <- tryCatch({
    if (is.null(pvals_after)) NA else suppressWarnings(min(pvals_after["primary_group",], na.rm=TRUE))
  }, error=function(e) NA)
  
  metrics <- data.frame(
    metric = c("pipeline", "lambda", "n_samples", "n_cpgs", "n_covariates_used", "covariates_used", "n_sv_used", "sv_used",
               "batch_sig_p_lt_0.05_before", "batch_min_p_before",
               "batch_sig_p_lt_0.05_after", "batch_min_p_after",
               "group_min_p_before", "group_min_p_after",
               "dropped_covariates",
               "perm_mean_sig", "perm_max_sig", "vp_primary_group_mean"),
    value = c(prefix, lambda_val, nrow(targets), nrow(betas), length(used_covariates),
              paste(used_covariates, collapse=";"),
              length(sv_cols), paste(sv_cols, collapse=";"),
              batch_before$sig_count, batch_before$min_p,
              batch_after$sig_count, batch_after$min_p,
              grp_min_before, grp_min_after,
              paste(unique(drop_log$Variable), collapse=";"),
              if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat=="mean_sig"] else NA,
              if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat=="max_sig"] else NA,
              vp_primary)
  )
  metrics_path <- file.path(out_dir, paste0(prefix, "_Metrics.csv"))
  write.csv(metrics, metrics_path, row.names = FALSE)
  
  return(res)
}

res_minfi <- run_pipeline(beta_minfi, "Minfi", anno_df)
res_sesame <- run_pipeline(beta_sesame, "Sesame", anno_df)
beta_intersect <- beta_minfi[common_cpgs, ]
res_intersect <- run_pipeline(beta_intersect, "Intersection", anno_df)

get_sig_counts <- function(res_df) {
  up <- sum(res_df$adj.P.Val < pval_thresh & res_df$logFC > lfc_thresh, na.rm=TRUE)
  down <- sum(res_df$adj.P.Val < pval_thresh & res_df$logFC < -lfc_thresh, na.rm=TRUE)
  return(list(up=up, down=down))
}

minfi_counts <- get_sig_counts(res_minfi)
sesame_counts <- get_sig_counts(res_sesame)
intersect_counts <- get_sig_counts(res_intersect)

summary_json <- sprintf(
'{
  "n_con": %d,
  "n_test": %d,
  "minfi_up": %d,
  "minfi_down": %d,
  "sesame_up": %d,
  "sesame_down": %d,
  "intersect_up": %d,
  "intersect_down": %d
}', 
  n_con, n_test,
  minfi_counts$up, minfi_counts$down,
  sesame_counts$up, sesame_counts$down,
  intersect_counts$up, intersect_counts$down
)

writeLines(summary_json, file.path(out_dir, "summary.json"))
message("Summary statistics saved to summary.json")

message("Saving session info...")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Analysis Complete.")
