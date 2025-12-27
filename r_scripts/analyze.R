#!/usr/bin/env Rscript

force_single_thread <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  options(mc.cores = 1)
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  }
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    RcppParallel::setThreadOptions(numThreads = 1)
  }
}

force_single_thread()

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
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(pvca))
suppressPackageStartupMessages(library(RefFreeEWAS))
suppressPackageStartupMessages(library(illuminaio))

# Workaround for variancePartition / lme4 compatibility issue
# variancePartition still calls lme4::findbars; in newer lme4 this may be a stub
# that errors and instructs to use reformulas. Patch lme4::findbars to point to
# reformulas::findbars when needed.
if (requireNamespace("lme4", quietly = TRUE) && requireNamespace("reformulas", quietly = TRUE)) {
  lme4_ns <- asNamespace("lme4")
  fb_lme4 <- tryCatch(get("findbars", envir = lme4_ns), error = function(e) NULL)
  fb_need_patch <- is.null(fb_lme4) || !identical(fb_lme4, reformulas::findbars)
  if (fb_need_patch) {
    tryCatch(unlockBinding("findbars", lme4_ns), error = function(e) {})
    assign("findbars", reformulas::findbars, envir = lme4_ns)
    tryCatch(lockBinding("findbars", lme4_ns), error = function(e) {})
    message("  [Patch] Set lme4::findbars <- reformulas::findbars to support variancePartition.")
  }
}

# variancePartition may import the stubbed findbars at load time; patch its namespace too.
if (requireNamespace("variancePartition", quietly = TRUE) && requireNamespace("reformulas", quietly = TRUE)) {
  vp_ns <- asNamespace("variancePartition")
  vp_imp <- parent.env(vp_ns) # imported bindings live in the namespace parent
  if (exists("findbars", envir = vp_imp, inherits = FALSE)) {
    fb_vp <- tryCatch(get("findbars", envir = vp_imp, inherits = FALSE), error = function(e) NULL)
    fb_need_patch <- !is.null(fb_vp) && !identical(fb_vp, reformulas::findbars)
    if (fb_need_patch) {
      tryCatch(unlockBinding("findbars", vp_imp), error = function(e) {})
      tryCatch(assign("findbars", reformulas::findbars, envir = vp_imp), error = function(e) {})
      tryCatch(lockBinding("findbars", vp_imp), error = function(e) {})
      message("  [Patch] Set variancePartition::findbars <- reformulas::findbars (compat fix).")
    }
  }
}

if (!exists("findbars")) {
  findbars <- reformulas::findbars
}

# Statistical Analysis Pipeline References:
# - Normalization: Triche et al. (2013) Nucleic Acids Res - Noob method
# - M-value transformation: Du et al. (2010) BMC Bioinformatics
# - Differential methylation: Ritchie et al. (2015) Nucleic Acids Res - limma
# - Batch effect correction: Leek et al. (2012) Nat Rev Genet - SVA
# - DMR analysis: Suderman et al. (2018) Bioinformatics - dmrff
# - Multiple testing: Benjamini & Hochberg (1995) JRSS-B - FDR

option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL, 
              help="Path to configure.tsv", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="Output directory for results", metavar="character"),
  make_option(c("--idat_dir"), type="character", default="",
              help="Path to IDAT directory (default: [config dir]/idat)", metavar="character"),
  make_option(c("--skip-sesame"), action="store_true", default=FALSE,
              help="Skip Sesame pipeline and intersection outputs (Minfi only)."),
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
  make_option(c("--include_clock_covariates"), action="store_true", default=FALSE,
              help="Compute epigenetic clocks and add them as metadata covariates for auto selection (writes *_configure_with_clocks.tsv)."),
  make_option(c("--permutations"), type="integer", default=0,
              help="Number of label permutations for null DMP counts (0 to skip)"),
  make_option(c("--vp_top"), type="integer", default=5000,
              help="Number of top-variable CpGs for variancePartition (default: 5000)"),
  make_option(c("--tissue"), type="character", default="Auto", 
              help="Tissue type for cell deconvolution (default: Auto = reference-free). Supported: Auto (RefFreeEWAS), Blood, CordBlood, DLPFC"),
  make_option(c("--positive_controls"), type="character", default=NULL,
              help="Comma-separated list of gene symbols (e.g. 'AHRR,CYP1A1') to verify as positive controls."),
  make_option(c("--id_column"), type="character", default="",
              help="Column in configure.tsv to treat as sample ID (for non-GEO datasets). If empty, tries GSM/geo_accession."),
  make_option(c("--min_total_size"), type="integer", default=6,
              help="Minimum total sample size required to proceed (default: 6)."),
  make_option(c("--qc_intensity_threshold"), type="double", default=10.5,
              help="Median M/U signal intensity threshold for sample QC (log2). Set <=0 to disable intensity-based sample drop (default: 10.5)"),
  make_option(c("--force_idat"), action="store_true", default=FALSE,
              help="Force reading IDATs when array sizes differ but types are similar (passes force=TRUE to read.metharray).")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Timestamped messages for all logging in this script
ts_message <- function(..., domain = NULL, appendLF = TRUE) {
  base::message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " ")),
                domain = domain, appendLF = appendLF)
}
message <- ts_message

# QC and filtering thresholds (kept as constants for transparency/reuse)
QC_MEDIAN_INTENSITY_THRESHOLD <- 10.5
QC_DETECTION_P_THRESHOLD <- 0.01
QC_SAMPLE_DET_FAIL_FRAC <- 0.05
SNP_MAF_THRESHOLD <- 0.01
LOGIT_OFFSET <- 1e-4
BETA_RANGE_MIN <- 0.05
AUTO_COVARIATE_ALPHA <- 0.01
MAX_PCS_FOR_COVARIATE_DETECTION <- 5
OVERCORRECTION_GUARD_RATIO <- 0.25
SV_GROUP_P_THRESHOLD <- 1e-6
SV_GROUP_ETA2_THRESHOLD <- 0.5
PVCA_MIN_SAMPLES <- 8
DMR_MAXGAP <- 500
DMR_P_CUTOFF <- 0.05
BATCH_EVAL_TOP_VAR <- 20000
MIN_GROUP_SIZE_WARN <- 3
MIN_TOTAL_SIZE_STOP <- NULL

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
include_clock_covariates <- opt$include_clock_covariates
cell_covariates <- character(0)
cell_counts_df <- NULL
id_col_override <- opt$id_column
MIN_TOTAL_SIZE_STOP <- max(2, opt$min_total_size)
force_idat <- opt$force_idat
QC_MEDIAN_INTENSITY_THRESHOLD <- ifelse(opt$qc_intensity_threshold <= 0, -Inf, opt$qc_intensity_threshold)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Set ExperimentHub/AnnotationHub Cache ---
sesame_cache_dir <- file.path(getwd(), "cache") 
cache_preexisting <- dir.exists(sesame_cache_dir)
if (!cache_preexisting) dir.create(sesame_cache_dir, recursive = TRUE)

# Force BiocFileCache/Sesame to use this directory
Sys.setenv(R_USER_CACHE_DIR = sesame_cache_dir)
Sys.setenv(XDG_CACHE_HOME = sesame_cache_dir)

options(ExperimentHub.cache = sesame_cache_dir)
options(AnnotationHub.cache = sesame_cache_dir)

message(paste("Setting ExperimentHub/AnnotationHub cache to:", sesame_cache_dir))

# Ensure essential Sesame data is cached
# Only warm the cache on first run to avoid redundant downloads
if (!cache_preexisting) {
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

save_static_plot <- function(p, filename, dir, width = 7, height = 4, dpi = 300) {
  tryCatch({
    ggsave(filename = file.path(dir, filename), plot = p, width = width, height = height, dpi = dpi)
  }, error = function(e) {
    message(sprintf("  Static plot save skipped (%s): %s", filename, e$message))
  })
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

drop_zero_variance_cols <- function(mat, label = "PCA") {
  if (is.null(mat) || ncol(mat) == 0) {
    return(list(mat = mat, dropped = 0L))
  }
  var_vals <- apply(mat, 2, var, na.rm = TRUE)
  keep <- is.finite(var_vals) & var_vals > 0
  dropped <- sum(!keep)
  if (dropped > 0) {
    message(sprintf("  %s: dropped %d zero-variance features before PCA.", label, dropped))
  }
  list(mat = mat[, keep, drop = FALSE], dropped = dropped)
}

safe_prcomp <- function(mat, label = "PCA", center = TRUE, scale. = TRUE) {
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
    message(sprintf("  %s: not enough data for PCA; skipping.", label))
    return(NULL)
  }
  if (scale.) {
    res <- drop_zero_variance_cols(mat, label = label)
    mat <- res$mat
    if (is.null(mat) || ncol(mat) < 2) {
      message(sprintf("  %s: fewer than 2 variable features; skipping PCA.", label))
      return(NULL)
    }
  }
  prcomp(mat, center = center, scale. = scale.)
}

# M-value transformation with a small offset to avoid log(0)
# Following Du et al. (2010) BMC Bioinformatics approach
logit_offset <- function(betas, offset = LOGIT_OFFSET) {
  log2((betas + offset) / (1 - betas + offset))
}

auto_detect_covariates <- function(targets, gsm_col, pca_scores, alpha=AUTO_COVARIATE_ALPHA, max_pcs=MAX_PCS_FOR_COVARIATE_DETECTION) {
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
    is_num <- is.numeric(val) || is.integer(val)
    if (uniq < 2) {
      drops <- rbind(drops, data.frame(Variable=cv, Reason="constant_or_single_level"))
      next
    }
    if (!is_num && uniq >= (n - 1)) {
      drops <- rbind(drops, data.frame(Variable=cv, Reason="high_cardinality_vs_samples"))
      next
    }
    if (!is_num && (is.character(val) || is.factor(val))) {
      lvl_counts <- table(val)
      if (any(lvl_counts < 2)) {
        drops <- rbind(drops, data.frame(Variable=cv, Reason="rare_level"))
        next
      }
    }
    # Check confounding with group (any level exclusive to one group or strong chi-square)
    if (!is_num && (is.character(val) || is.factor(val))) {
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
    } else if (is_num) {
      cc <- complete.cases(val, group_vals)
      if (sum(cc) >= 3) {
        gv <- group_vals[cc]
        vv <- val[cc]
        if (length(unique(gv)) >= 2) {
          p_grp <- tryCatch({
            if (length(unique(gv)) == 2) {
              t.test(vv ~ gv)$p.value
            } else {
              kruskal.test(vv ~ as.factor(gv))$p.value
            }
          }, error=function(e) NA)
          if (!is.na(p_grp) && p_grp < 1e-6) {
            drops <- rbind(drops, data.frame(Variable=cv, Reason="confounded_with_group"))
            next
          }
        }
      }
    }
    keep <- c(keep, cv)
  }
  list(keep = keep, dropped = drops)
}

rank_covariates <- function(targets, covariates, covar_log_df) {
  if (length(covariates) == 0) return(covariates)
  log_p <- numeric(0)
  if (!is.null(covar_log_df) && nrow(covar_log_df) > 0 && all(c("Variable", "MinP") %in% colnames(covar_log_df))) {
    log_p <- setNames(covar_log_df$MinP, covar_log_df$Variable)
  }
  priorities <- vapply(covariates, function(cv) {
    if (cv %in% names(log_p)) {
      return(-as.numeric(log_p[[cv]]))
    }
    val <- targets[[cv]]
    if (is.numeric(val) || is.integer(val)) {
      v <- stats::var(val, na.rm = TRUE)
      return(ifelse(is.finite(v), v, 0))
    }
    lvls <- length(unique(val))
    if (lvls == 0) return(0)
    1 / lvls
  }, numeric(1))
  covariates[order(priorities, decreasing = TRUE)]
}

cap_covariates_for_overcorrection <- function(targets, covariates, forced_covariates, covar_log_df, ratio = OVERCORRECTION_GUARD_RATIO) {
  max_terms <- max(1, floor(nrow(targets) * ratio))
  if (length(covariates) <= max_terms) {
    return(list(keep = covariates, dropped = character(0), cap = max_terms, forced_dropped = character(0)))
  }
  ranked_covars <- rank_covariates(targets, covariates, covar_log_df)
  forced_keep <- intersect(forced_covariates, covariates)
  forced_dropped <- character(0)
  if (length(forced_keep) > max_terms) {
    forced_dropped <- forced_keep[(max_terms + 1):length(forced_keep)]
    forced_keep <- forced_keep[1:max_terms]
  }
  remaining <- max_terms - length(forced_keep)
  if (remaining < 0) remaining <- 0
  keep_covars <- unique(c(forced_keep, head(setdiff(ranked_covars, forced_keep), remaining)))
  keep_covars <- keep_covars[1:min(length(keep_covars), max_terms)]
  drop_covars <- setdiff(covariates, keep_covars)
  list(keep = keep_covars, dropped = drop_covars, cap = max_terms, forced_dropped = forced_dropped)
}

filter_sv_by_group <- function(targets, sv_cols, group_col,
                               p_thresh = SV_GROUP_P_THRESHOLD,
                               eta2_thresh = SV_GROUP_ETA2_THRESHOLD) {
  drops <- data.frame(Variable=character(0), Reason=character(0), stringsAsFactors = FALSE)
  if (length(sv_cols) == 0) {
    return(list(keep = sv_cols, dropped = drops))
  }
  group_vals <- targets[[group_col]]
  for (sv in sv_cols) {
    if (!(sv %in% colnames(targets))) next
    val <- targets[[sv]]
    cc <- complete.cases(val, group_vals)
    if (sum(cc) < 3) next
    gv <- group_vals[cc]
    if (length(unique(gv)) < 2) next
    vv <- val[cc]
    if (sd(vv, na.rm = TRUE) == 0) {
      drops <- rbind(drops, data.frame(Variable = sv, Reason = "sv_constant"))
      next
    }
    p_grp <- tryCatch({
      if (length(unique(gv)) == 2) {
        t.test(vv ~ gv)$p.value
      } else {
        kruskal.test(vv ~ as.factor(gv))$p.value
      }
    }, error = function(e) NA)
    eta2 <- tryCatch({
      fit <- stats::lm(vv ~ as.factor(gv))
      an <- stats::anova(fit)
      if (nrow(an) >= 2) an[1, "Sum Sq"] / sum(an[, "Sum Sq"]) else NA
    }, error = function(e) NA)
    if ((!is.na(p_grp) && p_grp < p_thresh) || (!is.na(eta2) && eta2 > eta2_thresh)) {
      drops <- rbind(drops, data.frame(Variable = sv, Reason = "sv_confounded_with_group"))
    }
  }
  list(keep = setdiff(sv_cols, drops$Variable), dropped = drops)
}

filter_low_range <- function(betas, min_range = BETA_RANGE_MIN) {
  rng <- apply(betas, 1, function(x) {
    r <- range(x, na.rm = TRUE)
    if (any(!is.finite(r))) return(0)
    diff(r)
  })
  keep <- rng >= min_range
  list(mat = betas[keep, , drop = FALSE], removed = sum(!keep))
}

estimate_cell_counts_safe <- function(rgSet, tissue = "Blood", out_dir = NULL) {
  run_ref_free <- function() {
    message("Estimating cell composition using Reference-Free method (RefFreeEWAS)...")
    if (!requireNamespace("RefFreeEWAS", quietly = TRUE)) {
      message("  - RefFreeEWAS package not installed. Skipping.")
      return(NULL)
    }
    
    tryCatch({
      # RefFreeEWAS requires a beta matrix. We preprocess efficiently here.
      # using minfi's preprocessNoob for consistency with main pipeline
      mSet <- preprocessNoob(rgSet)
      betas <- getBeta(mSet)
      betas <- betas[rowSums(is.na(betas)) == 0, ] # Remove NAs
      
      # Select variable probes to speed up and improve signal
      # Top 10,000 variable probes is usually sufficient for deconvolution
      vars <- apply(betas, 1, var)
      top_idx <- head(order(vars, decreasing = TRUE), 10000)
      betas_sub <- betas[top_idx, ]
      
      # Run RefFreeCellMix
      # K=5 is a reasonable default for many tissues (e.g. blood has ~5-6 majors)
      # Estimating K can be unstable on small datasets, so fixed K is safer for automation.
      K_latent <- 5
      message(sprintf("  - Running RefFreeCellMix with K=%d latent components...", K_latent))
      
      # Initialize with K-means (standard RefFreeEWAS approach)
      res <- RefFreeEWAS::RefFreeCellMix(betas_sub, K = K_latent, verbose = FALSE)
      
      # Extract proportions (Omega)
      props <- res$Omega
      df <- as.data.frame(props)
      colnames(df) <- paste0("Cell_Latent", 1:ncol(df))
      df$SampleID <- rownames(df)
      
      if (!is.null(out_dir)) {
        write.csv(df, file.path(out_dir, "cell_counts_RefFree.csv"), row.names = FALSE)
      }
      
      message("  - Reference-Free estimation complete.")
      return(list(counts = df, reference = "RefFreeEWAS"))
      
    }, error = function(e) {
      message("  - Reference-Free estimation failed: ", e$message)
      return(NULL)
    })
  }

  # 1. Reference-Free Mode (Auto)
  if (tissue == "Auto") {
    return(run_ref_free())
  }

  # 2. Reference-Based Mode
  # Map tissue types to potential reference packages (prioritized order)
  ref_map <- list(
    "Blood" = c("FlowSorted.Blood.EPIC", "FlowSorted.Blood.450k"),
    "CordBlood" = c("FlowSorted.CordBlood.EPIC", "FlowSorted.CordBlood.450k", "FlowSorted.CordBloodCombined.450k"),
    "DLPFC" = c("FlowSorted.DLPFC.450k"), # Brain
    "Saliva" = c("FlowSorted.Saliva.450k") # Rare but exists
  )
  
  is_epicv2 <- FALSE
  if (requireNamespace("minfi", quietly = TRUE) &&
      exists(".isEPICv2", envir = asNamespace("minfi"), inherits = FALSE)) {
    is_epicv2 <- tryCatch(get(".isEPICv2", asNamespace("minfi"))(rgSet), error = function(e) FALSE)
  }
  if (is_epicv2) {
    message(sprintf("Cell composition references do not support EPIC v2; falling back to Auto (reference-free) for tissue '%s'.", tissue))
    return(run_ref_free())
  }

  refs <- ref_map[[tissue]]
  if (is.null(refs)) {
    message(sprintf("Reference-based cell composition skipped: No known references for tissue '%s'. (Supported: %s, or use 'Auto' for reference-free)", 
                    tissue, paste(names(ref_map), collapse=", ")))
    message("  - Falling back to reference-free (Auto) estimation.")
    return(run_ref_free())
  }

  minfi_est_fun <- NULL
  if (requireNamespace("minfi", quietly = TRUE)) {
    if (exists("estimateCellCounts2", envir = asNamespace("minfi"), inherits = FALSE)) {
      minfi_est_fun <- minfi::estimateCellCounts2
    } else if (exists("estimateCellCounts", envir = asNamespace("minfi"), inherits = FALSE)) {
      minfi_est_fun <- minfi::estimateCellCounts
    }
  }

  for (ref in refs) {
    if (!requireNamespace(ref, quietly = TRUE)) {
      # message(sprintf("  - Reference package '%s' not installed. Skipping...", ref))
      next
    }
    
    # Check if reference supports the query array type (roughly)
    # EPIC array usually works with 450k references via projection, but explicit EPIC refs are better.
    
    message(sprintf("Estimating cell composition for tissue '%s' using %s...", tissue, ref))
    
    # Determine reference platform string for minfi
    # Usually inferred by the package name, but estimateCellCounts wants 'referencePlatform' sometimes
    ref_platform <- if (grepl("EPIC", ref, ignore.case = TRUE)) "IlluminaHumanMethylationEPIC" else "IlluminaHumanMethylation450k"
    
    est_fun <- minfi_est_fun
    if (ref == "FlowSorted.Blood.EPIC" &&
        exists("estimateCellCounts2", envir = asNamespace("FlowSorted.Blood.EPIC"), inherits = FALSE)) {
      est_fun <- FlowSorted.Blood.EPIC::estimateCellCounts2
    }
    
    if (is.null(est_fun)) {
      message("Cell composition skipped (minfi::estimateCellCounts[2] not available).")
      return(NULL)
    }
    
    res <- tryCatch({
      arg_list <- list(
        rgSet = rgSet,
        compositeCellType = tissue,   # This must match the key in the FlowSorted package (usually "Blood", "DLPFC")
        referencePlatform = ref_platform,
        returnAll = FALSE
      )
      
      # Some specialized packages might use slightly different composite names
      # e.g. CordBlood might still expect "Blood" or specific subtypes. 
      # Standard minfi packages align compositeCellType with the tissue name in the package (e.g. "Blood", "DLPFC", "CordBlood")
      
      fn_args <- names(formals(est_fun))
      if ("processMethod" %in% fn_args) arg_list$processMethod <- "preprocessNoob"
      if ("normalizationMethod" %in% fn_args) arg_list$normalizationMethod <- "none"
      if ("meanPlot" %in% fn_args) arg_list$meanPlot <- FALSE
      
      do.call(est_fun, arg_list)
    }, error = function(e) {
      message("  - Cell composition estimation failed for ", ref, ": ", e$message)
      NULL
    })
    
    if (is.null(res)) next
    
    if (is.list(res) && !is.null(res$prop)) {
      res <- res$prop
    }
    df <- as.data.frame(res)
    colnames(df) <- paste0("Cell_", colnames(df))
    df$SampleID <- rownames(df)
    
    if (!is.null(out_dir)) {
      write.csv(df, file.path(out_dir, paste0("cell_counts_", ref, ".csv")), row.names = FALSE)
    }
    message(sprintf("  - Cell composition estimated using %s.", ref))
    return(list(counts = df, reference = ref))
  }
  
  message(sprintf("Reference-based cell composition skipped (all references for '%s' failed or missing).", tissue))
  message("  - Falling back to reference-free (Auto) estimation.")
  return(run_ref_free())
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

prepare_pvca_meta <- function(meta, factors) {
  if (length(factors) == 0) return(NULL)
  keep <- intersect(factors, colnames(meta))
  if (length(keep) == 0) return(NULL)
  out <- data.frame(row.names = rownames(meta))
  for (nm in keep) {
    val <- meta[[nm]]
    if (is.numeric(val)) {
      if (length(unique(val)) < 2) next
      # Robust binning for PVCA
      # 1. Try 5 bins
      breaks <- unique(quantile(val, probs = seq(0, 1, length.out = 5), na.rm = TRUE))
      # 2. If distinct values are few, quantiles might collapse. Ensure at least 2 breaks.
      if (length(breaks) < 2) {
         # Fallback: just min/max if they differ
         rng <- range(val, na.rm=TRUE)
         if (diff(rng) > 1e-8) breaks <- rng else next
      }
      if (length(breaks) == 1) next # still one point? skip
      
      # 3. Cut
      val <- tryCatch(
        cut(val, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE),
        error = function(e) NULL
      )
      if (is.null(val)) next
    } else {
      val <- as.factor(val)
    }
    val <- droplevels(val)
    if (length(levels(val)) < 2) next
    if (length(levels(val)) >= nrow(meta)) next
    out[[nm]] <- val
  }
  if (ncol(out) == 0) return(NULL)
  return(out)
}

run_pvca_assessment <- function(betas, meta, factors, prefix, out_dir, threshold = 0.6, max_probes = 5000, sample_ids = NULL) {
  pvca_meta <- prepare_pvca_meta(meta, factors)
  if (is.null(pvca_meta)) {
    message("  PVCA skipped: no usable factors or insufficient levels.")
    return(NULL)
  }
  if (nrow(pvca_meta) < PVCA_MIN_SAMPLES) {
    message(sprintf("  PVCA skipped: sample size too small for stable estimation (n=%d).", nrow(pvca_meta)))
    return(NULL)
  }
  if (!is.null(sample_ids) && nrow(pvca_meta) == length(sample_ids)) {
    rownames(pvca_meta) <- sample_ids
  }
  smp <- colnames(betas)
  if (!is.null(rownames(pvca_meta)) && all(smp %in% rownames(pvca_meta))) {
    pvca_meta <- pvca_meta[smp, , drop = FALSE]
  } else {
    rownames(pvca_meta) <- smp
  }
  betas <- pmin(pmax(betas, LOGIT_OFFSET), 1 - LOGIT_OFFSET)
  level_counts <- vapply(pvca_meta, function(x) length(levels(as.factor(x))), integer(1))
  # Drop highest-cardinality non-group factors until the crossed design fits the sample size
  while (length(level_counts) > 0 && prod(level_counts) >= nrow(pvca_meta) && length(level_counts) > 1) {
    drop_candidates <- setdiff(names(level_counts), "primary_group")
    if (length(drop_candidates) == 0) break
    drop_term <- drop_candidates[which.max(level_counts[drop_candidates])]
    pvca_meta[[drop_term]] <- NULL
    level_counts <- vapply(pvca_meta, function(x) length(levels(as.factor(x))), integer(1))
    message(sprintf("  PVCA: dropping '%s' to avoid over-specified crossing of factors.", drop_term))
  }
  if (length(level_counts) == 0) {
    message("  PVCA skipped: no factors left after dropping high-cardinality terms.")
    return(NULL)
  }
  if (prod(level_counts) >= nrow(pvca_meta)) {
    message("  PVCA skipped: crossed factor levels still exceed sample size.")
    return(NULL)
  }
  total_levels <- sum(vapply(pvca_meta, function(x) length(levels(as.factor(x))), integer(1)))
  if (total_levels >= nrow(pvca_meta)) {
    message("  PVCA skipped: total factor levels across terms exceed sample size.")
    return(NULL)
  }
  n_features <- nrow(betas)
  top_n <- min(max_probes, n_features)
  vars <- apply(betas, 1, var)
  top_idx <- head(order(vars, decreasing = TRUE), top_n)
  betas_use <- logit_offset(betas[top_idx, , drop=FALSE])
  pheno <- AnnotatedDataFrame(data = pvca_meta)
  eset <- ExpressionSet(assayData = as.matrix(betas_use), phenoData = pheno)
  
  pvca_res <- tryCatch(pvcaBatchAssess(eset, colnames(pvca_meta), threshold), error = function(e) e)
  if (inherits(pvca_res, "error") || !is.list(pvca_res) || is.null(pvca_res$dat) || is.null(pvca_res$label)) {
    message("  PVCA failed: ", if (inherits(pvca_res, "error")) pvca_res$message else "invalid PVCA result")
    return(NULL)
  }
  pvca_df <- data.frame(term = pvca_res$label, proportion = as.numeric(pvca_res$dat))
  if (nrow(pvca_df) == 0 || any(is.na(pvca_df$proportion))) {
    message("  PVCA failed: empty or invalid PVCA results.")
    return(NULL)
  }
  out_path <- file.path(out_dir, paste0(prefix, "_PVCA.csv"))
  write.csv(pvca_df, out_path, row.names = FALSE)
  message(paste("  PVCA saved to", out_path))
  # Try to render interactive bar chart; also emit a static PNG fallback for dashboard previews
  pvca_plot_name <- paste0(prefix, "_PVCA.html")
  pvca_png_name <- paste0(prefix, "_PVCA.png")
  tryCatch({
    p_pvca <- ggplot(pvca_df, aes(x = reorder(term, proportion), y = proportion, fill = term,
                                  text = paste0(term, ": ", scales::percent(proportion, accuracy = 0.1)))) +
      geom_col() +
      coord_flip() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Factor", y = "Proportion of Variance", title = paste0(prefix, " PVCA")) +
      theme_minimal() +
      theme(legend.position = "none")
    save_interactive_plot(p_pvca, pvca_plot_name, out_dir)
    tryCatch({
      ggsave(filename = file.path(out_dir, pvca_png_name), plot = p_pvca, width = 7, height = 4, dpi = 120)
    }, error = function(e) {
      message("  PVCA PNG save skipped: ", e$message)
    })
  }, error = function(e) {
    message("  PVCA plot failed: ", e$message)
  })
}

planet_predict_age_safe <- function(betas, type = "RPC") {
  ageCpGs <- NULL
  tryCatch({
    data("ageCpGs", package = "planet", envir = environment())
  }, error = function(e) {
    stop("planet ageCpGs dataset not available: ", e$message)
  })
  if (is.null(ageCpGs)) {
    stop("planet ageCpGs dataset not available.")
  }
  if (!(type %in% c("RPC", "CPC", "RRPC"))) {
    stop("Type must be one of \"CPC\", \"RPC\", or \"RRPC\"")
  }
  coef <- ageCpGs[ageCpGs[[type]] != 0, , drop = FALSE]
  cpgs <- intersect(coef$CpGs[coef$CpGs != "(Intercept)"], rownames(betas))
  if (length(cpgs)/(nrow(coef) - 1) < 0.5) {
    stop("Less than 50% of predictors were found.")
  }
  if (length(cpgs) < (nrow(coef) - 1)) {
    warning(paste("Only", length(cpgs), "out of", (nrow(coef) - 1), "predictors present."))
  } else {
    message(paste(length(cpgs), "of", (nrow(coef) - 1), "predictors present."))
  }
  coef_filtered <- coef[coef$CpGs %in% c("(Intercept)", cpgs), , drop = FALSE][[type]]
  betas_sub <- as.matrix(betas[cpgs, , drop = FALSE])
  betas_filtered <- cbind(1, t(betas_sub))
  as.vector(betas_filtered %*% coef_filtered)
}

compute_epigenetic_clocks <- function(betas, meta, id_col, prefix, out_dir, tissue = "Auto", return_covariates = FALSE) {
  # Safe ID extractor
  safe_id <- function(x) sub("\\.idat.*$", "", basename(as.character(x)))
  sample_ids <- safe_id(colnames(betas))
  clock_covariates <- NULL
  if (return_covariates) {
    clock_covariates <- data.frame(SampleID = sample_ids, stringsAsFactors = FALSE)
  }
  
  # --- 1. methylclock (Default / General Clocks) ---
  # Prioritize methylclock as it is more robust and supports newer clocks
  mc_success <- FALSE
  if (requireNamespace("methylclock", quietly = TRUE)) {
     message("  Computing epigenetic clocks using 'methylclock'...")
     tryCatch({
       mc_res <- methylclock::DNAmAge(betas, toBetas=FALSE, fastImp=TRUE, normalize=FALSE)
       mc_df <- as.data.frame(mc_res)
       
       # Align Sample IDs
       if (nrow(mc_df) == ncol(betas)) {
           mc_df$Sample <- colnames(betas)
       } else {
           mc_df$Sample <- sample_ids # Fallback
       }
       
       # Reorder columns
       mc_df <- mc_df[, c("Sample", setdiff(colnames(mc_df), "Sample"))]
       
       out_mc_path <- file.path(out_dir, paste0(prefix, "_Epigenetic_Age_methylclock.csv"))
       write.csv(mc_df, out_mc_path, row.names = FALSE)
       message("  Epigenetic clocks (methylclock) saved to ", out_mc_path)
       mc_success <- TRUE
       
       if (!is.null(clock_covariates)) {
         mc_cols <- setdiff(colnames(mc_df), "Sample")
         for (cc in mc_cols) {
           vals <- suppressWarnings(as.numeric(mc_df[[cc]]))
           if (all(is.na(vals))) next
           clock_covariates[[paste0("Clock_", cc)]] <- vals
         }
       }
       
       # Visualization
       meta_ids <- NULL
       if (!is.null(id_col) && nzchar(id_col) && id_col %in% colnames(meta)) {
           meta_ids <- safe_id(meta[[id_col]])
       } else if (!is.null(rownames(meta))) {
           meta_ids <- safe_id(rownames(meta))
       }
       idx <- if (!is.null(meta_ids)) match(sample_ids, meta_ids) else seq_along(sample_ids)
       age_candidates <- c("chronological_age", "Age", "age", "age_years", "Age_years")
       age_col <- NULL
       meta_cols_lower <- tolower(colnames(meta))
       for (cand in age_candidates) {
           hit <- which(meta_cols_lower == tolower(cand))
           if (length(hit) > 0) {
               age_col <- colnames(meta)[hit[1]]
               break
           }
       }
       
       if (!is.null(age_col) && "Horvath" %in% colnames(mc_df)) {
           raw_age <- meta[[age_col]][idx]
           # Force conversion to numeric (handles factors/characters like "45", "45 years" if clean)
           chron_age <- suppressWarnings(as.numeric(as.character(raw_age)))
           
           if (any(!is.na(chron_age)) && all(is.finite(chron_age[!is.na(chron_age)]))) {
               # Subset to valid samples only
               valid_mask <- !is.na(chron_age) & is.finite(chron_age)
               
               if (sum(valid_mask) >= 2) {
                   plot_df <- data.frame(
                     Sample = mc_df$Sample[valid_mask],
                     Chronological_Age = chron_age[valid_mask],
                     Horvath = mc_df$Horvath[valid_mask],
                     stringsAsFactors = FALSE
                   )
                   
                   p_mc <- ggplot(plot_df, aes(x = Chronological_Age, y = Horvath,
                                               text = paste0("Sample: ", Sample,
                                                             "<br>Horvath: ", round(Horvath, 2),
                                                             "<br>Chronological: ", round(Chronological_Age, 2)))) +
                     geom_point(size = 3, alpha = 0.7, color="purple") +
                     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
                     theme_minimal() +
                     labs(x = "Chronological Age", y = "Horvath DNAm Age", title = paste0(prefix, " Epigenetic Clock (methylclock)"))
                   
                   save_interactive_plot(p_mc, paste0(prefix, "_Epigenetic_Age_methylclock.html"), out_dir)
                   save_static_plot(p_mc, paste0(prefix, "_Epigenetic_Age_methylclock.png"), out_dir, width = 5, height = 4)
               }
           }
       }
     }, error = function(e) {
       message("  methylclock failed: ", e$message)
     })
  }

  # --- 2. wateRmelon (Legacy / Fallback) ---
  # Only run if methylclock failed or for comparison
  if (!mc_success && requireNamespace("wateRmelon", quietly = TRUE)) {
    betas_use <- pmin(pmax(betas, LOGIT_OFFSET), 1 - LOGIT_OFFSET)
    tryCatch({
      clock_res <- wateRmelon::agep(betas_use, method = "all")
      # ... (rest of wateRmelon logic omitted for brevity as methylclock is now default) ...
      # Saving logic would go here if we wanted strict fallback
      message("  wateRmelon fallback execution skipped (methylclock preferred).") 
    }, error = function(e) {
      # Silent fail or log
    })
  }
  
  # --- 3. planet (Placenta Only) ---
  # Only run if tissue is explicitly set to Placenta
  if (tolower(tissue) == "placenta" && requireNamespace("planet", quietly = TRUE)) {
     message("  Computing placental clocks using 'planet' (Tissue=Placenta)...")
     tryCatch({
       planet_use_fallback <- FALSE
       planet_predict <- function(type) {
         if (!planet_use_fallback) {
           res <- tryCatch(planet::predictAge(betas, type = type), error = function(e) e)
           if (!inherits(res, "error")) return(res)
           if (grepl("ageCpGs", res$message, fixed = TRUE)) {
             message("  planet::predictAge failed (", type, "): missing ageCpGs; switching to internal fallback.")
             planet_use_fallback <<- TRUE
           } else {
             message("  planet::predictAge failed (", type, "): ", res$message, " (trying fallback)")
           }
         }
         tryCatch(planet_predict_age_safe(betas, type = type), error = function(e2) {
           message("  planet fallback failed (", type, "): ", e2$message)
           rep(NA_real_, ncol(betas))
         })
       }
       age_rpc <- planet_predict("RPC")
       age_cpc <- planet_predict("CPC")
       
       planet_df <- data.frame(Sample = sample_ids)
       if (!all(is.na(age_rpc))) planet_df$RPC_Gestational_Age_Weeks <- as.numeric(age_rpc)
       if (!all(is.na(age_cpc))) planet_df$CPC_Gestational_Age_Weeks <- as.numeric(age_cpc)
       if (!is.null(clock_covariates)) {
         if ("RPC_Gestational_Age_Weeks" %in% colnames(planet_df)) {
           clock_covariates$PlacentaClock_RPC_GA_Weeks <- planet_df$RPC_Gestational_Age_Weeks
         }
         if ("CPC_Gestational_Age_Weeks" %in% colnames(planet_df)) {
           clock_covariates$PlacentaClock_CPC_GA_Weeks <- planet_df$CPC_Gestational_Age_Weeks
         }
       }
       
       out_planet_path <- file.path(out_dir, paste0(prefix, "_Placental_Age_planet.csv"))
       write.csv(planet_df, out_planet_path, row.names = FALSE)
       message("  Placental clocks (planet) saved to ", out_planet_path)
       
       if (all(c("RPC_Gestational_Age_Weeks", "CPC_Gestational_Age_Weeks") %in% colnames(planet_df))) {
         ok <- is.finite(planet_df$RPC_Gestational_Age_Weeks) & is.finite(planet_df$CPC_Gestational_Age_Weeks)
         if (sum(ok) >= 2) {
           p_planet <- ggplot(planet_df[ok, , drop = FALSE],
                              aes(x = RPC_Gestational_Age_Weeks, y = CPC_Gestational_Age_Weeks,
                                  text = paste0("Sample: ", Sample,
                                                "<br>RPC: ", round(RPC_Gestational_Age_Weeks, 2),
                                                "<br>CPC: ", round(CPC_Gestational_Age_Weeks, 2)))) +
             geom_point(size = 3, alpha = 0.8, color = "#2c3e50") +
             geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
             theme_minimal() +
             labs(
               title = paste0(prefix, " Placental clocks (planet)"),
               x = "RPC gestational age (weeks)",
               y = "CPC gestational age (weeks)"
             )
           save_interactive_plot(p_planet, paste0(prefix, "_Placental_Age_planet_scatter.html"), out_dir)
           save_static_plot(p_planet, paste0(prefix, "_Placental_Age_planet_scatter.png"), out_dir, width = 5, height = 4)
         }
       }
       
     }, error = function(e) {
       message("  planet clock failed: ", e$message)
     })
  }
  if (return_covariates) {
    if (is.null(clock_covariates) || ncol(clock_covariates) <= 1) return(NULL)
    return(clock_covariates)
  }
}

select_batch_factor <- function(meta, preferred = c("Sentrix_ID", "Sentrix_Position")) {
  for (nm in preferred) {
    if (nm %in% colnames(meta)) {
      vals <- meta[[nm]]
      if (length(unique(vals)) > 1 && length(unique(vals)) < nrow(meta)) {
        return(nm)
      }
    }
  }
  return(NULL)
}

is_batch_confounded <- function(meta, batch_col, group_col, p_thresh = 1e-6) {
  if (is.null(batch_col) || is.null(group_col)) return(FALSE)
  if (!(batch_col %in% colnames(meta)) || !(group_col %in% colnames(meta))) return(FALSE)
  tbl <- table(meta[[batch_col]], meta[[group_col]])
  if (any(tbl == 0)) return(TRUE)
  p_chisq <- tryCatch(suppressWarnings(chisq.test(tbl)$p.value), error = function(e) NA)
  if (!is.na(p_chisq) && p_chisq < p_thresh) return(TRUE)
  return(FALSE)
}

eval_batch_method <- function(M_mat, meta, group_col, batch_col, covariates, method) {
  meta_use <- meta
  sample_ids <- colnames(M_mat)
  if (!is.null(rownames(meta_use)) && all(sample_ids %in% rownames(meta_use))) {
    meta_use <- meta_use[sample_ids, , drop = FALSE]
  } else {
    match_col <- NA
    for (nm in colnames(meta_use)) {
      if (all(sample_ids %in% meta_use[[nm]])) {
        match_col <- nm
        break
      }
    }
    if (!is.na(match_col)) {
      idx <- match(sample_ids, meta_use[[match_col]])
      meta_use <- meta_use[idx, , drop = FALSE]
    } else if (nrow(meta_use) != length(sample_ids)) {
      stop("Sample alignment failed: metadata rows do not match expression columns.")
    }
  }
  rownames(meta_use) <- sample_ids
  meta_use[[batch_col]] <- as.factor(meta_use[[batch_col]])
  meta_use[[group_col]] <- as.factor(meta_use[[group_col]])
  cov_terms <- setdiff(covariates, batch_col)
  for (ct in cov_terms) {
    if (ct %in% colnames(meta_use) && !is.numeric(meta_use[[ct]]) && !is.integer(meta_use[[ct]])) {
      meta_use[[ct]] <- as.factor(meta_use[[ct]])
    }
  }
  batch <- meta_use[[batch_col]]
  group <- meta_use[[group_col]]
  
  # Speed: evaluate batch strategies on top-variable probes only
  if (nrow(M_mat) > BATCH_EVAL_TOP_VAR) {
    vars0 <- apply(M_mat, 1, var, na.rm = TRUE)
    eval_idx <- head(order(vars0, decreasing = TRUE), BATCH_EVAL_TOP_VAR)
    M_mat <- M_mat[eval_idx, , drop = FALSE]
  }
  
  # Helper to build model matrices safely
  mm_safe <- function(formula_str, data) {
    model.matrix(as.formula(formula_str), data = data)
  }
  
  M_corr <- M_mat
  if (method == "none") {
    M_corr <- M_mat
  } else if (method == "combat") {
    mod <- mm_safe(paste("~", paste(c(group_col, cov_terms), collapse = " + ")), meta_use)
    M_corr <- ComBat(dat = M_mat, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
  } else if (method == "limma") {
    design <- mm_safe(paste("~ 0 +", group_col), meta_use)
    cov_mat <- if (length(cov_terms) > 0) mm_safe(paste("~ 0 +", paste(cov_terms, collapse = " + ")), meta_use) else NULL
    M_corr <- removeBatchEffect(M_mat, batch = batch, covariates = cov_mat, design = design)
  } else if (method == "sva") {
    mod <- mm_safe(paste("~", paste(c(group_col, cov_terms), collapse = " + ")), meta_use)
    mod0_terms <- if (length(cov_terms) > 0) paste(cov_terms, collapse = " + ") else "1"
    mod0 <- mm_safe(paste("~", mod0_terms), meta_use)
    sva.obj <- sva(M_mat, mod, mod0)
    if (is.null(sva.obj$sv) || sva.obj$n.sv < 1) {
      M_corr <- M_mat
    } else {
      sv_mat <- sva.obj$sv
      if (nrow(sv_mat) != nrow(meta_use)) {
        message("    - SVA SV rows do not match samples; skipping SVA correction for evaluation.")
        M_corr <- M_mat
      } else {
        design_keep <- mm_safe(paste("~ 0 +", group_col), meta_use)
        M_corr <- removeBatchEffect(M_mat, covariates = sv_mat, design = design_keep)
      }
    }
  } else {
    stop("unknown method")
  }
  
  vars <- apply(M_corr, 1, var, na.rm = TRUE)
  top_idx <- head(order(vars, decreasing = TRUE), min(BATCH_EVAL_TOP_VAR, nrow(M_corr)))
  M_top_pca <- t(M_corr[top_idx, , drop = FALSE])
  pca <- safe_prcomp(M_top_pca, label = paste("Batch eval", method), center = TRUE, scale. = TRUE)
  n_pc_batch_sig <- 0
  if (!is.null(pca)) {
    pc_df <- as.data.frame(pca$x[, 1:min(10, ncol(pca$x))])
    pc_df$Batch <- batch
    pvals_pc_batch <- sapply(colnames(pc_df)[grep("^PC", colnames(pc_df))], function(pc) {
      fit <- aov(pc_df[[pc]] ~ pc_df$Batch)
      as.numeric(summary(fit)[[1]][["Pr(>F)"]][1])
    })
    fdr_pc_batch <- p.adjust(pvals_pc_batch, "BH")
    n_pc_batch_sig <- sum(fdr_pc_batch < 0.05, na.rm = TRUE)
  } else {
    message("    - Batch eval PCA skipped due to insufficient variance.")
  }
  
  form_terms <- c(sprintf("(1|%s)", batch_col))
  for (ct in cov_terms) {
    if (is.numeric(meta_use[[ct]]) || is.integer(meta_use[[ct]])) {
      form_terms <- c(form_terms, ct)
    } else {
      form_terms <- c(form_terms, sprintf("(1|%s)", ct))
    }
  }
  form_terms <- c(form_terms, sprintf("(1|%s)", group_col))
  form <- as.formula(paste("~", paste(form_terms, collapse = " + ")))
  M_top_varpart <- M_corr[top_idx, , drop = FALSE]
  varPart <- tryCatch(
    fitExtractVarPartModel(M_top_varpart, form, meta_use),
    error = function(e) {
      message("    variancePartition failed for ", method, ": ", e$message)
      NULL
    }
  )
  if (is.null(varPart)) {
    med_varFrac <- setNames(rep(0, length(c(batch_col, group_col))), c(batch_col, group_col))
  } else {
    med_varFrac <- apply(varPart, 2, median, na.rm = TRUE)
  }
  
  design_terms <- c(batch_col, group_col, cov_terms)
  design <- mm_safe(paste("~ 0 +", paste(design_terms, collapse = " + ")), meta_use)
  fit_batch <- lmFit(M_corr, design)
  fit_batch <- eBayes(fit_batch)
  batch_cols <- grep(paste0("^", batch_col), colnames(design))
  if (length(batch_cols) == 0) {
    prop_batch_sig <- NA_real_
  } else {
    contrast_batch <- diag(ncol(design))[, batch_cols, drop = FALSE]
    fit_con <- contrasts.fit(fit_batch, contrast_batch)
    fit_con <- eBayes(fit_con)
    p_batch <- fit_con$F.p.value
    fdr_batch <- p.adjust(p_batch, "BH")
    prop_batch_sig <- mean(fdr_batch < 0.05, na.rm = TRUE)
  }
  
  batch_var <- if (batch_col %in% names(med_varFrac)) med_varFrac[batch_col] else NA_real_
  group_var <- if (group_col %in% names(med_varFrac)) med_varFrac[group_col] else NA_real_
  if (is.na(batch_var)) batch_var <- 0
  if (is.na(group_var)) group_var <- 0
  if (is.na(prop_batch_sig)) prop_batch_sig <- 0
  score <- batch_var + prop_batch_sig + 0.1 * n_pc_batch_sig - 0.5 * group_var
  
  list(
    method = method,
    score = score,
    batch_var = batch_var,
    group_var = group_var,
    n_pc_batch_sig = n_pc_batch_sig,
    prop_batch_sig = prop_batch_sig,
    M_corr = M_corr
  )
}

apply_batch_correction <- function(betas, method, batch_col, covariates, targets, design) {
  # Apply the selected batch correction method to betas (logit space) while preserving group effects.
  if (is.null(method) || method == "none" || is.null(batch_col) || !(batch_col %in% colnames(targets))) {
    return(list(betas = betas, method = "none"))
  }
  batch_vals <- targets[[batch_col]]
  if (length(unique(batch_vals)) < 2) {
    return(list(betas = betas, method = "none"))
  }
  covariates <- covariates[covariates %in% colnames(targets)]
  cov_mat <- NULL
  if (length(covariates) > 0) {
    cov_formula <- as.formula(paste("~ 0 +", paste(covariates, collapse = " + ")))
    cov_mat <- model.matrix(cov_formula, data = targets)
  }
  design_keep <- model.matrix(~ primary_group, data = targets)
  M_mat <- logit_offset(betas)

  if (method == "combat") {
    mod <- design
    M_corr <- ComBat(dat = M_mat, batch = as.factor(batch_vals), mod = mod, par.prior = TRUE, prior.plots = FALSE)
  } else if (method == "limma") {
    M_corr <- removeBatchEffect(M_mat, batch = as.factor(batch_vals), covariates = cov_mat, design = design_keep)
  } else if (method == "sva") {
    # SVA handled via surrogate variables already included in the design.
    return(list(betas = betas, method = "sva"))
  } else {
    return(list(betas = betas, method = "none"))
  }

  betas_corr <- 2^M_corr / (1 + 2^M_corr)
  betas_corr <- pmin(pmax(betas_corr, LOGIT_OFFSET), 1 - LOGIT_OFFSET)
  list(betas = betas_corr, method = method)
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
idat_dir <- opt$idat_dir
if (!nzchar(idat_dir)) {
  idat_dir <- file.path(project_dir, "idat")
} else if (!grepl("^/", idat_dir)) {
  idat_dir <- file.path(project_dir, idat_dir)
}
if (!dir.exists(idat_dir)) {
  stop(paste("IDAT directory not found:", idat_dir))
}

files <- list.files(idat_dir, pattern = "_Grn.idat", full.names = TRUE)
basenames <- unique(sub("_Grn.idat.*", "", files))

# Resolve sample ID column (user override > GSM/geo_accession > Basename)
if (nzchar(id_col_override)) {
  if (!id_col_override %in% colnames(targets)) {
    stop(paste("Specified --id_column", id_col_override, "not found in configure.tsv"))
  }
  id_col <- id_col_override
} else {
  id_col <- grep("GSM|geo_accession", colnames(targets), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(id_col) && "Basename" %in% colnames(targets)) {
    id_col <- "Basename"
  }
}

if (is.na(id_col) || id_col == "") {
  stop("Could not identify a sample ID column. Provide --id_column or include GSM/geo_accession/Basename in configure.tsv.")
}

targets$SampleID <- trimws(targets[[id_col]])
if (any(is.na(targets$SampleID) | targets$SampleID == "")) {
  stop("Sample ID column contains missing/empty values. Please fill it before running analysis.")
}
if (any(duplicated(targets$SampleID))) {
  stop("Duplicate sample IDs detected. Ensure Sample IDs are unique or pass a different --id_column.")
}
gsm_col <- "SampleID"

if (!"Basename" %in% colnames(targets)) {
  targets$Basename <- NA_character_
}

find_basename_for_id <- function(sample_id, basenames) {
  exact <- basenames[basename(basenames) == sample_id]
  if (length(exact) > 0) return(exact[1])
  contains <- basenames[grepl(sample_id, basename(basenames), fixed = TRUE)]
  if (length(contains) > 0) return(contains[1])
  return(NA_character_)
}

resolve_basename <- function(candidate, sample_id, basenames, project_dir, idat_dir) {
  cand <- candidate
  if (is.na(cand) || cand == "") {
    cand <- find_basename_for_id(sample_id, basenames)
  } else {
    # If relative path, anchor to project directory (fallback to IDAT dir if provided)
    if (!grepl("^/", cand)) {
      cand_project <- file.path(project_dir, cand)
      cand_idat <- file.path(idat_dir, cand)
      cand <- cand_project
      test_paths <- c(paste0(cand_project, "_Grn.idat"), paste0(cand_project, "_Grn.idat.gz"))
      if (!any(file.exists(test_paths)) && idat_dir != project_dir) {
        test_paths_idat <- c(paste0(cand_idat, "_Grn.idat"), paste0(cand_idat, "_Grn.idat.gz"))
        if (any(file.exists(test_paths_idat))) cand <- cand_idat
      }
    }
    # If file is missing, try to fall back to discovered basenames
    test_paths <- c(paste0(cand, "_Grn.idat"), paste0(cand, "_Grn.idat.gz"))
    if (!any(file.exists(test_paths))) {
      alt <- find_basename_for_id(sample_id, basenames)
      if (!is.na(alt)) cand <- alt
    }
  }
  cand
}

for (i in seq_len(nrow(targets))) {
  targets$Basename[i] <- resolve_basename(targets$Basename[i], targets$SampleID[i], basenames, project_dir, idat_dir)
}

missing_base <- which(is.na(targets$Basename))
if (length(missing_base) > 0) {
  stop(paste0("Could not find IDAT basenames for: ", paste(targets$SampleID[missing_base], collapse = ", "),
              ". Ensure the IDAT filenames contain the sample IDs or set Basename explicitly."))
}

targets <- targets[!is.na(targets$Basename), ]
if (nrow(targets) == 0) stop("No matching IDAT files found for the defined samples.")

# Drop samples whose IDAT array size deviates from the modal value (prevents mixed-platform parsing errors)
idat_size <- vapply(targets$Basename, function(bn) {
  f <- paste0(bn, "_Grn.idat")
  if (!file.exists(f)) f <- paste0(bn, "_Grn.idat.gz")
  if (!file.exists(f)) return(NA_real_)
  tryCatch({
    readIDAT(f)[["nSNPsRead"]]
  }, error = function(e) NA_real_)
}, numeric(1))
mode_size <- NA_real_
if (sum(!is.na(idat_size)) > 0) {
  tbl <- sort(table(idat_size), decreasing = TRUE)
  mode_size <- as.numeric(names(tbl)[1])
  if (length(tbl) > 1) {
    mismatch_idx <- which(!is.na(idat_size) & idat_size != mode_size)
    if (length(mismatch_idx) > 0) {
      bad_ids <- targets[[gsm_col]][mismatch_idx]
      bad_sizes <- idat_size[mismatch_idx]
      message(sprintf("Dropping %d sample(s) with non-modal IDAT array size (mode=%s): %s",
                      length(mismatch_idx), format(mode_size, scientific = FALSE),
                      paste(paste0(bad_ids, " (", bad_sizes, ")"), collapse = ", ")))
      targets <- targets[-mismatch_idx, , drop = FALSE]
      idat_size <- idat_size[-mismatch_idx]
      if (nrow(targets) == 0) stop("All samples were dropped due to IDAT size mismatch.")
    }
  }
}

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

# Store counts before any QC filtering
n_samples_input <- nrow(targets)

# Count samples per group
n_con <- sum(tolower(targets$primary_group) == con_lower)
n_test <- sum(tolower(targets$primary_group) == test_lower)

message(paste("Found", nrow(targets), "samples for analysis."))
message(paste("  - Control Group:", group_con_in, paste0("(n=", n_con, ")")))
message(paste("  - Test Group:   ", group_test_in, paste0("(n=", n_test, ")")))

if (n_con < MIN_GROUP_SIZE_WARN || n_test < MIN_GROUP_SIZE_WARN) {
  message("WARNING: Very small sample size detected (n < 3 in one or both groups).")
  message("         Statistical power may be limited; interpret results with caution.")
}
if ((n_con + n_test) < MIN_TOTAL_SIZE_STOP) {
  stop(sprintf("ERROR: Total sample size too small (n < %d). Cannot proceed with reliable statistical analysis.", MIN_TOTAL_SIZE_STOP))
}

if (any(is.na(targets[[gsm_col]]))) {
  stop("Configuration has missing sample IDs; cannot align samples.")
}
if (any(duplicated(targets[[gsm_col]]))) {
  stop("Duplicate sample identifiers detected; please ensure unique sample identifiers.")
}
rownames(targets) <- targets[[gsm_col]]

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
if (force_idat) {
  message("  Forcing IDAT read despite differing array sizes (force=TRUE).")
}
rgSet <- read.metharray.exp(targets = targets, force = force_idat)

# A. Sample-Level QC
# Sample QC thresholds based on Illumina recommendations:
# - Median M/U signal intensity > 10.5 (log2 scale) indicates good quality
# - Poor quality samples can introduce bias in downstream analysis
detP <- detectionP(rgSet)
sample_fail <- colMeans(detP > QC_DETECTION_P_THRESHOLD, na.rm = TRUE) > QC_SAMPLE_DET_FAIL_FRAC
bad_samples_idx <- integer(0)
if (any(sample_fail)) {
    bad_samples <- colnames(detP)[sample_fail]
    message(sprintf("WARNING: Removing %d samples with > %.0f%% probes failing detection p>%.2g", length(bad_samples), QC_SAMPLE_DET_FAIL_FRAC * 100, QC_DETECTION_P_THRESHOLD))
    message(paste(bad_samples, collapse = ", "))
    keep <- !sample_fail
    rgSet <- rgSet[, keep]
    detP <- detP[, keep, drop = FALSE]
    targets <- targets[keep, , drop = FALSE]
}
if (is.finite(QC_MEDIAN_INTENSITY_THRESHOLD)) {
    message(sprintf("Performing Sample QC (Illumina recommendation: median M/U signal > %.1f log2)...", QC_MEDIAN_INTENSITY_THRESHOLD))
    qc <- getQC(preprocessRaw(rgSet))
    bad_samples_idx <- which(qc$mMed < QC_MEDIAN_INTENSITY_THRESHOLD | qc$uMed < QC_MEDIAN_INTENSITY_THRESHOLD)
    if (length(bad_samples_idx) > 0) {
        bad_samples <- rownames(qc)[bad_samples_idx]
        message(paste("WARNING: Removing", length(bad_samples), "samples due to low signal intensity (<", QC_MEDIAN_INTENSITY_THRESHOLD, "log2):"))
        message(paste(bad_samples, collapse=", "))
        targets <- targets[-bad_samples_idx, ]
        rgSet <- rgSet[, -bad_samples_idx]
        detP <- detP[, -bad_samples_idx]
        if (nrow(targets) < 2) stop("Too few samples remaining after QC.")
    }
} else {
    message("Performing Sample QC: intensity filter disabled (--qc_intensity_threshold<=0); skipping intensity-based removal.")
}
samples_failed_qc <- sum(sample_fail) + length(bad_samples_idx)

# Re-check group balance after QC filtering
n_con <- sum(targets$primary_group == clean_con)
n_test <- sum(targets$primary_group == clean_test)

if (n_con == 0 || n_test == 0) {
    stop(sprintf("ERROR: After QC filtering, at least one group has zero samples (Control: %d, Test: %d).", n_con, n_test))
}

if ((n_con + n_test) < MIN_TOTAL_SIZE_STOP) {
  stop(sprintf("ERROR: Total sample size too small after QC (n = %d, threshold = %d). Cannot proceed with reliable statistical analysis.", n_con + n_test, MIN_TOTAL_SIZE_STOP))
}

if (n_con < MIN_GROUP_SIZE_WARN || n_test < MIN_GROUP_SIZE_WARN) {
  message("WARNING: Very small sample size after QC (n < 3 in one or both groups).")
  message(sprintf("         Control: %d, Test: %d", n_con, n_test))
}

message(sprintf("Samples retained after QC - Control: %d, Test: %d (total: %d)", n_con, n_test, n_con + n_test))

# Save sample-level QC metrics and figures (publication-friendly)
tryCatch({
    qc_metrics <- data.frame(
      Sample = colnames(detP),
      Group = as.character(targets$primary_group[match(colnames(detP), rownames(targets))]),
      Detection_Fail_Fraction = colMeans(detP > QC_DETECTION_P_THRESHOLD, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (any(is.na(qc_metrics$Group))) qc_metrics$Group <- as.character(targets$primary_group)
    
    qc_for_plot <- NULL
    if (exists("qc", inherits = FALSE)) {
      qc_for_plot <- qc
    } else {
      qc_for_plot <- tryCatch(getQC(preprocessRaw(rgSet)), error = function(e) NULL)
    }
    if (!is.null(qc_for_plot) && nrow(qc_for_plot) > 0) {
      m_idx <- match(qc_metrics$Sample, rownames(qc_for_plot))
      qc_metrics$M_Median <- qc_for_plot$mMed[m_idx]
      qc_metrics$U_Median <- qc_for_plot$uMed[m_idx]
    }
    
    write.csv(qc_metrics, file.path(out_dir, "Sample_QC_Metrics.csv"), row.names = FALSE)
    
    p_det <- ggplot(qc_metrics, aes(x = reorder(Sample, Detection_Fail_Fraction),
                                   y = Detection_Fail_Fraction,
                                   fill = Group,
                                   text = paste0("Sample: ", Sample, "<br>Fail frac: ", round(Detection_Fail_Fraction, 4)))) +
      geom_col() +
      coord_flip() +
      geom_hline(yintercept = QC_SAMPLE_DET_FAIL_FRAC, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(
        title = "Sample QC: Detection P-value failure fraction",
        subtitle = sprintf("Threshold: > %.0f%% probes failing (P > %.3g)", QC_SAMPLE_DET_FAIL_FRAC * 100, QC_DETECTION_P_THRESHOLD),
        x = NULL,
        y = "Fraction of probes failing"
      ) +
      theme(legend.position = "bottom")
    save_interactive_plot(p_det, "Sample_QC_DetectionP_FailFraction.html", out_dir)
    save_static_plot(p_det, "Sample_QC_DetectionP_FailFraction.png", out_dir, width = 7, height = 4)
    
    if ("M_Median" %in% colnames(qc_metrics) && "U_Median" %in% colnames(qc_metrics) &&
        any(is.finite(qc_metrics$M_Median)) && any(is.finite(qc_metrics$U_Median))) {
      p_int <- ggplot(qc_metrics, aes(x = U_Median, y = M_Median, color = Group,
                                     text = paste0("Sample: ", Sample,
                                                   "<br>U median: ", round(U_Median, 2),
                                                   "<br>M median: ", round(M_Median, 2)))) +
        geom_point(size = 3, alpha = 0.8) +
        theme_minimal() +
        labs(
          title = "Sample QC: Median signal intensity (M/U)",
          subtitle = if (is.finite(QC_MEDIAN_INTENSITY_THRESHOLD)) {
            sprintf("Threshold: %.1f log2 (dashed lines)", QC_MEDIAN_INTENSITY_THRESHOLD)
          } else {
            "Intensity filter disabled"
          },
          x = "Unmethylated median (log2)",
          y = "Methylated median (log2)"
        )
      if (is.finite(QC_MEDIAN_INTENSITY_THRESHOLD)) {
        p_int <- p_int +
          geom_vline(xintercept = QC_MEDIAN_INTENSITY_THRESHOLD, linetype = "dashed", color = "red") +
          geom_hline(yintercept = QC_MEDIAN_INTENSITY_THRESHOLD, linetype = "dashed", color = "red")
      }
      save_interactive_plot(p_int, "Sample_QC_Intensity_Medians.html", out_dir)
      save_static_plot(p_int, "Sample_QC_Intensity_Medians.png", out_dir, width = 5, height = 4)
    }
}, error = function(e) {
    message("  QC plot export skipped: ", e$message)
})

# Cell composition estimation (before batch correction)
cell_est <- estimate_cell_counts_safe(rgSet, tissue = opt$tissue, out_dir = out_dir)
if (!is.null(cell_est)) {
  cell_df <- cell_est$counts
  cell_cols <- grep("^Cell_", colnames(cell_df), value = TRUE)
  clean_id <- function(x) sub("\\.idat.*$", "", basename(as.character(x)))
  match_idx <- match(clean_id(targets[[gsm_col]]), clean_id(cell_df$SampleID))
  if (all(is.na(match_idx)) && "Basename" %in% colnames(targets)) {
    match_idx <- match(clean_id(targets$Basename), clean_id(cell_df$SampleID))
  }
  matched <- !is.na(match_idx)
  if (any(matched)) {
    for (cc in cell_cols) {
      if (!cc %in% colnames(targets)) targets[[cc]] <- NA_real_
      targets[[cc]][matched] <- cell_df[[cc]][match_idx[matched]]
    }
    cell_covariates <- cell_cols
    write.csv(targets[, c(gsm_col, cell_cols), drop = FALSE], file.path(out_dir, "cell_counts_merged.csv"), row.names = FALSE)
    message(paste("  - Added cell composition covariates:", paste(cell_cols, collapse = ", ")))
    
    # [Over-correction Check] Evaluate association between Cell Types and Primary Group
    # High association indicates risk of over-adjustment (confounding)
    assoc_results <- data.frame(CellType = character(), P_Value = numeric(), Eta_Squared = numeric(), stringsAsFactors = FALSE)
    
    for (cc in cell_cols) {
        # ANOVA for numeric cell prop vs Factor group
        # Eta-squared ~ R-squared
        tryCatch({
            fit_aov <- aov(targets[[cc]] ~ targets$primary_group)
            summ <- summary(fit_aov)[[1]]
            ss_group <- summ["targets$primary_group", "Sum Sq"]
            ss_total <- sum(summ[, "Sum Sq"])
            eta_sq <- ss_group / ss_total
            p_val <- summ["targets$primary_group", "Pr(>F)"][1]
            
            assoc_results <- rbind(assoc_results, data.frame(CellType = cc, P_Value = p_val, Eta_Squared = eta_sq))
        }, error = function(e) {
            # Skip if error (e.g. constant values)
        })
    }
    
    if (nrow(assoc_results) > 0) {
        write.csv(assoc_results, file.path(out_dir, "Cell_Group_Association.csv"), row.names = FALSE)
        
        # Check for high confounding (Eta-squared > 0.5 is very strong association)
        high_confound <- assoc_results[assoc_results$Eta_Squared > 0.5, ]
        if (nrow(high_confound) > 0) {
            message("WARNING: Strong association detected between Cell Composition and Group. Risk of Over-correction!")
            for (i in 1:nrow(high_confound)) {
                message(sprintf("  - %s: Eta^2 = %.3f (P = %.3g)", high_confound$CellType[i], high_confound$Eta_Squared[i], high_confound$P_Value[i]))
            }
            message("  > Consider removing these covariates if they represent biological pathology rather than noise.")
        } else {
            message("  - Cell composition vs Group association checked: No strong confounding detected (all Eta^2 <= 0.5).")
        }
    }
  } else {
    message("  - Cell composition estimated but sample IDs did not match metadata; skipping merge.")
  }
}

# B. Normalization
message("Preprocessing (Noob)...")
mSet <- preprocessNoob(rgSet)
gmSet <- mapToGenome(mSet)

# C. Probe-Level QC
message("Performing Probe QC...")
nrow_raw <- nrow(gmSet)

# Fix: Filter detP to match gmSet probes first
common_probes <- intersect(rownames(detP), rownames(gmSet))
detP_sub <- detP[common_probes, ]
gmSet <- gmSet[common_probes, ]

# Identify probes passing detection p-value threshold
# We keep probes where P-value < threshold in ALL samples (strict) or fraction?
# Original logic: keep_detP <- rowSums(detP < QC_DETECTION_P_THRESHOLD) == ncol(gmSet)
# This implies ALL samples must pass.
pass_detP <- rowSums(detP_sub < QC_DETECTION_P_THRESHOLD) == ncol(gmSet)
gmSet <- gmSet[pass_detP, ]

probes_failed_detection <- length(common_probes) - sum(pass_detP)
message(sprintf("  - Removed %d probes with high detection P-value (threshold: %.3f).", probes_failed_detection, QC_DETECTION_P_THRESHOLD))

before_snps <- nrow(gmSet)
gmSet <- dropLociWithSnps(gmSet, snps=c("SBE","CpG"), maf=SNP_MAF_THRESHOLD)
n_snp_probes <- before_snps - nrow(gmSet)
message(paste("  - SNP filtering applied (MAF threshold:", SNP_MAF_THRESHOLD, "). Current probes:", nrow(gmSet)))

ann <- getAnnotation(gmSet)
keep_sex <- !(ann$chr %in% c("chrX", "chrY"))
before_sex <- nrow(gmSet)
gmSet <- gmSet[keep_sex, ]
n_sex_probes <- before_sex - nrow(gmSet)
message(paste("  - Removed Sex Chromosome probes. Final probes:", nrow(gmSet)))

qc_report <- data.frame(
  metric = c("Total_samples_input", "Samples_failed_QC", "Samples_passed_QC",
             "Total_probes_raw", "Probes_failed_detection", "Probes_with_SNPs",
             "Probes_sex_chromosomes", "Probes_final"),
  value = c(n_samples_input, samples_failed_qc, nrow(targets),
            nrow_raw, probes_failed_detection, n_snp_probes, n_sex_probes, nrow(gmSet))
)
write.csv(qc_report, file.path(out_dir, "QC_Summary.csv"), row.names = FALSE)

beta_minfi <- getBeta(gmSet)
if (ncol(beta_minfi) != nrow(targets)) {
  stop("Failed to align metadata to Minfi beta matrix: sample counts differ.")
}
# Preserve explicit sample IDs (no underscore trimming)
colnames(beta_minfi) <- targets[[gsm_col]]

anno_data <- getAnnotation(gmSet)
anno <- anno_data[, c("chr", "pos", "Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
anno_df <- as.data.frame(anno)
anno_df$CpG <- rownames(anno_df)

colnames(anno_df)[colnames(anno_df) == "UCSC_RefGene_Name"] <- "Gene"
colnames(anno_df)[colnames(anno_df) == "UCSC_RefGene_Group"] <- "Region"
colnames(anno_df)[colnames(anno_df) == "Relation_to_Island"] <- "Island_Context"


# --- 3. Sesame Analysis ---
message("--- Starting Sesame Analysis ---")
beta_sesame <- NULL
sesame_before_qc <- 0
if (opt$`skip-sesame`) {
  message("Sesame analysis skipped (--skip-sesame).")
} else {
  tryCatch({
    ssets <- lapply(targets$Basename, function(x) readIDATpair(x))
    sesame_preprocess <- function(x) {
      sdf <- noob(x)
      tryCatch(
        dyeBiasCorrTypeINorm(sdf),
        error = function(e) {
          message("  Sesame dyeBiasCorrTypeINorm failed; falling back to dyeBiasCorr: ", e$message)
          dyeBiasCorr(sdf)
        }
      )
    }
    betas_sesame_list <- lapply(ssets, function(x) getBetas(sesame_preprocess(x)))
    common_probes_sesame <- Reduce(intersect, lapply(betas_sesame_list, names))
    beta_sesame <- do.call(cbind, lapply(betas_sesame_list, function(x) x[common_probes_sesame]))
    colnames(beta_sesame) <- targets[[gsm_col]]
    beta_sesame <- na.omit(beta_sesame)
    sesame_before_qc <- nrow(beta_sesame)
    # Harmonize Sesame with Minfi QC/annotation footprint
    qc_probes <- rownames(beta_minfi)
    shared_qc <- intersect(rownames(beta_sesame), qc_probes)
    if (length(shared_qc) > 0) {
      beta_sesame <- beta_sesame[shared_qc, , drop = FALSE]
    }
    message(paste("Sesame processed", sesame_before_qc, "probes (after QC alignment:", nrow(beta_sesame), ")."))
  }, error = function(e) {
    message("Sesame analysis failed; continuing with Minfi only: ", e$message)
    beta_sesame <<- NULL
  })
}

# --- 4. Intersection ---
if (is.null(beta_sesame)) {
  common_cpgs <- character(0)
  message("Intersection: skipped (Sesame not available).")
} else {
  common_cpgs <- intersect(rownames(beta_minfi), rownames(beta_sesame))
  message(paste("Intersection:", length(common_cpgs), "CpGs common to Minfi and Sesame."))
}

# --- Pipeline Function ---

run_pipeline <- function(betas, prefix, annotation_df) {
  message(paste("Running pipeline for:", prefix))
  # Work on a local copy of targets to avoid altering the shared metadata across pipelines
  targets <- targets
  
  best_method <- "none"
  batch_evaluated <- FALSE
  common_ids <- intersect(rownames(betas), rownames(annotation_df))
  betas <- betas[common_ids, , drop = FALSE]
  curr_anno <- annotation_df[common_ids, , drop = FALSE]
  betas_full <- betas
  clock_computed <- FALSE
  
  range_filt <- filter_low_range(betas, min_range = BETA_RANGE_MIN)
  betas <- range_filt$mat
  curr_anno <- curr_anno[rownames(betas), , drop = FALSE]
  message(sprintf("  - Removed %d probes with beta range < %.2f", range_filt$removed, BETA_RANGE_MIN))
  
  # Drop probes with zero variance across samples to avoid PCA failures
  var_vals <- apply(betas, 1, var, na.rm = TRUE)
  keep_var <- var_vals > 0
  if (sum(keep_var) < 2) {
    stop(paste(prefix, "has fewer than 2 probes with non-zero variance; cannot proceed."))
  }
  betas <- betas[keep_var, , drop = FALSE]
  
  message(sprintf("  Computing PCA on %d probes across %d samples...", nrow(betas), ncol(betas)))
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
  save_static_plot(p_dist, paste0(prefix, "_Sample_Clustering_Distance.png"), out_dir, width = 7, height = 6)

  pca_res <- prcomp(t(betas), scale. = TRUE)
  pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = targets$primary_group)
  p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) + geom_point(size=3) + theme_minimal() + ggtitle(paste(prefix, "PCA (Raw)"))
  save_interactive_plot(p_pca, paste0(prefix, "_PCA_Before.html"), out_dir)
  save_static_plot(p_pca, paste0(prefix, "_PCA_Before.png"), out_dir, width = 5, height = 4)
  
  # Auto-detect covariates via PC association (top PCs) unless disabled
  covariates <- character(0)
  forced_covariates <- character(0)
  covar_log_path <- file.path(out_dir, paste0(prefix, "_AutoCovariates.csv"))
  drop_log <- data.frame(Variable=character(0), Reason=character(0))
  covar_log_df <- data.frame(Variable = character(0), MinP = numeric(0))

  if (include_clock_covariates) {
    clock_cov <- compute_epigenetic_clocks(betas_full, targets, id_col = gsm_col, prefix = prefix,
                                           out_dir = out_dir, tissue = opt$tissue, return_covariates = TRUE)
    clock_computed <- TRUE
    if (!is.null(clock_cov) && ncol(clock_cov) > 1) {
      match_idx <- match(targets[[gsm_col]], clock_cov$SampleID)
      added_cols <- setdiff(colnames(clock_cov), "SampleID")
      for (cc in added_cols) {
        targets[[cc]] <- clock_cov[[cc]][match_idx]
      }
      drop_clock <- character(0)
      for (cc in added_cols) {
        vals <- targets[[cc]]
        if (any(!is.finite(vals))) {
          drop_clock <- c(drop_clock, cc)
          drop_log <- rbind(drop_log, data.frame(Variable = cc, Reason = "clock_missing_values"))
          next
        }
        if (length(unique(vals)) < 2) {
          drop_clock <- c(drop_clock, cc)
          drop_log <- rbind(drop_log, data.frame(Variable = cc, Reason = "clock_constant"))
        }
      }
      if (length(drop_clock) > 0) {
        targets[, drop_clock] <- NULL
      }
      out_cfg <- file.path(out_dir, paste0(prefix, "_configure_with_clocks.tsv"))
      write.table(targets, out_cfg, sep = "\t", row.names = FALSE, quote = FALSE)
      message(paste("  Clock covariates merged; saved to", out_cfg))
    } else {
      message("  Clock covariates requested but no usable clock values were computed.")
    }
  }
  
  if (!disable_auto_cov) {
    covar_info <- auto_detect_covariates(targets, gsm_col, pca_res$x, alpha = AUTO_COVARIATE_ALPHA, max_pcs = MAX_PCS_FOR_COVARIATE_DETECTION)
    covar_log_df <- covar_info$log
    if (nrow(covar_log_df) > 0) write.csv(covar_log_df, covar_log_path, row.names = FALSE)
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
    forced_covariates <- unique(present_force)
    if (length(missing_force) > 0) {
      drop_log <- rbind(drop_log, data.frame(Variable = missing_force, Reason = "not_in_config"))
    }
  }
  if (length(cell_covariates) > 0) {
    present_cells <- intersect(cell_covariates, colnames(targets))
    varying_cells <- present_cells[vapply(present_cells, function(x) length(unique(targets[[x]])) > 1, logical(1))]
    if (length(varying_cells) > 0) {
      covariates <- unique(c(covariates, varying_cells))
    }
  }
  
  if (length(covariates) > 0) {
    message(paste("  Covariate candidates after filtering:", paste(covariates, collapse = ", ")))
  } else {
    message("  No covariate candidates retained before design drop.")
  }
  n_groups <- length(levels(targets$primary_group))
  max_covars <- nrow(targets) - n_groups - 1
  if (max_covars < 0) max_covars <- 0
  if (length(covariates) > max_covars) {
    ranked_covars <- rank_covariates(targets, covariates, covar_log_df)
    keep_covars <- head(ranked_covars, max_covars)
    drop_covars <- setdiff(covariates, keep_covars)
    if (length(drop_covars) > 0) {
      drop_log <- rbind(drop_log, data.frame(Variable = drop_covars, Reason = "covariate_limit_small_n"))
      covariates <- keep_covars
      message(paste("  Capping covariates to preserve residual df (max", max_covars, "): kept", paste(keep_covars, collapse = ", ")))
    }
  }

  guard_cov <- cap_covariates_for_overcorrection(targets, covariates, forced_covariates, covar_log_df)
  if (length(guard_cov$dropped) > 0) {
    drop_log <- rbind(drop_log, data.frame(Variable = guard_cov$dropped, Reason = "overcorrection_guard_covariates"))
    covariates <- guard_cov$keep
    message(paste("  Overcorrection guard: capped covariates to", guard_cov$cap, "(", paste(covariates, collapse = ", "), ")"))
    if (length(guard_cov$forced_dropped) > 0) {
      message(paste("  Overcorrection guard dropped forced covariates:", paste(guard_cov$forced_dropped, collapse = ", ")))
    }
  }

  # Drop samples with NA in any design variable to prevent model.matrix failures
  design_vars <- unique(c("primary_group", covariates))
  design_vars <- intersect(design_vars, colnames(targets))
  if (length(design_vars) > 0) {
    cc_idx <- complete.cases(targets[, design_vars, drop = FALSE])
    if (!all(cc_idx)) {
      dropped_ids <- rownames(targets)[!cc_idx]
      message(sprintf("  Removing %d samples with NA in design variables (%s): %s",
                      sum(!cc_idx), paste(design_vars, collapse = ", "),
                      paste(dropped_ids, collapse = ", ")))
      targets <- targets[cc_idx, , drop = FALSE]
      betas <- betas[, cc_idx, drop = FALSE]
    }
  }
  n_con_local <- sum(targets$primary_group == clean_con)
  n_test_local <- sum(targets$primary_group == clean_test)
  if (n_con_local == 0 || n_test_local == 0) {
    stop("After removing samples with NA covariates, one of the groups has zero samples.")
  }
  if (n_con_local < MIN_GROUP_SIZE_WARN || n_test_local < MIN_GROUP_SIZE_WARN) {
    message(sprintf("  WARNING: Very small sample size after NA filtering (Control: %d, Test: %d)", n_con_local, n_test_local))
  }

  guard_cov_post <- cap_covariates_for_overcorrection(targets, covariates, forced_covariates, covar_log_df)
  if (length(guard_cov_post$dropped) > 0) {
    drop_log <- rbind(drop_log, data.frame(Variable = guard_cov_post$dropped, Reason = "overcorrection_guard_covariates_post_na"))
    covariates <- guard_cov_post$keep
    message(paste("  Overcorrection guard (post-NA): capped covariates to", guard_cov_post$cap, "(", paste(covariates, collapse = ", "), ")"))
    if (length(guard_cov_post$forced_dropped) > 0) {
      message(paste("  Overcorrection guard (post-NA) dropped forced covariates:", paste(guard_cov_post$forced_dropped, collapse = ", ")))
    }
  }
  
  # PVCA before model fitting to quantify variance explained by group/batch factors
  pvca_factors <- unique(c("primary_group", covariates))
  run_pvca_assessment(betas, targets, pvca_factors, prefix, out_dir, threshold = 0.6, max_probes = 5000, sample_ids = colnames(betas))
  if (!clock_computed) {
    compute_epigenetic_clocks(betas, targets, id_col = gsm_col, prefix = prefix, out_dir = out_dir, tissue = opt$tissue)
  }

  # Compare batch correction strategies (none/ComBat/removeBatchEffect/SVA)
  batch_col <- select_batch_factor(targets, preferred = c("Sentrix_ID", "Sentrix_Position"))
  if (!is.null(batch_col) && is_batch_confounded(targets, batch_col, "primary_group")) {
    message(paste("  WARNING: Batch factor", batch_col, "is confounded with primary_group; skipping batch correction."))
    batch_col <- NULL
  }
  if (!is.null(batch_col) && nrow(targets) >= 4) {
    message(paste("  Evaluating batch correction methods using batch factor:", batch_col))
    M_mat <- logit_offset(betas)
    methods <- c("none", "combat", "limma", "sva")
    batch_results <- lapply(methods, function(m) {
      tryCatch(eval_batch_method(M_mat, targets, group_col = "primary_group", batch_col = batch_col, covariates = covariates, method = m),
               error = function(e) { message("    - ", m, " failed: ", e$message); NULL })
    })
    batch_results <- batch_results[!vapply(batch_results, is.null, logical(1))]
    if (length(batch_results) > 0) {
      batch_evaluated <- TRUE
      batch_df <- do.call(rbind, lapply(batch_results, function(x) {
        data.frame(method = x$method,
                   score = x$score,
                   batch_var = x$batch_var,
                   group_var = x$group_var,
                   n_pc_batch_sig = x$n_pc_batch_sig,
                   prop_batch_sig = x$prop_batch_sig)
      }))
      out_batch_path <- file.path(out_dir, paste0(prefix, "_BatchMethodComparison.csv"))
      write.csv(batch_df, out_batch_path, row.names = FALSE)
      best_idx <- which.min(batch_df$score)
      best_method <- batch_df$method[best_idx]
      message(paste("  Best batch method by score:", best_method, "(saved to", out_batch_path, ")"))
      
      # Visualization: Batch Method Comparison
      tryCatch({
        plot_df <- batch_df[, c("method", "batch_var", "group_var", "score")]
        # Normalize/Scale for visualization if needed, but raw values are informative enough for relative comparison
        # Reshape to long (avoid extra package dependencies)
        long_df <- rbind(
          data.frame(method = plot_df$method, variable = "batch_var", value = plot_df$batch_var),
          data.frame(method = plot_df$method, variable = "group_var", value = plot_df$group_var),
          data.frame(method = plot_df$method, variable = "score", value = plot_df$score)
        )
        
        # Add a flag for the best method
        long_df$is_best <- long_df$method == best_method
        
        p_comp <- ggplot(long_df, aes(x = method, y = value, fill = variable, alpha = is_best)) +
          geom_bar(stat = "identity", position = "dodge") +
          scale_alpha_manual(values = c("FALSE"=0.6, "TRUE"=1.0), guide="none") +
          scale_fill_manual(values = c("batch_var"="#e74c3c", "group_var"="#3498db", "score"="#95a5a6"),
                            labels = c("Batch Residual (Lower is better)", "Group Signal (Higher is better)", "Combined Score (Lower is better)")) +
          labs(title = paste(prefix, "Batch Correction Method Evaluation"),
               subtitle = paste("Selected Best Method:", best_method),
               y = "Metric Value", x = "Method", fill = "Metric") +
          theme_minimal()
          
        save_interactive_plot(p_comp, paste0(prefix, "_Batch_Method_Comparison.html"), out_dir)
        save_static_plot(p_comp, paste0(prefix, "_Batch_Method_Comparison.png"), out_dir, width = 7, height = 4)
      }, error = function(e) {
        message("  - Failed to generate batch comparison plot: ", e$message)
      })
    } else {
      message("  Batch method comparison skipped: no results.")
    }
  } else if (!is.null(batch_col)) {
    message("  Skipping batch method comparison (too few samples for stable evaluation).")
  } else {
    message("  No eligible batch factor found for batch method comparison.")
  }
  
  # Build Design Matrix
  formula_str <- "~ 0 + primary_group"
  if (length(covariates) > 0) formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = targets)
  colnames(design) <- gsub("primary_group", "", colnames(design))
  colnames(design) <- make.unique(make.names(colnames(design)))
  group_cols <- make.names(levels(targets$primary_group))
  
  sv_cols <- character(0)
  use_sva_in_model <- !disable_sva && (!batch_evaluated || best_method == "sva")
  
  # --- Surrogate Variable Analysis (SVA) ---
  if (!disable_sva && !use_sva_in_model) {
    message(paste("  SVA skipped for modeling to avoid double correction (selected method:", best_method, ")."))
  } else if (!disable_sva) {
    message("  Running SVA to detect hidden batch effects...")
    svobj <- NULL  # reset per pipeline to avoid stale state across runs
  
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
      # Speed: estimate SVs on top-variable probes only
      sva_mat <- logit_offset(betas)
      if (nrow(sva_mat) > BATCH_EVAL_TOP_VAR) {
        v0 <- apply(sva_mat, 1, var, na.rm = TRUE)
        top_idx <- head(order(v0, decreasing = TRUE), BATCH_EVAL_TOP_VAR)
        sva_mat <- sva_mat[top_idx, , drop = FALSE]
      }
      svobj <- sva(sva_mat, mod, mod0)
      return(svobj)
    }
    
    tryCatch({
        svobj <- run_sva(mod_sva, mod0_sva, "full")
        if (!is.null(svobj)) {
          n_sv <- svobj$n.sv
          if (n_sv == 0) message("    - No significant surrogate variables found.")
        }
    }, error = function(e) {
        message("    - SVA failed on full model: ", e$message)
        message("    - Attempting fallback: group-only SVA model...")
        mod_sva2 <- model.matrix(~ primary_group, data = targets)
        mod0_sva2 <- model.matrix(~ 1, data = targets)
        tryCatch({
          svobj <<- run_sva(mod_sva2, mod0_sva2, "group-only")
        }, error = function(e2) {
          message("    - SVA fallback also failed: ", e2$message)
          message("    - Proceeding without SVA. This may affect batch correction quality.")
          svobj <<- NULL
        })
    })
    
    if (!is.null(svobj)) {
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
            if (nrow(sv_mat) != nrow(design)) {
                message(paste("    - Skipping SVA integration: sv rows (", nrow(sv_mat),
                              ") != design rows (", nrow(design), ")", sep = ""))
            } else {
                design <- cbind(design, sv_mat)
                targets <- cbind(targets, as.data.frame(sv_mat))
                sv_cols <- colnames(sv_mat)
            }
        } else {
            message("    - No significant surrogate variables found.")
        }
    }
  } else {
    message("  SVA disabled by flag.")
  }
  if (length(sv_cols) > 0) {
    sv_filt <- filter_sv_by_group(targets, sv_cols, "primary_group")
    if (nrow(sv_filt$dropped) > 0) {
      drop_sv <- sv_filt$dropped$Variable
      message(paste("  Dropping SVs strongly associated with primary_group:", paste(drop_sv, collapse = ", ")))
      drop_log <- rbind(drop_log, sv_filt$dropped)
      drop_sv <- intersect(drop_sv, colnames(design))
      if (length(drop_sv) > 0) {
        design <- design[, setdiff(colnames(design), drop_sv), drop = FALSE]
      }
      targets <- targets[, setdiff(colnames(targets), drop_sv), drop = FALSE]
      sv_cols <- sv_filt$keep
    }
  }
  max_terms_total <- max(1, floor(nrow(targets) * OVERCORRECTION_GUARD_RATIO))
  if (length(sv_cols) > 0) {
    max_sv_allowed <- max(0, max_terms_total - length(covariates))
    if (max_sv_allowed < length(sv_cols)) {
      drop_sv <- sv_cols[(max_sv_allowed + 1):length(sv_cols)]
      if (length(drop_sv) > 0) {
        drop_log <- rbind(drop_log, data.frame(Variable = drop_sv, Reason = "overcorrection_guard_sv"))
        message(paste("  Overcorrection guard: capping SVs to", max_sv_allowed, "(dropping", paste(drop_sv, collapse = ", "), ")"))
        design <- design[, setdiff(colnames(design), drop_sv), drop = FALSE]
        targets <- targets[, setdiff(colnames(targets), drop_sv), drop = FALSE]
        sv_cols <- setdiff(sv_cols, drop_sv)
      }
    }
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
  covar_total <- length(used_covariates) + length(sv_cols)
  if (covar_total > floor(nrow(targets) / 3)) {
    message("WARNING: Number of covariates is large relative to sample size.")
    message(sprintf("  Covariates: %d, Samples: %d (ratio: %.2f)",
                    covar_total, nrow(targets),
                    covar_total / nrow(targets)))
    message("  This may reduce statistical power. Consider manual covariate selection.")
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
    clean_betas[!is.finite(clean_betas)] <- NA
    clean_betas <- pmin(pmax(clean_betas, LOGIT_OFFSET), 1 - LOGIT_OFFSET)
    
    pca_clean <- safe_prcomp(t(clean_betas), label = paste(prefix, "PCA corrected"), scale. = TRUE)
    if (!is.null(pca_clean)) {
      pca_clean_df <- data.frame(PC1 = pca_clean$x[,1], PC2 = pca_clean$x[,2], Group = targets$primary_group)
      cov_label <- if (length(sv_cols) > 0) "covariates + SVs" else "covariates"
      p_pca_clean <- ggplot(pca_clean_df, aes(x=PC1, y=PC2, color=Group)) + 
        geom_point(size=3) + theme_minimal() + ggtitle(paste(prefix, "PCA (Corrected:", cov_label, ")"))
      save_interactive_plot(p_pca_clean, paste0(prefix, "_PCA_After_Correction.html"), out_dir)
      save_static_plot(p_pca_clean, paste0(prefix, "_PCA_After_Correction.png"), out_dir, width = 5, height = 4)
    } else {
      message("  PCA after correction skipped due to insufficient variance.")
    }
    
    # PVCA after correction to quantify residual batch/group contribution
    run_pvca_assessment(clean_betas, targets, pvca_factors, paste0(prefix, "_AfterCorrection"), out_dir, threshold = 0.6, max_probes = 5000, sample_ids = colnames(clean_betas))
  }

  # Apply the selected batch correction method to the modeling matrix
  betas_for_model <- betas
  applied_batch_method <- "none"
  if (!is.null(batch_col) && (batch_col %in% colnames(targets)) && length(unique(targets[[batch_col]])) > 1) {
    bc_res <- apply_batch_correction(betas, best_method, batch_col, all_covariates, targets, design)
    betas_for_model <- bc_res$betas
    applied_batch_method <- bc_res$method
    if (applied_batch_method != "none") {
      message(paste("  Batch correction applied for modeling using:", applied_batch_method))
    } else {
      message("  Batch correction not applied to modeling matrix (method=none or unsupported).")
    }
  }
  betas <- betas_for_model

  # Effect sizes on the beta scale (paper-friendly interpretability)
  con_mask <- targets$primary_group == clean_con
  test_mask <- targets$primary_group == clean_test
  mean_beta_con <- rowMeans(betas[, con_mask, drop = FALSE], na.rm = TRUE)
  mean_beta_test <- rowMeans(betas[, test_mask, drop = FALSE], na.rm = TRUE)
  delta_beta <- mean_beta_test - mean_beta_con
  
  groups <- levels(targets$primary_group)
  if (length(groups) < 2) {
      message("Skipping stats: Not enough groups.")
      return(list(res = NULL, n_con = n_con_local, n_test = n_test_local, n_samples = nrow(targets)))
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

  res$Mean_Beta_Con <- mean_beta_con[res$CpG]
  res$Mean_Beta_Test <- mean_beta_test[res$CpG]
  res$Delta_Beta <- delta_beta[res$CpG]
  
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
  save_static_plot(p_qq, paste0(prefix, "_QQPlot.png"), out_dir, width = 5, height = 4)

  plot_res <- res[order(res$P.Value), ]
  if (nrow(plot_res) > max_points) {
      message(paste("Subsampling top", max_points, "probes for interactive plots..."))
      plot_res <- plot_res[1:max_points, ]
  }
  
  plot_res$diffexpressed <- "NO"
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC > lfc_thresh] <- "UP"
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC < -lfc_thresh] <- "DOWN"
  
  subtitle_str <- paste0("Control: ", group_con_in, " (n=", n_con_local, ") vs Test: ", group_test_in, " (n=", n_test_local, ")")
  
  p_vol <- ggplot(plot_res, aes(x=logFC, y=-log10(P.Value), color=diffexpressed, 
                  text=paste("CpG:", CpG,
                             "<br>Gene:", Gene,
                             "<br>DeltaBeta:", round(Delta_Beta, 4),
                             "<br>Region:", Region,
                             "<br>Island:", Island_Context))) +
    geom_point(alpha=0.5) + theme_minimal() +
    scale_color_manual(values=c("DOWN"="blue", "NO"="grey", "UP"="red")) +
    labs(title = paste(prefix, "Volcano Plot"), 
         subtitle = subtitle_str,
         color = "Status")
  save_interactive_plot(p_vol, paste0(prefix, "_Volcano.html"), out_dir)
  save_static_plot(p_vol, paste0(prefix, "_Volcano.png"), out_dir, width = 6, height = 5)
  
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
                        text=paste("CpG:", CpG,
                                   "<br>Gene:", Gene,
                                   "<br>DeltaBeta:", round(Delta_Beta, 4),
                                   "<br>Region:", Region,
                                   "<br>Island:", Island_Context))) +
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
      save_static_plot(p_man, paste0(prefix, "_Manhattan.png"), out_dir, width = 8, height = 4)
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
          save_static_plot(p_heat_top, paste0(prefix, "_Top100_Heatmap.png"), out_dir, width = 8, height = 10)
      }
  }

  col_order <- c("CpG", "Gene", setdiff(colnames(res), c("CpG", "Gene")))
  res <- res[, col_order]
  res <- res[order(res$P.Value, na.last = TRUE), ]

  write.csv(res, file.path(out_dir, paste0(prefix, "_DMPs_full.csv")), row.names = FALSE)
  top_for_table <- head(res, min(10000, nrow(res)))
  save_datatable(top_for_table, paste0(prefix, "_Top_DMPs.html"), out_dir)
  
  # --- D. DMR Analysis (dmrff) ---
  if (nrow(targets) < 4) {
    message("  Skipping DMR analysis (dmrff) due to very small sample size (<4).")
    dmr_res <- data.frame()
  } else {
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
  
  if (length(dmr_est) < 100) {
    message("    - Too few CpGs for reliable DMR analysis (< 100). Skipping DMR detection.")
    dmr_res <- data.frame()
  } else {
    annotate_dmr <- function(df, anno_tbl) {
      if (nrow(df) == 0) return(df)
      anno_tbl$Gene <- ifelse(is.na(anno_tbl$Gene), "", anno_tbl$Gene)
      anno_tbl$Region <- ifelse(is.na(anno_tbl$Region), "", anno_tbl$Region)
      gene_region <- function(chr, start, end, field) {
        idx <- which(anno_tbl$chr == chr & anno_tbl$pos >= start & anno_tbl$pos <= end)
        if (length(idx) == 0) return("")
        vals <- unique(unlist(strsplit(anno_tbl[[field]][idx], ";")))
        vals <- vals[nzchar(vals)]
        if (length(vals) == 0) return("")
        paste(head(vals, 5), collapse = ";")
      }
      df$Genes <- mapply(gene_region, df$chr, df$start, df$end, MoreArgs = list(field = "Gene"))
      df$Regions <- mapply(gene_region, df$chr, df$start, df$end, MoreArgs = list(field = "Region"))
      df
    }

    tryCatch({
      dmr_res <- dmrff(
        estimate = dmr_est,
        se = dmr_se,
        p.value = dmr_p,
        methylation = betas_dmr,
        chr = dmr_anno$chr,
        pos = dmr_anno$pos,
        maxgap = DMR_MAXGAP,
        p.cutoff = DMR_P_CUTOFF
      )
      
      if (nrow(dmr_res) > 0) {
          dmr_res <- annotate_dmr(dmr_res, curr_anno)
          dmr_res <- dmr_res[order(dmr_res$p.adjust, dmr_res$p.value, na.last = TRUE), ]
          message(paste("    - Found", nrow(dmr_res), "DMRs."))
          
          # Save Table
          dmr_filename <- paste0(prefix, "_DMRs_Table.html")
          save_datatable(head(dmr_res, 3000), dmr_filename, out_dir)
          write.csv(dmr_res, file.path(out_dir, paste0(prefix, "_DMRs.csv")))
          
          # 1. DMR Volcano Plot
          # dmrff returns 'estimate' (avg LogFC) and 'p.adjust'
          dmr_res$diffexpressed <- "NO"
          dmr_res$diffexpressed[dmr_res$p.adjust < DMR_P_CUTOFF & dmr_res$estimate > 0] <- "UP"
          dmr_res$diffexpressed[dmr_res$p.adjust < DMR_P_CUTOFF & dmr_res$estimate < 0] <- "DOWN"
          
          p_vol_dmr <- ggplot(dmr_res, aes(x=estimate, y=-log10(p.adjust), color=diffexpressed,
                          text=paste("Region:", paste0(chr, ":", start, "-", end), "<br>nCpG:", n, "<br>P-adj:", p.adjust))) +
              geom_point(alpha=0.6, size=2) + theme_minimal() +
              scale_color_manual(values=c("DOWN"="blue", "NO"="grey", "UP"="red")) +
              labs(title = paste(prefix, "DMR Volcano Plot"), 
                   subtitle = "dmrff analysis",
                   x = "Mean Methylation Difference (Estimate)", y = "-log10 Adjusted P-value",
                   color = "Status")
          save_interactive_plot(p_vol_dmr, paste0(prefix, "_DMR_Volcano.html"), out_dir)
          save_static_plot(p_vol_dmr, paste0(prefix, "_DMR_Volcano.png"), out_dir, width = 6, height = 5)
          
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
              save_static_plot(p_man_dmr, paste0(prefix, "_DMR_Manhattan.png"), out_dir, width = 8, height = 4)
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
               save_static_plot(p_heat_dmr, paste0(prefix, "_Top_DMRs_Heatmap.png"), out_dir, width = 8, height = 8)
          }
          
      } else {
          message("    - No DMRs found.")
      }
    }, error = function(e) {
        message("    - DMR analysis failed: ", e$message)
    })
  }
  }

  evaluate_batch <- function(data_betas, meta, label, file_suffix, allowed_vars) {
     allowed <- unique(c("primary_group", allowed_vars))
     allowed <- intersect(allowed, colnames(meta))
     keep_cols <- vapply(allowed, function(nm) {
         if (nm == "primary_group") return(TRUE)
         if (grepl("^SV", nm)) return(length(unique(meta[[nm]])) > 1)
         uniq <- length(unique(meta[[nm]]))
         return(uniq > 1 & uniq < nrow(meta))
     }, logical(1))
     test_meta <- meta[, allowed[keep_cols], drop=FALSE]
     
     pca <- safe_prcomp(t(data_betas), label = paste(label, "PCA"), scale. = TRUE)
     if (is.null(pca)) {
         empty_pvals <- matrix(NA_real_, nrow = ncol(test_meta), ncol = 1)
         rownames(empty_pvals) <- colnames(test_meta)
         colnames(empty_pvals) <- "PC1"
         csv_name <- paste0(prefix, "_Batch_Evaluation_", file_suffix, "_Table.csv")
         write.csv(empty_pvals, file.path(out_dir, csv_name))
         p_heat <- ggplot() + theme_void() + ggtitle(paste(label, "- PCA skipped (insufficient variance)"))
         return(list(plot = p_heat, pvals = empty_pvals))
     }
     
     eig <- (pca$sdev)^2
     var_expl <- eig/sum(eig)
     
     n_pcs <- min(10, ncol(pca$x))
     pcs <- pca$x[, 1:n_pcs]
     
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
  save_static_plot(p_eval_before$plot, paste0(prefix, "_Batch_Evaluation_Before.png"), out_dir, width = 7, height = 5)
  pvals_before <- p_eval_before$pvals
  
  p_eval_after <- NULL
  pvals_after <- NULL
  if (length(all_covariates) > 0) {
      cov_formula <- as.formula(paste("~ 0 +", paste(all_covariates, collapse = " + ")))
      cov_mat <- model.matrix(cov_formula, data = targets)
      
      design_keep <- model.matrix(~ primary_group, data=targets)
      clean_betas <- removeBatchEffect(betas, covariates=cov_mat, design=design_keep)
      clean_betas[!is.finite(clean_betas)] <- NA
      clean_betas <- pmin(pmax(clean_betas, LOGIT_OFFSET), 1 - LOGIT_OFFSET)
      
      p_eval_after <- evaluate_batch(clean_betas, targets, paste(prefix, "After"), "After", allowed_vars)
      pvals_after <- p_eval_after$pvals
      save_interactive_plot(p_eval_after$plot, paste0(prefix, "_Batch_Evaluation_After.html"), out_dir)
      save_static_plot(p_eval_after$plot, paste0(prefix, "_Batch_Evaluation_After.png"), out_dir, width = 7, height = 5)
  }
  
  # Permutation test on top-variable CpGs (optional)
  perm_summary <- NULL
  if (perm_n > 0) {
      message(paste("  Running", perm_n, "permutations for null DMP counts (top", vp_top, "CpGs)..."))
      top_perm <- head(order(apply(betas, 1, var), decreasing = TRUE), min(vp_top, nrow(betas)))
      m_perm <- logit_offset(betas[top_perm, , drop=FALSE])
      perm_results <- data.frame(run = integer(0), sig_count = integer(0), min_p = numeric(0))
      perm_group_cols <- make.names(levels(targets$primary_group))
      perm_start <- Sys.time()
      
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
          
          # Progress every 10 permutations (and at the end)
          if (i %% 10 == 0 || i == perm_n) {
              elapsed <- as.numeric(difftime(Sys.time(), perm_start, units = "mins"))
              eta <- (elapsed / i) * (perm_n - i)
              message(sprintf("    [perm %d/%d] elapsed: %.1f min, ETA: %.1f min", i, perm_n, elapsed, eta))
          }
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
    metric = c("pipeline", "lambda", "batch_method_applied", "n_samples", "n_cpgs", "n_covariates_used", "covariates_used", "n_sv_used", "sv_used",
               "batch_sig_p_lt_0.05_before", "batch_min_p_before",
               "batch_sig_p_lt_0.05_after", "batch_min_p_after",
               "group_min_p_before", "group_min_p_after",
               "dropped_covariates",
               "perm_mean_sig", "perm_max_sig", "vp_primary_group_mean"),
    value = c(prefix, lambda_val, applied_batch_method, nrow(targets), nrow(betas), length(used_covariates),
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
  
  # Positive Control Check
  if (!is.null(opt$positive_controls)) {
      pc_genes <- trimws(strsplit(opt$positive_controls, ",")[[1]])
      # res$Gene contains "GeneA;GeneB" format. We search loosely.
      pc_hits <- data.frame()
      for (pg in pc_genes) {
          # Regex search for exact gene symbol (surrounded by ; or start/end)
          # Simplified: just grep the symbol
          idx <- grep(paste0("(^|;)", pg, "(;|$)"), res$Gene)
          if (length(idx) > 0) {
              hits <- res[idx, c("CpG", "Gene", "logFC", "P.Value", "adj.P.Val")]
              hits$Target <- pg
              # Take top hit by P-value
              top_hit <- hits[which.min(hits$P.Value), ]
              pc_hits <- rbind(pc_hits, top_hit)
          } else {
              # Not found row
              pc_hits <- rbind(pc_hits, data.frame(CpG="Not Found", Gene=NA, logFC=NA, P.Value=NA, adj.P.Val=NA, Target=pg))
          }
      }
      
      if (nrow(pc_hits) > 0) {
          pc_out <- file.path(out_dir, paste0(prefix, "_Positive_Controls.csv"))
          write.csv(pc_hits, pc_out, row.names=FALSE)
          message(paste("  - Positive Control Check:", paste(nrow(pc_hits), "genes checked.")))
          # Print brief summary of found genes
          found_pcs <- pc_hits[!is.na(pc_hits$P.Value), ]
          if (nrow(found_pcs) > 0) {
              message("    > Found:")
              print(found_pcs[, c("Target", "logFC", "adj.P.Val")])
          }
      }
  }
  
  return(list(res = res, n_con = n_con_local, n_test = n_test_local, n_samples = nrow(targets)))
}

minfi_out <- run_pipeline(beta_minfi, "Minfi", anno_df)
sesame_out <- NULL
if (!is.null(beta_sesame)) {
  sesame_out <- run_pipeline(beta_sesame, "Sesame", anno_df)
} else {
  message("Sesame pipeline skipped (no Sesame betas available).")
}
res_minfi <- if (!is.null(minfi_out)) minfi_out$res else NULL
res_sesame <- if (!is.null(sesame_out)) sesame_out$res else NULL

# --- Consensus (Intersection) between Minfi and Sesame ---
# For publication-ready "double validation", we define consensus DMPs as probes that are significant
# in BOTH pipelines with the same direction (and passing the same thresholds).
consensus_counts <- list(up = 0, down = 0)
if (is.null(res_minfi) || is.null(res_sesame)) {
  message("Warning: Consensus (Intersection) skipped because one or more pipelines produced no results.")
} else {
  tryCatch({
    keep_cols <- intersect(
      c("CpG", "Gene", "chr", "pos", "Region", "Island_Context", "logFC", "P.Value", "adj.P.Val"),
      colnames(res_minfi)
    )
    a <- res_minfi[, keep_cols, drop = FALSE]
    b <- res_sesame[, intersect(c("CpG", "logFC", "P.Value", "adj.P.Val"), colnames(res_sesame)), drop = FALSE]
    concord <- merge(a, b, by = "CpG", suffixes = c(".Minfi", ".Sesame"))
    concord <- concord[is.finite(concord$logFC.Minfi) & is.finite(concord$logFC.Sesame), , drop = FALSE]
    
    is_up <- concord$adj.P.Val.Minfi < pval_thresh &
      concord$adj.P.Val.Sesame < pval_thresh &
      concord$logFC.Minfi > lfc_thresh &
      concord$logFC.Sesame > lfc_thresh
    
    is_down <- concord$adj.P.Val.Minfi < pval_thresh &
      concord$adj.P.Val.Sesame < pval_thresh &
      concord$logFC.Minfi < -lfc_thresh &
      concord$logFC.Sesame < -lfc_thresh
    
    consensus_counts$up <- sum(is_up, na.rm = TRUE)
    consensus_counts$down <- sum(is_down, na.rm = TRUE)
    
    consensus_df <- concord[is_up | is_down, , drop = FALSE]
    if (nrow(consensus_df) > 0) {
      consensus_df$logFC_mean <- rowMeans(consensus_df[, c("logFC.Minfi", "logFC.Sesame")], na.rm = TRUE)
      consensus_df$P.Value <- pmax(consensus_df$P.Value.Minfi, consensus_df$P.Value.Sesame, na.rm = TRUE)
      consensus_df$adj.P.Val <- pmax(consensus_df$adj.P.Val.Minfi, consensus_df$adj.P.Val.Sesame, na.rm = TRUE)
      consensus_df <- consensus_df[order(consensus_df$P.Value, consensus_df$adj.P.Val), , drop = FALSE]
      
      out_cons_csv <- file.path(out_dir, "Intersection_Consensus_DMPs.csv")
      write.csv(consensus_df, out_cons_csv, row.names = FALSE)
      save_datatable(head(consensus_df, min(10000, nrow(consensus_df))), "Intersection_Consensus_DMPs.html", out_dir)
      message(paste("Consensus DMPs saved to", out_cons_csv))
    }
    
    # Concordance plot (logFC Minfi vs Sesame)
    concord$Consensus <- (is_up | is_down)
    plot_df <- concord
    if (nrow(plot_df) > max_points) {
      set.seed(12345)
      idx <- sample(seq_len(nrow(plot_df)), max_points)
      plot_df <- plot_df[idx, , drop = FALSE]
    }
    r_val <- suppressWarnings(cor(concord$logFC.Minfi, concord$logFC.Sesame, use = "complete.obs"))
    p_conc <- ggplot(plot_df, aes(x = logFC.Minfi, y = logFC.Sesame, color = Consensus)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#e74c3c")) +
      theme_minimal() +
      labs(
        title = "Minfi vs Sesame logFC concordance",
        subtitle = sprintf("Pearson r = %.3f (points subsampled to max_plots=%d for rendering)", r_val, max_points),
        x = "log2FC (Minfi)",
        y = "log2FC (Sesame)",
        color = "Consensus"
      )
    save_interactive_plot(p_conc, "Intersection_LogFC_Concordance.html", out_dir)
    save_static_plot(p_conc, "Intersection_LogFC_Concordance.png", out_dir, width = 6, height = 6)
    
    # Overlap summary plot (counts)
    sig_minfi <- concord$adj.P.Val.Minfi < pval_thresh & abs(concord$logFC.Minfi) > lfc_thresh
    sig_sesame <- concord$adj.P.Val.Sesame < pval_thresh & abs(concord$logFC.Sesame) > lfc_thresh
    n_both <- sum(is_up | is_down, na.rm = TRUE)
    overlap_df <- data.frame(
      Category = c("Minfi only", "Sesame only", "Both (consensus)"),
      Count = c(sum(sig_minfi, na.rm = TRUE) - n_both, sum(sig_sesame, na.rm = TRUE) - n_both, n_both)
    )
    p_ov <- ggplot(overlap_df, aes(x = Category, y = Count, fill = Category, text = Count)) +
      geom_col() +
      scale_fill_manual(values = c("Minfi only" = "#3498db", "Sesame only" = "#2ecc71", "Both (consensus)" = "#e74c3c")) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = "Significant DMP overlap (Minfi vs Sesame)", y = "Count", x = NULL)
    save_interactive_plot(p_ov, "Intersection_Significant_Overlap.html", out_dir)
    save_static_plot(p_ov, "Intersection_Significant_Overlap.png", out_dir, width = 7, height = 4)
    
  }, error = function(e) {
    message("Warning: Consensus (Intersection) outputs failed: ", e$message)
  })
}

get_sig_counts <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(list(up = 0, down = 0))
  up <- sum(res_df$adj.P.Val < pval_thresh & res_df$logFC > lfc_thresh, na.rm=TRUE)
  down <- sum(res_df$adj.P.Val < pval_thresh & res_df$logFC < -lfc_thresh, na.rm=TRUE)
  return(list(up=up, down=down))
}

minfi_counts <- get_sig_counts(res_minfi)
sesame_counts <- get_sig_counts(res_sesame)
intersect_counts <- consensus_counts
summary_n_con <- if (!is.null(minfi_out) && !is.null(minfi_out$n_con)) {
  minfi_out$n_con
} else if (!is.null(sesame_out) && !is.null(sesame_out$n_con)) {
  sesame_out$n_con
} else {
  n_con
}
summary_n_test <- if (!is.null(minfi_out) && !is.null(minfi_out$n_test)) {
  minfi_out$n_test
} else if (!is.null(sesame_out) && !is.null(sesame_out$n_test)) {
  sesame_out$n_test
} else {
  n_test
}

tryCatch({
  summary_json <- sprintf(
    '{"n_con": %d, "n_test": %d, "minfi_up": %d, "minfi_down": %d, "sesame_up": %d, "sesame_down": %d, "intersect_up": %d, "intersect_down": %d}', 
    summary_n_con, summary_n_test,
    minfi_counts$up, minfi_counts$down,
    sesame_counts$up, sesame_counts$down,
    intersect_counts$up, intersect_counts$down
  )
  
  writeLines(summary_json, file.path(out_dir, "summary.json"))
  message("Summary statistics saved to summary.json")
  
  sva_rule <- ifelse(disable_sva, "disabled", "best_method_or_no_batch")
  params_json <- sprintf(
    '{"pval_threshold": %.3f, "lfc_threshold": %.2f, "snp_maf": %.3f, "qc_intensity_threshold": %.1f, "detection_p_threshold": %.3f, "auto_covariate_alpha": %.3f, "auto_covariate_max_pcs": %d, "overcorrection_guard_ratio": %.2f, "sva_enabled": %s, "sva_inclusion_rule": "%s", "clock_covariates_enabled": %s, "sv_group_p_threshold": %.1e, "sv_group_eta2_threshold": %.2f, "pvca_min_samples": %d, "permutations": %d, "vp_top": %d, "dmr_maxgap": %d, "dmr_p_cutoff": %.3f, "logit_offset": %.6f, "seed": %d}',
    pval_thresh, lfc_thresh, SNP_MAF_THRESHOLD, QC_MEDIAN_INTENSITY_THRESHOLD,
    QC_DETECTION_P_THRESHOLD, AUTO_COVARIATE_ALPHA, MAX_PCS_FOR_COVARIATE_DETECTION,
    OVERCORRECTION_GUARD_RATIO, ifelse(disable_sva, "false", "true"), sva_rule, ifelse(include_clock_covariates, "true", "false"),
    SV_GROUP_P_THRESHOLD, SV_GROUP_ETA2_THRESHOLD, PVCA_MIN_SAMPLES, perm_n, vp_top, DMR_MAXGAP, DMR_P_CUTOFF, LOGIT_OFFSET, 12345
  )
  writeLines(params_json, file.path(out_dir, "analysis_parameters.json"))
  message("Analysis parameters saved to analysis_parameters.json")
  
  git_hash <- tryCatch(system("git rev-parse HEAD 2>/dev/null", intern=TRUE), error = function(e) character(0))
  if (length(git_hash) > 0) {
    writeLines(git_hash, file.path(out_dir, "code_version.txt"))
    message("Code version saved to code_version.txt")
  }
}, error = function(e) {
  message("Warning: failed to write summary/parameters: ", e$message)
})

message("Saving session info...")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

tryCatch({
  qc_path <- file.path(out_dir, "QC_Summary.csv")
  qc_tbl <- NULL
  if (file.exists(qc_path)) {
    qc_tbl <- read.csv(qc_path, stringsAsFactors = FALSE)
  }
  qc_val <- function(metric, default = NA) {
    if (is.null(qc_tbl) || nrow(qc_tbl) == 0) return(default)
    v <- qc_tbl$value[qc_tbl$metric == metric]
    if (length(v) == 0) return(default)
    v[1]
  }
  qc_intensity_str <- if (is.finite(QC_MEDIAN_INTENSITY_THRESHOLD)) {
    sprintf("median methylated/unmethylated intensity < %.1f (log2)", QC_MEDIAN_INTENSITY_THRESHOLD)
  } else {
    "disabled"
  }
  methods_lines <- c(
    "# IlluMeta analysis methods (auto-generated)",
    "",
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("- Config: `%s`", config_file),
    sprintf("- Output: `%s`", out_dir),
    sprintf("- Groups: control=`%s` (n=%s), test=`%s` (n=%s)", group_con_in, summary_n_con, group_test_in, summary_n_test),
    sprintf("- Tissue: `%s`", opt$tissue),
    sprintf("- Significance thresholds: BH FDR < %.3f and |log2FC| > %.2f", pval_thresh, lfc_thresh),
    "",
    "## Overview",
    "IlluMeta performs an end-to-end DNA methylation analysis from raw Illumina IDAT files, running two independent normalization pipelines (minfi and sesame) and reporting both per-pipeline results and a consensus (intersection) call set.",
    "",
    "## Data input",
    "- Raw IDATs are read from the project `idat/` directory.",
    sprintf("- Samples are defined by `%s` and `primary_group` in `configure.tsv`.", gsm_col),
    "",
    "## Sample-level QC",
    sprintf("- Detection P-value: samples with > %.0f%% probes failing (P > %.3g) are excluded.", QC_SAMPLE_DET_FAIL_FRAC * 100, QC_DETECTION_P_THRESHOLD),
    sprintf("- Signal intensity QC: %s.", qc_intensity_str),
    "- Mixed-array safeguard: samples whose IDAT array size deviates from the modal array size are excluded (unless `--force_idat`).",
    "",
    "## Probe-level QC (minfi-derived probe set)",
    sprintf("- Detection P-value filter: probes are retained only if P < %.3g in all retained samples.", QC_DETECTION_P_THRESHOLD),
    sprintf("- SNP filtering: probes overlapping common SNPs (SBE/CpG; MAF ≥ %.2f) are removed using `minfi::dropLociWithSnps()`.", SNP_MAF_THRESHOLD),
    "- Sex chromosome probes (chrX/chrY) are removed.",
    "",
    "## Normalization (two pipelines)",
    "- **minfi**: `preprocessNoob()` followed by `mapToGenome()`; beta values are extracted with `getBeta()`.",
    "- **sesame**: `noob()` + `dyeBiasCorrTypeINorm()`; beta values are extracted with `getBetas()`.",
    "",
    "## Covariates and batch control",
    sprintf("- Automatic covariate discovery: metadata variables associated with the top PCs (alpha=%.3f; up to %d PCs) are considered, then filtered for stability/confounding.", AUTO_COVARIATE_ALPHA, MAX_PCS_FOR_COVARIATE_DETECTION),
    "- Small-n safeguard: if covariate count would eliminate residual degrees of freedom, covariates are capped using PC-association/variance ranking to preserve model stability.",
    sprintf("- Overcorrection guard: total covariates + SVs are capped at %.0f%% of sample size (excess terms are dropped to preserve power).", OVERCORRECTION_GUARD_RATIO * 100),
    "- Cell composition: if `--tissue Auto`, reference-free deconvolution is performed via RefFreeEWAS (K=5 latent components); these components are optionally added as covariates after basic sanity checks.",
    sprintf("- Surrogate variable analysis (SVA): enabled unless `--disable_sva`; SVs are estimated on top-variable probes and included in the model only when selected as the best batch strategy or when no batch factor is evaluated (to avoid double correction). SVs strongly associated with the group (P < %.1e or Eta^2 > %.2f) are excluded to avoid over-correction.", SV_GROUP_P_THRESHOLD, SV_GROUP_ETA2_THRESHOLD),
    "- Epigenetic clock covariates: when `--include_clock_covariates` is enabled, clock outputs are merged into metadata and considered by auto covariate selection (clocks with missing/constant values are excluded).",
    "- Batch method comparison: if a suitable unconfounded batch factor exists, multiple correction strategies are evaluated and a best method is chosen for modeling.",
    "",
    "## Differential methylation (DMP)",
    sprintf("- Beta values are transformed to M-values with a logit offset of %.6f.", LOGIT_OFFSET),
    "- Differential methylation is tested using limma (`lmFit`/`eBayes`) with a group contrast (test − control).",
    "- Multiple testing is controlled by Benjamini–Hochberg FDR.",
    "",
    "## Differentially methylated regions (DMR)",
    sprintf("- DMRs are called with `dmrff` (maxgap=%d; p.cutoff=%.3f).", DMR_MAXGAP, DMR_P_CUTOFF),
    "",
    "## Consensus (intersection) call set",
    "- Consensus DMPs are defined as CpGs significant in **both** minfi and sesame with the **same direction** under the same thresholds.",
    "- Consensus outputs: `Intersection_Consensus_DMPs.csv`, `Intersection_Consensus_DMPs.html`, and concordance/overlap plots.",
    "",
    "## Reproducibility artifacts",
    "- `analysis_parameters.json`: run parameters and thresholds.",
    "- `sessionInfo.txt`: full R session/package versions.",
    "- `code_version.txt`: git commit hash (when available).",
    "",
    "## QC summary (from QC_Summary.csv)",
    sprintf("- Total_samples_input: %s", qc_val("Total_samples_input", NA)),
    sprintf("- Samples_failed_QC: %s", qc_val("Samples_failed_QC", NA)),
    sprintf("- Samples_passed_QC: %s", qc_val("Samples_passed_QC", NA)),
    sprintf("- Total_probes_raw: %s", qc_val("Total_probes_raw", NA)),
    sprintf("- Probes_final: %s", qc_val("Probes_final", NA))
  )
  writeLines(methods_lines, file.path(out_dir, "methods.md"))
  message("Methods summary saved to methods.md")
}, error = function(e) {
  message("Warning: failed to write methods.md: ", e$message)
})

message("Analysis Complete.")
