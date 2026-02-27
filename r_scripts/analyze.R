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

# Respect ILLUMETA_SESAME_SINGLE_THREAD env var (default "1" = single-thread)
if (Sys.getenv("ILLUMETA_SESAME_SINGLE_THREAD", "1") == "1") {
  force_single_thread()
}

# ── Headless rendering ──────────────────────────────────────────────────────
# ggplotly() internally opens a temporary graphics device.  On headless
# systems (Docker, HPC, SSH) this fails with "X11 is not available" (Linux)
# or "unable to open connection to QuartzCore" (macOS).  We detect whether a
# working raster device is available and, if not, fall back to cairo or pdf.
.illumeta_ensure_graphics_device <- function() {
  # Quick probe: can the current default png() actually open?
  ok <- tryCatch({
    tf <- tempfile(fileext = ".png")
    on.exit(unlink(tf), add = TRUE)
    grDevices::png(tf, width = 1, height = 1)
    grDevices::dev.off()
    TRUE
  }, error = function(e) FALSE)

  if (ok) return(invisible(NULL))  # default device works — nothing to do

  # Default device failed.  Try cairo, then pdf as last resort.
  if (capabilities("cairo")) {
    options(bitmapType = "cairo")
    options(device = function(...) grDevices::png(type = "cairo", ...))
  } else {
    options(device = grDevices::pdf)
  }
  invisible(NULL)
}
.illumeta_ensure_graphics_device()

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
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(reformulas))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(pvca))
suppressPackageStartupMessages(library(illuminaio))

if (Sys.getenv("ILLUMETA_TRACEBACK", unset = "") == "1") {
  options(error = function() {
    tb <- capture.output(traceback(2))
    message("TRACEBACK:")
    message(paste(tb, collapse = "\n"))
    quit(status = 1)
  })
}

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

# =============================================================================
# IlluMeta Statistical Analysis Pipeline (analyze.R)
# =============================================================================
#
# DESCRIPTION:
#   Core statistical analysis engine for IlluMeta. Implements dual-pipeline
#   methylation analysis using minfi (Noob normalization) and sesame, with
#   consensus intersection for high-confidence CpG calls.
#
# MAIN FEATURES:
#   - Dual-pipeline normalization (minfi Noob + sesame)
#   - Batch effect correction (SVA, ComBat, limma::removeBatchEffect)
#   - Variance partitioning (PVCA) for batch assessment
#   - Differential methylation analysis (DMPs via limma, DMRs via dmrff)
#   - Cell composition deconvolution (EpiDISH, RefFreeEWAS)
#   - Sample-size adaptive robustness (CRF)
#   - Epigenetic clock estimation (methylclock, planet)
#
# STATISTICAL METHODS AND REFERENCES:
#   - Normalization: Triche et al. (2013) Nucleic Acids Res - Noob method
#   - M-value transformation: Du et al. (2010) BMC Bioinformatics
#   - Differential methylation: Ritchie et al. (2015) Nucleic Acids Res - limma
#   - Batch effect correction: Leek et al. (2012) Nat Rev Genet - SVA
#   - DMR analysis: Suderman et al. (2018) bioRxiv doi:10.1101/508556 - dmrff
#   - Multiple testing: Benjamini & Hochberg (1995) JRSS-B - FDR
#   - Cell deconvolution: Teschendorff et al. (2017) BMC Bioinformatics - EpiDISH
#   - Variance partition: Hoffman & Schadt (2016) BMC Bioinformatics
#
# KNOWN LIMITATIONS:
#   - Sesame pthread errors may occur on some systems; single-thread mode is
#     enforced by default (see force_single_thread())
#   - Small sample sizes (n < 12) trigger "minimal" CRF tier with limited
#     statistical power warnings
#   - EPIC v2 requires R 4.4+ and Bioconductor 3.19+
#
# AUTHOR: Keunsoo Kang
# LICENSE: Apache-2.0
# =============================================================================

option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL, 
              help="Path to configure.tsv", metavar="character"),
  make_option(c("--config_yaml"), type="character", default="",
              help="Optional YAML config (default: config.yaml next to configure.tsv)", metavar="character"),
  make_option(c("--preset"), type="character", default="conservative",
              help="Optimization preset (conservative or aggressive)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="Output directory for results", metavar="character"),
  make_option(c("--idat_dir"), type="character", default="",
              help="Path to IDAT directory (default: [config dir]/idat)", metavar="character"),
  make_option(c("--skip-sesame"), action="store_true", default=FALSE,
              help="Skip Sesame pipeline and intersection outputs (Minfi only)."),
  make_option(c("--sesame_typeinorm"), action="store_true", default=FALSE,
              help="Enable sesame dyeBiasCorrTypeINorm (default: disabled for stability)."),
  make_option(c("-m", "--max_plots"), type="integer", default=10000, 
              help="Max points for interactive plots (default: 10000)", metavar="integer"),
  make_option(c("-p", "--pval"), type="double", default=0.05, 
              help="Adjusted P-value threshold (default: 0.05)", metavar="double"),
  make_option(c("-l", "--lfc"), type="double", default=0.5, 
              help="LogFC threshold (default: 0.5)", metavar="double"),
  make_option(c("--delta_beta"), type="double", default=0, 
              help="Absolute Delta Beta threshold (default: 0; disabled)", metavar="double"),
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
  make_option(c("--beginner_safe"), action="store_true", default=FALSE,
              help="Enable beginner-safe mode (stricter group checks and conservative thresholds)."),
  make_option(c("--marker_list"), type="character", default="",
              help="Optional CpG marker list (TSV/CSV with CpG column or one CpG per line) for signal preservation checks"),
  make_option(c("--cross_reactive_list"), type="character", default="",
              help="Optional cross-reactive probe list (TSV/CSV with CpG column or one CpG per line) for filtering"),
  make_option(c("--unsafe-skip-cross-reactive"), action="store_true", default=FALSE,
              help="UNSAFE: skip mandatory cross-reactive probe filtering (requires explicit flag)."),
  make_option(c("--sex-mismatch-action"), type="character", default="",
              help="Sex mismatch check action: stop|drop|ignore (default: stop)."),
  make_option(c("--sex_check_column"), type="character", default="",
              help="Metadata column to use for sex mismatch check (default: auto-detect)."),
  make_option(c("--permutations"), type="integer", default=20,
              help="Number of label permutations for null DMP counts (default: 20; set 0 to use config minimum; set calibration.permutations_min=0 to disable)"),
  make_option(c("--vp_top"), type="integer", default=5000,
              help="Number of top-variable CpGs for variancePartition (default: 5000)"),
  make_option(c("--tissue"), type="character", default="Auto", 
              help="Tissue type for cell deconvolution (default: Auto = reference-free). Supported: Auto (RefFreeEWAS), Blood, CordBlood, DLPFC, Saliva, Placenta"),
  make_option(c("--cell_reference"), type="character", default="",
              help="Custom cell reference (package name, package::object, or .rds/.rda path)"),
  make_option(c("--cell_reference_platform"), type="character", default="",
              help="Reference platform string for custom cell reference (e.g. IlluminaHumanMethylationEPIC)"),
  make_option(c("--positive_controls"), type="character", default=NULL,
              help="Comma-separated list of gene symbols (e.g. 'AHRR,CYP1A1') to verify as positive controls."),
  make_option(c("--id_column"), type="character", default="",
              help="Column in configure.tsv to treat as sample ID (for non-GEO datasets). If empty, tries GSM/geo_accession."),
  make_option(c("--min_total_size"), type="integer", default=6,
              help="Minimum total sample size required to proceed (default: 6)."),
  make_option(c("--qc_intensity_threshold"), type="double", default=9.0,
              help="Median M/U signal intensity threshold for sample QC (log2). Set <=0 to disable intensity-based sample drop (default: 9.0)"),
  make_option(c("--qc_detection_p_threshold"), type="double", default=NA_real_,
              help="Detection P-value threshold for sample/probe QC (default: 0.05)."),
  make_option(c("--qc_sample_fail_frac"), type="double", default=NA_real_,
              help="Sample fail fraction for detection P-value QC (default: 0.20)."),
  make_option(c("--dmr_min_cpgs"), type="integer", default=NA_integer_,
              help="Minimum CpGs per DMR to retain (default: 2)."),
  make_option(c("--dmr_maxgap"), type="integer", default=NA_integer_,
              help="Maximum gap between CpGs for DMRs (default: 500 or config)."),
  make_option(c("--beginner_safe_delta_beta"), type="double", default=NA_real_,
              help="Minimum |DeltaBeta| enforced in beginner-safe mode (default: 0.05 or config)."),
  make_option(c("--logit_offset"), type="double", default=NA_real_,
              help="Logit offset for M-value transform (default: 1e-4 or config)."),
  make_option(c("--batch_column"), type="character", default="",
              help="Override batch column name for correction (skips auto selection)."),
  make_option(c("--batch_method"), type="character", default="",
              help="Override batch correction method: none|combat|limma|sva (skips auto selection)."),
  make_option(c("--tier3-min-total-n"), type="integer", default=NA_integer_,
              help="Tier3 minimum total N required to run stratified/meta (default: 20)."),
  make_option(c("--tier3-min-per-group-per-stratum"), type="integer", default=NA_integer_,
              help="Tier3 minimum samples per group per stratum (default: 5)."),
  make_option(c("--tier3-on-fail"), type="character", default="",
              help="Tier3 action when ineligible: stop|skip (default: stop)."),
  make_option(c("--auto-covariates-enabled"), type="character", default="",
              help="Enable/disable auto covariates (true/false); overrides config."),
  make_option(c("--auto-covariates-exclude-group-associated"),
              type="character", default="", help="Exclude covariates strongly associated with group (true/false)."),
  make_option(c("--auto-covariates-group-assoc-p-threshold"),
              type="double", default=NA_real_, help="P-value threshold for group association exclusion."),
  make_option(c("--auto-covariates-max-cor"),
              type="double", default=NA_real_, help="Max allowed absolute correlation among auto covariates."),
  make_option(c("--cell-adjustment-on-high-eta2"),
              type="character", default="", help="Action when cell vs group Eta^2 exceeds threshold: warn|stop."),
  make_option(c("--force_idat"), action="store_true", default=FALSE,
              help="Force reading IDATs when array sizes differ but types are similar (passes force=TRUE to read.metharray).")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

# Timestamped messages for all logging in this script
ts_message <- function(..., domain = NULL, appendLF = TRUE) {
  base::message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " ")),
                domain = domain, appendLF = appendLF)
}
message <- ts_message

# QC and filtering thresholds (kept as constants for transparency/reuse)
QC_MEDIAN_INTENSITY_THRESHOLD <- 9.0
QC_DETECTION_P_THRESHOLD <- 0.05
QC_SAMPLE_DET_FAIL_FRAC <- 0.20
SNP_MAF_THRESHOLD <- 0.01
LOGIT_OFFSET <- 1e-4
BETA_RANGE_MIN <- 0.05
SESAME_NATIVE_NA_MAX_FRAC <- 0.1
SESAME_NATIVE_KNN_K <- 10
SESAME_NATIVE_KNN_MAX_ROWS <- 50000
SESAME_NATIVE_KNN_REF_ROWS <- 20000
SESAME_NATIVE_IMPUTE_METHOD <- "knn"
AUTO_COVARIATE_ALPHA <- 0.01
MAX_PCS_FOR_COVARIATE_DETECTION <- 5
OVERCORRECTION_GUARD_RATIO <- 0.25
UNDERCORRECTION_GUARD_MIN_SV <- 1
SV_GROUP_P_THRESHOLD <- 1e-6
SV_GROUP_ETA2_THRESHOLD <- 0.5
PVCA_MIN_SAMPLES <- 8
DMR_MAXGAP <- 500
DMR_P_CUTOFF_DEFAULT <- 0.05
DMR_LABEL_TOP_N <- 15
DMR_MIN_CPGS <- 2
BATCH_EVAL_TOP_VAR <- 20000
MIN_GROUP_SIZE_WARN <- 3
MIN_TOTAL_SIZE_STOP <- NULL

CONFIG_DEFAULTS <- list(
  preset = "conservative",
  min_overlap_per_cell = 2,
  logit_offset = 1e-4,
  qc = list(
    detection_p_threshold = 0.05,
    sample_fail_frac = 0.20
  ),
  dmr = list(
    min_cpgs = 2,
    maxgap = 500
  ),
  beginner_safe_delta_beta = 0.05,
  confounding = list(
    tier1_v = 0.2,
    tier2_v = 0.4,
    tier1_r2 = 0.1,
    tier2_r2 = 0.25
  ),
  batch_candidate_patterns = c("sentrix", "slide", "array", "plate", "chip", "batch"),
  batch_candidates = NULL,
  batch_override = list(
    column = "",
    method = ""
  ),
  forced_adjust = character(0),
  allowed_adjust = character(0),
  do_not_adjust = character(0),
  missingness = list(max_frac = 0.2, impute_max_frac = 0.1),
  collinearity = list(cor_threshold = 0.9),
  calibration = list(ks_p = 0.05, permutations_min = 20),
  unsafe = list(
    skip_cross_reactive = FALSE,
    skip_sex_check = FALSE
  ),
  caf = list(
    enabled = TRUE,
    profile = "auto",
    target_fpr = 0.05,
    weights = list(
      conservative = list(calibration = 0.40, preservation = 0.35, batch = 0.25, ncs = 0.00),
      discovery = list(calibration = 0.25, preservation = 0.25, batch = 0.50, ncs = 0.00)
    )
  ),
  crf = list(
    enabled = TRUE,
    max_probes_mmc = 20000,
    max_probes_rss = 20000,
    housekeeping_path = "",
    ncs_lambda_ci = list(
      enabled = TRUE,
      n_boot = 200,
      sample_frac = 0.8
    ),
    sample_tiers = list(
      minimal = list(
        threshold = 12,
        min_per_group = 3,
        mmc = list(methods = c("none", "limma"), concordance_levels = c("core", "discordant")),
        ncs = list(types = c("snp"), lambda_range = c(0.7, 1.3), interpret = "reference_only"),
        rss = list(mode = "bootstrap", n_iterations = 100, top_k = c(20, 50), overlap_threshold = 0.20),
        warnings = c("Very small sample size: interpret all results with extreme caution.",
                     "Statistical power severely limited.",
                     "Consider this analysis exploratory only.")
      ),
      small = list(
        threshold = 24,
        min_per_group = 6,
        mmc = list(methods = c("none", "combat_np", "limma"), concordance_levels = c("core", "consensus", "weak")),
        ncs = list(types = c("snp"), lambda_range = c(0.8, 1.2), interpret = "cautious"),
        rss = list(mode = "leave_pair_out", n_iterations = 50, top_k = c(30, 100, 200), overlap_threshold = 0.30),
        warnings = c("Small sample size: results should be validated independently.",
                     "SVA excluded due to instability.")
      ),
      moderate = list(
        threshold = 50,
        min_per_group = 12,
        mmc = list(methods = c("none", "combat", "limma", "sva"), concordance_levels = c("core", "consensus", "weak"),
                   sva_note = "SVA may be slightly unstable."),
        ncs = list(types = c("snp", "housekeeping"), lambda_range = c(0.9, 1.1), interpret = "standard"),
        rss = list(mode = "split", n_iterations = 10, top_k = c(50, 200, 500), overlap_threshold = 0.40),
        warnings = character(0)
      ),
      large = list(
        threshold = NA,
        min_per_group = 25,
        mmc = list(methods = c("none", "combat", "limma", "sva"), concordance_levels = c("core", "consensus", "weak")),
        ncs = list(types = c("snp", "housekeeping"), lambda_range = c(0.9, 1.1), interpret = "standard"),
        rss = list(mode = "split", n_iterations = 10, top_k = c(100, 500, 1000), overlap_threshold = 0.50),
        warnings = character(0)
      )
    )
  ),
  auto_covariates = list(
    enabled = TRUE,
    exclude_if_group_associated = TRUE,
    group_assoc_p_threshold = 1e-3,
    max_cor = 0.9
  ),
  lambda_guard = list(
    enabled = TRUE,
    threshold = 1.5,
    min_samples = 8,
    action = "warn"
  ),
  sva = list(
    group_p_threshold = 1e-6,
    group_eta2_threshold = 0.5
  ),
  variance_partition = list(
    autoscale_numeric = TRUE,
    autoscale_on_fail = TRUE
  ),
  tier3_meta = list(
    method = "auto",
    i2_threshold = 0.5,
    min_total_n = 20,
    min_per_group_per_stratum = 5,
    on_fail = "stop",
    min_total_warn = 20,
    min_stratum_warn = 6
  ),
  cell_adjustment = list(
    on_high_eta2 = "warn"
  ),
  cell_deconv = list(
    ref_free_k = 5,
    ref_free_max_probes = 10000,
    confound_eta2_threshold = 0.5,
    confound_action = "warn" # warn | drop
  ),
  markers = list(
    path = ""
  ),
  cross_reactive = list(
    enabled = TRUE,
    path = "",
    local_dir = "references/probe_blacklists",
    use_maxprobes = TRUE
  ),
  sex_check = list(
    enabled = TRUE,
    action = "stop",
    column = ""
  ),
  scoring_presets = list(
    conservative = list(
      weights = list(batch = 0.20, bio = 0.35, cal = 0.25, stab = 0.20),
      guards = list(batch_reduction_min = 0.20, pc_batch_r2_reduction = 0.30, bio_min = 0.85),
      top_k = 5
    ),
    aggressive = list(
      weights = list(batch = 0.35, bio = 0.25, cal = 0.15, stab = 0.25),
      guards = list(batch_reduction_min = 0.25, mix_improve = 0.05, bio_min = 0.80),
      top_k = 5
    )
  )
)

merge_config_defaults <- function(base, override) {
  if (is.null(override) || length(override) == 0) return(base)
  out <- base
  for (nm in names(override)) {
    if (is.list(base[[nm]]) && is.list(override[[nm]])) {
      out[[nm]] <- merge_config_defaults(base[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

load_yaml_config <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    message("Warning: yaml package not available; config.yaml ignored.")
    return(NULL)
  }
  tryCatch(yaml::read_yaml(path), error = function(e) {
    message("Warning: failed to parse config.yaml: ", e$message)
    NULL
  })
}

write_yaml_safe <- function(obj, path) {
  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(obj, path)
    return(invisible(TRUE))
  }
  lines <- c()
  for (nm in names(obj)) {
    val <- obj[[nm]]
    if (is.list(val)) {
      lines <- c(lines, paste0(nm, ":"))
      sublines <- paste0("  ", names(val), ": ", unlist(val))
      lines <- c(lines, sublines)
    } else if (length(val) > 1) {
      lines <- c(lines, paste0(nm, ": [", paste(val, collapse = ", "), "]"))
    } else {
      lines <- c(lines, paste0(nm, ": ", val))
    }
  }
  writeLines(lines, path)
  invisible(TRUE)
}

normalize_preset <- function(preset) {
  if (is.null(preset) || !nzchar(preset)) return("conservative")
  preset <- tolower(trimws(preset))
  if (!preset %in% c("conservative", "aggressive")) preset <- "conservative"
  preset
}

parse_bool_flag <- function(x, default = NA) {
  if (is.null(x)) return(default)
  if (is.logical(x) && length(x) == 1) return(x)
  if (!nzchar(as.character(x))) return(default)
  val <- tolower(trimws(as.character(x)))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

decision_log <- data.frame(
  timestamp = character(0),
  stage = character(0),
  decision = character(0),
  value = character(0),
  reason = character(0),
  metrics = character(0),
  stringsAsFactors = FALSE
)

log_decision <- function(stage, decision, value, reason = "", metrics = NULL) {
  metrics_str <- ""
  if (!is.null(metrics) && length(metrics) > 0) {
    metrics_str <- paste(sprintf("%s=%s", names(metrics), metrics), collapse = ";")
  }
  decision_log <<- rbind(
    decision_log,
    data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stage = stage,
      decision = decision,
      value = value,
      reason = reason,
      metrics = metrics_str,
      stringsAsFactors = FALSE
    )
  )
}

safe_ebayes <- function(fit, label = "", robust = TRUE) {
  if (is.null(fit)) return(NULL)
  out <- tryCatch(eBayes(fit, robust = robust), error = function(e) {
    if (robust) {
      # Retry without robust if it fails (e.g., too few residuals)
      out2 <- tryCatch(eBayes(fit, robust = FALSE), error = function(e2) NULL)
      if (!is.null(out2)) {
        msg <- if (nzchar(label)) paste0(label, ": ") else ""
        message("  eBayes robust=TRUE failed; falling back to robust=FALSE: ", msg, conditionMessage(e))
        return(out2)
      }
    }
    msg <- if (nzchar(label)) paste0(label, ": ") else ""
    message("  eBayes skipped: ", msg, conditionMessage(e))
    return(NULL)
  })
  if (is.null(out)) return(NULL)
  if (!is.null(out$sigma) && all(!is.finite(out$sigma))) {
    msg <- if (nzchar(label)) paste0(label, ": ") else ""
    message("  eBayes skipped: ", msg, "no finite residual standard deviations.")
    return(NULL)
  }
  out
}

if (is.null(opt$config) || is.null(opt$group_con) || is.null(opt$group_test)){
  print_help(opt_parser)
  stop("Configuration file, group_con, and group_test must be supplied", call.=FALSE)
}

config_file <- normalizePath(opt$config, winslash = "/", mustWork = FALSE)
project_dir <- dirname(config_file)
config_yaml_path <- opt$config_yaml
if (nzchar(config_yaml_path) && !grepl("^(/|[A-Za-z]:)", config_yaml_path)) {
  # Relative config_yaml: resolve relative to configure.tsv directory, not CWD
  config_yaml_path <- normalizePath(file.path(project_dir, config_yaml_path), winslash = "/", mustWork = FALSE)
} else if (nzchar(config_yaml_path)) {
  config_yaml_path <- normalizePath(config_yaml_path, winslash = "/", mustWork = FALSE)
}
if (!nzchar(config_yaml_path)) {
  config_yaml_path <- file.path(dirname(config_file), "config.yaml")
}
config_raw <- load_yaml_config(config_yaml_path)
config_settings <- merge_config_defaults(CONFIG_DEFAULTS, config_raw)
preset_name <- normalize_preset(config_settings$preset)
cli_preset <- normalize_preset(opt$preset)
if (!identical(cli_preset, "conservative")) {
  preset_name <- cli_preset
}
config_settings$preset <- preset_name
scoring_preset <- config_settings$scoring_presets[[preset_name]]
log_decision("config", "preset", preset_name,
             reason = if (file.exists(config_yaml_path)) "config_yaml_or_cli" else "default")
out_dir <- opt$out
max_points <- opt$max_plots
pval_thresh <- opt$pval
lfc_thresh <- opt$lfc
beginner_safe <- isTRUE(opt$beginner_safe)
delta_beta_thresh <- suppressWarnings(as.numeric(opt$delta_beta))
if (!is.finite(delta_beta_thresh)) delta_beta_thresh <- 0
delta_beta_thresh <- abs(delta_beta_thresh)
beginner_safe_delta_beta <- suppressWarnings(as.numeric(if (!is.null(config_settings$beginner_safe_delta_beta))
  config_settings$beginner_safe_delta_beta else 0.05))
if (!is.na(opt$beginner_safe_delta_beta)) beginner_safe_delta_beta <- opt$beginner_safe_delta_beta
if (!is.finite(beginner_safe_delta_beta) || beginner_safe_delta_beta < 0) {
  beginner_safe_delta_beta <- 0.05
}
if (beginner_safe && (is.na(delta_beta_thresh) || delta_beta_thresh < beginner_safe_delta_beta)) {
  delta_beta_thresh <- beginner_safe_delta_beta
}
group_con_in <- opt$group_con
group_test_in <- opt$group_test
perm_n <- opt$permutations
if (!is.finite(perm_n) || perm_n < 0) perm_n <- 0
if (perm_n == 0) perm_n <- config_settings$calibration$permutations_min
vp_top <- opt$vp_top
disable_auto_cov <- opt$disable_auto_covariates
disable_sva <- opt$disable_sva
include_cov <- ifelse(opt$include_covariates == "", character(0), trimws(strsplit(opt$include_covariates, ",")[[1]]))
include_clock_covariates <- opt$include_clock_covariates
marker_list_path <- opt$marker_list
if (!nzchar(marker_list_path) && !is.null(config_settings$markers$path)) {
  marker_list_path <- as.character(config_settings$markers$path)
}
if (nzchar(marker_list_path) && !grepl("^/", marker_list_path)) {
  marker_list_path <- file.path(project_dir, marker_list_path)
}
if (nzchar(marker_list_path) && !file.exists(marker_list_path)) {
  warning(sprintf("Marker list not found: %s (skipping marker preservation checks)", marker_list_path))
  marker_list_path <- ""
}
cross_reactive_path <- opt$cross_reactive_list
if (!nzchar(cross_reactive_path) && !is.null(config_settings$cross_reactive$path)) {
  cross_reactive_path <- as.character(config_settings$cross_reactive$path)
}
if (nzchar(cross_reactive_path) && !grepl("^/", cross_reactive_path)) {
  cross_reactive_path <- file.path(project_dir, cross_reactive_path)
}
if (nzchar(cross_reactive_path) && !file.exists(cross_reactive_path)) {
  warning(sprintf("Cross-reactive probe list not found: %s", cross_reactive_path))
  cross_reactive_path <- ""
}
cross_reactive_enabled <- isTRUE(config_settings$cross_reactive$enabled)
cross_reactive_use_maxprobes <- isTRUE(config_settings$cross_reactive$use_maxprobes)
cross_reactive_local_dir <- ""
if (!is.null(config_settings$cross_reactive$local_dir)) {
  cross_reactive_local_dir <- as.character(config_settings$cross_reactive$local_dir)
}
if (nzchar(cross_reactive_local_dir) && !grepl("^/", cross_reactive_local_dir)) {
  cross_reactive_local_dir <- file.path(project_dir, cross_reactive_local_dir)
}
unsafe_skip_cross_reactive <- isTRUE(config_settings$unsafe$skip_cross_reactive)
if (isTRUE(opt$unsafe_skip_cross_reactive)) {
  unsafe_skip_cross_reactive <- TRUE
}
if (!cross_reactive_enabled && !unsafe_skip_cross_reactive) {
  message("  Cross-reactive filtering is mandatory; overriding config to enable it (use --unsafe-skip-cross-reactive to bypass).")
  cross_reactive_enabled <- TRUE
}
if (nzchar(cross_reactive_path)) cross_reactive_enabled <- TRUE

sex_check_enabled <- isTRUE(config_settings$sex_check$enabled)
unsafe_skip_sex_check <- isTRUE(config_settings$unsafe$skip_sex_check)
sex_check_action <- opt$sex_mismatch_action
if (!nzchar(sex_check_action) && !is.null(config_settings$sex_check$action)) {
  sex_check_action <- as.character(config_settings$sex_check$action)
}
sex_check_action <- tolower(trimws(sex_check_action))
if (!nzchar(sex_check_action)) sex_check_action <- "stop"
if (!sex_check_action %in% c("stop", "drop", "ignore")) {
  warning(sprintf("Unknown sex mismatch action '%s'; defaulting to 'stop'.", sex_check_action))
  sex_check_action <- "stop"
}
if (unsafe_skip_sex_check) sex_check_action <- "ignore"
if (!sex_check_enabled && !unsafe_skip_sex_check) {
  message("  Sex check is mandatory for safety; overriding config to enable it (use --sex-mismatch-action ignore to bypass).")
  sex_check_enabled <- TRUE
}
sex_check_column <- opt$sex_check_column
if (is.null(sex_check_column) || !nzchar(sex_check_column)) {
  sex_check_column <- ""
}
if (!nzchar(sex_check_column) && !is.null(config_settings$sex_check$column)) {
  sex_check_column <- as.character(config_settings$sex_check$column)
}

auto_cov_enabled <- isTRUE(config_settings$auto_covariates$enabled)
auto_cov_enabled_cli <- parse_bool_flag(opt$auto_covariates_enabled, default = NA)
if (!is.na(auto_cov_enabled_cli)) auto_cov_enabled <- auto_cov_enabled_cli
if (isTRUE(disable_auto_cov)) auto_cov_enabled <- FALSE
auto_cov_exclude_group <- isTRUE(config_settings$auto_covariates$exclude_if_group_associated)
auto_cov_exclude_cli <- parse_bool_flag(opt$auto_covariates_exclude_group_associated, default = NA)
if (!is.na(auto_cov_exclude_cli)) auto_cov_exclude_group <- auto_cov_exclude_cli
auto_cov_group_p <- suppressWarnings(as.numeric(config_settings$auto_covariates$group_assoc_p_threshold))
if (!is.na(opt$auto_covariates_group_assoc_p_threshold)) {
  auto_cov_group_p <- opt$auto_covariates_group_assoc_p_threshold
}
if (!is.finite(auto_cov_group_p) || auto_cov_group_p <= 0) auto_cov_group_p <- 1e-3
auto_cov_max_cor <- suppressWarnings(as.numeric(config_settings$auto_covariates$max_cor))
if (!is.na(opt$auto_covariates_max_cor)) auto_cov_max_cor <- opt$auto_covariates_max_cor
if (!is.finite(auto_cov_max_cor)) auto_cov_max_cor <- NA_real_
batch_col_override <- opt$batch_column
if (!nzchar(batch_col_override) && !is.null(config_settings$batch_override$column)) {
  batch_col_override <- as.character(config_settings$batch_override$column)
}
batch_method_override <- opt$batch_method
if (!nzchar(batch_method_override) && !is.null(config_settings$batch_override$method)) {
  batch_method_override <- as.character(config_settings$batch_override$method)
}
batch_method_override <- tolower(trimws(batch_method_override))
if (nzchar(batch_method_override) && !batch_method_override %in% c("none", "combat", "limma", "sva")) {
  warning(sprintf("Unknown batch_method override '%s'; ignoring.", batch_method_override))
  batch_method_override <- ""
}
cell_reference <- opt$cell_reference
cell_reference_platform <- opt$cell_reference_platform
tissue_use <- opt$tissue
tissue_source <- "arg"
cell_covariates <- character(0)
cell_counts_df <- NULL
sesame_native_impute_method <- "none"
sesame_cell_est <- NULL
sesame_cell_reference <- "minfi-derived"
id_col_override <- opt$id_column
MIN_TOTAL_SIZE_STOP <- max(2, opt$min_total_size)
if (beginner_safe && MIN_TOTAL_SIZE_STOP < 8) {
  MIN_TOTAL_SIZE_STOP <- 8
}
force_idat <- opt$force_idat
QC_MEDIAN_INTENSITY_THRESHOLD <- ifelse(opt$qc_intensity_threshold <= 0, -Inf, opt$qc_intensity_threshold)
qc_det_p <- config_settings$qc$detection_p_threshold
if (!is.na(opt$qc_detection_p_threshold)) qc_det_p <- opt$qc_detection_p_threshold
QC_DETECTION_P_THRESHOLD <- qc_det_p
qc_fail_frac <- config_settings$qc$sample_fail_frac
if (!is.na(opt$qc_sample_fail_frac)) qc_fail_frac <- opt$qc_sample_fail_frac
QC_SAMPLE_DET_FAIL_FRAC <- qc_fail_frac
logit_offset_cfg <- if (!is.null(config_settings$logit_offset)) as.numeric(config_settings$logit_offset) else NA_real_
if (!is.na(opt$logit_offset)) logit_offset_cfg <- opt$logit_offset
if (is.finite(logit_offset_cfg) && logit_offset_cfg > 0 && logit_offset_cfg < 0.5) {
  LOGIT_OFFSET <- logit_offset_cfg
} else if (!is.na(logit_offset_cfg)) {
  warning(sprintf("Invalid logit_offset %.6f; using default %.6f.", logit_offset_cfg, LOGIT_OFFSET))
}
sva_cfg <- config_settings$sva
if (!is.null(sva_cfg$group_p_threshold)) {
  sva_p <- suppressWarnings(as.numeric(sva_cfg$group_p_threshold))
  if (is.finite(sva_p) && sva_p > 0) SV_GROUP_P_THRESHOLD <- sva_p
}
if (!is.null(sva_cfg$group_eta2_threshold)) {
  sva_eta2 <- suppressWarnings(as.numeric(sva_cfg$group_eta2_threshold))
  if (is.finite(sva_eta2) && sva_eta2 > 0) SV_GROUP_ETA2_THRESHOLD <- sva_eta2
}
cell_deconv_cfg <- config_settings$cell_deconv
cell_ref_free_k <- suppressWarnings(as.integer(if (!is.null(cell_deconv_cfg$ref_free_k)) cell_deconv_cfg$ref_free_k else 5))
if (!is.finite(cell_ref_free_k) || cell_ref_free_k < 2) cell_ref_free_k <- 5
cell_ref_free_max_probes <- suppressWarnings(as.integer(if (!is.null(cell_deconv_cfg$ref_free_max_probes)) cell_deconv_cfg$ref_free_max_probes else 10000))
if (!is.finite(cell_ref_free_max_probes) || cell_ref_free_max_probes < 1000) cell_ref_free_max_probes <- 10000
cell_confound_eta2_threshold <- suppressWarnings(as.numeric(if (!is.null(cell_deconv_cfg$confound_eta2_threshold))
  cell_deconv_cfg$confound_eta2_threshold else 0.5))
if (!is.finite(cell_confound_eta2_threshold) || cell_confound_eta2_threshold <= 0) cell_confound_eta2_threshold <- 0.5
cell_confound_action <- if (!is.null(cell_deconv_cfg$confound_action)) as.character(cell_deconv_cfg$confound_action) else "warn"
cell_confound_action <- tolower(trimws(cell_confound_action))
if (!cell_confound_action %in% c("warn", "drop")) cell_confound_action <- "warn"
cell_adjustment_action <- "warn"
if (!is.null(config_settings$cell_adjustment$on_high_eta2)) {
  cell_adjustment_action <- as.character(config_settings$cell_adjustment$on_high_eta2)
}
if (nzchar(opt$cell_adjustment_on_high_eta2)) {
  cell_adjustment_action <- opt$cell_adjustment_on_high_eta2
}
cell_adjustment_action <- tolower(trimws(cell_adjustment_action))
if (!cell_adjustment_action %in% c("warn", "stop")) cell_adjustment_action <- "warn"

caf_cfg <- config_settings$caf
caf_enabled <- isTRUE(caf_cfg$enabled)
caf_target_fpr <- suppressWarnings(as.numeric(if (!is.null(caf_cfg$target_fpr)) caf_cfg$target_fpr else 0.05))
if (!is.finite(caf_target_fpr) || caf_target_fpr <= 0) caf_target_fpr <- 0.05

tier3_min_total_n <- suppressWarnings(as.integer(if (!is.null(config_settings$tier3_meta$min_total_n))
  config_settings$tier3_meta$min_total_n else 20))
if (!is.na(opt$tier3_min_total_n)) tier3_min_total_n <- opt$tier3_min_total_n
if (!is.finite(tier3_min_total_n) || tier3_min_total_n < 2) tier3_min_total_n <- 20
tier3_min_per_group_per_stratum <- suppressWarnings(as.integer(if (!is.null(config_settings$tier3_meta$min_per_group_per_stratum))
  config_settings$tier3_meta$min_per_group_per_stratum else 5))
if (!is.na(opt$tier3_min_per_group_per_stratum)) tier3_min_per_group_per_stratum <- opt$tier3_min_per_group_per_stratum
if (!is.finite(tier3_min_per_group_per_stratum) || tier3_min_per_group_per_stratum < 2) tier3_min_per_group_per_stratum <- 5
tier3_on_fail <- if (!is.null(config_settings$tier3_meta$on_fail)) as.character(config_settings$tier3_meta$on_fail) else "stop"
if (nzchar(opt$tier3_on_fail)) tier3_on_fail <- opt$tier3_on_fail
tier3_on_fail <- tolower(trimws(tier3_on_fail))
if (!tier3_on_fail %in% c("stop", "skip")) tier3_on_fail <- "stop"
if (beginner_safe) {
  message(sprintf("Beginner-safe mode enabled: min_total_size >= %d; |DeltaBeta| >= %.3f", MIN_TOTAL_SIZE_STOP, delta_beta_thresh))
}
dmr_min_cpgs <- config_settings$dmr$min_cpgs
if (!is.na(opt$dmr_min_cpgs)) dmr_min_cpgs <- opt$dmr_min_cpgs
if (!is.finite(dmr_min_cpgs) || dmr_min_cpgs < 1) dmr_min_cpgs <- 1
DMR_MIN_CPGS <- as.integer(dmr_min_cpgs)
dmr_maxgap <- if (!is.null(config_settings$dmr$maxgap)) config_settings$dmr$maxgap else DMR_MAXGAP
if (!is.na(opt$dmr_maxgap)) dmr_maxgap <- opt$dmr_maxgap
if (!is.finite(dmr_maxgap) || dmr_maxgap < 1) dmr_maxgap <- DMR_MAXGAP
DMR_MAXGAP <- as.integer(dmr_maxgap)
dmr_p_cutoff <- pval_thresh
if (!is.finite(dmr_p_cutoff) || dmr_p_cutoff <= 0) {
  dmr_p_cutoff <- DMR_P_CUTOFF_DEFAULT
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(out_dir)) stop(paste("Cannot create output directory:", out_dir))
results_root <- file.path(out_dir, "results")
results_dirs <- list(
  minfi = file.path(results_root, "minfi"),
  sesame = file.path(results_root, "sesame"),
  sesame_native = file.path(results_root, "sesame_native"),
  consensus = file.path(results_root, "consensus"),
  reports = file.path(results_root, "reports"),
  logs = file.path(results_root, "logs")
)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(results_dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

# --- Set ExperimentHub/AnnotationHub Cache ---
sesame_cache_dir <- file.path(project_dir, "cache")
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
  tryCatch({
    # Prefer explicit `text` tooltips when available (prevents hover picking a
    # layer without tooltip and makes overlapping points usable).
    has_text <- FALSE
    tryCatch({
      if (!is.null(p$mapping) && !is.null(p$mapping$text)) {
        has_text <- TRUE
      } else if (!is.null(p$layers) && length(p$layers) > 0) {
        for (ly in p$layers) {
          if (!is.null(ly$mapping) && !is.null(ly$mapping$text)) {
            has_text <- TRUE
            break
          }
        }
      }
    }, error = function(e) {
      has_text <<- FALSE
    })
    pp <- if (isTRUE(has_text)) ggplotly(p, tooltip = "text") else ggplotly(p)
    pp <- plotly::layout(pp, hovermode = "closest")
    saveWidget(pp, file = file.path(dir, filename), selfcontained = TRUE)
  }, error = function(e) {
    message(sprintf("  Interactive plot save skipped (%s): %s", filename, e$message))
  })
}

STATIC_PLOT_DPI <- 300
STATIC_PLOT_FONT_BASE <- 13
STATIC_PLOT_TITLE_SIZE <- 14
STATIC_PLOT_SUBTITLE_SIZE <- 12
STATIC_PLOT_AXIS_TEXT_SIZE <- 11
STATIC_PLOT_AXIS_TITLE_SIZE <- 12
STATIC_PLOT_LEGEND_TEXT_SIZE <- 11
STATIC_PLOT_LEGEND_TITLE_SIZE <- 12

apply_static_theme <- function(p) {
  p + theme(
    text = element_text(size = STATIC_PLOT_FONT_BASE),
    plot.title = element_text(size = STATIC_PLOT_TITLE_SIZE),
    plot.subtitle = element_text(size = STATIC_PLOT_SUBTITLE_SIZE),
    axis.text = element_text(size = STATIC_PLOT_AXIS_TEXT_SIZE),
    axis.title = element_text(size = STATIC_PLOT_AXIS_TITLE_SIZE),
    legend.text = element_text(size = STATIC_PLOT_LEGEND_TEXT_SIZE),
    legend.title = element_text(size = STATIC_PLOT_LEGEND_TITLE_SIZE)
  )
}

save_static_plot <- function(p, filename, dir, width = 7, height = 4, dpi = STATIC_PLOT_DPI) {
  p_use <- apply_static_theme(p)
  tryCatch({
    ggsave(filename = file.path(dir, filename), plot = p_use, width = width, height = height, dpi = dpi)
  }, error = function(e) {
    message(sprintf("  Static plot save skipped (%s): %s", filename, e$message))
  })
  # Also save PDF version for publication-ready figures
  pdf_name <- if (grepl("\\.png$", filename, ignore.case = TRUE)) {
    sub("\\.png$", ".pdf", filename, ignore.case = TRUE)
  } else {
    paste0(filename, ".pdf")
  }
  tryCatch({
    ggsave(filename = file.path(dir, pdf_name), plot = p_use, width = width, height = height, device = "pdf")
  }, error = function(e) {
    message(sprintf("  PDF plot save skipped (%s): %s", pdf_name, e$message))
  })
}

save_datatable <- function(df, filename, dir, sort_col = NULL) {
  # Identify P-value column index (0-based for JS)
  if (is.null(sort_col) || !(sort_col %in% colnames(df))) {
    sort_candidates <- c("P.Value", "p.value", "p.adjust", "adj.P.Val")
    sort_col <- sort_candidates[sort_candidates %in% colnames(df)][1]
  }
  p_idx <- if (!is.null(sort_col)) match(sort_col, colnames(df)) - 1 else 0
  if (length(p_idx) == 0 || is.na(p_idx)) p_idx <- 0 # Fallback
  
  dt <- datatable(df, rownames = FALSE, extensions = 'Buttons', options = list(
    dom = 'Bfrtip',
    pageLength = 25,
    buttons = c('csv', 'excel'),
    order = list(list(p_idx, 'asc'))
  ))
  saveWidget(dt, file = file.path(dir, filename), selfcontained = TRUE)
}

write_matrix_tsv_gz <- function(mat, path, chunk_size = 50000L) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(invisible(FALSE))
  ok <- tryCatch({
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    header <- c("CpG", colnames(mat))
    writeLines(paste(header, collapse = "\t"), con)
    n <- nrow(mat)
    for (i in seq(1, n, by = chunk_size)) {
      idx <- i:min(i + chunk_size - 1L, n)
      chunk <- mat[idx, , drop = FALSE]
      out <- cbind(CpG = rownames(chunk), chunk)
      utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    TRUE
  }, error = function(e) {
    message("  WARNING: Failed to write matrix file: ", path, " (", e$message, ")")
    FALSE
  })
  invisible(ok)
}

with_muffled_warnings <- function(expr, label = NULL, patterns = NULL) {
  seen <- character(0)
  handler <- function(w) {
    msg <- conditionMessage(w)
    if (is.null(patterns) || any(grepl(patterns, msg, ignore.case = TRUE))) {
      seen <<- unique(c(seen, msg))
      invokeRestart("muffleWarning")
    }
  }
  res <- withCallingHandlers(expr, warning = handler)
  if (!is.null(label) && length(seen) > 0) {
    if (any(grepl("boundary \\(singular\\) fit", seen, ignore.case = TRUE))) {
      message(sprintf("  %s: singular fit warning suppressed; consider reducing covariates.", label))
    } else if (length(seen) == 1) {
      message(sprintf("  %s: warning suppressed (%s).", label, seen[1]))
    } else {
      message(sprintf("  %s: %d warnings suppressed.", label, length(seen)))
    }
  }
  res
}

with_muffled_conditions <- function(expr, label = NULL, patterns = NULL) {
  seen <- character(0)
  warn_handler <- function(w) {
    msg <- conditionMessage(w)
    if (is.null(patterns) || any(grepl(patterns, msg, ignore.case = TRUE))) {
      seen <<- unique(c(seen, msg))
      invokeRestart("muffleWarning")
    }
  }
  msg_handler <- function(m) {
    msg <- conditionMessage(m)
    if (is.null(patterns) || any(grepl(patterns, msg, ignore.case = TRUE))) {
      seen <<- unique(c(seen, msg))
      invokeRestart("muffleMessage")
    }
  }
  res <- withCallingHandlers(expr, warning = warn_handler, message = msg_handler)
  if (!is.null(label) && length(seen) > 0) {
    if (any(grepl("boundary \\(singular\\) fit", seen, ignore.case = TRUE))) {
      message(sprintf("  %s: singular fit warning suppressed; consider reducing covariates.", label))
    } else if (length(seen) == 1) {
      message(sprintf("  %s: warning/message suppressed (%s).", label, seen[1]))
    } else {
      message(sprintf("  %s: %d warning/message(s) suppressed.", label, length(seen)))
    }
  }
  res
}

delta_beta_pass <- function(x) {
  if (!is.finite(delta_beta_thresh) || delta_beta_thresh <= 0) {
    return(rep(TRUE, length(x)))
  }
  is.finite(x) & abs(x) >= delta_beta_thresh
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

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

reduction_score <- function(baseline, current) {
  if (!is.finite(baseline) || baseline <= 0) return(0)
  clamp01((baseline - current) / baseline)
}

compute_pc_batch_r2 <- function(pca_scores, batch, max_pcs = 5) {
  if (is.null(pca_scores) || ncol(pca_scores) < 1) return(NA_real_)
  batch <- normalize_meta_vals(batch)
  if (length(unique(batch[!is.na(batch)])) < 2) return(NA_real_)
  n_pcs <- min(max_pcs, ncol(pca_scores))
  r2s <- numeric(0)
  for (k in seq_len(n_pcs)) {
    df <- data.frame(pc = pca_scores[, k], batch = batch)
    fit <- tryCatch(lm(pc ~ batch, data = df), error = function(e) NULL)
    if (is.null(fit)) next
    r2s <- c(r2s, summary(fit)$r.squared)
  }
  if (length(r2s) == 0) return(NA_real_)
  median(r2s, na.rm = TRUE)
}

compute_knn_mixing <- function(pca_scores, batch, k = 10, max_samples = 500) {
  if (is.null(pca_scores) || nrow(pca_scores) < (k + 2)) return(NA_real_)
  n <- nrow(pca_scores)
  idx <- seq_len(n)
  if (n > max_samples) {
    set.seed(seed_value)
    idx <- sort(sample(idx, max_samples))
  }
  pcs <- pca_scores[idx, , drop = FALSE]
  batch_sub <- batch[idx]
  d <- as.matrix(dist(pcs))
  mix_scores <- numeric(length(idx))
  for (i in seq_along(idx)) {
    ord <- order(d[i, ], na.last = NA)
    ord <- ord[ord != i]
    nn <- head(ord, k)
    same <- sum(batch_sub[nn] == batch_sub[i])
    mix_scores[i] <- 1 - (same / k)
  }
  mean(mix_scores, na.rm = TRUE)
}

detect_bio_controls <- function(targets) {
  cols <- colnames(targets)
  cols <- setdiff(cols, c("primary_group", "Basename", "filenames", "SampleID"))
  sex_cols <- cols[grepl("sex|gender", tolower(cols))]
  age_cols <- cols[grepl("age", tolower(cols))]
  cell_cols <- cols[grepl("^cell_|cell", tolower(cols))]
  unique(c(sex_cols, age_cols, cell_cols))
}

compute_genomic_lambda <- function(pvals) {
  pvals <- suppressWarnings(as.numeric(pvals))
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals < 1]
  if (length(pvals) < 2) return(NA_real_)
  chisq_vals <- suppressWarnings(qchisq(1 - pvals, 1))
  chisq_vals <- chisq_vals[is.finite(chisq_vals)]
  if (length(chisq_vals) < 2) return(NA_real_)
  median(chisq_vals, na.rm = TRUE) / qchisq(0.5, 1)
}

compute_bio_preservation <- function(pca_before, pca_after, targets, bio_vars, max_pcs = 5) {
  if (is.null(pca_before) || is.null(pca_after)) return(list(score = 1, detail = NULL))
  if (length(bio_vars) == 0) return(list(score = 1, detail = NULL))
  n_pcs <- min(max_pcs, ncol(pca_before$x), ncol(pca_after$x))
  scores <- numeric(0)
  detail <- data.frame(Variable = character(0), Before = numeric(0), After = numeric(0), Ratio = numeric(0))
  for (var in bio_vars) {
    if (!var %in% colnames(targets)) next
    vals <- normalize_meta_vals(targets[[var]])
    non_missing <- vals[!is.na(vals)]
    if (length(unique(non_missing)) < 2) next
    before_vals <- numeric(0)
    after_vals <- numeric(0)
    for (k in seq_len(n_pcs)) {
      if (is.numeric(vals) || is.integer(vals)) {
        cc <- complete.cases(pca_before$x[, k], vals)
        if (sum(cc) >= 3) {
          fit_b <- tryCatch(lm(pca_before$x[cc, k] ~ vals[cc]), error = function(e) NULL)
          if (!is.null(fit_b)) before_vals <- c(before_vals, summary(fit_b)$r.squared)
        }
        cc2 <- complete.cases(pca_after$x[, k], vals)
        if (sum(cc2) >= 3) {
          fit_a <- tryCatch(lm(pca_after$x[cc2, k] ~ vals[cc2]), error = function(e) NULL)
          if (!is.null(fit_a)) after_vals <- c(after_vals, summary(fit_a)$r.squared)
        }
      } else {
        cc <- complete.cases(pca_before$x[, k], vals)
        if (sum(cc) >= 3 && length(unique(vals[cc])) >= 2) {
          fit_b <- tryCatch(aov(pca_before$x[cc, k] ~ as.factor(vals[cc])), error = function(e) NULL)
          if (!is.null(fit_b)) {
            ss <- summary(fit_b)[[1]]
            before_vals <- c(before_vals, ss$`Sum Sq`[1] / sum(ss$`Sum Sq`))
          }
        }
        cc2 <- complete.cases(pca_after$x[, k], vals)
        if (sum(cc2) >= 3 && length(unique(vals[cc2])) >= 2) {
          fit_a <- tryCatch(aov(pca_after$x[cc2, k] ~ as.factor(vals[cc2])), error = function(e) NULL)
          if (!is.null(fit_a)) {
            ss <- summary(fit_a)[[1]]
            after_vals <- c(after_vals, ss$`Sum Sq`[1] / sum(ss$`Sum Sq`))
          }
        }
      }
    }
    if (length(before_vals) == 0 || length(after_vals) == 0) next
    before_med <- median(before_vals, na.rm = TRUE)
    after_med <- median(after_vals, na.rm = TRUE)
    if (!is.finite(before_med) || before_med <= 0) next
    ratio <- clamp01(after_med / before_med)
    detail <- rbind(detail, data.frame(Variable = var, Before = before_med, After = after_med, Ratio = ratio))
    scores <- c(scores, ratio)
  }
  if (length(scores) == 0) return(list(score = 1, detail = detail))
  list(score = mean(scores, na.rm = TRUE), detail = detail)
}

compute_batch_stability <- function(betas, targets, batch_col, covariates, group_col, max_probes = 5000) {
  if (is.null(batch_col) || !(batch_col %in% colnames(targets))) return(NA_real_)
  overlap_batches <- identify_overlap_batches(targets, batch_col, group_col)
  if (length(overlap_batches) < 2) return(NA_real_)
  vars <- apply(betas, 1, var, na.rm = TRUE)
  top_idx <- head(order(vars, decreasing = TRUE), min(max_probes, nrow(betas)))
  eff_list <- list()
  for (b in overlap_batches) {
    idx <- targets[[batch_col]] == b
    sub_targets <- targets[idx, , drop = FALSE]
    sub_betas <- betas[top_idx, idx, drop = FALSE]
    covars <- intersect(covariates, colnames(sub_targets))
    if (length(covars) > 0) {
      cov_drop <- drop_single_level_covariates(sub_targets, covars)
      covars <- cov_drop$keep
    }
    formula_str <- "~ 0 + primary_group"
    if (length(covars) > 0) formula_str <- paste(formula_str, "+", paste(covars, collapse = " + "))
    design <- tryCatch(model.matrix(as.formula(formula_str), data = sub_targets), error = function(e) NULL)
    if (is.null(design)) next
    colnames(design) <- gsub("primary_group", "", colnames(design))
    colnames(design) <- make.names(colnames(design), unique = TRUE)
    group_cols <- make.names(levels(sub_targets$primary_group))
    group_cols <- intersect(group_cols, colnames(design))
    design_ld <- drop_linear_dependencies(design, group_cols = group_cols)
    design <- design_ld$mat
    if (ncol(design) < 2) next
    m_vals <- logit_offset(sub_betas)
    fit <- tryCatch(lmFit(m_vals, design), error = function(e) NULL)
    if (is.null(fit)) next
    if (length(group_cols) < 2) next
    contrast_str <- paste0(group_cols[2], "-", group_cols[1])
    cm <- tryCatch(makeContrasts(contrasts = contrast_str, levels = design), error = function(e) NULL)
    if (is.null(cm)) next
    fit2 <- tryCatch(contrasts.fit(fit, cm), error = function(e) NULL)
    if (is.null(fit2)) next
    fit2 <- safe_ebayes(fit2, "batch_stability")
    if (is.null(fit2)) next
    eff_list[[b]] <- fit2$coefficients[, 1]
  }
  if (length(eff_list) < 2) return(NA_real_)
  batches <- names(eff_list)
  cors <- numeric(0)
  for (i in seq_len(length(batches) - 1)) {
    for (j in (i + 1):length(batches)) {
      a <- eff_list[[batches[i]]]
      b <- eff_list[[batches[j]]]
      cors <- c(cors, suppressWarnings(cor(a, b, use = "pairwise.complete.obs")))
    }
  }
  if (length(cors) == 0) return(NA_real_)
  clamp01((median(cors, na.rm = TRUE) + 1) / 2)
}

#' Run Permutation-Based Null Calibration Test
#'
#' Performs label permutation to assess type I error calibration. Shuffles
#' group labels within batches (if batch_col provided) to maintain batch
#' structure, then runs differential analysis to generate null p-value
#' distributions.
#'
#' @param betas Numeric matrix of beta values (CpGs x samples)
#' @param targets Data frame with sample metadata
#' @param group_col Character; column name for group labels
#' @param covariates Character vector of covariate column names
#' @param batch_col Character; batch column for stratified permutation (or NULL)
#' @param perm_n Integer; number of permutations to run
#' @param vp_top Integer; number of top-variable CpGs to use
#' @return Data frame with columns: run, ks_p (KS test p-value vs uniform),
#'   lambda (genomic inflation factor), n_sig (number of significant at 0.05)
#' @details
#' KS p-value < 0.05 suggests deviation from uniform null distribution.
#' Lambda near 1.0 indicates well-calibrated test statistics.
run_permutation_uniformity <- function(betas, targets, group_col, covariates, batch_col, perm_n, vp_top) {
  if (perm_n <= 0) return(NULL)
  top_perm <- head(order(apply(betas, 1, var), decreasing = TRUE), min(vp_top, nrow(betas)))
  m_perm <- logit_offset(betas[top_perm, , drop = FALSE])
  perm_results <- data.frame(run = integer(0), ks_p = numeric(0), lambda = numeric(0))
  perm_group_cols <- make.names(levels(targets[[group_col]]))
  for (i in seq_len(perm_n)) {
    perm_labels <- targets[[group_col]]
    if (!is.null(batch_col) && batch_col %in% colnames(targets)) {
      idx_list <- split(seq_len(nrow(targets)), targets[[batch_col]])
      perm_labels <- perm_labels
      for (idx in idx_list) {
        perm_labels[idx] <- sample(perm_labels[idx])
      }
    } else {
      perm_labels <- sample(perm_labels)
    }
    perm_df <- data.frame(perm_labels = perm_labels)
    perm_design <- model.matrix(~ 0 + perm_labels, data = perm_df)
    colnames(perm_design) <- gsub("perm_labels", "", colnames(perm_design))
    if (length(covariates) > 0) {
      cov_formula <- as.formula(paste("~ 0 +", paste(covariates, collapse = " + ")))
      cov_mat <- model.matrix(cov_formula, data = targets)
      perm_design <- cbind(perm_design, cov_mat)
    }
    perm_ld <- drop_linear_dependencies(perm_design, group_cols = perm_group_cols)
    perm_design <- perm_ld$mat
    # Find group columns by name rather than assuming fixed positions
    grp_idx <- match(perm_group_cols, colnames(perm_design))
    grp_idx <- grp_idx[!is.na(grp_idx)]
    if (length(grp_idx) < 2) next
    cont_vec <- rep(0, ncol(perm_design))
    cont_vec[grp_idx[2]] <- 1; cont_vec[grp_idx[1]] <- -1
    cm_perm <- matrix(cont_vec, ncol = 1)
    rownames(cm_perm) <- colnames(perm_design)
    fitp <- lmFit(m_perm, perm_design)
    fitp2 <- contrasts.fit(fitp, cm_perm)
    fitp2 <- safe_ebayes(fitp2, "perm_uniformity")
    if (is.null(fitp2)) next
    resp <- topTable(fitp2, coef = 1, number = Inf, adjust.method = "BH")
    pvals <- suppressWarnings(as.numeric(resp$P.Value))
    pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
    if (length(pvals) < 2) next
    ks_p <- tryCatch(ks.test(pvals, "punif")$p.value, error = function(e) NA_real_)
    lambda_val <- compute_genomic_lambda(pvals)
    perm_results <- rbind(perm_results, data.frame(run = i, ks_p = ks_p, lambda = lambda_val))
  }
  perm_results
}

#' Select Optimal Batch Correction Strategy
#'
#' Evaluates multiple batch correction methods (none, ComBat, limma, SVA) and
#' covariate combinations to select the optimal strategy based on a weighted
#' scoring system (Correction Adequacy Framework, CAF).
#'
#' @param betas Numeric matrix of beta values (CpGs x samples)
#' @param targets Data frame with sample metadata
#' @param batch_col Character; column name for batch factor
#' @param batch_tier Integer; confounding tier (0-3, where 3 = non-identifiable)
#' @param covariate_sets Named list of covariate vectors to evaluate
#' @param group_col Character; column name for group variable
#' @param config_settings List; configuration from config.yaml
#' @param scoring_preset List; weights and guards from scoring_presets
#' @param perm_n Integer; number of permutations for calibration
#' @param vp_top Integer; number of top-variable CpGs for evaluation
#' @param prefix Character; output file prefix
#' @param out_dir Character; output directory path
#' @param method_pool Character vector; methods to evaluate (default: all)
#' @return List with components:
#'   - best_method: Selected correction method
#'   - best_covariates: Selected covariate set
#'   - candidates: Data frame of all evaluated candidates with scores
#' @details
#' Scoring weights (from config.yaml scoring_presets):
#'   - batch: Batch effect reduction (median variance proportion)
#'   - bio: Biological signal preservation (R2 with bio variables)
#'   - cal: Calibration (permutation uniformity)
#'   - stab: Stability (cross-batch effect correlation)
select_batch_strategy <- function(betas, targets, batch_col, batch_tier, covariate_sets, group_col,
                                  config_settings, scoring_preset, perm_n, vp_top, prefix, out_dir,
                                  method_pool = NULL) {
  methods <- if (!is.null(method_pool) && length(method_pool) > 0) method_pool else c("none", "combat", "limma", "sva")
  if (is.null(batch_col) || nrow(targets) < 4) {
    return(list(best_method = "none", best_covariates = covariate_sets[[1]], candidates = NULL))
  }
  if (!is.null(batch_tier) && batch_tier == 3) {
    return(list(best_method = "none", best_covariates = covariate_sets[[1]], candidates = NULL))
  }
  batch_vals <- targets[[batch_col]]
  M_mat <- logit_offset(betas)
  bio_vars <- detect_bio_controls(targets)
  cand_rows <- data.frame()
  for (set_name in names(covariate_sets)) {
    cov_set <- covariate_sets[[set_name]]
    base_res <- eval_batch_method(M_mat, targets, group_col = group_col, batch_col = batch_col, covariates = cov_set, method = "none")
    base_pca <- safe_prcomp(t(base_res$M_corr), label = paste("Batch eval", set_name, "baseline"), scale. = TRUE)
    base_pc_r2 <- if (!is.null(base_pca)) compute_pc_batch_r2(base_pca$x, batch_vals) else NA_real_
    base_mix <- if (!is.null(base_pca)) compute_knn_mixing(base_pca$x, batch_vals) else NA_real_
    for (m in methods) {
      if (!is.null(batch_tier) && batch_tier == 2 && m == "limma") {
        # keep limma but mark as restricted in notes
      }
      res <- eval_batch_method(M_mat, targets, group_col = group_col, batch_col = batch_col, covariates = cov_set, method = m)
      if (isTRUE(res$failed)) {
        note_msg <- res$fail_reason
        if (!is.null(res$note) && nzchar(res$note)) {
          note_msg <- paste(note_msg, res$note, sep = " | ")
        }
        cand_rows <- rbind(cand_rows, data.frame(
          cov_set = set_name,
          method = m,
          batch_score = NA_real_,
          batch_reduction = NA_real_,
          pc_r2_reduction = NA_real_,
          mix_improve = NA_real_,
          bio_score = NA_real_,
          cal_score = NA_real_,
          stab_score = NA_real_,
          total_score = -Inf,
          status = "FAILED",
          note = note_msg,
          stringsAsFactors = FALSE
        ))
        next
      }
      pca_after <- safe_prcomp(t(res$M_corr), label = paste("Batch eval", set_name, m), scale. = TRUE)
      pc_r2_after <- if (!is.null(pca_after)) compute_pc_batch_r2(pca_after$x, batch_vals) else NA_real_
      mix_after <- if (!is.null(pca_after)) compute_knn_mixing(pca_after$x, batch_vals) else NA_real_
      bio_pres <- compute_bio_preservation(base_pca, pca_after, targets, bio_vars)
      batch_red <- mean(c(
        reduction_score(base_res$batch_var, res$batch_var),
        reduction_score(base_res$prop_batch_sig, res$prop_batch_sig),
        reduction_score(base_res$n_pc_batch_sig, res$n_pc_batch_sig)
      ), na.rm = TRUE)
      pc_r2_red <- reduction_score(base_pc_r2, pc_r2_after)
      mix_score <- NA_real_
      if (is.finite(base_mix) && is.finite(mix_after)) {
        mix_score <- clamp01((mix_after - base_mix) / max(1e-6, 1 - base_mix))
      }
      batch_score <- mean(c(batch_red, pc_r2_red, mix_score), na.rm = TRUE)
      bio_score <- bio_pres$score
      cal_score <- NA_real_
      stab_score <- NA_real_
      total_pre <- scoring_preset$weights$batch * ifelse(is.finite(batch_score), batch_score, 0.5) +
        scoring_preset$weights$bio * ifelse(is.finite(bio_score), bio_score, 0.5) +
        scoring_preset$weights$cal * 0.5 +
        scoring_preset$weights$stab * 0.5
      cand_rows <- rbind(cand_rows, data.frame(
        cov_set = set_name,
        method = m,
        batch_score = batch_score,
        batch_reduction = batch_red,
        pc_r2_reduction = pc_r2_red,
        mix_improve = mix_score,
        bio_score = bio_score,
        cal_score = cal_score,
        stab_score = stab_score,
        total_score = total_pre,
        status = "PENDING",
        note = if (!is.null(res$note)) res$note else "",
        stringsAsFactors = FALSE
      ))
    }
  }
  if (nrow(cand_rows) == 0) {
    return(list(best_method = "none", best_covariates = covariate_sets[[1]], candidates = NULL))
  }
  top_k <- scoring_preset$top_k
  cand_rows <- cand_rows[order(-cand_rows$total_score), ]
  eval_idx <- seq_len(min(top_k, nrow(cand_rows)))
  for (i in eval_idx) {
    row <- cand_rows[i, ]
    if (row$status == "FAILED") next
    cov_set <- covariate_sets[[row$cov_set]]
    formula_str <- "~ 0 + primary_group"
    if (length(cov_set) > 0) formula_str <- paste(formula_str, "+", paste(cov_set, collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = targets)
    colnames(design) <- gsub("primary_group", "", colnames(design))
    group_cols <- make.names(levels(targets[[group_col]]))
    design_ld <- drop_linear_dependencies(design, group_cols = group_cols)
    design <- design_ld$mat
    betas_corr <- betas
    if (row$method %in% c("combat", "limma")) {
      bc_res <- apply_batch_correction(betas, row$method, batch_col, cov_set, targets, design)
      betas_corr <- bc_res$betas
    }
    stab <- compute_batch_stability(betas_corr, targets, batch_col, cov_set, group_col)
    perm_res <- run_permutation_uniformity(betas_corr, targets, group_col, cov_set, batch_col, perm_n, vp_top)
    ks_p_med <- if (!is.null(perm_res) && nrow(perm_res) > 0) median(perm_res$ks_p, na.rm = TRUE) else NA_real_
    cal_score <- ifelse(is.finite(ks_p_med), clamp01(ks_p_med / config_settings$calibration$ks_p), NA_real_)
    total_score <- scoring_preset$weights$batch * ifelse(is.finite(row$batch_score), row$batch_score, 0.5) +
      scoring_preset$weights$bio * ifelse(is.finite(row$bio_score), row$bio_score, 0.5) +
      scoring_preset$weights$cal * ifelse(is.finite(cal_score), cal_score, 0.5) +
      scoring_preset$weights$stab * ifelse(is.finite(stab), stab, 0.5)
    guard_batch <- is.finite(row$batch_score) && row$batch_score >= scoring_preset$guards$batch_reduction_min
    guard_bio <- is.finite(row$bio_score) && row$bio_score >= scoring_preset$guards$bio_min
    guard_cal <- is.finite(ks_p_med) && ks_p_med >= config_settings$calibration$ks_p
    if (guard_batch && guard_bio && guard_cal) {
      borderline <- 0
      if (row$batch_score < (scoring_preset$guards$batch_reduction_min + 0.05)) borderline <- borderline + 1
      if (row$bio_score < (scoring_preset$guards$bio_min + 0.05)) borderline <- borderline + 1
      if (ks_p_med < (config_settings$calibration$ks_p * 1.5)) borderline <- borderline + 1
      status <- if (borderline >= 2) "YELLOW" else "GREEN"
    } else {
      status <- "RED"
    }
    cand_rows$stab_score[i] <- stab
    cand_rows$cal_score[i] <- cal_score
    cand_rows$total_score[i] <- total_score
    cand_rows$status[i] <- status
    if (!is.null(perm_res) && nrow(perm_res) > 0) {
      write.csv(perm_res, file.path(out_dir, paste0(prefix, "_PermUniformity_", row$cov_set, "_", row$method, ".csv")),
                row.names = FALSE)
    }
  }
  cand_rows <- cand_rows[order(-cand_rows$total_score), ]
  best_rows <- cand_rows[!(cand_rows$status %in% c("RED", "FAILED")), , drop = FALSE]
  if (nrow(best_rows) == 0) {
    best_rows <- cand_rows[cand_rows$status != "FAILED", , drop = FALSE]
  }
  if (nrow(best_rows) == 0) {
    return(list(best_method = "none", best_covariates = covariate_sets[[1]], candidates = cand_rows))
  }
  best <- best_rows[1, ]
  list(best_method = best$method, best_covariates = covariate_sets[[best$cov_set]], candidates = cand_rows)
}

# M-value transformation with a small offset to avoid log(0)
# Following Du et al. (2010) BMC Bioinformatics approach
logit_offset <- function(betas, offset = LOGIT_OFFSET) {
  log2((betas + offset) / (1 - betas + offset))
}

load_marker_list <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  sep <- ""
  if (grepl("\t", first_line, fixed = TRUE)) sep <- "\t"
  if (sep == "" && grepl(",", first_line, fixed = TRUE)) sep <- ","
  if (sep == "") {
    vals <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
    vals <- trimws(vals)
    vals <- vals[vals != ""]
    return(unique(vals))
  }
  df <- tryCatch(read.delim(path, sep = sep, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  col <- names(df)[1]
  for (cand in c("CpG", "cpg", "probe", "Probe", "marker", "Marker")) {
    if (cand %in% names(df)) {
      col <- cand
      break
    }
  }
  vals <- trimws(as.character(df[[col]]))
  vals <- vals[!grepl("^#", vals)]
  vals <- vals[vals != ""]
  unique(vals)
}

load_probe_list_from_file <- function(path, col_candidates = character(0)) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  sep <- ""
  if (grepl("\t", first_line, fixed = TRUE)) sep <- "\t"
  if (sep == "" && grepl(",", first_line, fixed = TRUE)) sep <- ","
  if (sep == "") {
    vals <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
    vals <- trimws(vals)
    vals <- vals[!grepl("^#", vals)]
    vals <- vals[vals != ""]
    return(unique(vals))
  }
  df <- tryCatch(read.delim(path, sep = sep, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  col <- names(df)[1]
  if (length(col_candidates) > 0) {
    for (cand in col_candidates) {
      if (cand %in% names(df)) {
        col <- cand
        break
      }
    }
  }
  vals <- trimws(as.character(df[[col]]))
  vals <- vals[vals != ""]
  unique(vals)
}

looks_like_cpg_ids <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return(FALSE)
  mean(grepl("^(cg|ch)", ids, ignore.case = TRUE)) >= 0.5
}

extract_cpgs_from_object <- function(obj, col_candidates = character(0)) {
  if (is.null(obj) || is.function(obj)) return(NULL)
  if (is.factor(obj)) obj <- as.character(obj)
  if (is.character(obj)) return(unique(obj))
  if (is.data.frame(obj) || is.matrix(obj)) {
    if (!is.null(colnames(obj)) && length(col_candidates) > 0) {
      for (cand in col_candidates) {
        if (cand %in% colnames(obj)) {
          vals <- trimws(as.character(obj[, cand]))
          vals <- vals[vals != ""]
          return(unique(vals))
        }
      }
    }
    rn <- rownames(obj)
    if (!is.null(rn) && looks_like_cpg_ids(rn)) return(unique(rn))
  }
  if (is.list(obj)) {
    for (nm in c("cpgs", "CpG", "probes", "probe_ids", "Probe", "probe", "targets", "TargetID")) {
      if (!is.null(obj[[nm]])) {
        out <- extract_cpgs_from_object(obj[[nm]], col_candidates = col_candidates)
        if (!is.null(out) && length(out) > 0) return(out)
      }
    }
  }
  NULL
}

normalize_array_key <- function(array_type) {
  if (is.na(array_type) || !nzchar(array_type)) return("")
  gsub("[^a-z0-9]", "", tolower(array_type))
}

normalize_array_label <- function(array_type) {
  key <- normalize_array_key(array_type)
  if (grepl("epicv2", key) || grepl("950", key)) return("EPICv2")
  if (grepl("epic", key) || grepl("850", key)) return("EPIC")
  if (grepl("450", key)) return("450K")
  ""
}

select_local_cross_reactive_file <- function(array_type, local_dir) {
  if (!nzchar(local_dir)) return("")
  if (!dir.exists(local_dir)) return("")
  key <- normalize_array_key(array_type)
  candidates <- c()
  if (grepl("epicv2", key) || grepl("950", key)) {
    candidates <- c("EPICv2", "epicv2")
  } else if (grepl("epic", key) || grepl("850", key)) {
    candidates <- c("EPIC", "EPICv1", "epic")
  } else if (grepl("450", key)) {
    candidates <- c("450K", "450k", "HM450", "hm450")
  }
  prefixes <- c("cross_reactive", "xreactive", "xreactive_probes")
  exts <- c("tsv", "csv", "txt")
  for (cand in candidates) {
    for (pref in prefixes) {
      for (ext in exts) {
        path <- file.path(local_dir, paste0(pref, "_", cand, ".", ext))
        if (file.exists(path)) return(path)
      }
    }
  }
  ""
}

get_cross_reactive_from_maxprobes <- function(array_type) {
  if (!requireNamespace("maxprobes", quietly = TRUE)) return(NULL)
  array_label <- normalize_array_label(array_type)
  if (!nzchar(array_label)) return(NULL)
  res <- tryCatch(maxprobes::xreactive_probes(array_type = array_label), error = function(e) NULL)
  if (is.null(res) || length(res) == 0) return(NULL)
  list(probes = unique(as.character(res)), source = paste0("maxprobes::xreactive_probes(", array_label, ")"),
       version = as.character(utils::packageVersion("maxprobes")))
}

resolve_cross_reactive_probes <- function(path, use_maxprobes = TRUE, array_type = NA_character_,
                                          local_dir = "", allow_missing = FALSE) {
  col_candidates <- c("CpG", "cpg", "probe", "Probe", "probe_id", "Probe_ID",
                      "TargetID", "IlmnID", "Name", "ID")
  if (nzchar(path) && file.exists(path)) {
    vals <- load_probe_list_from_file(path, col_candidates = col_candidates)
    if (!is.null(vals) && length(vals) > 0) {
      return(list(probes = unique(vals), source = path, version = "local"))
    }
  }
  if (isTRUE(use_maxprobes)) {
    max_res <- get_cross_reactive_from_maxprobes(array_type)
    if (!is.null(max_res) && !is.null(max_res$probes) && length(max_res$probes) > 0) {
      return(max_res)
    }
  }
  local_path <- select_local_cross_reactive_file(array_type, local_dir)
  if (nzchar(local_path) && file.exists(local_path)) {
    vals <- load_probe_list_from_file(local_path, col_candidates = col_candidates)
    if (!is.null(vals) && length(vals) > 0) {
      return(list(probes = unique(vals), source = local_path, version = "local"))
    }
  }
  if (allow_missing) return(NULL)
  NULL
}

# Backward-compatible alias
load_cross_reactive_list <- function(path, use_maxprobes = TRUE, array_type = NA_character_, local_dir = "") {
  resolve_cross_reactive_probes(path = path, use_maxprobes = use_maxprobes,
                                array_type = array_type, local_dir = local_dir, allow_missing = TRUE)
}

apply_cross_reactive_to_gmSet <- function(gmSet, cross_reactive_probes, use_maxprobes = TRUE) {
  removed <- character(0)
  source <- ""
  before_ids <- rownames(gmSet)
  if (isTRUE(use_maxprobes) && requireNamespace("maxprobes", quietly = TRUE) &&
      exists("dropXreactiveLoci", envir = asNamespace("maxprobes"), inherits = FALSE)) {
    gm_try <- tryCatch(maxprobes::dropXreactiveLoci(gmSet), error = function(e) NULL)
    if (!is.null(gm_try)) {
      gmSet <- gm_try
      removed <- setdiff(before_ids, rownames(gmSet))
      source <- "maxprobes::dropXreactiveLoci"
      return(list(gmSet = gmSet, removed = removed, source = source))
    }
  }
  if (!is.null(cross_reactive_probes) && length(cross_reactive_probes) > 0) {
    keep <- !(rownames(gmSet) %in% cross_reactive_probes)
    removed <- rownames(gmSet)[!keep]
    gmSet <- gmSet[keep, ]
    source <- "list"
  }
  list(gmSet = gmSet, removed = removed, source = source)
}

apply_cross_reactive_to_matrix <- function(mat, cross_reactive_probes) {
  if (is.null(mat) || is.null(cross_reactive_probes) || length(cross_reactive_probes) == 0) {
    return(list(mat = mat, removed = character(0)))
  }
  keep <- !(rownames(mat) %in% cross_reactive_probes)
  removed <- rownames(mat)[!keep]
  mat <- mat[keep, , drop = FALSE]
  list(mat = mat, removed = removed)
}

normalize_sex_vector <- function(x) {
  if (is.null(x)) return(character(0))
  v <- as.character(x)
  v <- trimws(tolower(v))
  v_clean <- gsub("[^a-z0-9]", "", v)
  out <- rep(NA_character_, length(v_clean))
  if (length(v_clean) > 0) {
    num_vals <- suppressWarnings(as.numeric(v_clean))
    if (sum(!is.na(num_vals)) == sum(v_clean != "" & !is.na(v_clean))) {
      uniq <- sort(unique(num_vals[is.finite(num_vals)]))
      if (length(uniq) > 0 && all(uniq %in% c(0, 1))) {
        out[num_vals == 1] <- "M"
        out[num_vals == 0] <- "F"
        return(out)
      }
      if (length(uniq) > 0 && all(uniq %in% c(1, 2))) {
        out[num_vals == 1] <- "M"
        out[num_vals == 2] <- "F"
        return(out)
      }
    }
  }
  male_tokens <- c("m", "male", "man", "boy", "xy", "males")
  female_tokens <- c("f", "female", "woman", "girl", "xx", "females")
  out[v_clean %in% male_tokens] <- "M"
  out[v_clean %in% female_tokens] <- "F"
  out
}

resolve_sex_column <- function(targets, preferred = "") {
  if (nzchar(preferred)) {
    if (preferred %in% colnames(targets)) return(preferred)
    message(sprintf("  Sex check: specified column '%s' not found; falling back to auto-detect.", preferred))
  }
  cols <- colnames(targets)
  hits <- cols[grepl("sex|gender", tolower(cols))]
  if (length(hits) == 0) return("")
  hits[1]
}

extract_predicted_sex <- function(sex_obj) {
  if (is.null(sex_obj)) return(NULL)
  if (is.data.frame(sex_obj) || inherits(sex_obj, "DataFrame")) {
    for (col in c("predictedSex", "sex", "Sex", "PredictedSex", "predicted_sex")) {
      if (col %in% colnames(sex_obj)) {
        return(as.character(sex_obj[[col]]))
      }
    }
    return(NULL)
  }
  if (is.character(sex_obj) || is.factor(sex_obj)) return(as.character(sex_obj))
  NULL
}

#' Detect Sex Mismatches Between Metadata and Methylation Prediction
#'
#' Compares reported biological sex from metadata against predicted sex
#' from X/Y chromosome probe intensities using minfi::getSex(). Generates
#' detailed reports and diagnostic plots for quality control.
#'
#' @param rgSet RGChannelSet object containing raw IDAT data
#' @param targets Data frame with sample metadata including sex/gender column
#' @param gsm_col Character string specifying sample ID column in targets
#' @param out_dir Output directory for reports and plots
#' @param action Action on mismatch: "stop" (error), "warn", or "flag"
#' @param column Optional specific column name for sex metadata
#'
#' @return List containing:
#'   - mismatch_samples: Vector of sample IDs with sex mismatches
#'   - mismatch_count: Number of mismatched samples
#'   - column: Name of sex column used
#'   - skipped: TRUE if check was skipped (with reason)
#'
#' @details
#' Sex prediction uses X and Y chromosome median intensities (xMed, yMed).
#' Metadata values are normalized (M/Male/1 -> M, F/Female/0 -> F).
#' Outputs Sex_Check_Summary.csv with detailed per-sample results.
#'
#' @seealso minfi::getSex for prediction algorithm details
run_sex_mismatch_check <- function(rgSet, targets, gsm_col, out_dir, action = "stop", column = "") {
  if (is.null(rgSet) || ncol(rgSet) == 0) return(list(skipped = TRUE, reason = "no_data"))
  sex_col <- resolve_sex_column(targets, column)
  if (!nzchar(sex_col)) {
    message("  Sex check skipped: no metadata sex/gender column found.")
    return(list(skipped = TRUE, reason = "no_column"))
  }
  sex_obj <- tryCatch(minfi::getSex(rgSet), error = function(e) {
    message("  Sex check failed: minfi::getSex() error: ", e$message)
    NULL
  })
  if (is.null(sex_obj)) return(list(skipped = TRUE, reason = "getsex_failed", column = sex_col))
  pred_raw <- extract_predicted_sex(sex_obj)
  if (is.null(pred_raw)) {
    message("  Sex check skipped: predicted sex column not found in getSex() output.")
    return(list(skipped = TRUE, reason = "predicted_missing", column = sex_col))
  }
  sample_ids <- colnames(rgSet)
  pred_vec <- rep(NA_character_, length(sample_ids))
  names(pred_vec) <- sample_ids
  if (!is.null(rownames(sex_obj)) && any(rownames(sex_obj) %in% sample_ids)) {
    pred_vec[rownames(sex_obj)] <- pred_raw
  } else if (!is.null(names(pred_raw)) && any(names(pred_raw) %in% sample_ids)) {
    pred_vec[names(pred_raw)] <- pred_raw
  } else if (length(pred_raw) == length(sample_ids)) {
    pred_vec <- pred_raw
    names(pred_vec) <- sample_ids
  }
  xMed <- NULL
  yMed <- NULL
  if (is.data.frame(sex_obj) || inherits(sex_obj, "DataFrame")) {
    if ("xMed" %in% colnames(sex_obj)) xMed <- as.numeric(sex_obj$xMed)
    if ("yMed" %in% colnames(sex_obj)) yMed <- as.numeric(sex_obj$yMed)
    if (!is.null(rownames(sex_obj)) && any(rownames(sex_obj) %in% sample_ids)) {
      xMed <- xMed[match(sample_ids, rownames(sex_obj))]
      yMed <- yMed[match(sample_ids, rownames(sex_obj))]
    }
  }
  meta_idx <- match(sample_ids, rownames(targets))
  if (all(is.na(meta_idx)) && gsm_col %in% colnames(targets)) {
    meta_idx <- match(sample_ids, targets[[gsm_col]])
  }
  meta_vals <- targets[[sex_col]][meta_idx]
  meta_norm <- normalize_sex_vector(meta_vals)
  pred_norm <- normalize_sex_vector(pred_vec)
  mismatch <- !is.na(meta_norm) & !is.na(pred_norm) & meta_norm != pred_norm
  mismatch_samples <- sample_ids[mismatch]
  status <- ifelse(is.na(meta_norm) | is.na(pred_norm), "ambiguous", ifelse(mismatch, "mismatch", "match"))
  xMed_vec <- if (is.null(xMed)) rep(NA_real_, length(sample_ids)) else xMed
  yMed_vec <- if (is.null(yMed)) rep(NA_real_, length(sample_ids)) else yMed
  out_df <- data.frame(
    SampleID = sample_ids,
    Sex_Metadata = meta_vals,
    Sex_Metadata_Norm = meta_norm,
    Sex_Predicted = pred_vec,
    Sex_Predicted_Norm = pred_norm,
    xMed = xMed_vec,
    yMed = yMed_vec,
    Status = status,
    stringsAsFactors = FALSE
  )
  write.csv(out_df, file.path(out_dir, "Sex_Check_Summary.csv"), row.names = FALSE)
  write.csv(out_df, file.path(out_dir, "Sex_Check.csv"), row.names = FALSE)

  if (!is.null(xMed) && !is.null(yMed) && any(is.finite(xMed)) && any(is.finite(yMed))) {
    plot_df <- out_df[is.finite(out_df$xMed) & is.finite(out_df$yMed), , drop = FALSE]
    if (nrow(plot_df) > 0) {
      p_sex <- ggplot(plot_df, aes(x = xMed, y = yMed, color = Sex_Predicted_Norm, shape = Status, label = SampleID)) +
        geom_point(size = 3, alpha = 0.8) +
        theme_minimal() +
        ggtitle("Sex Check: X vs Y Median Intensity") +
        xlab("X median (xMed)") + ylab("Y median (yMed)")
      if (any(plot_df$Status == "mismatch")) {
        p_sex <- p_sex + ggrepel::geom_text_repel(data = plot_df[plot_df$Status == "mismatch", ],
                                                  size = 3, max.overlaps = 20)
      }
      save_interactive_plot(p_sex, "Sex_Check_Plot.html", out_dir)
      save_static_plot(p_sex, "Sex_Check_Plot.png", out_dir, width = 5, height = 4)
    }
  }

  mismatch_count <- sum(mismatch, na.rm = TRUE)
  if (mismatch_count > 0) {
    message(sprintf("WARNING: Sex mismatch detected in %d sample(s).", mismatch_count))
    message(sprintf("  - Column: %s | Action: %s", sex_col, action))
  } else {
    message("  Sex check: no mismatches detected.")
  }
  list(mismatch_samples = mismatch_samples, mismatch_count = mismatch_count, column = sex_col)
}

safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, use = "complete.obs", method = method))
}

compute_effect_size_consistency <- function(raw_res, corr_res, raw_delta_beta, corr_delta_beta,
                                            pval_thresh, lfc_thresh, top_n = 1000) {
  if (is.null(raw_res) || is.null(corr_res)) return(NULL)
  if (nrow(raw_res) == 0 || nrow(corr_res) == 0) return(NULL)
  if (!"CpG" %in% colnames(raw_res)) raw_res$CpG <- rownames(raw_res)
  if (!"CpG" %in% colnames(corr_res)) corr_res$CpG <- rownames(corr_res)
  raw_df <- data.frame(
    CpG = raw_res$CpG,
    logFC_raw = raw_res$logFC,
    adj_raw = raw_res$adj.P.Val,
    stringsAsFactors = FALSE
  )
  if (!is.null(raw_delta_beta)) {
    raw_df$delta_beta_raw <- raw_delta_beta[raw_df$CpG]
  }
  corr_df <- data.frame(
    CpG = corr_res$CpG,
    logFC_corr = corr_res$logFC,
    adj_corr = corr_res$adj.P.Val,
    stringsAsFactors = FALSE
  )
  if (!is.null(corr_delta_beta)) {
    corr_df$delta_beta_corr <- corr_delta_beta[corr_df$CpG]
  }
  merged <- merge(raw_df, corr_df, by = "CpG")
  if (nrow(merged) == 0) return(NULL)
  logfc_corr <- safe_cor(merged$logFC_raw, merged$logFC_corr, method = "pearson")
  logfc_spear <- safe_cor(merged$logFC_raw, merged$logFC_corr, method = "spearman")
  delta_corr <- NA_real_
  if ("delta_beta_raw" %in% colnames(merged) && "delta_beta_corr" %in% colnames(merged)) {
    delta_corr <- safe_cor(merged$delta_beta_raw, merged$delta_beta_corr, method = "pearson")
  }
  sign_conc <- mean(sign(merged$logFC_raw) == sign(merged$logFC_corr), na.rm = TRUE)

  sig_raw <- merged$adj_raw < pval_thresh & abs(merged$logFC_raw) > lfc_thresh
  sig_corr <- merged$adj_corr < pval_thresh & abs(merged$logFC_corr) > lfc_thresh
  sig_union <- sig_raw | sig_corr
  sig_conc <- if (any(sig_union, na.rm = TRUE)) {
    mean(sign(merged$logFC_raw[sig_union]) == sign(merged$logFC_corr[sig_union]), na.rm = TRUE)
  } else {
    NA_real_
  }
  sig_overlap <- sum(sig_raw & sig_corr, na.rm = TRUE)
  sig_jaccard <- if (sum(sig_raw | sig_corr, na.rm = TRUE) > 0) {
    sig_overlap / sum(sig_raw | sig_corr, na.rm = TRUE)
  } else {
    NA_real_
  }

  top_n <- min(top_n, nrow(merged))
  top_corr <- merged[order(merged$adj_corr), ][seq_len(top_n), , drop = FALSE]
  top_raw <- merged[order(merged$adj_raw), ][seq_len(top_n), , drop = FALSE]
  top_overlap <- length(intersect(top_corr$CpG, top_raw$CpG))
  top_union <- length(unique(c(top_corr$CpG, top_raw$CpG)))
  top_jaccard <- if (top_union > 0) top_overlap / top_union else NA_real_
  top_logfc_corr <- safe_cor(top_corr$logFC_raw, top_corr$logFC_corr, method = "pearson")

  list(
    logFC_corr_all = logfc_corr,
    logFC_spearman_all = logfc_spear,
    delta_beta_corr_all = delta_corr,
    sign_concordance_all = sign_conc,
    sign_concordance_sig = sig_conc,
    sig_jaccard = sig_jaccard,
    top_n = top_n,
    top_overlap = top_overlap,
    top_jaccard = top_jaccard,
    top_logFC_corr = top_logfc_corr
  )
}

evaluate_marker_retention <- function(corr_res, markers, pval_thresh, lfc_thresh, raw_res = NULL) {
  if (is.null(markers) || length(markers) == 0) return(NULL)
  if (is.null(corr_res) || nrow(corr_res) == 0) return(NULL)
  if (!"CpG" %in% colnames(corr_res)) corr_res$CpG <- rownames(corr_res)
  corr_hits <- corr_res[corr_res$CpG %in% markers, , drop = FALSE]
  markers_total <- length(unique(markers))
  markers_in_data <- length(unique(corr_hits$CpG))
  sig_mask <- corr_hits$adj.P.Val < pval_thresh & abs(corr_hits$logFC) > lfc_thresh
  markers_sig <- sum(sig_mask, na.rm = TRUE)
  sig_frac <- if (markers_in_data > 0) markers_sig / markers_in_data else NA_real_

  sign_conc <- NA_real_
  if (!is.null(raw_res) && nrow(corr_hits) > 0) {
    if (!"CpG" %in% colnames(raw_res)) raw_res$CpG <- rownames(raw_res)
    raw_hits <- raw_res[raw_res$CpG %in% corr_hits$CpG, c("CpG", "logFC"), drop = FALSE]
    merged <- merge(raw_hits, corr_hits[, c("CpG", "logFC")], by = "CpG", suffixes = c("_raw", "_corr"))
    if (nrow(merged) > 0) {
      sign_conc <- mean(sign(merged$logFC_raw) == sign(merged$logFC_corr), na.rm = TRUE)
    }
  }

  list(
    markers_total = markers_total,
    markers_in_data = markers_in_data,
    markers_sig = markers_sig,
    markers_sig_fraction = sig_frac,
    markers_sign_concordance = sign_conc
  )
}

score_grade <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x >= 0.85) return("EXCELLENT")
  if (x >= 0.70) return("GOOD")
  if (x >= 0.50) return("FAIR")
  "POOR"
}

fmt_val <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  sprintf(paste0("%.", digits, "f"), x)
}

resolve_caf_weights <- function(caf_cfg, preset_name = "conservative") {
  # Conservative default keeps the original 3-term CAI definition.
  # NCS is optional and contributes only when explicitly given a non-zero weight.
  # If any component is unavailable, finite components are renormalized in compute_caf_scores().
  default_weights <- c(calibration = 0.40, preservation = 0.35, batch = 0.25, ncs = 0.00)
  if (is.null(caf_cfg)) return(default_weights)
  profile <- ifelse(is.null(caf_cfg$profile), "auto", tolower(as.character(caf_cfg$profile)))
  if (profile == "auto") {
    profile <- ifelse(preset_name == "aggressive", "discovery", "conservative")
  }
  w_in <- NULL
  if (!is.null(caf_cfg$weights$calibration)) {
    w_in <- caf_cfg$weights
  } else if (!is.null(caf_cfg$weights[[profile]])) {
    w_in <- caf_cfg$weights[[profile]]
  } else {
    w_in <- default_weights
  }
  w_in <- suppressWarnings(unlist(w_in))
  # Overlay user-provided weights onto defaults to keep backward compatibility when new keys are missing.
  w <- default_weights
  for (nm in names(w_in)) {
    if (nm %in% names(w)) w[[nm]] <- suppressWarnings(as.numeric(w_in[[nm]]))
  }
  if (any(!is.finite(w))) w <- default_weights
  w_sum <- sum(w)
  if (!is.finite(w_sum) || w_sum <= 0) return(default_weights)
  w / w_sum
}

caf_weights <- resolve_caf_weights(caf_cfg, preset_name)

compute_caf_scores <- function(perm_summary, perm_n_probes,
                               effect_metrics, marker_metrics,
                               batch_before, batch_after,
                               weights, target_fpr = 0.05,
                               observed_lambda_all = NA_real_,
                               observed_lambda_vp_top = NA_real_,
                               ncs_score = NA_real_) {
  perm_mean_sig <- NA_real_
  perm_lambda <- NA_real_
  perm_ks <- NA_real_
  if (!is.null(perm_summary) && nrow(perm_summary) > 0) {
    perm_mean_sig <- perm_summary$value[perm_summary$stat == "mean_sig"]
    perm_lambda <- perm_summary$value[perm_summary$stat == "lambda_median"]
    perm_ks <- perm_summary$value[perm_summary$stat == "ks_p_median"]
  }
  perm_mean_sig <- suppressWarnings(as.numeric(perm_mean_sig))
  perm_lambda <- suppressWarnings(as.numeric(perm_lambda))
  perm_ks <- suppressWarnings(as.numeric(perm_ks))
  perm_n_probes <- suppressWarnings(as.numeric(perm_n_probes))
  null_fpr <- if (is.finite(perm_mean_sig) && is.finite(perm_n_probes) && perm_n_probes > 0) {
    perm_mean_sig / perm_n_probes
  } else {
    NA_real_
  }
  cal_score <- if (is.finite(null_fpr) && is.finite(target_fpr) && target_fpr > 0) {
    clamp01(1 - abs(null_fpr - target_fpr) / target_fpr)
  } else {
    NA_real_
  }

  marker_ret <- if (!is.null(marker_metrics)) suppressWarnings(as.numeric(marker_metrics$markers_sig_fraction)) else NA_real_
  effect_corr <- if (!is.null(effect_metrics)) suppressWarnings(as.numeric(effect_metrics$logFC_corr_all)) else NA_real_
  pres_score <- if (is.finite(marker_ret) && is.finite(effect_corr)) {
    clamp01(marker_ret * effect_corr)
  } else if (is.finite(effect_corr)) {
    clamp01(effect_corr)
  } else {
    NA_real_
  }

  before_sig <- if (!is.null(batch_before)) suppressWarnings(as.numeric(batch_before$sig_count)) else NA_real_
  after_sig <- if (!is.null(batch_after)) suppressWarnings(as.numeric(batch_after$sig_count)) else NA_real_
  batch_reduction <- if (is.finite(before_sig) && before_sig > 0 && is.finite(after_sig)) {
    clamp01((before_sig - after_sig) / before_sig)
  } else {
    NA_real_
  }

  ncs_score <- suppressWarnings(as.numeric(ncs_score))
  ncs_score <- if (is.finite(ncs_score)) clamp01(ncs_score) else NA_real_

  scores <- c(calibration = cal_score, preservation = pres_score, batch = batch_reduction, ncs = ncs_score)
  w <- weights
  keep <- is.finite(scores)
  cai <- if (any(keep)) sum(scores[keep] * w[keep]) / sum(w[keep]) else NA_real_

  observed_lambda_all <- suppressWarnings(as.numeric(observed_lambda_all))
  observed_lambda_vp_top <- suppressWarnings(as.numeric(observed_lambda_vp_top))
  # Heuristic context metric: only compute ratio when the *overall* observed lambda is inflated.
  # Note: this does not separate biology from group-linked confounding; it is only a rough
  # enrichment-vs-null indicator and must be interpreted with PVCA/covariate context.
  lambda_ratio <- if (is.finite(observed_lambda_all) && observed_lambda_all > 1.2 &&
                      is.finite(observed_lambda_vp_top) &&
                      is.finite(perm_lambda) && perm_lambda > 0) {
    observed_lambda_vp_top / perm_lambda
  } else {
    NA_real_
  }

  list(
    null_fpr = null_fpr,
    null_lambda = perm_lambda,
    null_ks_p = perm_ks,
    calibration_score = cal_score,
    marker_retention = marker_ret,
    effect_concordance = effect_corr,
    preservation_score = pres_score,
    batch_reduction = batch_reduction,
    ncs_score = ncs_score,
    cai = cai,
    weights = w,
    observed_lambda_all = observed_lambda_all,
    observed_lambda_vp_top = observed_lambda_vp_top,
    lambda_ratio = lambda_ratio
  )
}

flag_ok <- function(cond) {
  if (isTRUE(cond)) "OK" else "WARN"
}

build_caf_report_lines <- function(prefix, applied_batch_method, n_con, n_test,
                                   perm_n, perm_n_probes, caf_res, marker_metrics, target_fpr,
                                   observed_lambda_all = NA_real_,
                                   observed_lambda_vp_top = NA_real_) {
  if (is.null(caf_res)) return(NULL)
  total_n <- n_con + n_test
  null_fpr <- caf_res$null_fpr
  cal_ok <- is.finite(null_fpr) && abs(null_fpr - target_fpr) <= target_fpr
  lambda_ok <- is.finite(caf_res$null_lambda) && caf_res$null_lambda >= 0.8 && caf_res$null_lambda <= 1.2
  ks_ok <- is.finite(caf_res$null_ks_p) && caf_res$null_ks_p > 0.05
  marker_count <- if (!is.null(marker_metrics)) marker_metrics$markers_in_data else NA
  marker_ret <- caf_res$marker_retention
  effect_corr <- caf_res$effect_concordance
  observed_lambda_all <- suppressWarnings(as.numeric(observed_lambda_all))
  observed_lambda_vp_top <- suppressWarnings(as.numeric(observed_lambda_vp_top))
  lines <- c(
    "CORRECTION ADEQUACY REPORT",
    "===============================================",
    sprintf("Pipeline: %s", prefix),
    sprintf("Batch method: %s", applied_batch_method),
    sprintf("Samples: %d (Control=%d, Test=%d)", total_n, n_con, n_test),
    "",
    "1. CALIBRATION (Permutation-based)",
    sprintf("  Permutations run: %s", ifelse(is.finite(perm_n), perm_n, "NA")),
    sprintf("  Null FP rate: %s (target %.2f) [%s]", fmt_val(null_fpr), target_fpr, flag_ok(cal_ok)),
    sprintf("  Null Lambda: %s [%s]", fmt_val(caf_res$null_lambda), flag_ok(lambda_ok)),
    sprintf("  Observed Lambda (all CpGs): %s", fmt_val(observed_lambda_all)),
    sprintf("  Observed Lambda (vp_top CpGs): %s (n=%s)", fmt_val(observed_lambda_vp_top),
            ifelse(is.finite(perm_n_probes), as.integer(perm_n_probes), "NA")),
    sprintf("  Lambda Ratio (vp_top obs/null): %s [%s]", fmt_val(caf_res$lambda_ratio),
            if (is.finite(caf_res$lambda_ratio)) {
              if (caf_res$lambda_ratio > 2) "strong deviation vs null (heuristic; review PVCA/covariates)"
              else if (caf_res$lambda_ratio > 1.2) "moderate deviation vs null (heuristic; review PVCA/covariates)"
              else "close to null (heuristic; weak enrichment vs null)"
            } else {
              if (is.finite(observed_lambda_all) && observed_lambda_all <= 1.2) "N/A (observed lambda not inflated)"
              else "N/A"
            }),
    sprintf("  KS test p-value: %s [%s]", fmt_val(caf_res$null_ks_p), flag_ok(ks_ok)),
    sprintf("  CALIBRATION SCORE: %s [%s]", fmt_val(caf_res$calibration_score), score_grade(caf_res$calibration_score)),
    "",
    "2. SIGNAL PRESERVATION",
    sprintf("  Known markers tested: %s", ifelse(is.finite(marker_count), marker_count, "NA")),
    sprintf("  Marker retention: %s", fmt_val(marker_ret)),
    sprintf("  Effect concordance (logFC r): %s", fmt_val(effect_corr)),
    sprintf("  PRESERVATION SCORE: %s [%s]", fmt_val(caf_res$preservation_score), score_grade(caf_res$preservation_score)),
    "",
    "3. BATCH REMOVAL",
    sprintf("  Batch reduction ratio: %s", fmt_val(caf_res$batch_reduction)),
    sprintf("  BATCH REMOVAL SCORE: %s [%s]", fmt_val(caf_res$batch_reduction), score_grade(caf_res$batch_reduction)),
    "",
    "4. NEGATIVE CONTROL STABILITY (NCS)",
    sprintf("  NCS score: %s [%s]", fmt_val(caf_res$ncs_score), score_grade(caf_res$ncs_score)),
    "",
    "OVERALL CAI",
    sprintf("  CAI: %s [%s]", fmt_val(caf_res$cai), score_grade(caf_res$cai))
  )
  if (is.null(marker_metrics) || !is.finite(marker_ret)) {
    lines <- c(lines, "", "NOTE: Marker list not provided or no markers found; preservation uses effect concordance only.")
  }
  if (!is.finite(perm_n_probes)) {
    lines <- c(lines, "", "NOTE: Permutation probe count unavailable; calibration may be incomplete.")
  }
  lines
}

resolve_crf_tier <- function(n_con, n_test, crf_cfg) {
  total_n <- n_con + n_test
  min_pg <- min(n_con, n_test)
  tier_cfg <- if (!is.null(crf_cfg$sample_tiers)) crf_cfg$sample_tiers else list()
  get_tier <- function(name, fallback) {
    if (is.null(tier_cfg[[name]])) return(fallback)
    merge_config_defaults(fallback, tier_cfg[[name]])
  }
  minimal_def <- list(threshold = 12, min_per_group = 3,
                      mmc = list(methods = c("none", "limma"), concordance_levels = c("core", "discordant")),
                      ncs = list(types = c("snp"), lambda_range = c(0.7, 1.3), interpret = "reference_only"),
                      rss = list(mode = "bootstrap", n_iterations = 100, top_k = c(20, 50), overlap_threshold = 0.20),
                      warnings = c("CAUTION: Very small sample size (n<12). Statistical power is severely limited.",
                                   "Treat all results as exploratory hypotheses requiring independent validation.",
                                   "False discovery rate control may be unreliable at this sample size."))
  small_def <- list(threshold = 24, min_per_group = 6,
                    mmc = list(methods = c("none", "combat_np", "limma"), concordance_levels = c("core", "consensus", "weak")),
                    ncs = list(types = c("snp"), lambda_range = c(0.8, 1.2), interpret = "cautious"),
                    rss = list(mode = "leave_pair_out", n_iterations = 50, top_k = c(30, 100, 200), overlap_threshold = 0.30),
                    warnings = c("NOTE: Small sample size (n<24). Results should be validated in an independent cohort.",
                                 "SVA is excluded due to instability at this sample size.",
                                 "Effect size estimates may be imprecise; focus on direction rather than magnitude."))
  moderate_def <- list(threshold = 50, min_per_group = 12,
                       mmc = list(methods = c("none", "combat", "limma", "sva"), concordance_levels = c("core", "consensus", "weak")),
                       ncs = list(types = c("snp", "housekeeping"), lambda_range = c(0.9, 1.1), interpret = "standard"),
                       rss = list(mode = "split", n_iterations = 10, top_k = c(50, 200, 500), overlap_threshold = 0.40),
                       warnings = character(0))
  large_def <- list(threshold = NA, min_per_group = 25,
                    mmc = list(methods = c("none", "combat", "limma", "sva"), concordance_levels = c("core", "consensus", "weak")),
                    ncs = list(types = c("snp", "housekeeping"), lambda_range = c(0.9, 1.1), interpret = "standard"),
                    rss = list(mode = "split", n_iterations = 10, top_k = c(100, 500, 1000), overlap_threshold = 0.50),
                    warnings = character(0))

  minimal_def <- get_tier("minimal", minimal_def)
  small_def <- get_tier("small", small_def)
  moderate_def <- get_tier("moderate", moderate_def)
  large_def <- get_tier("large", large_def)

  thr_min <- suppressWarnings(as.numeric(minimal_def$threshold))
  thr_small <- suppressWarnings(as.numeric(small_def$threshold))
  thr_mod <- suppressWarnings(as.numeric(moderate_def$threshold))
  min_pg_min <- suppressWarnings(as.numeric(minimal_def$min_per_group))
  min_pg_small <- suppressWarnings(as.numeric(small_def$min_per_group))
  min_pg_mod <- suppressWarnings(as.numeric(moderate_def$min_per_group))
  min_pg_large <- suppressWarnings(as.numeric(large_def$min_per_group))
  if (!is.finite(thr_min)) thr_min <- 12
  if (!is.finite(thr_small)) thr_small <- 24
  if (!is.finite(thr_mod)) thr_mod <- 50
  if (!is.finite(min_pg_min)) min_pg_min <- 3
  if (!is.finite(min_pg_small)) min_pg_small <- 6
  if (!is.finite(min_pg_mod)) min_pg_mod <- 12
  if (!is.finite(min_pg_large)) min_pg_large <- 25

  tier_name <- "large"
  tier_def <- large_def
  if (total_n < thr_min || min_pg < min_pg_small) {
    tier_name <- "minimal"
    tier_def <- minimal_def
  } else if (total_n < thr_small || min_pg < min_pg_mod) {
    tier_name <- "small"
    tier_def <- small_def
  } else if (total_n < thr_mod || min_pg < min_pg_large) {
    tier_name <- "moderate"
    tier_def <- moderate_def
  }

  mmc_methods <- if (!is.null(tier_def$mmc$methods)) as.character(tier_def$mmc$methods) else character(0)
  combat_par_prior <- TRUE
  combat_allowed <- any(grepl("^combat", mmc_methods))
  if (any(mmc_methods %in% c("combat_np", "combat-np"))) {
    combat_par_prior <- FALSE
    mmc_methods <- gsub("combat_np", "combat", mmc_methods)
  }
  mmc_methods <- unique(mmc_methods)
  sva_allowed <- "sva" %in% mmc_methods

  warnings <- tier_def$warnings
  if (min_pg < min_pg_min) {
    warnings <- unique(c(warnings, "Sample size below minimum threshold; robustness assessment may be unreliable."))
  }
  imbalance_ratio <- if (max(n_con, n_test) > 0) min_pg / max(n_con, n_test) else 0
  if (is.finite(imbalance_ratio) && imbalance_ratio < 0.4) {
    warnings <- unique(c(warnings, "Highly unbalanced groups may bias stability assessments."))
  }

  list(
    n_con = n_con,
    n_test = n_test,
    tier = tier_name,
    total_n = total_n,
    min_per_group = min_pg,
    min_required = min_pg_min,
    eligible = min_pg >= min_pg_min,
    mmc_methods = mmc_methods,
    mmc_levels = tier_def$mmc$concordance_levels,
    combat_par_prior = combat_par_prior,
    combat_allowed = combat_allowed,
    sva_allowed = sva_allowed,
    ncs_types = tier_def$ncs$types,
    ncs_lambda_range = tier_def$ncs$lambda_range,
    ncs_interpret = tier_def$ncs$interpret,
    rss_mode = tier_def$rss$mode,
    rss_iterations = tier_def$rss$n_iterations,
    rss_top_k = tier_def$rss$top_k,
    rss_overlap_threshold = tier_def$rss$overlap_threshold,
    warnings = warnings
  )
}

subset_top_variable <- function(betas, max_probes) {
  if (!is.finite(max_probes) || max_probes <= 0 || nrow(betas) <= max_probes) return(betas)
  vars <- apply(betas, 1, var, na.rm = TRUE)
  idx <- head(order(vars, decreasing = TRUE), max_probes)
  betas[idx, , drop = FALSE]
}

compute_mmc_summary <- function(res_list, top_k = 100, levels = c("core", "consensus", "weak"),
                                pval_thresh = 0.05, lfc_thresh = 0) {
  if (is.null(res_list) || length(res_list) < 2) return(NULL)
  methods <- names(res_list)
  sig_sets <- lapply(res_list, function(df) {
    if (!"CpG" %in% colnames(df)) df$CpG <- rownames(df)
    df <- df[order(df$adj.P.Val), ]
    df$CpG
  })
  effect_maps <- lapply(res_list, function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!"CpG" %in% colnames(df)) df$CpG <- rownames(df)
    eff_col <- NULL
    if ("logFC" %in% colnames(df)) eff_col <- "logFC"
    if (is.null(eff_col) && "Effect" %in% colnames(df)) eff_col <- "Effect"
    if (is.null(eff_col)) return(NULL)
    stats <- df[, c("CpG", eff_col), drop = FALSE]
    setNames(stats[[eff_col]], stats$CpG)
  })
  top_k <- as.integer(top_k)
  top_k <- top_k[is.finite(top_k) & top_k > 0]
  if (length(top_k) == 0) return(NULL)
  out <- data.frame()
  m_count <- length(methods)
  for (k in top_k) {
    top_sets <- lapply(sig_sets, function(cpgs) head(cpgs, min(k, length(cpgs))))
    all_cpgs <- unique(unlist(top_sets))
    freq <- table(unlist(top_sets))
    core_n <- sum(freq >= m_count)
    core_pct <- if (length(all_cpgs) > 0) core_n / length(all_cpgs) else NA_real_
    consensus_n <- if (m_count >= 3) sum(freq >= (m_count - 1)) else NA_real_
    weak_n <- if (m_count >= 3) sum(freq >= 2) else NA_real_
    discordant_n <- if (m_count == 2) length(all_cpgs) - core_n else NA_real_
    spearman_vals <- numeric(0)
    direction_vals <- numeric(0)
    if (length(effect_maps) >= 2) {
      for (i in seq_len(length(effect_maps) - 1)) {
        for (j in (i + 1):length(effect_maps)) {
          eff_i <- effect_maps[[i]]
          eff_j <- effect_maps[[j]]
          if (is.null(eff_i) || is.null(eff_j)) next
          common_ids <- intersect(names(eff_i), names(eff_j))
          if (length(common_ids) < 3) next
          common_ids <- intersect(common_ids, all_cpgs)
          if (length(common_ids) < 3) next
          v_i <- eff_i[common_ids]
          v_j <- eff_j[common_ids]
          ok <- is.finite(v_i) & is.finite(v_j)
          if (sum(ok) < 3) next
          rho <- suppressWarnings(cor(v_i[ok], v_j[ok], method = "spearman"))
          if (is.finite(rho)) spearman_vals <- c(spearman_vals, rho)
          dir_cons <- mean(sign(v_i[ok]) == sign(v_j[ok]))
          if (is.finite(dir_cons)) direction_vals <- c(direction_vals, dir_cons)
        }
      }
    }
    spearman_mean <- if (length(spearman_vals) > 0) mean(spearman_vals) else NA_real_
    spearman_min <- if (length(spearman_vals) > 0) min(spearman_vals) else NA_real_
    spearman_max <- if (length(spearman_vals) > 0) max(spearman_vals) else NA_real_
    dir_mean <- if (length(direction_vals) > 0) mean(direction_vals) else NA_real_
    dir_min <- if (length(direction_vals) > 0) min(direction_vals) else NA_real_
    dir_max <- if (length(direction_vals) > 0) max(direction_vals) else NA_real_

    # Composite (heuristic) MMC score combining identity + effect size + direction concordance.
    # Keep this conservative and interpretable; all components are in [0, 1] after scaling.
    identity_score <- core_pct
    effect_score <- if (is.finite(spearman_mean)) (spearman_mean + 1) / 2 else NA_real_
    direction_score <- dir_mean
    mmc_weights <- c(identity = 0.5, effect = 0.3, direction = 0.2)
    mmc_vals <- c(identity = identity_score, effect = effect_score, direction = direction_score)
    mmc_ok <- is.finite(mmc_vals)
    mmc_composite <- if (any(mmc_ok)) {
      sum(mmc_weights[mmc_ok] * mmc_vals[mmc_ok]) / sum(mmc_weights[mmc_ok])
    } else {
      NA_real_
    }
    out <- rbind(out, data.frame(
      top_k = k,
      methods = m_count,
      core = core_n,
      core_pct = core_pct,
      consensus = consensus_n,
      weak = weak_n,
      discordant = discordant_n,
      total_unique = length(all_cpgs),
      spearman_mean = spearman_mean,
      spearman_min = spearman_min,
      spearman_max = spearman_max,
      spearman_pairs = length(spearman_vals),
      direction_mean = dir_mean,
      direction_min = dir_min,
      direction_max = dir_max,
      direction_pairs = length(direction_vals),
      mmc_composite = mmc_composite,
      stringsAsFactors = FALSE
    ))
  }
  out
}

compute_negative_control_stats <- function(betas_nc, design, group_cols,
                                           lambda_ci = FALSE, n_boot = 0L, sample_frac = 0.8) {
  if (is.null(betas_nc) || nrow(betas_nc) < 2) return(NULL)
  res <- tryCatch(run_limma_dmp(betas_nc, design, group_cols), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  pvals <- suppressWarnings(as.numeric(res$P.Value))
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  if (length(pvals) < 2) return(NULL)
  sig_rate <- mean(pvals < 0.05, na.rm = TRUE)
  chisq_vals <- suppressWarnings(qchisq(1 - pvals, 1))
  chisq_vals <- chisq_vals[is.finite(chisq_vals)]
  lambda_val <- compute_genomic_lambda(pvals)
  ci_low <- NA_real_
  ci_high <- NA_real_
  boot_n <- 0L
  if (isTRUE(lambda_ci) && is.finite(lambda_val)) {
    n_boot <- suppressWarnings(as.integer(n_boot))
    if (!is.finite(n_boot)) n_boot <- 0L
    sample_frac <- suppressWarnings(as.numeric(sample_frac))
    if (!is.finite(sample_frac) || sample_frac <= 0 || sample_frac > 1) sample_frac <- 0.8
    if (n_boot >= 20 && length(chisq_vals) >= 10) {
      boot_n <- n_boot
      boot_l <- numeric(n_boot)
      boot_size <- max(5L, as.integer(round(length(chisq_vals) * sample_frac)))
      for (i in seq_len(n_boot)) {
        idx <- sample.int(length(chisq_vals), boot_size, replace = TRUE)
        boot_l[i] <- median(chisq_vals[idx]) / qchisq(0.5, 1)
      }
      qs <- suppressWarnings(stats::quantile(boot_l, probs = c(0.025, 0.975), na.rm = TRUE))
      if (length(qs) == 2) {
        ci_low <- qs[[1]]
        ci_high <- qs[[2]]
      }
    }
  }
  list(lambda = lambda_val, sig_rate = sig_rate, n = nrow(betas_nc),
       lambda_ci_low = ci_low, lambda_ci_high = ci_high, lambda_boot_n = boot_n)
}

compute_rss_overlap <- function(betas, targets, group_col, covariates, mode, n_iter, top_k, max_probes = 0) {
  if (is.null(betas) || nrow(betas) == 0) return(NULL)
  betas_use <- subset_top_variable(betas, max_probes)
  covariates <- covariates[covariates %in% colnames(targets)]
  formula_str <- "~ 0 + primary_group"
  if (length(covariates) > 0) {
    formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  }
  design_full <- tryCatch(model.matrix(as.formula(formula_str), data = targets), error = function(e) NULL)
  if (is.null(design_full) || ncol(design_full) < 2) return(NULL)
  colnames(design_full) <- make.names(gsub("primary_group", "", colnames(design_full)))
  group_levels <- make.names(levels(targets[[group_col]]))
  ref_res <- tryCatch(run_limma_dmp(betas_use, design_full, group_levels), error = function(e) NULL)
  if (is.null(ref_res) || nrow(ref_res) == 0) return(NULL)
  if (!"CpG" %in% colnames(ref_res)) ref_res$CpG <- rownames(ref_res)
  ref_res <- ref_res[order(ref_res$adj.P.Val), ]
  ref_effect <- NULL
  if ("logFC" %in% colnames(ref_res)) {
    ref_effect <- setNames(ref_res$logFC, ref_res$CpG)
  }
  ref_sets <- lapply(top_k, function(k) head(ref_res$CpG, min(k, nrow(ref_res))))
  names(ref_sets) <- as.character(top_k)

  # Rank-Biased Overlap (RBO) for ranked lists (Webber et al., 2010).
  # p controls top-weighting; p=0.9 gives substantial weight to the top ranks.
  rbo_p <- 0.9
  rbo_score <- function(a, b, p = rbo_p) {
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    if (!is.finite(p) || p <= 0 || p >= 1) p <- 0.9
    k <- min(length(a), length(b))
    if (k < 1) return(NA_real_)
    sa <- character(0)
    sb <- character(0)
    A <- numeric(k)
    for (d in seq_len(k)) {
      sa <- c(sa, a[d])
      sb <- c(sb, b[d])
      A[d] <- length(intersect(sa, sb)) / d
    }
    (1 - p) * sum(p^(0:(k - 1)) * A) + (p^k) * A[k]
  }

  iter_overlaps <- lapply(top_k, function(k) numeric(0))
  names(iter_overlaps) <- as.character(top_k)
  iter_jaccards <- lapply(top_k, function(k) numeric(0))
  names(iter_jaccards) <- as.character(top_k)
  iter_rbos <- lapply(top_k, function(k) numeric(0))
  names(iter_rbos) <- as.character(top_k)
  iter_rss <- lapply(top_k, function(k) numeric(0))
  names(iter_rss) <- as.character(top_k)
  iter_signs <- lapply(top_k, function(k) numeric(0))
  names(iter_signs) <- as.character(top_k)

  for (i in seq_len(n_iter)) {
    idx_keep <- rep(TRUE, nrow(targets))
    sample_idx <- NULL
    if (mode == "split") {
      idx_keep <- rep(FALSE, nrow(targets))
      for (lvl in levels(targets[[group_col]])) {
        grp_idx <- which(targets[[group_col]] == lvl)
        n_take <- floor(length(grp_idx) / 2)
        take_idx <- sample(grp_idx, n_take)
        idx_keep[take_idx] <- TRUE
      }
    } else if (mode == "leave_pair_out") {
      idx_keep <- rep(TRUE, nrow(targets))
      for (lvl in levels(targets[[group_col]])) {
        grp_idx <- which(targets[[group_col]] == lvl)
        if (length(grp_idx) >= 2) {
          drop_idx <- sample(grp_idx, 2)
          idx_keep[drop_idx] <- FALSE
        }
      }
    } else if (mode == "bootstrap") {
      keep_idx <- c()
      for (lvl in levels(targets[[group_col]])) {
        grp_idx <- which(targets[[group_col]] == lvl)
        if (length(grp_idx) > 0) {
          keep_idx <- c(keep_idx, sample(grp_idx, length(grp_idx), replace = TRUE))
        }
      }
      sample_idx <- keep_idx
    }
    if (!is.null(sample_idx)) {
      if (length(sample_idx) < 4) next
      betas_iter <- betas_use[, sample_idx, drop = FALSE]
      targets_iter <- targets[sample_idx, , drop = FALSE]
    } else {
      if (sum(idx_keep) < 4) next
      betas_iter <- betas_use[, idx_keep, drop = FALSE]
      targets_iter <- targets[idx_keep, , drop = FALSE]
    }
    if (length(unique(targets_iter[[group_col]])) < 2) next
    design_iter <- tryCatch(model.matrix(as.formula(formula_str), data = targets_iter), error = function(e) NULL)
    if (is.null(design_iter) || ncol(design_iter) < 2) next
    colnames(design_iter) <- make.names(gsub("primary_group", "", colnames(design_iter)))
    res_iter <- tryCatch(run_limma_dmp(betas_iter, design_iter, make.names(levels(targets_iter[[group_col]]))), error = function(e) NULL)
    if (is.null(res_iter) || nrow(res_iter) == 0) next
    if (!"CpG" %in% colnames(res_iter)) res_iter$CpG <- rownames(res_iter)
    res_iter <- res_iter[order(res_iter$adj.P.Val), ]
    iter_effect <- NULL
    if ("logFC" %in% colnames(res_iter)) {
      iter_effect <- setNames(res_iter$logFC, res_iter$CpG)
    }
    for (k in top_k) {
      ref_set <- ref_sets[[as.character(k)]]
      iter_set <- head(res_iter$CpG, min(k, nrow(res_iter)))
      inter_n <- length(intersect(ref_set, iter_set))
      union_n <- length(unique(c(ref_set, iter_set)))
      overlap <- inter_n / max(1, k)
      jaccard <- if (union_n > 0) inter_n / union_n else NA_real_
      rbo <- rbo_score(ref_set, iter_set, p = rbo_p)
      rss_val <- NA_real_
      if (is.finite(jaccard) && is.finite(rbo)) {
        rss_val <- 0.5 * jaccard + 0.5 * rbo
      } else if (is.finite(jaccard)) {
        rss_val <- jaccard
      } else if (is.finite(rbo)) {
        rss_val <- rbo
      }
      iter_overlaps[[as.character(k)]] <- c(iter_overlaps[[as.character(k)]], overlap)
      iter_jaccards[[as.character(k)]] <- c(iter_jaccards[[as.character(k)]], jaccard)
      iter_rbos[[as.character(k)]] <- c(iter_rbos[[as.character(k)]], rbo)
      iter_rss[[as.character(k)]] <- c(iter_rss[[as.character(k)]], rss_val)
      if (!is.null(ref_effect) && !is.null(iter_effect)) {
        common_ids <- intersect(ref_set, iter_set)
        if (length(common_ids) >= 3) {
          ref_vals <- ref_effect[common_ids]
          iter_vals <- iter_effect[common_ids]
          ok <- is.finite(ref_vals) & is.finite(iter_vals)
          if (sum(ok) >= 3) {
            sign_cons <- mean(sign(ref_vals[ok]) == sign(iter_vals[ok]))
            iter_signs[[as.character(k)]] <- c(iter_signs[[as.character(k)]], sign_cons)
          }
        }
      }
    }
  }
  out <- data.frame()
  for (k in top_k) {
    vals <- iter_overlaps[[as.character(k)]]
    j_vals <- iter_jaccards[[as.character(k)]]
    rbo_vals <- iter_rbos[[as.character(k)]]
    rss_vals <- iter_rss[[as.character(k)]]
    sign_vals <- iter_signs[[as.character(k)]]
    out <- rbind(out, data.frame(top_k = k,
                                 overlap_mean = if (length(vals) > 0) mean(vals, na.rm = TRUE) else NA_real_,
                                 overlap_sd = if (length(vals) > 1) sd(vals, na.rm = TRUE) else NA_real_,
                                 jaccard_mean = if (length(j_vals) > 0) mean(j_vals, na.rm = TRUE) else NA_real_,
                                 jaccard_sd = if (length(j_vals) > 1) sd(j_vals, na.rm = TRUE) else NA_real_,
                                 rbo_mean = if (length(rbo_vals) > 0) mean(rbo_vals, na.rm = TRUE) else NA_real_,
                                 rbo_sd = if (length(rbo_vals) > 1) sd(rbo_vals, na.rm = TRUE) else NA_real_,
                                 rss_mean = if (length(rss_vals) > 0) mean(rss_vals, na.rm = TRUE) else NA_real_,
                                 rss_sd = if (length(rss_vals) > 1) sd(rss_vals, na.rm = TRUE) else NA_real_,
                                 rbo_p = rbo_p,
                                 sign_mean = if (length(sign_vals) > 0) mean(sign_vals, na.rm = TRUE) else NA_real_,
                                 sign_sd = if (length(sign_vals) > 1) sd(sign_vals, na.rm = TRUE) else NA_real_,
                                 n_iterations = length(vals)))
  }
  out
}

build_crf_report_lines <- function(tier_info, mmc_summary, ncs_summary, rss_summary, prefix,
                                   pvca_ci_before = NA_real_, pvca_ci_after = NA_real_) {
  tier_name <- toupper(tier_info$tier)
  lines <- c(
    "CORRECTION ROBUSTNESS REPORT (CRF)",
    "",
    sprintf("Pipeline: %s", prefix),
    sprintf("Samples: %d (Control=%d, Test=%d)", tier_info$total_n, tier_info$n_con, tier_info$n_test),
    sprintf("Sample Tier: %s", tier_name)
  )
  if (length(tier_info$warnings) > 0) {
    lines <- c(lines, "", "IMPORTANT LIMITATIONS:", paste0("- ", tier_info$warnings))
  }
  lines <- c(lines, "", "SECTION 1: MULTI-METHOD CONCORDANCE")
  if (is.null(mmc_summary) || nrow(mmc_summary) == 0) {
    lines <- c(lines, "  MMC was not run (insufficient methods or data).")
  } else {
    lines <- c(lines, sprintf("  Methods compared: %s", paste(tier_info$mmc_methods, collapse = ", ")))
    for (i in seq_len(nrow(mmc_summary))) {
      row <- mmc_summary[i, ]
      core_pct_disp <- if (!is.null(row$core_pct) && is.finite(row$core_pct)) sprintf(" (%.0f%%)", 100 * row$core_pct) else ""
      rho_disp <- if (is.finite(row$spearman_mean)) sprintf(", rho=%.3f", row$spearman_mean) else ""
      dir_disp <- if (!is.null(row$direction_mean) && is.finite(row$direction_mean)) sprintf(", dir=%.3f", row$direction_mean) else ""
      mmc_disp <- if (!is.null(row$mmc_composite) && is.finite(row$mmc_composite)) sprintf(", MMC=%.3f", row$mmc_composite) else ""
      lines <- c(lines, sprintf("  Top K=%d: core=%s%s, consensus=%s, weak=%s, discordant=%s (unique=%s)%s%s%s",
                                row$top_k,
                                fmt_val(row$core), core_pct_disp,
                                fmt_val(row$consensus), fmt_val(row$weak),
                                fmt_val(row$discordant), fmt_val(row$total_unique),
                                rho_disp, dir_disp, mmc_disp))
    }
  }
  lines <- c(lines, "", "SECTION 2: NEGATIVE CONTROL STABILITY")
  if (is.null(ncs_summary) || nrow(ncs_summary) == 0) {
    lines <- c(lines, "  Negative control analysis not available.")
  } else {
    for (i in seq_len(nrow(ncs_summary))) {
      row <- ncs_summary[i, ]
      lambda_ci_disp <- ""
      if (!is.null(row$lambda_ci_low) && is.finite(row$lambda_ci_low) &&
          !is.null(row$lambda_ci_high) && is.finite(row$lambda_ci_high)) {
        lambda_ci_disp <- sprintf(", CI=[%.3f, %.3f]", row$lambda_ci_low, row$lambda_ci_high)
      }
      ncs_score_disp <- ""
      if (!is.null(row$ncs_score) && is.finite(row$ncs_score)) {
        ncs_score_disp <- sprintf(", score=%.2f", row$ncs_score)
      }
      lines <- c(lines, sprintf("  %s (%s): lambda=%.3f%s%s, sig_rate=%.3f (n=%d)",
                                row$type, row$stage, row$lambda, lambda_ci_disp,
                                ncs_score_disp, row$sig_rate, row$n))
    }
  }
  lines <- c(lines, "", "SECTION 3: RESAMPLING STABILITY (RSS)")
  if (is.null(rss_summary) || nrow(rss_summary) == 0) {
    lines <- c(lines, "  Stability analysis not available.")
  } else {
    for (i in seq_len(nrow(rss_summary))) {
      row <- rss_summary[i, ]
      sign_disp <- ""
      if (!is.null(row$sign_mean_raw) && is.finite(row$sign_mean_raw) &&
          !is.null(row$sign_mean_corr) && is.finite(row$sign_mean_corr)) {
        sign_disp <- sprintf(", sign(before)=%.3f, sign(after)=%.3f", row$sign_mean_raw, row$sign_mean_corr)
      }
      rss_disp <- ""
      if (!is.null(row$rss_mean_corr) && is.finite(row$rss_mean_corr)) {
        j <- if (!is.null(row$jaccard_mean_corr) && is.finite(row$jaccard_mean_corr)) sprintf("%.3f", row$jaccard_mean_corr) else "NA"
        r <- if (!is.null(row$rbo_mean_corr) && is.finite(row$rbo_mean_corr)) sprintf("%.3f", row$rbo_mean_corr) else "NA"
        p <- if (!is.null(row$rbo_p) && is.finite(row$rbo_p)) sprintf("%.2f", row$rbo_p) else "NA"
        rss_disp <- sprintf(", RSS(after)=%.3f (J=%s, RBO=%s, p=%s)", row$rss_mean_corr, j, r, p)
      }
      lines <- c(lines, sprintf("  Top K=%d: overlap before=%.3f (sd=%.3f), after=%.3f (sd=%.3f)%s%s",
                                row$top_k, row$overlap_mean_raw, row$overlap_sd_raw,
                                row$overlap_mean_corr, row$overlap_sd_corr, sign_disp, rss_disp))
    }
  }

  pvca_ci_before <- suppressWarnings(as.numeric(pvca_ci_before))
  pvca_ci_after <- suppressWarnings(as.numeric(pvca_ci_after))
  if (is.finite(pvca_ci_before) || is.finite(pvca_ci_after)) {
    lines <- c(lines, "", "SECTION 4: CONFOUNDING VARIANCE DECOMPOSITION (PVCA)")
    if (is.finite(pvca_ci_before)) {
      lines <- c(lines, sprintf("  PVCA confounding index (batch/group) before correction: %.3f", pvca_ci_before))
    }
    if (is.finite(pvca_ci_after)) {
      lines <- c(lines, sprintf("  PVCA confounding index (batch/group) after correction: %.3f", pvca_ci_after))
      cvd_score <- 1 / (1 + pvca_ci_after)
      lines <- c(lines, sprintf("  CVD score (1/(1+CI)): %.3f (higher = safer)", cvd_score))
    }
  }
  lines
}

#' Correction Robustness Framework (CRF) Assessment
#'
#' Evaluates batch correction robustness through multiple complementary metrics:
#' negative control stability, multi-method concordance, resampling stability,
#' and PVCA-based confounding variance decomposition.
#' Produces a comprehensive CRF report with tier-appropriate assessments.
#'
#' @param betas_raw Raw beta matrix before batch correction
#' @param betas_corr Corrected beta matrix after batch correction
#' @param targets Sample metadata data frame
#' @param covariates Character vector of covariate column names
#' @param tier_info List from get_crf_tier() containing sample tier info
#' @param batch_col Character name of batch variable column
#' @param applied_batch_method String identifying which correction was applied
#' @param out_dir Output directory for CRF reports
#' @param prefix Pipeline prefix (e.g., "Minfi_Noob", "Sesame_pOOBAH")
#' @param pval_thresh P-value threshold for significance
#' @param lfc_thresh Log-fold change threshold for effect size
#' @param config_settings Full config list including crf settings
#' @param rgSet Optional RGChannelSet for SNP probe analysis
#' @param mmc_res Optional multi-method concordance results
#'
#' @return List containing CRF assessment results and pass/fail status
#'
#' @details
#' CRF includes four assessment axes:
#' 1. Negative Control Stability (NCS): Lambda on control probes (SNP, OOB)
#' 2. Multi-Method Concordance (MMC): Agreement across correction methods
#' 3. Resampling Stability Score (RSS): Top-K hit overlap across resampled replicates
#' 4. Confounding Variance Decomposition (CVD): PVCA-derived confounding index
#'
#' @seealso get_crf_tier for sample tier classification
run_crf_assessment <- function(betas_raw, betas_corr, targets, covariates, tier_info,
                               batch_col, applied_batch_method, out_dir, prefix,
                               pval_thresh, lfc_thresh,
                               config_settings, rgSet = NULL, mmc_res = NULL,
                               pvca_ci_before = NA_real_, pvca_ci_after = NA_real_) {
  if (is.null(tier_info) || !isTRUE(config_settings$crf$enabled)) return(NULL)
  crf_cfg <- config_settings$crf
  if (!isTRUE(tier_info$eligible)) {
    warn_lines <- c("CRF assessment skipped: sample size below minimum threshold.")
    lines <- c("CORRECTION ROBUSTNESS REPORT (CRF)",
               "", sprintf("Pipeline: %s", prefix),
               sprintf("Samples: %d (Control=%d, Test=%d)", tier_info$total_n, tier_info$n_con, tier_info$n_test),
               sprintf("Sample Tier: %s", toupper(tier_info$tier)), "",
               warn_lines)
    writeLines(lines, file.path(out_dir, paste0(prefix, "_CRF_Report.txt")))
    return(list(tier = tier_info$tier))
  }

  mmc_summary <- NULL
  if (!is.null(mmc_res) && !is.null(mmc_res$res_list) && length(mmc_res$res_list) >= 2) {
    mmc_summary <- compute_mmc_summary(mmc_res$res_list,
                                       top_k = tier_info$rss_top_k,
                                       levels = tier_info$mmc_levels,
                                       pval_thresh = pval_thresh,
                                       lfc_thresh = lfc_thresh)
    if (!is.null(mmc_summary)) {
      write.csv(mmc_summary, file.path(out_dir, paste0(prefix, "_CRF_MMC_Summary.csv")), row.names = FALSE)
    }
  }

  ncs_summary <- data.frame()
  if (!is.null(betas_raw) && !is.null(betas_corr)) {
    ncs_types <- tier_info$ncs_types
    if (!is.null(ncs_types) && length(ncs_types) > 0) {
      ncs_ci_cfg <- crf_cfg$ncs_lambda_ci
      ncs_ci_enabled <- !is.null(ncs_ci_cfg) && isTRUE(ncs_ci_cfg$enabled)
      ncs_ci_boot <- if (!is.null(ncs_ci_cfg$n_boot)) ncs_ci_cfg$n_boot else 0L
      ncs_ci_frac <- if (!is.null(ncs_ci_cfg$sample_frac)) ncs_ci_cfg$sample_frac else 0.8
      # SNP probes from minfi
      if ("snp" %in% ncs_types && !is.null(rgSet) && requireNamespace("minfi", quietly = TRUE)) {
        snp_beta <- tryCatch(minfi::getSnpBeta(rgSet), error = function(e) NULL)
        if (!is.null(snp_beta)) {
          common_samples <- intersect(colnames(snp_beta), colnames(betas_raw))
          snp_beta <- snp_beta[, common_samples, drop = FALSE]
          if (ncol(snp_beta) >= 4) {
            design_nc <- tryCatch(model.matrix(~ 0 + primary_group, data = targets[common_samples, , drop = FALSE]),
                                  error = function(e) NULL)
              if (!is.null(design_nc) && ncol(design_nc) >= 2) {
                colnames(design_nc) <- make.names(gsub("primary_group", "", colnames(design_nc)))
                group_cols <- make.names(levels(targets$primary_group))
                stats_raw <- compute_negative_control_stats(snp_beta, design_nc, group_cols,
                                                            lambda_ci = ncs_ci_enabled,
                                                            n_boot = ncs_ci_boot, sample_frac = ncs_ci_frac)
                snp_corr <- snp_beta
                if (nzchar(applied_batch_method) && applied_batch_method != "none" &&
                    !is.null(batch_col) && batch_col %in% colnames(targets)) {
                  bc_snp <- apply_batch_correction(snp_beta, applied_batch_method, batch_col, covariates,
                                                 targets[common_samples, , drop = FALSE], design_nc)
                  snp_corr <- bc_snp$betas
                }
                stats_corr <- compute_negative_control_stats(snp_corr, design_nc, group_cols,
                                                             lambda_ci = ncs_ci_enabled,
                                                             n_boot = ncs_ci_boot, sample_frac = ncs_ci_frac)
                if (!is.null(stats_raw)) {
                  ncs_summary <- rbind(ncs_summary, data.frame(type = "snp", stage = "raw",
                                                              lambda = stats_raw$lambda, sig_rate = stats_raw$sig_rate,
                                                              n = stats_raw$n,
                                                              lambda_ci_low = stats_raw$lambda_ci_low,
                                                              lambda_ci_high = stats_raw$lambda_ci_high,
                                                              lambda_boot_n = stats_raw$lambda_boot_n))
                }
                if (!is.null(stats_corr)) {
                  ncs_summary <- rbind(ncs_summary, data.frame(type = "snp", stage = "corrected",
                                                              lambda = stats_corr$lambda, sig_rate = stats_corr$sig_rate,
                                                              n = stats_corr$n,
                                                              lambda_ci_low = stats_corr$lambda_ci_low,
                                                              lambda_ci_high = stats_corr$lambda_ci_high,
                                                              lambda_boot_n = stats_corr$lambda_boot_n))
                }
              }
            }
        }
      }
      # Housekeeping probes if provided
      if ("housekeeping" %in% ncs_types) {
        hk_path <- if (!is.null(crf_cfg$housekeeping_path)) as.character(crf_cfg$housekeeping_path) else ""
        if (nzchar(hk_path) && !grepl("^/", hk_path)) hk_path <- file.path(project_dir, hk_path)
        if (nzchar(hk_path) && file.exists(hk_path)) {
          hk_ids <- load_probe_list_from_file(hk_path, col_candidates = c("CpG", "cpg", "probe"))
          hk_ids <- intersect(hk_ids, rownames(betas_raw))
          if (length(hk_ids) >= 10) {
            hk_raw <- betas_raw[hk_ids, , drop = FALSE]
            hk_corr <- betas_corr[hk_ids, , drop = FALSE]
            design_nc <- tryCatch(model.matrix(~ 0 + primary_group, data = targets), error = function(e) NULL)
            if (!is.null(design_nc) && ncol(design_nc) >= 2) {
              colnames(design_nc) <- make.names(gsub("primary_group", "", colnames(design_nc)))
              group_cols <- make.names(levels(targets$primary_group))
              stats_raw <- compute_negative_control_stats(hk_raw, design_nc, group_cols,
                                                          lambda_ci = ncs_ci_enabled,
                                                          n_boot = ncs_ci_boot, sample_frac = ncs_ci_frac)
              stats_corr <- compute_negative_control_stats(hk_corr, design_nc, group_cols,
                                                           lambda_ci = ncs_ci_enabled,
                                                           n_boot = ncs_ci_boot, sample_frac = ncs_ci_frac)
              if (!is.null(stats_raw)) {
                ncs_summary <- rbind(ncs_summary, data.frame(type = "housekeeping", stage = "raw",
                                                            lambda = stats_raw$lambda, sig_rate = stats_raw$sig_rate,
                                                            n = stats_raw$n,
                                                            lambda_ci_low = stats_raw$lambda_ci_low,
                                                            lambda_ci_high = stats_raw$lambda_ci_high,
                                                            lambda_boot_n = stats_raw$lambda_boot_n))
              }
              if (!is.null(stats_corr)) {
                ncs_summary <- rbind(ncs_summary, data.frame(type = "housekeeping", stage = "corrected",
                                                            lambda = stats_corr$lambda, sig_rate = stats_corr$sig_rate,
                                                            n = stats_corr$n,
                                                            lambda_ci_low = stats_corr$lambda_ci_low,
                                                            lambda_ci_high = stats_corr$lambda_ci_high,
                                                            lambda_boot_n = stats_corr$lambda_boot_n))
              }
            }
          }
        }
      }
    }
  }
  if (nrow(ncs_summary) > 0) {
    # Convert NCS lambda to a continuous score in [0, 1] to make it usable for CAF/Verdict.
    # Score is tier-adaptive: larger tolerance for minimal tiers, stricter for larger cohorts.
    ncs_tol <- 0.2
    if (!is.null(tier_info$tier)) {
      if (tier_info$tier == "minimal") ncs_tol <- 0.3
      if (tier_info$tier == "small") ncs_tol <- 0.2
      if (tier_info$tier == "moderate") ncs_tol <- 0.15
      if (tier_info$tier == "large") ncs_tol <- 0.1
    }
    ncs_summary$ncs_score_tol <- ncs_tol
    ncs_summary$ncs_score <- ifelse(is.finite(ncs_summary$lambda),
                                    pmax(0, pmin(1, 1 - abs(ncs_summary$lambda - 1.0) / ncs_tol)),
                                    NA_real_)
    write.csv(ncs_summary, file.path(out_dir, paste0(prefix, "_CRF_NCS_Summary.csv")), row.names = FALSE)
  }

  rss_summary <- data.frame()
  if (!is.null(betas_raw) && !is.null(betas_corr)) {
    top_k <- tier_info$rss_top_k
    mode <- tier_info$rss_mode
    n_iter <- tier_info$rss_iterations
    max_probes <- suppressWarnings(as.integer(crf_cfg$max_probes_rss))
    if (!is.finite(max_probes)) max_probes <- 0
    covariates_use <- covariates
    rss_raw <- compute_rss_overlap(betas_raw, targets, "primary_group", covariates_use,
                                   mode = mode, n_iter = n_iter, top_k = top_k, max_probes = max_probes)
    rss_corr <- compute_rss_overlap(betas_corr, targets, "primary_group", covariates_use,
                                    mode = mode, n_iter = n_iter, top_k = top_k, max_probes = max_probes)
    if (!is.null(rss_raw) || !is.null(rss_corr)) {
      rss_summary <- data.frame(top_k = top_k,
                                overlap_mean_raw = if (!is.null(rss_raw)) rss_raw$overlap_mean else NA_real_,
                                overlap_sd_raw = if (!is.null(rss_raw)) rss_raw$overlap_sd else NA_real_,
                                jaccard_mean_raw = if (!is.null(rss_raw)) rss_raw$jaccard_mean else NA_real_,
                                jaccard_sd_raw = if (!is.null(rss_raw)) rss_raw$jaccard_sd else NA_real_,
                                rbo_mean_raw = if (!is.null(rss_raw)) rss_raw$rbo_mean else NA_real_,
                                rbo_sd_raw = if (!is.null(rss_raw)) rss_raw$rbo_sd else NA_real_,
                                rss_mean_raw = if (!is.null(rss_raw)) rss_raw$rss_mean else NA_real_,
                                rss_sd_raw = if (!is.null(rss_raw)) rss_raw$rss_sd else NA_real_,
                                overlap_mean_corr = if (!is.null(rss_corr)) rss_corr$overlap_mean else NA_real_,
                                overlap_sd_corr = if (!is.null(rss_corr)) rss_corr$overlap_sd else NA_real_,
                                jaccard_mean_corr = if (!is.null(rss_corr)) rss_corr$jaccard_mean else NA_real_,
                                jaccard_sd_corr = if (!is.null(rss_corr)) rss_corr$jaccard_sd else NA_real_,
                                rbo_mean_corr = if (!is.null(rss_corr)) rss_corr$rbo_mean else NA_real_,
                                rbo_sd_corr = if (!is.null(rss_corr)) rss_corr$rbo_sd else NA_real_,
                                rss_mean_corr = if (!is.null(rss_corr)) rss_corr$rss_mean else NA_real_,
                                rss_sd_corr = if (!is.null(rss_corr)) rss_corr$rss_sd else NA_real_,
                                rbo_p = if (!is.null(rss_corr)) rss_corr$rbo_p else if (!is.null(rss_raw)) rss_raw$rbo_p else NA_real_,
                                sign_mean_raw = if (!is.null(rss_raw)) rss_raw$sign_mean else NA_real_,
                                sign_sd_raw = if (!is.null(rss_raw)) rss_raw$sign_sd else NA_real_,
                                sign_mean_corr = if (!is.null(rss_corr)) rss_corr$sign_mean else NA_real_,
                                sign_sd_corr = if (!is.null(rss_corr)) rss_corr$sign_sd else NA_real_)
      write.csv(rss_summary, file.path(out_dir, paste0(prefix, "_CRF_RSS_Summary.csv")), row.names = FALSE)
    }
  }

  report_lines <- build_crf_report_lines(tier_info, mmc_summary, ncs_summary, rss_summary, prefix,
                                         pvca_ci_before = pvca_ci_before,
                                         pvca_ci_after = pvca_ci_after)
  report_path <- file.path(out_dir, paste0(prefix, "_CRF_Report.txt"))
  writeLines(report_lines, report_path)

  ncs_best_score <- NA_real_
  if (nrow(ncs_summary) > 0 && "ncs_score" %in% colnames(ncs_summary)) {
    pick <- ncs_summary
    if ("stage" %in% colnames(pick)) {
      corr <- pick[tolower(pick$stage) == "corrected", , drop = FALSE]
      if (nrow(corr) > 0) pick <- corr
    }
    if ("type" %in% colnames(pick)) {
      snp <- pick[tolower(pick$type) == "snp", , drop = FALSE]
      if (nrow(snp) > 0) pick <- snp
    }
    if ("n" %in% colnames(pick)) {
      pick <- pick[order(-suppressWarnings(as.numeric(pick$n))), , drop = FALSE]
    }
    ncs_best_score <- suppressWarnings(as.numeric(pick$ncs_score[1]))
  }

  list(
    report_path = report_path,
    mmc_summary = mmc_summary,
    ncs_summary = ncs_summary,
    rss_summary = rss_summary,
    ncs_best_score = ncs_best_score,
    pvca_ci_before = pvca_ci_before,
    pvca_ci_after = pvca_ci_after
  )
}

marker_list <- load_marker_list(marker_list_path)
if (!is.null(marker_list)) {
  message(sprintf("Loaded marker list with %d CpGs for preservation checks.", length(marker_list)))
}

auto_detect_covariates <- function(targets, gsm_col, pca_scores, alpha=AUTO_COVARIATE_ALPHA, max_pcs=MAX_PCS_FOR_COVARIATE_DETECTION) {
  # Use PC-covariate association to decide covariates to include
  meta_cols <- setdiff(colnames(targets), c("primary_group", "Basename", "filenames", gsm_col))
  keep <- c()
  log_df <- data.frame(Variable=character(0), MinP=numeric(0), MinPC=integer(0), Type=character(0),
                       stringsAsFactors = FALSE)
  
  if (length(meta_cols) == 0 || ncol(pca_scores) == 0) {
    return(list(selected = keep, log = log_df))
  }
  
  # Limit PCs checked
  pca_use <- pca_scores[, 1:min(max_pcs, ncol(pca_scores)), drop=FALSE]
  
  for (col in meta_cols) {
    val <- normalize_meta_vals(targets[[col]])
    non_missing <- val[!is.na(val)]
    if (length(unique(non_missing)) < 2) next
    cc <- !is.na(val)
    if (sum(cc) < 3) next
    val_cc <- val[cc]
    if (length(unique(val_cc)) < 2) next
    pca_cc <- pca_use[cc, , drop=FALSE]

    vtype <- if (is.numeric(val_cc)) "numeric" else "categorical"
    
    # Skip if almost perfectly collinear with group (chi-square check for factors)
    if (is.character(val_cc) || is.factor(val_cc)) {
      tbl <- table(val_cc, targets$primary_group[cc])
      p_chisq <- safe_chisq_p(tbl)
      if (!is.na(p_chisq) && p_chisq < 1e-6) next
    }
    
    p_vec <- rep(NA_real_, ncol(pca_cc))
    for (k in seq_len(ncol(pca_cc))) {
      if (is.numeric(val_cc)) {
        fit <- lm(pca_cc[,k] ~ val_cc)
        coefs <- summary(fit)$coefficients
        if (nrow(coefs) >= 2) p_vec[k] <- coefs[2,4]
      } else {
        p_vec[k] <- safe_kruskal_p(pca_cc[,k], val_cc)
      }
    }
    
    min_p <- suppressWarnings(min(p_vec, na.rm = TRUE))
    if (!is.finite(min_p)) next
    min_pc <- if (all(is.na(p_vec))) NA_integer_ else which.min(p_vec)[1]

    log_df <- rbind(log_df, data.frame(Variable = col, MinP = min_p, MinPC = min_pc, Type = vtype))
    # Bonferroni-correct alpha over the number of PCs tested to control
    # the per-variable false selection rate (min-p across PCs is anti-conservative).
    n_pcs_tested <- sum(is.finite(p_vec))
    adj_alpha <- if (n_pcs_tested > 1) alpha / n_pcs_tested else alpha
    if (min_p < adj_alpha) keep <- c(keep, col)
  }
  
  return(list(selected = keep, log = log_df))
}

safe_chisq_p <- function(tbl) {
  if (is.null(tbl) || any(dim(tbl) == 0)) return(NA_real_)
  n <- suppressWarnings(sum(tbl))
  if (!is.finite(n) || n <= 0) return(NA_real_)
  tryCatch(suppressWarnings(chisq.test(tbl)$p.value), error = function(e) NA_real_)
}

normalize_meta_vals <- function(val) {
  if (is.character(val)) {
    val <- trimws(val)
    val[val == ""] <- NA
  }
  val
}

safe_kruskal_p <- function(x, g, min_total = 3) {
  g <- normalize_meta_vals(g)
  cc <- complete.cases(x, g)
  if (sum(cc) < min_total) return(NA_real_)
  x <- x[cc]
  g <- g[cc]
  if (length(unique(g)) < 2) return(NA_real_)
  tryCatch(kruskal.test(x ~ as.factor(g))$p.value, error = function(e) NA_real_)
}

safe_ttest_p <- function(x, g, min_per_group = 2) {
  g <- normalize_meta_vals(g)
  cc <- complete.cases(x, g)
  if (sum(cc) < (min_per_group * 2)) return(NA_real_)
  x <- x[cc]
  g <- g[cc]
  if (length(unique(g)) != 2) return(NA_real_)
  tab <- table(g)
  if (any(tab < min_per_group)) return(NA_real_)
  tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
}

cramers_v <- function(tbl) {
  if (any(dim(tbl) == 0)) return(NA_real_)
  n <- suppressWarnings(sum(tbl))
  if (!is.finite(n) || n <= 0) return(NA_real_)
  chi2 <- tryCatch(suppressWarnings(chisq.test(tbl)$statistic), error = function(e) NA_real_)
  if (!is.finite(chi2)) return(NA_real_)
  k <- min(nrow(tbl) - 1, ncol(tbl) - 1)
  if (k <= 0) return(NA_real_)
  sqrt(chi2 / (n * k))
}

group_assoc_stats <- function(val, group_vec) {
  if (length(val) != length(group_vec)) return(list(p = NA_real_, eta2 = NA_real_, method = ""))
  if (is.numeric(val)) {
    df <- data.frame(val = val, group = as.factor(group_vec))
    df <- df[is.finite(df$val) & !is.na(df$group), , drop = FALSE]
    if (nrow(df) < 3 || length(unique(df$group)) < 2) return(list(p = NA_real_, eta2 = NA_real_, method = ""))
    fit <- tryCatch(lm(val ~ group, data = df), error = function(e) NULL)
    if (is.null(fit)) return(list(p = NA_real_, eta2 = NA_real_, method = ""))
    an <- tryCatch(anova(fit), error = function(e) NULL)
    if (is.null(an) || nrow(an) < 1) return(list(p = NA_real_, eta2 = NA_real_, method = ""))
    ss_total <- sum(an$`Sum Sq`, na.rm = TRUE)
    ss_group <- an$`Sum Sq`[1]
    eta2 <- if (is.finite(ss_total) && ss_total > 0) ss_group / ss_total else NA_real_
    pval <- an$`Pr(>F)`[1]
    return(list(p = pval, eta2 = eta2, method = "anova"))
  }
  tbl <- table(as.factor(val), as.factor(group_vec))
  pval <- safe_chisq_p(tbl)
  v <- cramers_v(tbl)
  list(p = pval, eta2 = v, method = "chisq")
}

map_mm_columns_to_vars <- function(mm_cols, vars) {
  vars_clean <- make.names(vars)
  ord <- order(nchar(vars_clean), decreasing = TRUE)
  vars_clean <- vars_clean[ord]
  vars_orig <- vars[ord]
  col_map <- setNames(rep("", length(mm_cols)), mm_cols)
  for (i in seq_along(mm_cols)) {
    col <- mm_cols[i]
    matched <- ""
    for (j in seq_along(vars_clean)) {
      if (startsWith(col, vars_clean[j])) {
        matched <- vars_orig[j]
        break
      }
    }
    col_map[i] <- matched
  }
  col_map
}

compute_covariate_collinearity <- function(targets, vars) {
  vars <- intersect(vars, colnames(targets))
  if (length(vars) < 2) {
    return(data.frame(Variable = vars, Max_Correlation = NA_real_, Max_Correlation_With = "", stringsAsFactors = FALSE))
  }
  meta <- targets[, vars, drop = FALSE]
  for (v in vars) {
    if (is.character(meta[[v]])) meta[[v]] <- as.factor(meta[[v]])
  }
  mm <- tryCatch(model.matrix(~ . -1, data = meta), error = function(e) NULL)
  if (is.null(mm) || ncol(mm) < 2) {
    return(data.frame(Variable = vars, Max_Correlation = NA_real_, Max_Correlation_With = "", stringsAsFactors = FALSE))
  }
  cor_mat <- suppressWarnings(cor(mm, use = "pairwise.complete.obs"))
  mm_cols <- colnames(mm)
  col_to_var <- map_mm_columns_to_vars(mm_cols, vars)
  out <- data.frame(Variable = vars, Max_Correlation = NA_real_, Max_Correlation_With = "",
                    stringsAsFactors = FALSE)
  for (i in seq_along(vars)) {
    v <- vars[i]
    cols_v <- mm_cols[col_to_var == v]
    cols_other <- mm_cols[col_to_var != v & nzchar(col_to_var)]
    if (length(cols_v) == 0 || length(cols_other) == 0) next
    sub_cor <- cor_mat[cols_v, cols_other, drop = FALSE]
    if (length(sub_cor) == 0) next
    idx <- arrayInd(which.max(abs(sub_cor)), dim(sub_cor))
    max_val <- sub_cor[idx[1], idx[2]]
    out$Max_Correlation[i] <- max_val
    col_other <- cols_other[idx[2]]
    if (col_other %in% names(col_to_var)) {
      out$Max_Correlation_With[i] <- col_to_var[[col_other]]
    }
  }
  out
}

build_covariate_candidate_report <- function(targets, covar_log_df, group_col,
                                             group_assoc_p_threshold = 1e-3) {
  if (is.null(covar_log_df) || nrow(covar_log_df) == 0) return(covar_log_df)
  vars <- covar_log_df$Variable
  group_vec <- targets[[group_col]]
  assoc_p <- numeric(length(vars))
  assoc_eta2 <- numeric(length(vars))
  assoc_method <- character(length(vars))
  for (i in seq_along(vars)) {
    v <- vars[i]
    stats <- group_assoc_stats(targets[[v]], group_vec)
    assoc_p[i] <- stats$p
    assoc_eta2[i] <- stats$eta2
    assoc_method[i] <- stats$method
  }
  collin <- compute_covariate_collinearity(targets, vars)
  report <- covar_log_df
  report$Group_Assoc_P <- assoc_p
  report$Group_Assoc_Eta2 <- assoc_eta2
  report$Group_Assoc_Method <- assoc_method
  report$Group_Assoc_Flag <- ifelse(is.finite(assoc_p) & assoc_p < group_assoc_p_threshold, TRUE, FALSE)
  if (!is.null(collin) && nrow(collin) > 0) {
    idx <- match(report$Variable, collin$Variable)
    report$Max_Correlation <- collin$Max_Correlation[idx]
    report$Max_Correlation_With <- collin$Max_Correlation_With[idx]
  } else {
    report$Max_Correlation <- NA_real_
    report$Max_Correlation_With <- ""
  }
  report
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
    val <- normalize_meta_vals(targets[[cv]])
    non_missing <- val[!is.na(val)]
    uniq <- length(unique(non_missing))
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
      lvl_counts <- table(non_missing)
      if (any(lvl_counts < 2)) {
        drops <- rbind(drops, data.frame(Variable=cv, Reason="rare_level"))
        next
      }
    }
    # Check confounding with group (any level exclusive to one group or strong chi-square)
    if (!is_num && (is.character(val) || is.factor(val))) {
      cc <- complete.cases(val, group_vals)
      if (sum(cc) < 3) next
      tbl <- table(val[cc], group_vals[cc])
      if (any(tbl == 0)) {
        drops <- rbind(drops, data.frame(Variable=cv, Reason="confounded_with_group"))
        next
      }
      p_chisq <- safe_chisq_p(tbl)
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
          p_grp <- if (length(unique(gv)) == 2) safe_ttest_p(vv, gv) else safe_kruskal_p(vv, gv)
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

drop_single_level_covariates <- function(targets, covariates) {
  keep <- character(0)
  drops <- data.frame(Variable=character(0), Reason=character(0), stringsAsFactors = FALSE)
  for (cv in covariates) {
    if (!(cv %in% colnames(targets))) next
    vals <- targets[[cv]]
    if (all(is.na(vals)) || (is.character(vals) && all(trimws(vals) == ""))) {
      drops <- rbind(drops, data.frame(Variable = cv, Reason = "all_missing"))
      next
    }
    if (is.numeric(vals) || is.integer(vals)) {
      vv <- vals[is.finite(vals)]
      if (length(vv) < 2 || length(unique(vv)) < 2) {
        drops <- rbind(drops, data.frame(Variable = cv, Reason = "constant_or_single_level"))
        next
      }
    } else {
      f <- droplevels(as.factor(vals))
      if (length(levels(f)) < 2) {
        drops <- rbind(drops, data.frame(Variable = cv, Reason = "single_level"))
        next
      }
    }
    keep <- c(keep, cv)
  }
  list(keep = keep, dropped = drops)
}

filter_covariates_with_batch <- function(targets, covariates, batch_col, p_thresh = 1e-6) {
  keep <- c()
  drops <- data.frame(Variable=character(0), Reason=character(0), stringsAsFactors = FALSE)
  covariates <- covariates[covariates %in% colnames(targets)]
  if (is.null(batch_col) || !(batch_col %in% colnames(targets))) {
    return(list(keep = covariates, dropped = drops))
  }
  batch_vals <- targets[[batch_col]]
  for (cv in covariates) {
    val <- targets[[cv]]
    if (is.numeric(val) || is.integer(val)) {
      df <- data.frame(val = val, batch = batch_vals)
      df <- df[complete.cases(df), , drop = FALSE]
      if (nrow(df) < 3) {
        keep <- c(keep, cv)
        next
      }
      by_batch <- split(df$val, df$batch)
      within_var <- sapply(by_batch, function(x) {
        if (length(unique(x)) < 2) 0 else var(x)
      })
      if (length(within_var) > 0 && all(within_var == 0) &&
          length(unique(df$batch)) > 1 && length(unique(df$val)) > 1) {
        drops <- rbind(drops, data.frame(Variable = cv, Reason = "confounded_with_batch"))
        next
      }
      fit <- tryCatch(lm(val ~ batch, data = df), error = function(e) NULL)
      if (!is.null(fit)) {
        r2 <- summary(fit)$r.squared
        if (is.finite(r2) && r2 > 0.999) {
          drops <- rbind(drops, data.frame(Variable = cv, Reason = "confounded_with_batch"))
          next
        }
      }
    } else {
      if (is_batch_confounded(targets, batch_col, cv, p_thresh = p_thresh)) {
        drops <- rbind(drops, data.frame(Variable = cv, Reason = "confounded_with_batch"))
        next
      }
    }
    keep <- c(keep, cv)
  }
  list(keep = keep, dropped = drops)
}

scale_numeric_covariates <- function(meta, cols = NULL, exclude = character(0)) {
  out <- meta
  scaled <- character(0)
  if (is.null(cols)) cols <- colnames(out)
  cols <- intersect(cols, colnames(out))
  cols <- setdiff(cols, exclude)
  for (nm in cols) {
    vals <- out[[nm]]
    if (!is.numeric(vals) && !is.integer(vals)) next
    vv <- vals[is.finite(vals)]
    if (length(vv) < 2) next
    sd_val <- sd(vv, na.rm = TRUE)
    if (!is.finite(sd_val) || sd_val == 0) next
    out[[nm]] <- as.numeric(scale(vals))
    scaled <- c(scaled, nm)
  }
  list(meta = out, scaled = scaled)
}

normalize_vp_config <- function(vp_cfg) {
  if (is.null(vp_cfg) || length(vp_cfg) == 0) {
    return(list(autoscale_numeric = FALSE, autoscale_on_fail = FALSE))
  }
  if (is.null(vp_cfg$autoscale_numeric)) vp_cfg$autoscale_numeric <- FALSE
  if (is.null(vp_cfg$autoscale_on_fail)) vp_cfg$autoscale_on_fail <- FALSE
  vp_cfg
}

fit_variance_partition_safe <- function(M_mat, form, meta, vp_cfg = NULL, label = "variancePartition") {
  vp_cfg <- normalize_vp_config(vp_cfg)
  vars <- all.vars(form)
  meta_use <- meta
  scaled_cols <- character(0)
  if (isTRUE(vp_cfg$autoscale_numeric)) {
    scale_res <- scale_numeric_covariates(meta_use, cols = vars)
    meta_use <- scale_res$meta
    scaled_cols <- scale_res$scaled
    if (length(scaled_cols) > 0) {
      message(sprintf("    - %s: scaled numeric covariates (%s)", label, paste(scaled_cols, collapse = ", ")))
    }
  }
  run_fit <- function(meta_in) {
    with_muffled_conditions(
      fitExtractVarPartModel(M_mat, form, meta_in),
      label = label,
      patterns = "boundary \\(singular\\) fit"
    )
  }
  vp_res <- tryCatch(run_fit(meta_use), error = function(e) {
    message(sprintf("    - %s failed: %s", label, e$message))
    if (isTRUE(vp_cfg$autoscale_on_fail) && !isTRUE(vp_cfg$autoscale_numeric)) {
      scale_res <- scale_numeric_covariates(meta, cols = vars)
      if (length(scale_res$scaled) > 0) {
        message(sprintf("    - %s retry with scaled covariates: %s", label, paste(scale_res$scaled, collapse = ", ")))
        return(tryCatch(run_fit(scale_res$meta), error = function(e2) {
          message(sprintf("    - %s failed after autoscale: %s", label, e2$message))
          NULL
        }))
      }
    }
    NULL
  })
  list(result = vp_res, scaled_cols = scaled_cols)
}

build_combat_design <- function(targets, group_col, covariates, batch_col, p_thresh = 1e-6) {
  covariates <- covariates[covariates %in% colnames(targets)]
  drop_log <- data.frame(Variable=character(0), Reason=character(0), stringsAsFactors = FALSE)
  if (length(covariates) > 0) {
    single_filt <- drop_single_level_covariates(targets, covariates)
    covariates <- single_filt$keep
    if (nrow(single_filt$dropped) > 0) drop_log <- rbind(drop_log, single_filt$dropped)
  }
  if (length(covariates) > 0 && !is.null(batch_col) && batch_col %in% colnames(targets)) {
    batch_filter <- filter_covariates_with_batch(targets, covariates, batch_col, p_thresh = p_thresh)
    covariates <- batch_filter$keep
    if (nrow(batch_filter$dropped) > 0) drop_log <- rbind(drop_log, batch_filter$dropped)
  }
  formula_str <- paste("~ 0 +", group_col)
  if (length(covariates) > 0) {
    formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  }
  mod <- model.matrix(as.formula(formula_str), data = targets)
  colnames(mod) <- gsub(group_col, "", colnames(mod), fixed = TRUE)
  group_cols <- make.names(levels(targets[[group_col]]))
  mod_ld <- drop_linear_dependencies(mod, group_cols = group_cols)
  mod <- mod_ld$mat
  if (length(mod_ld$dropped) > 0) {
    drop_log <- rbind(drop_log, data.frame(Variable = mod_ld$dropped, Reason = "linear_dependency"))
  }
  list(mod = mod, covariates = covariates, dropped = drop_log)
}

apply_covariate_governance <- function(targets, covariates, forced_covariates, config_settings, drop_log) {
  forced_cfg <- config_settings$forced_adjust
  allowed_cfg <- config_settings$allowed_adjust
  do_not_cfg <- config_settings$do_not_adjust
  if (length(forced_cfg) > 0) {
    present_force <- intersect(forced_cfg, colnames(targets))
    missing_force <- setdiff(forced_cfg, present_force)
    if (length(missing_force) > 0) {
      drop_log <- rbind(drop_log, data.frame(Variable = missing_force, Reason = "forced_not_in_config"))
    }
    covariates <- unique(c(covariates, present_force))
    forced_covariates <- unique(c(forced_covariates, present_force))
  }
  if (length(allowed_cfg) > 0) {
    covariates <- intersect(covariates, allowed_cfg)
  }
  if (length(do_not_cfg) > 0) {
    blocked <- intersect(covariates, do_not_cfg)
    if (length(blocked) > 0) {
      drop_log <- rbind(drop_log, data.frame(Variable = blocked, Reason = "do_not_adjust"))
      covariates <- setdiff(covariates, blocked)
    }
  }
  max_frac <- config_settings$missingness$max_frac
  impute_frac <- config_settings$missingness$impute_max_frac
  for (cv in covariates) {
    vals <- targets[[cv]]
    miss_mask <- is.na(vals) | (is.character(vals) & trimws(vals) == "")
    miss_frac <- mean(miss_mask)
    if (!is.finite(miss_frac)) miss_frac <- 0
    if (miss_frac > max_frac) {
      drop_log <- rbind(drop_log, data.frame(Variable = cv, Reason = "missingness_high"))
      covariates <- setdiff(covariates, cv)
      forced_covariates <- setdiff(forced_covariates, cv)
      next
    }
    if (miss_frac > 0 && miss_frac <= impute_frac) {
      if (is.numeric(vals) || is.integer(vals)) {
        imp <- median(vals, na.rm = TRUE)
        vals[miss_mask] <- imp
      } else {
        vals <- as.character(vals)
        vals[miss_mask] <- "Unknown"
        vals <- as.factor(vals)
      }
      targets[[cv]] <- vals
      drop_log <- rbind(drop_log, data.frame(Variable = cv, Reason = "imputed_missing"))
    } else if (miss_frac > 0) {
      drop_log <- rbind(drop_log, data.frame(Variable = cv, Reason = "missingness_not_imputed"))
      covariates <- setdiff(covariates, cv)
      forced_covariates <- setdiff(forced_covariates, cv)
    }
  }
  num_covs <- covariates[vapply(covariates, function(x) {
    vals <- targets[[x]]
    is.numeric(vals) || is.integer(vals)
  }, logical(1))]
  cor_thresh <- config_settings$collinearity$cor_threshold
  if (length(num_covs) >= 2) {
    num_mat <- targets[, num_covs, drop = FALSE]
    cor_mat <- suppressWarnings(cor(num_mat, use = "pairwise.complete.obs"))
    if (is.matrix(cor_mat)) {
      drop_set <- character(0)
      for (i in seq_len(ncol(cor_mat) - 1)) {
        for (j in (i + 1):ncol(cor_mat)) {
          if (!is.finite(cor_mat[i, j])) next
          if (abs(cor_mat[i, j]) >= cor_thresh) {
            a <- colnames(cor_mat)[i]
            b <- colnames(cor_mat)[j]
            drop_choice <- if (a %in% forced_covariates && !(b %in% forced_covariates)) b else a
            drop_set <- unique(c(drop_set, drop_choice))
          }
        }
      }
      if (length(drop_set) > 0) {
        drop_log <- rbind(drop_log, data.frame(Variable = drop_set, Reason = "collinearity"))
        covariates <- setdiff(covariates, drop_set)
        forced_covariates <- setdiff(forced_covariates, drop_set)
      }
    }
  }
  list(targets = targets, covariates = covariates, forced_covariates = forced_covariates, drop_log = drop_log)
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

generate_covariate_sets <- function(targets, covariates, forced_covariates, covar_log_df) {
  sets <- list(full = covariates)
  if (length(covariates) >= 4) {
    ranked <- rank_covariates(targets = targets, covariates = covariates, covar_log_df = covar_log_df)
    keep_n <- max(1, floor(length(ranked) / 2))
    reduced <- unique(c(forced_covariates, head(ranked, keep_n)))
    sets$reduced <- reduced
  }
  if (length(forced_covariates) > 0) {
    sets$forced <- forced_covariates
  }
  sets
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
    val <- normalize_meta_vals(targets[[sv]])
    cc <- complete.cases(val, group_vals)
    if (sum(cc) < 3) next
    gv <- group_vals[cc]
    if (length(unique(gv)) < 2) next
    vv <- val[cc]
    if (sd(vv, na.rm = TRUE) == 0) {
      drops <- rbind(drops, data.frame(Variable = sv, Reason = "sv_constant"))
      next
    }
    p_grp <- if (length(unique(gv)) == 2) safe_ttest_p(vv, gv) else safe_kruskal_p(vv, gv)
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

normalize_epicv2_ids <- function(betas_in) {
  betas_in <- as.matrix(betas_in)
  storage.mode(betas_in) <- "numeric"
  ids <- rownames(betas_in)
  if (is.null(ids)) {
    return(list(betas = betas_in, normalized = FALSE, collapsed = FALSE, dup_count = 0))
  }
  base_ids <- sub("_.+$", "", ids)
  valid_base <- grepl("^(cg|ch)[0-9]+$", base_ids) | grepl("^ch\\.[0-9]+$", base_ids) | grepl("^rs[0-9]+$", base_ids)
  has_suffix <- ids != base_ids & valid_base
  if (!any(has_suffix)) {
    return(list(betas = betas_in, normalized = FALSE, collapsed = FALSE, dup_count = 0))
  }
  base_ids[!has_suffix] <- ids[!has_suffix]
  dup_count <- sum(duplicated(base_ids))
  if (dup_count == 0) {
    rownames(betas_in) <- base_ids
    return(list(betas = betas_in, normalized = TRUE, collapsed = FALSE, dup_count = 0))
  }
  collapse_res <- tryCatch({
    na_mask <- is.na(betas_in)
    betas_zero <- betas_in
    betas_zero[na_mask] <- 0
    sums <- rowsum(betas_zero, group = base_ids, reorder = FALSE)
    counts <- rowsum((!na_mask) * 1, group = base_ids, reorder = FALSE)
    betas_out <- sums / counts
    betas_out[counts == 0] <- NA_real_
    list(betas = betas_out)
  }, error = function(e) e)
  if (inherits(collapse_res, "error")) {
    keep_idx <- !duplicated(base_ids)
    betas_out <- betas_in[keep_idx, , drop = FALSE]
    rownames(betas_out) <- base_ids[keep_idx]
    return(list(betas = betas_out, normalized = TRUE, collapsed = FALSE, dup_count = dup_count))
  }
  list(betas = collapse_res$betas, normalized = TRUE, collapsed = TRUE, dup_count = dup_count)
}

impute_row_means <- function(mat) {
  if (!anyNA(mat)) return(list(mat = mat, imputed = 0L))
  row_means <- rowMeans(mat, na.rm = TRUE)
  na_idx <- which(is.na(mat), arr.ind = TRUE)
  mat[na_idx] <- row_means[na_idx[, 1]]
  list(mat = mat, imputed = nrow(na_idx))
}

impute_knn <- function(mat, na_max_frac = SESAME_NATIVE_NA_MAX_FRAC,
                       k = SESAME_NATIVE_KNN_K,
                       max_rows = SESAME_NATIVE_KNN_MAX_ROWS,
                       ref_rows = SESAME_NATIVE_KNN_REF_ROWS) {
  if (!anyNA(mat)) return(list(mat = mat, imputed = 0L, method = "none"))
  if (!requireNamespace("impute", quietly = TRUE)) {
    return(list(mat = mat, imputed = 0L, method = "missing_impute_pkg"))
  }

  na_row_counts <- rowSums(is.na(mat))
  missing_rows <- which(na_row_counts > 0)
  if (length(missing_rows) == 0) return(list(mat = mat, imputed = 0L, method = "none"))
  if (length(missing_rows) > max_rows) {
    return(list(mat = mat, imputed = 0L, method = "too_many_missing_rows"))
  }

  complete_rows <- which(na_row_counts == 0)
  if (length(complete_rows) == 0) {
    return(list(mat = mat, imputed = 0L, method = "no_complete_rows"))
  }
  ref_take <- min(length(complete_rows), ref_rows)
  ref_rows <- head(complete_rows, ref_take)
  work_rows <- c(missing_rows, ref_rows)
  work_mat <- mat[work_rows, , drop = FALSE]
  k_use <- min(k, nrow(work_mat) - 1L)
  if (k_use < 1) {
    return(list(mat = mat, imputed = 0L, method = "k_too_small"))
  }

  imp <- tryCatch(
    with_muffled_conditions(
      impute::impute.knn(work_mat, k = k_use, rowmax = na_max_frac, colmax = 1),
      label = NULL,
      patterns = "no finite values|missing values|all\\s+NA"
    ),
    error = function(e) e
  )
  if (inherits(imp, "error")) {
    return(list(mat = mat, imputed = 0L, method = "knn_failed"))
  }
  out <- mat
  out[missing_rows, ] <- imp$data[seq_along(missing_rows), , drop = FALSE]
  before_na <- sum(is.na(mat[missing_rows, , drop = FALSE]))
  after_na <- sum(is.na(out[missing_rows, , drop = FALSE]))
  list(mat = out, imputed = max(0L, before_na - after_na), method = "knn")
}

prepare_sesame_native <- function(betas, na_max_frac = SESAME_NATIVE_NA_MAX_FRAC) {
  if (is.null(betas) || nrow(betas) == 0) {
    return(list(mat = betas, dropped_all_na = 0L, dropped_na = 0L, imputed = 0L, impute_method = "none"))
  }
  na_frac <- rowMeans(is.na(betas))
  drop_all <- is.na(na_frac) | na_frac >= 1
  dropped_all_na <- sum(drop_all)
  betas <- betas[!drop_all, , drop = FALSE]
  if (nrow(betas) == 0) {
    return(list(mat = betas, dropped_all_na = dropped_all_na, dropped_na = 0L, imputed = 0L, impute_method = "none"))
  }
  na_frac <- na_frac[!drop_all]
  keep <- na_frac <= na_max_frac
  dropped_na <- sum(!keep)
  betas <- betas[keep, , drop = FALSE]
  if (nrow(betas) == 0) {
    return(list(mat = betas, dropped_all_na = dropped_all_na, dropped_na = dropped_na, imputed = 0L, impute_method = "none"))
  }
  imputed <- 0L
  impute_method <- "none"
  if (anyNA(betas)) {
    knn_res <- impute_knn(betas, na_max_frac = na_max_frac)
    if (knn_res$method == "knn") {
      betas <- knn_res$mat
      imputed <- knn_res$imputed
      impute_method <- "knn"
    } else {
      impute_res <- impute_row_means(betas)
      betas <- impute_res$mat
      imputed <- impute_res$imputed
      impute_method <- "row_mean_fallback"
    }
  }
  list(mat = betas, dropped_all_na = dropped_all_na, dropped_na = dropped_na,
       imputed = imputed, impute_method = impute_method)
}

extract_sesame_cc_vector <- function(res) {
  if (is.null(res)) return(NULL)
  if (is.list(res) && !is.null(res$prop)) res <- res$prop
  if (is.data.frame(res) || is.matrix(res)) {
    if (nrow(res) == 1 && ncol(res) >= 1) {
      vec <- as.numeric(res[1, ])
      names(vec) <- colnames(res)
      return(vec)
    }
    if (ncol(res) == 1 && nrow(res) >= 1) {
      vec <- as.numeric(res[, 1])
      names(vec) <- rownames(res)
      return(vec)
    }
    return(NULL)
  }
  if (!is.numeric(res)) return(NULL)
  vec <- as.numeric(res)
  if (is.null(names(res)) || any(!nzchar(names(res)))) {
    names(vec) <- paste0("CellType", seq_along(vec))
  } else {
    names(vec) <- names(res)
  }
  vec
}

estimate_sesame_cell_counts <- function(sdf_list, sample_ids, out_dir = NULL) {
  if (length(sdf_list) == 0) return(NULL)
  if (!requireNamespace("sesame", quietly = TRUE)) return(NULL)
  fun_candidates <- c("estimateCellComposition", "estimateCellComposition2")
  for (fn in fun_candidates) {
    if (!exists(fn, envir = asNamespace("sesame"), inherits = FALSE)) next
    fun <- get(fn, envir = asNamespace("sesame"))
    probe <- tryCatch(fun(sdf_list[[1]]), error = function(e) e)
    if (inherits(probe, "error")) next
    vec <- extract_sesame_cc_vector(probe)
    if (is.null(vec)) next
    cell_types <- names(vec)
    mat <- vapply(sdf_list, function(sdf) {
      res <- tryCatch(fun(sdf), error = function(e) NULL)
      vec2 <- extract_sesame_cc_vector(res)
      if (is.null(vec2)) return(rep(NA_real_, length(cell_types)))
      vec2 <- vec2[cell_types]
      as.numeric(vec2)
    }, FUN.VALUE = numeric(length(cell_types)))
    df <- as.data.frame(t(mat))
    colnames(df) <- paste0("Cell_", cell_types)
    df$SampleID <- sample_ids
    if (!is.null(out_dir)) {
      write.csv(df, file.path(out_dir, "cell_counts_sesame.csv"), row.names = FALSE)
    }
    return(list(counts = df, reference = paste0("sesame::", fn)))
  }
  NULL
}

load_epidish_reference <- function(ref_name) {
  if (!requireNamespace("EpiDISH", quietly = TRUE)) return(NULL)
  if (exists(ref_name, envir = asNamespace("EpiDISH"), inherits = FALSE)) {
    ref <- get(ref_name, envir = asNamespace("EpiDISH"))
    if (!is.null(ref)) return(ref)
  }
  ref_env <- new.env(parent = emptyenv())
  tryCatch(data(list = ref_name, package = "EpiDISH", envir = ref_env), error = function(e) NULL)
  if (exists(ref_name, envir = ref_env, inherits = FALSE)) {
    return(ref_env[[ref_name]])
  }
  NULL
}

extract_epidish_fractions <- function(res, sample_ids) {
  mat <- NULL
  if (is.list(res)) {
    if (!is.null(res$estF)) {
      mat <- res$estF
    } else if (!is.null(res$cellProp)) {
      mat <- res$cellProp
    }
  }
  if (is.null(mat) && (is.matrix(res) || is.data.frame(res))) {
    mat <- res
  }
  if (is.null(mat)) return(NULL)
  mat <- as.data.frame(mat)
  n_samples <- length(sample_ids)
  if (nrow(mat) == n_samples) {
    if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
      rownames(mat) <- sample_ids
    }
  } else if (ncol(mat) == n_samples) {
    mat <- as.data.frame(t(mat))
    rownames(mat) <- sample_ids
  } else {
    return(NULL)
  }
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    colnames(mat) <- paste0("CellType", seq_len(ncol(mat)))
  }
  mat
}

estimate_cell_counts_epidish <- function(betas, tissue = "Blood", out_dir = NULL) {
  if (!requireNamespace("EpiDISH", quietly = TRUE)) return(NULL)
  if (is.null(betas) || nrow(betas) == 0) return(NULL)
  if (!tolower(tissue) %in% c("blood")) return(NULL)

  norm_res <- normalize_epicv2_ids(betas)
  betas_use <- norm_res$betas
  if (isTRUE(norm_res$collapsed)) {
    message("  - EPICv2-style CpG IDs detected; collapsed for EpiDISH reference alignment.")
  } else if (isTRUE(norm_res$normalized)) {
    message("  - EPICv2-style CpG IDs normalized for EpiDISH reference alignment.")
  }
  betas_use <- pmin(pmax(betas_use, 0), 1)

  ref_candidates <- c("centDHSbloodDMC.m", "IDOL", "centDHS")
  ref <- NULL
  ref_name <- NULL
  for (cand in ref_candidates) {
    ref_try <- load_epidish_reference(cand)
    if (!is.null(ref_try)) {
      ref <- ref_try
      ref_name <- cand
      break
    }
  }
  if (is.null(ref)) return(NULL)

  common <- intersect(rownames(betas_use), rownames(ref))
  if (length(common) < min(50, nrow(ref))) return(NULL)
  betas_use <- betas_use[common, , drop = FALSE]
  ref_use <- ref[common, , drop = FALSE]

  epifun <- NULL
  epifun_name <- NULL
  if (exists("epidish", envir = asNamespace("EpiDISH"), inherits = FALSE)) {
    epifun <- get("epidish", envir = asNamespace("EpiDISH"))
    epifun_name <- "epidish"
  } else if (exists("RPC", envir = asNamespace("EpiDISH"), inherits = FALSE)) {
    epifun <- get("RPC", envir = asNamespace("EpiDISH"))
    epifun_name <- "RPC"
  }
  if (is.null(epifun)) return(NULL)

  fn_args <- names(formals(epifun))
  arg_list <- list()
  if ("beta.m" %in% fn_args) {
    arg_list$beta.m <- betas_use
  } else if ("dat" %in% fn_args) {
    arg_list$dat <- betas_use
  } else {
    arg_list[[fn_args[1]]] <- betas_use
  }
  if ("ref.m" %in% fn_args) {
    arg_list$ref.m <- ref_use
  } else if ("ref" %in% fn_args) {
    arg_list$ref <- ref_use
  }
  if ("method" %in% fn_args) arg_list$method <- "RPC"
  if ("robust" %in% fn_args) arg_list$robust <- TRUE

  res <- tryCatch(
    with_muffled_conditions(do.call(epifun, arg_list), label = NULL),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  frac_mat <- extract_epidish_fractions(res, colnames(betas_use))
  if (is.null(frac_mat)) return(NULL)
  df <- as.data.frame(frac_mat)
  colnames(df) <- paste0("Cell_", colnames(df))
  df$SampleID <- rownames(df)
  if (!is.null(out_dir)) {
    write.csv(df, file.path(out_dir, "cell_counts_epidish.csv"), row.names = FALSE)
  }
  list(counts = df, reference = sprintf("EpiDISH::%s(%s)", epifun_name, ref_name))
}

estimate_cell_counts_planet <- function(betas, out_dir = NULL) {
  if (!requireNamespace("planet", quietly = TRUE)) return(NULL)
  if (!requireNamespace("minfi", quietly = TRUE)) return(NULL)
  if (is.null(betas) || nrow(betas) == 0 || ncol(betas) == 0) return(NULL)

  ref_env <- new.env(parent = emptyenv())
  tryCatch(data(list = "plCellCpGsThird", package = "planet", envir = ref_env), error = function(e) NULL)
  if (!exists("plCellCpGsThird", envir = ref_env, inherits = FALSE)) return(NULL)
  ref <- ref_env$plCellCpGsThird
  if (is.null(ref) || !(is.matrix(ref) || is.data.frame(ref))) return(NULL)
  ref <- as.matrix(ref)

  norm_res <- normalize_epicv2_ids(betas)
  betas_use <- norm_res$betas
  if (isTRUE(norm_res$collapsed)) {
    message("  - EPICv2-style CpG IDs detected; collapsed for planet reference alignment.")
  } else if (isTRUE(norm_res$normalized)) {
    message("  - EPICv2-style CpG IDs normalized for planet reference alignment.")
  }

  common <- intersect(rownames(betas_use), rownames(ref))
  if (length(common) < min(50, nrow(ref))) return(NULL)
  betas_use <- betas_use[common, , drop = FALSE]
  ref_use <- ref[common, , drop = FALSE]

  proj_fun <- NULL
  if (exists("projectCellType", envir = asNamespace("minfi"), inherits = FALSE)) {
    proj_fun <- get("projectCellType", envir = asNamespace("minfi"))
  }
  if (is.null(proj_fun)) return(NULL)

  res <- tryCatch(proj_fun(betas_use, ref_use, lessThanOne = FALSE), error = function(e) NULL)
  if (is.null(res)) return(NULL)

  mat <- as.data.frame(res)
  if (nrow(mat) != ncol(betas_use) && ncol(mat) == ncol(betas_use)) {
    mat <- as.data.frame(t(mat))
  }
  if (nrow(mat) != ncol(betas_use)) return(NULL)
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    rownames(mat) <- colnames(betas_use)
  }
  colnames(mat) <- paste0("Cell_", colnames(mat))
  mat$SampleID <- rownames(mat)
  if (!is.null(out_dir)) {
    write.csv(mat, file.path(out_dir, "cell_counts_planet.csv"), row.names = FALSE)
  }
  list(counts = mat, reference = "planet::plCellCpGsThird (projectCellType)")
}

estimate_cell_counts_reffree <- function(betas, out_dir = NULL, label = "RefFreeEWAS (Sesame)",
                                         K_latent = 5, max_probes = 10000) {
  if (!requireNamespace("RefFreeEWAS", quietly = TRUE)) return(NULL)
  if (is.null(betas) || nrow(betas) < 10 || ncol(betas) < 3) return(NULL)

  betas_use <- as.matrix(betas)
  storage.mode(betas_use) <- "numeric"
  betas_use <- betas_use[rowSums(is.na(betas_use)) == 0, , drop = FALSE]
  if (nrow(betas_use) < 10) return(NULL)

  vars <- apply(betas_use, 1, var)
  vars[!is.finite(vars)] <- 0
  top_n <- min(max_probes, length(vars))
  top_idx <- order(vars, decreasing = TRUE)[seq_len(top_n)]
  betas_sub <- betas_use[top_idx, , drop = FALSE]

  res <- tryCatch({
    suppressMessages(with_muffled_warnings(
      RefFreeEWAS::RefFreeCellMix(betas_sub, K = K_latent, verbose = FALSE),
      patterns = "ward",
      label = NULL
    ))
  }, error = function(e) NULL)
  if (is.null(res) || is.null(res$Omega)) return(NULL)

  df <- as.data.frame(res$Omega)
  colnames(df) <- paste0("Cell_Latent", seq_len(ncol(df)))
  if (is.null(rownames(df)) || any(!nzchar(rownames(df)))) {
    rownames(df) <- colnames(betas_sub)
  }
  df$SampleID <- rownames(df)
  if (!is.null(out_dir)) {
    write.csv(df, file.path(out_dir, "cell_counts_RefFree_Sesame.csv"), row.names = FALSE)
  }
  list(counts = df, reference = label)
}

append_cell_deconv_summary <- function(out_dir, summary_row) {
  if (is.null(out_dir) || !nzchar(out_dir) || is.null(summary_row)) return(invisible(NULL))
  path <- file.path(out_dir, "Cell_Deconvolution_Summary.csv")
  existing <- NULL
  if (file.exists(path)) {
    existing <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  if (is.null(existing)) existing <- data.frame()
  new_df <- rbind(existing, summary_row)
  write.csv(new_df, path, row.names = FALSE)
  invisible(new_df)
}

compute_dmr_delta_beta <- function(dmr_res, dmr_anno, mean_con, mean_test) {
  if (nrow(dmr_res) == 0 || nrow(dmr_anno) == 0) return(dmr_res)
  dmr_res$Mean_Beta_Con <- NA_real_
  dmr_res$Mean_Beta_Test <- NA_real_
  dmr_res$Delta_Beta <- NA_real_

  anno_chr <- as.character(dmr_anno$chr)
  res_chr <- as.character(dmr_res$chr)
  common_chr <- intersect(unique(res_chr), unique(anno_chr))
  if (length(common_chr) == 0) return(dmr_res)

  for (chr in common_chr) {
    idx_region <- which(res_chr == chr)
    idx_cpg <- which(anno_chr == chr)
    if (length(idx_cpg) == 0) next

    pos <- dmr_anno$pos[idx_cpg]
    keep_cpg <- is.finite(pos) & is.finite(mean_con[idx_cpg]) & is.finite(mean_test[idx_cpg])
    if (!any(keep_cpg)) next
    idx_cpg <- idx_cpg[keep_cpg]
    pos <- pos[keep_cpg]
    ord <- order(pos)
    pos_ord <- pos[ord]
    con_ord <- mean_con[idx_cpg][ord]
    test_ord <- mean_test[idx_cpg][ord]
    con_cum <- cumsum(con_ord)
    test_cum <- cumsum(test_ord)

    starts <- dmr_res$start[idx_region]
    ends <- dmr_res$end[idx_region]
    lo <- findInterval(starts - 1, pos_ord) + 1
    hi <- findInterval(ends, pos_ord)
    valid <- hi >= lo & hi > 0
    if (!any(valid)) next

    base_con <- numeric(length(lo))
    base_test <- numeric(length(lo))
    idx_base <- which(valid & lo > 1)
    if (length(idx_base) > 0) {
      base_con[idx_base] <- con_cum[lo[idx_base] - 1]
      base_test[idx_base] <- test_cum[lo[idx_base] - 1]
    }

    con_sum <- rep(NA_real_, length(lo))
    test_sum <- rep(NA_real_, length(lo))
    con_sum[valid] <- con_cum[hi[valid]] - base_con[valid]
    test_sum[valid] <- test_cum[hi[valid]] - base_test[valid]
    n_cpg <- hi - lo + 1

    con_mean <- con_sum / n_cpg
    test_mean <- test_sum / n_cpg
    valid_idx <- idx_region[valid]
    dmr_res$Mean_Beta_Con[valid_idx] <- con_mean[valid]
    dmr_res$Mean_Beta_Test[valid_idx] <- test_mean[valid]
    dmr_res$Delta_Beta[valid_idx] <- test_mean[valid] - con_mean[valid]
  }

  dmr_res
}

normalize_tissue <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  val <- trimws(as.character(x)[1])
  if (!nzchar(val)) return(NA_character_)
  key <- gsub("[^a-z0-9]", "", tolower(val))
  map <- c(
    "blood" = "Blood",
    "wholeblood" = "Blood",
    "peripheralblood" = "Blood",
    "cordblood" = "CordBlood",
    "umbilicalcordblood" = "CordBlood",
    "dlpfc" = "DLPFC",
    "dorsolateralprefrontalcortex" = "DLPFC",
    "placenta" = "Placenta",
    "saliva" = "Saliva"
  )
  if (key %in% names(map)) return(map[[key]])
  NA_character_
}

infer_tissue_from_config <- function(targets) {
  cand_cols <- c("tissue", "Tissue", "sample_tissue", "Sample_Tissue", "sample_type", "Sample_Type")
  for (col in cand_cols) {
    if (!col %in% colnames(targets)) next
    vals <- unique(trimws(as.character(targets[[col]])))
    vals <- vals[vals != ""]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    norm_vals <- unique(na.omit(vapply(vals, normalize_tissue, character(1))))
    if (length(norm_vals) == 1) {
      return(list(tissue = norm_vals[1], source = col, raw = vals[1]))
    }
  }
  hint <- infer_tissue_from_text(targets)
  if (!is.null(hint)) return(hint)
  NULL
}

infer_tissue_from_text <- function(targets) {
  if (is.null(targets) || nrow(targets) == 0) return(NULL)
  char_cols <- colnames(targets)[vapply(targets, function(x) is.character(x) || is.factor(x), logical(1))]
  if (length(char_cols) == 0) return(NULL)
  name_hits <- char_cols[grepl("tissue|source|character|cell|organ|specimen|biosample", tolower(char_cols))]
  if (length(name_hits) == 0) name_hits <- char_cols
  text_mat <- apply(targets[, name_hits, drop = FALSE], 1, function(row) {
    vals <- na.omit(as.character(row))
    vals <- vals[vals != ""]
    paste(vals, collapse = ";")
  })
  text_mat <- tolower(text_mat)
  patterns <- list(
    Blood = c("blood", "pbmc", "leuk", "lymph", "bcell", "tcell", "monocyte"),
    CordBlood = c("cord", "umbilical"),
    Placenta = c("placenta", "chorion", "trophoblast"),
    Saliva = c("saliva"),
    DLPFC = c("dlpfc", "prefrontal", "cortex")
  )
  counts <- vapply(patterns, function(pats) {
    sum(grepl(paste(pats, collapse = "|"), text_mat))
  }, integer(1))
  if (all(counts == 0)) return(NULL)
  best <- names(counts)[which.max(counts)][1]
  if (!is.null(best) && counts[[best]] >= ceiling(0.6 * length(text_mat))) {
    return(list(tissue = best, source = "text_heuristic", raw = best))
  }
  NULL
}

detect_array_type <- function(rgSet) {
  if (!requireNamespace("minfi", quietly = TRUE)) return(NA_character_)
  try_check <- function(fn) {
    if (exists(fn, envir = asNamespace("minfi"), inherits = FALSE)) {
      return(tryCatch(get(fn, asNamespace("minfi"))(rgSet), error = function(e) FALSE))
    }
    FALSE
  }
  if (try_check(".isEPICv2")) return("EPICv2")
  if (try_check(".isEPIC")) return("EPIC")
  if (try_check(".is450k")) return("450k")
  ann <- tryCatch(minfi::annotation(rgSet), error = function(e) NULL)
  if (!is.null(ann) && "array" %in% names(ann)) {
    return(as.character(ann["array"]))
  }
  NA_character_
}

reference_platform_from_array <- function(array_type) {
  if (is.na(array_type) || !nzchar(array_type)) return("")
  key <- gsub("[^a-z0-9]", "", tolower(array_type))
  if (grepl("epicv2", key) || grepl("950", key)) return("IlluminaHumanMethylationEPICv2")
  if (grepl("epic", key) || grepl("850", key)) return("IlluminaHumanMethylationEPIC")
  if (grepl("450k", key) || grepl("hm450", key) || grepl("450", key)) return("IlluminaHumanMethylation450k")
  ""
}

load_cell_reference <- function(ref_spec) {
  if (!nzchar(ref_spec)) return(NULL)
  if (file.exists(ref_spec)) {
    if (grepl("\\.rds$", ref_spec, ignore.case = TRUE)) {
      ref <- tryCatch(readRDS(ref_spec), error = function(e) NULL)
      if (!is.null(ref)) {
        return(list(reference = ref, label = basename(ref_spec), pkg = NA_character_))
      }
    }
    if (grepl("\\.(rda|rdata)$", ref_spec, ignore.case = TRUE)) {
      ref_env <- new.env(parent = emptyenv())
      tryCatch(load(ref_spec, envir = ref_env), error = function(e) NULL)
      objs <- ls(ref_env)
      if (length(objs) > 0) {
        if (length(objs) > 1) {
          message("  - Multiple objects found in reference file; using the first: ", objs[1])
        }
        return(list(reference = ref_env[[objs[1]]], label = basename(ref_spec), pkg = NA_character_))
      }
    }
    return(NULL)
  }
  if (grepl("::", ref_spec, fixed = TRUE)) {
    parts <- strsplit(ref_spec, "::", fixed = TRUE)[[1]]
    pkg <- parts[1]
    obj <- parts[2]
    if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
    ref <- tryCatch(get(obj, envir = asNamespace(pkg)), error = function(e) NULL)
    if (!is.null(ref)) {
      return(list(reference = ref, label = ref_spec, pkg = pkg))
    }
  }
  if (requireNamespace(ref_spec, quietly = TRUE)) {
    ref_env <- new.env(parent = emptyenv())
    tryCatch(data(list = ref_spec, package = ref_spec, envir = ref_env), error = function(e) NULL)
    if (exists(ref_spec, envir = ref_env, inherits = FALSE)) {
      return(list(reference = ref_env[[ref_spec]], label = ref_spec, pkg = ref_spec))
    }
  }
  NULL
}

estimate_cell_counts_safe <- function(rgSet, tissue = "Blood", out_dir = NULL,
                                      cell_reference = "", cell_reference_platform = "",
                                      array_type = NA_character_) {
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
      # Top variable probes is usually sufficient for deconvolution
      vars <- apply(betas, 1, var)
      top_idx <- head(order(vars, decreasing = TRUE), min(cell_ref_free_max_probes, nrow(betas)))
      betas_sub <- betas[top_idx, ]
      
      # Run RefFreeCellMix
      # Fixed K is safer for automation; configurable via cell_deconv.ref_free_k.
      K_latent <- cell_ref_free_k
      message(sprintf("  - Running RefFreeCellMix with K=%d latent components...", K_latent))
      
      # Initialize with K-means (standard RefFreeEWAS approach)
      res <- suppressMessages(with_muffled_warnings(
        RefFreeEWAS::RefFreeCellMix(betas_sub, K = K_latent, verbose = FALSE),
        patterns = "ward",
        label = NULL
      ))
      
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

  if (nzchar(cell_reference)) {
    ref_info <- load_cell_reference(cell_reference)
    if (!is.null(ref_info)) {
      ref_tissue <- tissue
      if (ref_tissue == "Auto") {
        message("Custom reference provided with tissue=Auto; using compositeCellType='Blood'.")
        ref_tissue <- "Blood"
      }
      message(sprintf("Estimating cell composition using custom reference: %s...", ref_info$label))
      est_fun <- NULL
      if (!is.na(ref_info$pkg) && requireNamespace(ref_info$pkg, quietly = TRUE)) {
        if (exists("estimateCellCounts2", envir = asNamespace(ref_info$pkg), inherits = FALSE)) {
          est_fun <- get("estimateCellCounts2", envir = asNamespace(ref_info$pkg))
        }
      }
      if (is.null(est_fun)) {
        if (requireNamespace("minfi", quietly = TRUE)) {
          if (exists("estimateCellCounts2", envir = asNamespace("minfi"), inherits = FALSE)) {
            est_fun <- minfi::estimateCellCounts2
          } else if (exists("estimateCellCounts", envir = asNamespace("minfi"), inherits = FALSE)) {
            est_fun <- minfi::estimateCellCounts
          }
        }
      }
      if (!is.null(est_fun)) {
        ref_platform <- cell_reference_platform
        if (!nzchar(ref_platform)) {
          ref_platform <- reference_platform_from_array(array_type)
        }
        res <- tryCatch({
          arg_list <- list(
            rgSet = rgSet,
            compositeCellType = ref_tissue,
            returnAll = FALSE
          )
          fn_args <- names(formals(est_fun))
          if ("reference" %in% fn_args) arg_list$reference <- ref_info$reference
          if ("referencePlatform" %in% fn_args && nzchar(ref_platform)) arg_list$referencePlatform <- ref_platform
          if ("processMethod" %in% fn_args) arg_list$processMethod <- "preprocessNoob"
          if ("normalizationMethod" %in% fn_args) arg_list$normalizationMethod <- "none"
          if ("meanPlot" %in% fn_args) arg_list$meanPlot <- FALSE
          suppressWarnings(do.call(est_fun, arg_list))
        }, error = function(e) {
          message("  - Custom reference cell composition failed: ", e$message)
          NULL
        })
        if (!is.null(res)) {
          if (is.list(res) && !is.null(res$prop)) {
            res <- res$prop
          }
          df <- as.data.frame(res)
          colnames(df) <- paste0("Cell_", colnames(df))
          df$SampleID <- rownames(df)
          if (!is.null(out_dir)) {
            write.csv(df, file.path(out_dir, "cell_counts_custom_reference.csv"), row.names = FALSE)
          }
          message("  - Cell composition estimated using custom reference.")
          return(list(counts = df, reference = ref_info$label))
        }
      } else {
        message("  - Custom reference skipped: estimateCellCounts[2] not available.")
      }
    } else {
      message("  - Unable to load custom cell reference; falling back to defaults.")
    }
  }

  if (tolower(tissue) == "placenta") {
    message("Estimating cell composition for tissue 'Placenta' using planet reference...")
    betas <- tryCatch({
      mSet <- preprocessNoob(rgSet)
      getBeta(mSet)
    }, error = function(e) NULL)
    if (!is.null(betas)) {
      res <- estimate_cell_counts_planet(betas, out_dir = out_dir)
      if (!is.null(res)) {
        message("  - Cell composition estimated using planet (Placenta reference).")
        return(res)
      }
    }
    message("  - Placenta reference-based estimation failed; falling back to reference-free.")
    return(run_ref_free())
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
      
      suppressWarnings(do.call(est_fun, arg_list))
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

#' Principal Variance Component Analysis (PVCA) for Batch Effect Assessment
#'
#' Quantifies variance explained by technical and biological factors using
#' a mixed-effects model approach. Identifies batch effects by examining the
#' proportion of variance attributable to each experimental factor.
#'
#' @param betas Beta value matrix (probes x samples)
#' @param meta Sample metadata data frame
#' @param factors Character vector of factor column names to assess
#' @param prefix Pipeline prefix for output file naming
#' @param out_dir Output directory for PVCA results
#' @param threshold Variance threshold for factor significance (default: 0.6)
#' @param max_probes Maximum probes to use for estimation (default: 5000)
#' @param sample_ids Optional vector of sample identifiers
#'
#' @return Data frame of variance components by factor, or NULL if skipped
#'
#' @details
#' PVCA is skipped when: insufficient samples (< PVCA_MIN_SAMPLES),
#' no usable factors, or crossed factor levels exceed sample size.
#' Factors are automatically dropped if their cardinality would cause
#' over-specification of the mixed model.
#'
#' @references
#' Bushel et al. (2007) Simultaneous clustering of gene expression data
#' with clinical chemistry and pathological evaluations reveals phenotypic
#' prototypes. BMC Systems Biology 1:15
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
  
  pvca_res <- tryCatch(
    with_muffled_conditions(
      pvcaBatchAssess(eset, colnames(pvca_meta), threshold),
      label = "PVCA",
      patterns = "boundary \\(singular\\) fit"
    ),
    error = function(e) e
  )
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
  pvca_df
}

compute_pvca_confounding_index <- function(pvca_df, batch_term, group_term = "primary_group") {
  if (is.null(pvca_df) || nrow(pvca_df) == 0) return(NA_real_)
  batch_term <- as.character(batch_term)
  group_term <- as.character(group_term)
  b <- suppressWarnings(as.numeric(pvca_df$proportion[pvca_df$term == batch_term]))
  g <- suppressWarnings(as.numeric(pvca_df$proportion[pvca_df$term == group_term]))
  if (length(b) == 0 || length(g) == 0) return(NA_real_)
  b <- b[1]
  g <- g[1]
  if (!is.finite(b) || !is.finite(g) || g <= 0) return(NA_real_)
  b / g
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

  normalize_clock_betas <- function(betas_in) {
    betas_in <- as.matrix(betas_in)
    storage.mode(betas_in) <- "numeric"
    ids <- rownames(betas_in)
    if (is.null(ids)) return(list(betas = betas_in, normalized = FALSE))
    has_suffix <- grepl("^cg[0-9]+_", ids)
    if (!any(has_suffix)) return(list(betas = betas_in, normalized = FALSE))

    base_ids <- ids
    base_ids[has_suffix] <- sub("_.+$", "", base_ids[has_suffix])
    dup_count <- sum(duplicated(base_ids))
    if (dup_count == 0) {
      rownames(betas_in) <- base_ids
      return(list(betas = betas_in, normalized = TRUE, collapsed = FALSE, dup_count = 0))
    }

    collapse_res <- tryCatch({
      na_mask <- is.na(betas_in)
      betas_zero <- betas_in
      betas_zero[na_mask] <- 0
      sums <- rowsum(betas_zero, group = base_ids, reorder = FALSE)
      counts <- rowsum((!na_mask) * 1, group = base_ids, reorder = FALSE)
      betas_out <- sums / counts
      betas_out[counts == 0] <- NA_real_
      list(betas = betas_out)
    }, error = function(e) e)
    if (inherits(collapse_res, "error")) {
      message("  EPICv2 CpG collapse failed; using first occurrence for duplicate IDs: ", collapse_res$message)
      keep_idx <- !duplicated(base_ids)
      betas_out <- betas_in[keep_idx, , drop = FALSE]
      rownames(betas_out) <- base_ids[keep_idx]
      return(list(betas = betas_out, normalized = TRUE, collapsed = FALSE, dup_count = dup_count))
    }
    return(list(betas = collapse_res$betas, normalized = TRUE, collapsed = TRUE, dup_count = dup_count))
  }
  
  # --- 1. methylclock (Default / General Clocks) ---
  # Prioritize methylclock as it is more robust and supports newer clocks
  mc_success <- FALSE
  mc_available <- suppressWarnings(suppressMessages(requireNamespace("methylclock", quietly = TRUE)))
  if (isTRUE(mc_available)) {
     message("  Computing epigenetic clocks using 'methylclock'...")
     tryCatch({
       clock_betas <- betas
       norm_res <- normalize_clock_betas(clock_betas)
       clock_betas <- norm_res$betas
       if (isTRUE(norm_res$normalized)) {
         if (isTRUE(norm_res$collapsed)) {
           message("  EPICv2-style CpG IDs detected; collapsed to base cg IDs for clock computation (duplicates: ", norm_res$dup_count, ").")
         } else {
           message("  EPICv2-style CpG IDs detected; normalized CpG IDs for clock computation.")
         }
       }
       mc_attempts <- list(
         list(label = "fastImp=TRUE, cell.count=TRUE", args = list(fastImp = TRUE, cell.count = TRUE)),
         list(label = "fastImp=FALSE, cell.count=TRUE", args = list(fastImp = FALSE, cell.count = TRUE)),
         list(label = "fastImp=FALSE, cell.count=FALSE", args = list(fastImp = FALSE, cell.count = FALSE))
       )
       mc_res <- NULL
       mc_err <- NULL
       for (att in mc_attempts) {
         mc_try <- tryCatch(
           suppressMessages(suppressWarnings(do.call(methylclock::DNAmAge, c(list(x = clock_betas, toBetas = FALSE, normalize = FALSE), att$args)))),
           error = function(e) e
         )
         if (!inherits(mc_try, "error")) {
           mc_res <- mc_try
           message("  methylclock succeeded (", att$label, ").")
           break
         } else {
           mc_err <- mc_try
           message("  methylclock attempt failed (", att$label, "): ", mc_try$message)
         }
       }
       if (is.null(mc_res)) {
         stop(mc_err$message)
       }
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

detect_batch_candidates <- function(meta, group_col, config_settings) {
  cand <- character(0)
  if (!is.null(config_settings$batch_candidates) && length(config_settings$batch_candidates) > 0) {
    cand <- intersect(config_settings$batch_candidates, colnames(meta))
  } else {
    patterns <- config_settings$batch_candidate_patterns
    if (length(patterns) > 0) {
      pattern <- paste(patterns, collapse = "|")
      cand <- colnames(meta)[grepl(pattern, tolower(colnames(meta)))]
    }
    cand <- unique(c(c("Sentrix_ID", "Sentrix_Position"), cand))
    cand <- intersect(cand, colnames(meta))
  }
  cand <- setdiff(cand, group_col)
  cand <- cand[vapply(cand, function(x) {
    vals <- meta[[x]]
    if (all(is.na(vals))) return(FALSE)
    uniq <- length(unique(vals[!is.na(vals)]))
    uniq > 1 && uniq < nrow(meta)
  }, logical(1))]
  cand
}

compute_cramers_v <- function(tbl) {
  if (is.null(tbl) || any(dim(tbl) < 2)) return(NA_real_)
  n <- sum(tbl)
  if (!is.finite(n) || n <= 0) return(NA_real_)
  chi <- tryCatch(suppressWarnings(chisq.test(tbl, correct = FALSE)), error = function(e) NULL)
  if (is.null(chi)) return(NA_real_)
  if (is.null(chi$statistic) || !is.finite(chi$statistic)) return(NA_real_)
  r <- nrow(tbl)
  c <- ncol(tbl)
  denom <- n * (min(r - 1, c - 1))
  if (denom <= 0) return(NA_real_)
  sqrt(as.numeric(chi$statistic) / denom)
}

assess_batch_confounding <- function(meta, batch_cols, group_col, config_settings) {
  out <- data.frame(
    batch = character(0),
    voi_type = character(0),
    n_samples = integer(0),
    n_batch_levels = integer(0),
    n_voi_levels = integer(0),
    min_cell_count = integer(0),
    pct_empty = numeric(0),
    min_balance = numeric(0),
    chi_p = numeric(0),
    cramer_v = numeric(0),
    r2 = numeric(0),
    tier = integer(0),
    identifiability_warning = character(0),
    stringsAsFactors = FALSE
  )
  if (length(batch_cols) == 0) return(out)
  min_overlap <- config_settings$min_overlap_per_cell
  t1_v <- config_settings$confounding$tier1_v
  t2_v <- config_settings$confounding$tier2_v
  t1_r2 <- config_settings$confounding$tier1_r2
  t2_r2 <- config_settings$confounding$tier2_r2
  for (batch_col in batch_cols) {
    if (!(batch_col %in% colnames(meta)) || !(group_col %in% colnames(meta))) next
    vals <- meta[[batch_col]]
    voi <- meta[[group_col]]
    cc <- complete.cases(vals, voi)
    if (sum(cc) < 2) next
    vals <- vals[cc]
    voi <- voi[cc]
    n <- length(voi)
    voi_is_num <- (is.numeric(voi) || is.integer(voi)) && length(unique(voi)) > 2
    if (voi_is_num) {
      batch_f <- as.factor(vals)
      counts <- table(batch_f)
      min_cell <- suppressWarnings(min(counts))
      r2 <- tryCatch(summary(lm(voi ~ batch_f))$r.squared, error = function(e) NA_real_)
      tier <- 0
      if (!is.na(r2) && r2 >= t2_r2) tier <- 2 else if (!is.na(r2) && r2 >= t1_r2) tier <- 1
      if (!is.na(min_cell) && min_cell < min_overlap) tier <- 3
      out <- rbind(out, data.frame(
        batch = batch_col,
        voi_type = "continuous",
        n_samples = n,
        n_batch_levels = length(unique(batch_f)),
        n_voi_levels = length(unique(voi)),
        min_cell_count = ifelse(is.na(min_cell), 0L, as.integer(min_cell)),
        pct_empty = 0,
        min_balance = NA_real_,
        chi_p = NA_real_,
        cramer_v = NA_real_,
        r2 = r2,
        tier = tier,
        identifiability_warning = ifelse(tier == 3, "non_identifiable", ""),
        stringsAsFactors = FALSE
      ))
    } else {
      batch_f <- as.factor(vals)
      voi_f <- as.factor(voi)
      tbl <- table(batch_f, voi_f)
      min_cell <- suppressWarnings(min(tbl))
      pct_empty <- sum(tbl == 0) / length(tbl)
      min_balance <- suppressWarnings(min(apply(prop.table(tbl, 1), 1, function(x) min(x))))
      chi_p <- safe_chisq_p(tbl)
      cramer_v <- compute_cramers_v(tbl)
      voi_in_batch <- apply(tbl, 2, function(x) sum(x > 0))
      batch_in_voi <- apply(tbl, 1, function(x) sum(x > 0))
      tier3 <- any(batch_in_voi == 1) || any(voi_in_batch == 1) ||
        (!is.na(min_cell) && min_cell < min_overlap)
      tier <- 0
      if (!is.na(cramer_v) && cramer_v >= t2_v) tier <- 2 else if (!is.na(cramer_v) && cramer_v >= t1_v) tier <- 1
      if (tier3) tier <- 3
      out <- rbind(out, data.frame(
        batch = batch_col,
        voi_type = "categorical",
        n_samples = n,
        n_batch_levels = length(unique(batch_f)),
        n_voi_levels = length(unique(voi_f)),
        min_cell_count = ifelse(is.na(min_cell), 0L, as.integer(min_cell)),
        pct_empty = pct_empty,
        min_balance = min_balance,
        chi_p = chi_p,
        cramer_v = cramer_v,
        r2 = NA_real_,
        tier = tier,
        identifiability_warning = ifelse(tier == 3, "non_identifiable", ""),
        stringsAsFactors = FALSE
      ))
    }
  }
  out
}

plot_confounding_heatmap <- function(meta, batch_col, group_col, prefix, out_dir) {
  if (!(batch_col %in% colnames(meta)) || !(group_col %in% colnames(meta))) return(NULL)
  tbl <- table(meta[[batch_col]], meta[[group_col]])
  df <- as.data.frame(tbl)
  colnames(df) <- c("Batch", "VOI", "Count")
  p <- ggplot(df, aes(x = VOI, y = Batch, fill = Count, text = Count)) +
    geom_tile() +
    geom_text(aes(label = Count), size = 3, color = "black") +
    scale_fill_gradient(low = "white", high = "#2c7fb8") +
    theme_minimal() +
    labs(title = paste(prefix, "Batch-VOI Overlap:", batch_col), x = group_col, y = batch_col)
  save_interactive_plot(p, paste0(prefix, "_BatchVOI_Heatmap_", batch_col, ".html"), out_dir)
  save_static_plot(p, paste0(prefix, "_BatchVOI_Heatmap_", batch_col, ".png"), out_dir, width = 6, height = 4)
  p
}

plot_confounding_summary <- function(conf_df, prefix, out_dir) {
  if (is.null(conf_df) || nrow(conf_df) == 0) return(NULL)
  df <- conf_df
  df$metric <- ifelse(df$voi_type == "continuous", df$r2, df$cramer_v)
  df$metric_label <- ifelse(df$voi_type == "continuous", "R2", "CramersV")
  p <- ggplot(df, aes(x = batch, y = metric, fill = as.factor(tier), text = metric)) +
    geom_col() +
    theme_minimal() +
    labs(title = paste(prefix, "Confounding Summary"), x = "Batch candidate", y = "Association strength", fill = "Tier")
  save_interactive_plot(p, paste0(prefix, "_Confounding_Summary.html"), out_dir)
  save_static_plot(p, paste0(prefix, "_Confounding_Summary.png"), out_dir, width = 6, height = 4)
  p
}

identify_overlap_batches <- function(meta, batch_col, group_col) {
  if (is.null(batch_col) || !(batch_col %in% colnames(meta))) return(character(0))
  tbl <- table(meta[[batch_col]], meta[[group_col]])
  if (nrow(tbl) == 0 || ncol(tbl) == 0) return(character(0))
  keep <- rownames(tbl)[apply(tbl, 1, function(x) sum(x > 0) >= 2)]
  keep
}

is_batch_confounded <- function(meta, batch_col, group_col, p_thresh = 1e-6) {
  if (is.null(batch_col) || is.null(group_col)) return(FALSE)
  if (!(batch_col %in% colnames(meta)) || !(group_col %in% colnames(meta))) return(FALSE)
  tbl <- table(meta[[batch_col]], meta[[group_col]])
  if (any(tbl == 0)) return(TRUE)
  p_chisq <- safe_chisq_p(tbl)
  if (!is.na(p_chisq) && p_chisq < p_thresh) return(TRUE)
  return(FALSE)
}

eval_batch_method <- function(M_mat, meta, group_col, batch_col, covariates, method) {
  combat_par_prior <- if (exists("combat_par_prior", inherits = TRUE)) combat_par_prior else TRUE
  combat_allowed <- if (exists("combat_allowed", inherits = TRUE)) combat_allowed else TRUE
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
  note <- ""
  cov_terms <- setdiff(covariates, batch_col)
  cov_terms_use <- cov_terms
  if (length(cov_terms_use) > 0 && !is.null(batch_col)) {
    batch_filter <- filter_covariates_with_batch(meta_use, cov_terms_use, batch_col)
    cov_terms_use <- batch_filter$keep
    if (nrow(batch_filter$dropped) > 0) {
      note <- paste0("batch_confounded_covariates=", paste(unique(batch_filter$dropped$Variable), collapse = ";"))
    }
  }
  for (ct in cov_terms_use) {
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
    tryCatch(model.matrix(as.formula(formula_str), data = data), error = function(e) NULL)
  }
  
  M_corr <- M_mat
  if (method == "none") {
    M_corr <- M_mat
  } else if (method == "combat") {
    if (!isTRUE(combat_allowed)) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = "ComBat disabled for small sample size.",
                  note = note))
    }
    combat_setup <- build_combat_design(meta_use, group_col, cov_terms_use, batch_col)
    if (nrow(combat_setup$dropped) > 0) {
      drop_note <- paste0("combat_drop=", paste(unique(combat_setup$dropped$Variable), collapse = ";"))
      if (nzchar(note)) {
        note <- paste(note, drop_note, sep = " | ")
      } else {
        note <- drop_note
      }
    }
    mod <- combat_setup$mod
    if (ncol(mod) < 2) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = "ComBat design has <2 columns after filtering.",
                  note = note))
    }
    combat_try <- tryCatch(
      ComBat(dat = M_mat, batch = batch, mod = mod, par.prior = combat_par_prior, prior.plots = FALSE),
      error = function(e) e
    )
    if (inherits(combat_try, "error")) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = combat_try$message, note = note))
    }
    M_corr <- combat_try
  } else if (method == "limma") {
    design <- mm_safe(paste("~ 0 +", group_col), meta_use)
    if (is.null(design) || ncol(design) < 1) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = "Limma design failed.", note = note))
    }
    cov_mat <- if (length(cov_terms_use) > 0) mm_safe(paste("~ 0 +", paste(cov_terms_use, collapse = " + ")), meta_use) else NULL
    if (length(cov_terms_use) > 0 && (is.null(cov_mat) || ncol(cov_mat) < 1)) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = "Limma covariate design failed.", note = note))
    }
    M_corr <- tryCatch(
      removeBatchEffect(M_mat, batch = batch, covariates = cov_mat, design = design),
      error = function(e) {
        message(sprintf("    - removeBatchEffect (limma) failed: %s", conditionMessage(e)))
        M_mat
      }
    )
  } else if (method == "sva") {
    mod <- mm_safe(paste("~", paste(c(group_col, cov_terms_use), collapse = " + ")), meta_use)
    mod0_terms <- if (length(cov_terms_use) > 0) paste(cov_terms_use, collapse = " + ") else "1"
    mod0 <- mm_safe(paste("~", mod0_terms), meta_use)
    if (is.null(mod) || is.null(mod0) || ncol(mod) < 2 || ncol(mod0) < 1) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = "SVA design has <2 columns after filtering.",
                  note = note))
    }
    sva_try <- tryCatch(sva(M_mat, mod, mod0), error = function(e) e)
    if (inherits(sva_try, "error")) {
      return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                  n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_mat,
                  failed = TRUE, fail_reason = sva_try$message, note = note))
    }
    sva.obj <- sva_try
    if (is.null(sva.obj$sv) || sva.obj$n.sv < 1) {
      M_corr <- M_mat
    } else {
      sv_mat <- sva.obj$sv
      if (nrow(sv_mat) != nrow(meta_use)) {
        message("    - SVA SV rows do not match samples; skipping SVA correction for evaluation.")
        M_corr <- M_mat
      } else {
        design_keep <- mm_safe(paste("~ 0 +", group_col), meta_use)
        M_corr <- tryCatch(
          removeBatchEffect(M_mat, covariates = sv_mat, design = design_keep),
          error = function(e) {
            message(sprintf("    - removeBatchEffect (SVA) failed: %s", conditionMessage(e)))
            M_mat
          }
        )
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
  for (ct in cov_terms_use) {
    if (is.numeric(meta_use[[ct]]) || is.integer(meta_use[[ct]])) {
      form_terms <- c(form_terms, ct)
    } else {
      form_terms <- c(form_terms, sprintf("(1|%s)", ct))
    }
  }
  form_terms <- c(form_terms, sprintf("(1|%s)", group_col))
  form <- as.formula(paste("~", paste(form_terms, collapse = " + ")))
  M_top_varpart <- M_corr[top_idx, , drop = FALSE]
  vp_cfg <- if (exists("config_settings")) config_settings$variance_partition else NULL
  vp_label <- paste("variancePartition", method)
  vp_fit <- fit_variance_partition_safe(M_top_varpart, form, meta_use, vp_cfg = vp_cfg, label = vp_label)
  varPart <- vp_fit$result
  if (is.null(varPart)) {
    med_varFrac <- setNames(rep(0, length(c(batch_col, group_col))), c(batch_col, group_col))
  } else {
    med_varFrac <- apply(varPart, 2, median, na.rm = TRUE)
  }
  
  design_terms <- c(batch_col, group_col, cov_terms_use)
  design <- mm_safe(paste("~ 0 +", paste(design_terms, collapse = " + ")), meta_use)
  if (is.null(design) || ncol(design) < 2) {
    return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_corr,
                failed = TRUE, fail_reason = "Batch design failed.", note = note))
  }
  fit_batch <- tryCatch(lmFit(M_corr, design), error = function(e) NULL)
  fit_batch <- safe_ebayes(fit_batch, "batch_eval")
  if (is.null(fit_batch)) {
    return(list(method = method, score = -Inf, batch_var = NA_real_, group_var = NA_real_,
                n_pc_batch_sig = NA_real_, prop_batch_sig = NA_real_, M_corr = M_corr,
                failed = TRUE, fail_reason = "Batch eBayes failed.", note = note))
  }
  batch_cols <- grep(paste0("^", batch_col), colnames(design))
  if (length(batch_cols) == 0) {
    prop_batch_sig <- NA_real_
  } else {
    contrast_batch <- diag(ncol(design))[, batch_cols, drop = FALSE]
    fit_con <- tryCatch(contrasts.fit(fit_batch, contrast_batch), error = function(e) NULL)
    fit_con <- safe_ebayes(fit_con, "batch_eval_contrast")
    if (is.null(fit_con)) {
      prop_batch_sig <- NA_real_
    } else {
      p_batch <- fit_con$F.p.value
      fdr_batch <- p.adjust(p_batch, "BH")
      prop_batch_sig <- mean(fdr_batch < 0.05, na.rm = TRUE)
    }
  }
  
  batch_var <- if (batch_col %in% names(med_varFrac)) med_varFrac[batch_col] else NA_real_
  group_var <- if (group_col %in% names(med_varFrac)) med_varFrac[group_col] else NA_real_
  if (is.na(batch_var)) batch_var <- 0
  if (is.na(group_var)) group_var <- 0
  if (is.na(prop_batch_sig)) prop_batch_sig <- 0
  score <- batch_var + prop_batch_sig + 0.1 * n_pc_batch_sig
  
  list(
    method = method,
    score = score,
    batch_var = batch_var,
    group_var = group_var,
    n_pc_batch_sig = n_pc_batch_sig,
    prop_batch_sig = prop_batch_sig,
    M_corr = M_corr,
    failed = FALSE,
    fail_reason = "",
    note = note
  )
}

apply_batch_correction <- function(betas, method, batch_col, covariates, targets, design,
                                   group_col = "primary_group") {
  combat_par_prior <- if (exists("combat_par_prior", inherits = TRUE)) combat_par_prior else TRUE
  combat_allowed <- if (exists("combat_allowed", inherits = TRUE)) combat_allowed else TRUE
  # Apply the selected batch correction method to betas (logit space) while preserving group effects.
  if (is.null(method) || method == "none" || is.null(batch_col) || !(batch_col %in% colnames(targets))) {
    return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
  }
  batch_vals <- targets[[batch_col]]
  if (length(unique(batch_vals)) < 2) {
    return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
  }
  covariates <- covariates[covariates %in% colnames(targets)]
  batch_filter <- filter_covariates_with_batch(targets, covariates, batch_col)
  covariates_use <- batch_filter$keep
  if (nrow(batch_filter$dropped) > 0) {
    message(paste("  Batch-confounded covariates removed:", paste(unique(batch_filter$dropped$Variable), collapse = ", ")))
  }
  cov_mat <- NULL
  if (length(covariates_use) > 0) {
    cov_formula <- as.formula(paste("~ 0 +", paste(covariates_use, collapse = " + ")))
    cov_mat <- model.matrix(cov_formula, data = targets)
  }
  design_keep <- model.matrix(as.formula(paste("~", group_col)), data = targets)
  M_mat <- logit_offset(betas)
  method_used <- method

  if (method == "combat") {
    if (!isTRUE(combat_allowed)) {
      message("  ComBat disabled for small sample size; skipping.")
      return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
    }
    combat_setup <- build_combat_design(targets, group_col, covariates_use, batch_col)
    if (nrow(combat_setup$dropped) > 0) {
      message(paste("  ComBat: dropped covariates/terms:", paste(unique(combat_setup$dropped$Variable), collapse = ", ")))
    }
    mod <- combat_setup$mod
    if (ncol(mod) < 2) {
      message("  ComBat skipped: design has <2 columns after filtering.")
      return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
    }
    combat_try <- tryCatch(
      ComBat(dat = M_mat, batch = as.factor(batch_vals), mod = mod, par.prior = combat_par_prior, prior.plots = FALSE),
      error = function(e) e
    )
    if (inherits(combat_try, "error")) {
      message(paste("  ComBat failed:", combat_try$message, "-> falling back to limma."))
      limma_try <- tryCatch(
        removeBatchEffect(M_mat, batch = as.factor(batch_vals), covariates = cov_mat, design = design_keep),
        error = function(e) e
      )
      if (inherits(limma_try, "error")) {
        message(paste("  limma fallback failed:", limma_try$message, "-> skipping batch correction."))
        return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
      }
      M_corr <- limma_try
      method_used <- "limma"
    } else {
      M_corr <- combat_try
    }
  } else if (method == "limma") {
    limma_try <- tryCatch(
      removeBatchEffect(M_mat, batch = as.factor(batch_vals), covariates = cov_mat, design = design_keep),
      error = function(e) e
    )
    if (inherits(limma_try, "error")) {
      message(paste("  limma batch correction failed:", limma_try$message, "-> skipping batch correction."))
      return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
    }
    M_corr <- limma_try
  } else if (method == "sva") {
    # SVA handled via surrogate variables already included in the design.
    return(list(betas = betas, mvals = logit_offset(betas), method = "sva"))
  } else {
    return(list(betas = betas, mvals = logit_offset(betas), method = "none"))
  }

  betas_corr <- 2^M_corr / (1 + 2^M_corr)
  betas_corr <- pmin(pmax(betas_corr, LOGIT_OFFSET), 1 - LOGIT_OFFSET)
  list(betas = betas_corr, mvals = M_corr, method = method_used)
}

meta_analysis_fixed <- function(effects, ses) {
  w <- 1 / (ses ^ 2)
  w[!is.finite(w)] <- 0
  effects[!is.finite(effects)] <- NA
  w[is.na(effects)] <- 0
  w_sum <- rowSums(w)
  k_eff <- rowSums(w > 0)
  meta_beta <- rowSums(w * effects, na.rm = TRUE) / w_sum
  meta_se <- sqrt(1 / w_sum)
  z <- meta_beta / meta_se
  p <- 2 * pnorm(-abs(z))
  q <- rowSums(w * (effects - meta_beta) ^ 2, na.rm = TRUE)
  i2 <- ifelse(q > 0, pmax(0, (q - (k_eff - 1)) / q), NA_real_)
  meta_beta[k_eff < 2] <- NA_real_
  meta_se[k_eff < 2] <- NA_real_
  p[k_eff < 2] <- NA_real_
  q[k_eff < 2] <- NA_real_
  i2[k_eff < 2] <- NA_real_
  list(beta = meta_beta, se = meta_se, p = p, q = q, i2 = i2, k = k_eff)
}

meta_analysis_random <- function(effects, ses) {
  w <- 1 / (ses ^ 2)
  w[!is.finite(w)] <- 0
  effects[!is.finite(effects)] <- NA
  w[is.na(effects)] <- 0
  w_sum <- rowSums(w)
  k_eff <- rowSums(w > 0)
  beta_fixed <- rowSums(w * effects, na.rm = TRUE) / w_sum
  q <- rowSums(w * (effects - beta_fixed) ^ 2, na.rm = TRUE)
  c_val <- w_sum - rowSums(w ^ 2) / w_sum
  tau2 <- (q - (k_eff - 1)) / c_val
  tau2[!is.finite(tau2)] <- NA_real_
  tau2 <- pmax(0, tau2)
  w_star <- 1 / (ses ^ 2 + tau2)
  w_star[!is.finite(w_star)] <- 0
  w_star[is.na(effects)] <- 0
  w_star_sum <- rowSums(w_star)
  beta_random <- rowSums(w_star * effects, na.rm = TRUE) / w_star_sum
  se_random <- sqrt(1 / w_star_sum)
  z <- beta_random / se_random
  p <- 2 * pnorm(-abs(z))
  i2 <- ifelse(q > 0, pmax(0, (q - (k_eff - 1)) / q), NA_real_)
  beta_random[k_eff < 2] <- NA_real_
  se_random[k_eff < 2] <- NA_real_
  p[k_eff < 2] <- NA_real_
  q[k_eff < 2] <- NA_real_
  i2[k_eff < 2] <- NA_real_
  tau2[k_eff < 2] <- NA_real_
  list(beta = beta_random, se = se_random, p = p, q = q, i2 = i2, k = k_eff, tau2 = tau2)
}

run_stratified_meta_analysis <- function(betas, targets, batch_col, group_col, covariates,
                                         forced_covariates, prefix, out_dir, pval_thresh, lfc_thresh, delta_beta_thresh) {
  if (is.null(batch_col) || !(batch_col %in% colnames(targets))) return(NULL)
  overlap_batches <- identify_overlap_batches(targets, batch_col, group_col)
  if (length(overlap_batches) < 2) {
    message("  Stratified meta-analysis skipped: insufficient overlap strata.")
    return(NULL)
  }
  message(paste("  Tier3 fallback: stratified EWAS + meta-analysis across", length(overlap_batches), "batches."))
  strata <- overlap_batches
  effects <- list()
  ses <- list()
  for (b in strata) {
    idx <- targets[[batch_col]] == b
    sub_targets <- targets[idx, , drop = FALSE]
    sub_betas <- betas[, idx, drop = FALSE]
    group_levels <- levels(sub_targets[[group_col]])
    if (length(group_levels) < 2) next
    if (sum(sub_targets[[group_col]] == group_levels[1]) < 1 ||
        sum(sub_targets[[group_col]] == group_levels[2]) < 1) {
      next
    }
    sub_covars <- intersect(covariates, colnames(sub_targets))
    if (length(sub_covars) > 0) {
      single_filt <- drop_single_level_covariates(sub_targets, sub_covars)
      if (nrow(single_filt$dropped) > 0) {
        message(sprintf("  - Dropping single-level covariates in stratum %s: %s",
                        b, paste(single_filt$dropped$Variable, collapse = ", ")))
      }
      sub_covars <- single_filt$keep
    }
    build_design <- function(covars) {
      formula_str <- "~ 0 + primary_group"
      if (length(covars) > 0) formula_str <- paste(formula_str, "+", paste(covars, collapse = " + "))
      design <- model.matrix(as.formula(formula_str), data = sub_targets)
      colnames(design) <- gsub("primary_group", "", colnames(design))
      design <- design[, unique(colnames(design)), drop = FALSE]
      group_cols <- make.names(levels(sub_targets$primary_group))
      design_ld <- drop_linear_dependencies(design, group_cols = group_cols)
      list(mat = design_ld$mat)
    }
    covars_use <- sub_covars
    dropped_for_df <- character(0)
    forced_keep <- intersect(forced_covariates, covars_use)
    if (length(covars_use) > 0) {
      ranked <- rank_covariates(sub_targets, covars_use, covar_log_df = NULL)
      drop_non_forced <- setdiff(ranked, forced_keep)
      drop_forced <- intersect(ranked, forced_keep)
      drop_queue <- c(drop_non_forced, drop_forced)
    } else {
      drop_queue <- character(0)
    }
    design_obj <- build_design(covars_use)
    while (!is.null(design_obj$mat) && nrow(design_obj$mat) <= ncol(design_obj$mat) && length(covars_use) > 0) {
      if (length(drop_queue) == 0) break
      drop_cv <- drop_queue[1]
      drop_queue <- drop_queue[-1]
      covars_use <- setdiff(covars_use, drop_cv)
      dropped_for_df <- c(dropped_for_df, drop_cv)
      design_obj <- build_design(covars_use)
    }
    if (length(dropped_for_df) > 0) {
      message(sprintf("  - Reduced covariates in stratum %s to preserve df: dropped %s",
                      b, paste(dropped_for_df, collapse = ", ")))
      forced_dropped <- intersect(dropped_for_df, forced_keep)
      if (length(forced_dropped) > 0) {
        message(sprintf("  - Stratum %s: forced covariates dropped due to df (%s)",
                        b, paste(forced_dropped, collapse = ", ")))
      }
    }
    design <- design_obj$mat
    if (ncol(design) < 2) next
    if (nrow(design) <= ncol(design)) {
      message(sprintf("  - Skipping stratum %s: insufficient degrees of freedom after covariate reduction (n=%d, p=%d).", b, nrow(design), ncol(design)))
      next
    }
    group_cols <- intersect(make.names(levels(sub_targets$primary_group)), colnames(design))
    if (length(group_cols) < 2) {
      message(sprintf("  - Skipping stratum %s: missing group columns after design reduction.", b))
      next
    }
    m_vals <- logit_offset(sub_betas)
    fit <- lmFit(m_vals, design)
    contrast_str <- paste0(group_cols[2], "-", group_cols[1])
    cm <- makeContrasts(contrasts = contrast_str, levels = design)
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- safe_ebayes(fit2, label = paste0("stratum ", b))
    if (is.null(fit2)) next
    coef_idx <- 1
    eff <- fit2$coefficients[, coef_idx]
    se <- fit2$stdev.unscaled[, coef_idx] * fit2$sigma
    effects[[b]] <- eff
    ses[[b]] <- se
    out_df <- data.frame(CpG = rownames(fit2$coefficients), logFC = eff, SE = se)
    write.csv(out_df, file.path(out_dir, paste0(prefix, "_Stratum_", b, "_EWAS.csv")), row.names = FALSE)
  }
  if (length(effects) > 0) {
    keep <- vapply(names(effects), function(nm) {
      eff <- effects[[nm]]
      se <- ses[[nm]]
      if (is.null(eff) || is.null(se)) return(FALSE)
      if (length(eff) == 0 || length(se) == 0) return(FALSE)
      if (length(eff) != length(se)) return(FALSE)
      any(is.finite(eff) & is.finite(se))
    }, logical(1))
    if (any(!keep)) {
      dropped <- names(effects)[!keep]
      message(sprintf("  - Dropping stratum(s) with no finite effects/SE: %s", paste(dropped, collapse = ", ")))
    }
    effects <- effects[keep]
    ses <- ses[keep]
  }
  if (length(effects) < 2) {
    message("  Stratified meta-analysis skipped: fewer than 2 valid strata after QC.")
    log_decision("tier3", "meta_status", "skipped",
                 reason = "insufficient_valid_strata",
                 metrics = list(valid_strata = length(effects), total_strata = length(strata)))
    return(NULL)
  }
  eff_mat <- do.call(cbind, effects)
  se_mat <- do.call(cbind, ses)
  if (is.null(eff_mat) || is.null(se_mat) || ncol(eff_mat) < 2 || ncol(se_mat) < 2 || nrow(eff_mat) == 0) {
    message("  Stratified meta-analysis skipped: invalid effect/SE matrix after QC.")
    log_decision("tier3", "meta_status", "skipped",
                 reason = "invalid_effect_matrix",
                 metrics = list(valid_strata = length(effects), total_strata = length(strata)))
    return(NULL)
  }
  if (!any(rowSums(is.finite(eff_mat) & is.finite(se_mat)) >= 2)) {
    message("  Stratified meta-analysis skipped: no CpGs with >=2 valid strata.")
    log_decision("tier3", "meta_status", "skipped",
                 reason = "no_valid_cpgs",
                 metrics = list(valid_strata = length(effects), total_strata = length(strata)))
    return(NULL)
  }
  meta_fixed <- meta_analysis_fixed(eff_mat, se_mat)
  meta_random <- meta_analysis_random(eff_mat, se_mat)
  tier3_cfg <- if (exists("config_settings")) config_settings$tier3_meta else NULL
  meta_method_cfg <- if (!is.null(tier3_cfg$method)) as.character(tier3_cfg$method) else "auto"
  meta_method_cfg <- tolower(meta_method_cfg)
  if (!meta_method_cfg %in% c("auto", "fixed", "random")) meta_method_cfg <- "auto"
  i2_threshold <- suppressWarnings(as.numeric(if (!is.null(tier3_cfg$i2_threshold)) tier3_cfg$i2_threshold else 0.5))
  if (!is.finite(i2_threshold)) i2_threshold <- 0.5
  i2_median <- suppressWarnings(median(meta_random$i2, na.rm = TRUE))
  i2_mean <- suppressWarnings(mean(meta_random$i2, na.rm = TRUE))
  meta_method <- "fixed"
  if (meta_method_cfg == "random") {
    meta_method <- "random"
  } else if (meta_method_cfg == "fixed") {
    meta_method <- "fixed"
  } else {
    if (is.finite(i2_median) && i2_median >= i2_threshold) meta_method <- "random"
  }
  meta_use <- if (meta_method == "random") meta_random else meta_fixed
  log_decision("tier3", "meta_method", meta_method,
               reason = ifelse(meta_method_cfg == "auto", "auto_i2", meta_method_cfg),
               metrics = list(i2_median = i2_median, i2_mean = i2_mean, i2_threshold = i2_threshold))
  meta_df <- data.frame(
    CpG = rownames(betas),
    logFC = meta_use$beta,
    SE = meta_use$se,
    P.Value = meta_use$p,
    Q = meta_use$q,
    I2 = meta_use$i2,
    k = meta_use$k,
    tau2 = if (!is.null(meta_use$tau2)) meta_use$tau2 else NA_real_,
    meta_method = meta_method
  )
  meta_df$adj.P.Val <- p.adjust(meta_df$P.Value, "BH")
  meta_df <- meta_df[order(meta_df$P.Value), ]
  write.csv(meta_df, file.path(out_dir, paste0(prefix, "_Stratified_Meta_DMPs.csv")), row.names = FALSE)
  save_datatable(head(meta_df, min(10000, nrow(meta_df))), paste0(prefix, "_Stratified_Meta_DMPs.html"), out_dir)
  top_n <- min(10, nrow(meta_df))
  if (top_n > 0) {
    top_cpgs <- meta_df$CpG[1:top_n]
    forest_df <- data.frame()
    for (b in names(effects)) {
      eff_vec <- effects[[b]][top_cpgs]
      se_vec <- ses[[b]][top_cpgs]
      forest_df <- rbind(forest_df, data.frame(CpG = top_cpgs, Stratum = b, logFC = eff_vec, SE = se_vec))
    }
    forest_df <- rbind(forest_df, data.frame(CpG = top_cpgs, Stratum = "META", logFC = meta_df$logFC[1:top_n], SE = meta_df$SE[1:top_n]))
    p_forest <- ggplot(forest_df, aes(x = logFC, y = Stratum)) +
      geom_point() +
      geom_errorbarh(aes(xmin = logFC - 1.96 * SE, xmax = logFC + 1.96 * SE), height = 0.2) +
      facet_wrap(~ CpG, scales = "free_x") +
      theme_minimal() +
      labs(title = paste(prefix, "Meta-analysis Forest Plot"), x = "Effect (logFC)", y = "Stratum")
    save_interactive_plot(p_forest, paste0(prefix, "_Meta_Forest.html"), out_dir)
    save_static_plot(p_forest, paste0(prefix, "_Meta_Forest.png"), out_dir, width = 8, height = 5)
  }
  if (any(is.finite(meta_df$I2))) {
    i2_df <- data.frame(I2 = meta_df$I2)
    p_i2 <- ggplot(i2_df, aes(x = I2)) +
      geom_histogram(bins = 30, fill = "#4c78a8", color = "white") +
      theme_minimal() +
      labs(title = paste(prefix, "Heterogeneity (I2)"), x = "I2", y = "Count")
    save_interactive_plot(p_i2, paste0(prefix, "_Meta_I2.html"), out_dir)
    save_static_plot(p_i2, paste0(prefix, "_Meta_I2.png"), out_dir, width = 6, height = 4)
  }
  strata_used <- names(effects)
  stratum_sizes <- if (length(strata_used) > 0) {
    vapply(strata_used, function(b) sum(targets[[batch_col]] == b), integer(1))
  } else {
    integer(0)
  }
  list(
    res = meta_df,
    method = meta_method,
    i2_median = i2_median,
    i2_mean = i2_mean,
    n_strata = length(strata_used),
    min_stratum_n = ifelse(length(stratum_sizes) > 0, min(stratum_sizes), NA_integer_),
    median_stratum_n = ifelse(length(stratum_sizes) > 0, median(stratum_sizes), NA_integer_),
    total_n = ifelse(length(stratum_sizes) > 0, sum(stratum_sizes), NA_integer_),
    strata_used = strata_used
  )
}

check_tier3_eligibility <- function(targets, batch_col, group_col,
                                    min_total_n = 20,
                                    min_per_group_per_stratum = 5,
                                    out_dir = NULL) {
  if (is.null(batch_col) || !(batch_col %in% colnames(targets))) {
    return(list(eligible = FALSE, reason = "missing_batch_column"))
  }
  if (!(group_col %in% colnames(targets))) {
    return(list(eligible = FALSE, reason = "missing_group_column"))
  }
  counts <- as.data.frame(table(targets[[batch_col]], targets[[group_col]]), stringsAsFactors = FALSE)
  colnames(counts) <- c("Batch", "Group", "N")
  counts$N <- as.integer(counts$N)
  counts$Pass_Group_Min <- counts$N >= min_per_group_per_stratum
  total_n <- nrow(targets)
  # Only check overlap batches (batches containing both groups) for eligibility;
  # non-overlap batches (one group only) are excluded from stratified analysis anyway.
  groups_in_data <- unique(as.character(counts$Group[counts$N > 0]))
  n_groups <- length(groups_in_data)
  overlap_batches <- if (n_groups >= 2) {
    grp_per_batch <- tapply(counts$N > 0, counts$Batch, sum)
    names(grp_per_batch)[grp_per_batch >= n_groups]
  } else character(0)
  overlap_counts <- counts[counts$Batch %in% overlap_batches, , drop = FALSE]
  min_stratum_n <- if (nrow(overlap_counts) > 0) min(overlap_counts$N) else 0
  batch_ok <- tapply(overlap_counts$Pass_Group_Min, overlap_counts$Batch, all)
  pass_total <- total_n >= min_total_n
  pass_strata <- length(batch_ok) > 0 && all(batch_ok)
  eligible <- isTRUE(pass_total && pass_strata)
  counts$Batch_Pass <- batch_ok[match(counts$Batch, names(batch_ok))]
  counts$Total_N <- total_n
  counts$Min_Group_N <- min_stratum_n
  counts$Eligible <- eligible
  counts$Fail_Reason <- ""
  if (!pass_total) counts$Fail_Reason <- "total_n_lt_min"
  if (!pass_strata) {
    counts$Fail_Reason <- paste0(counts$Fail_Reason, ifelse(nzchar(counts$Fail_Reason), ";", ""), "stratum_group_lt_min")
  }
  if (!is.null(out_dir)) {
    write.csv(counts, file.path(out_dir, "Tier3_Eligibility.csv"), row.names = FALSE)
  }
  list(eligible = eligible, total_n = total_n, min_stratum_n = min_stratum_n,
       counts = counts, reason = ifelse(eligible, "eligible",
                                        ifelse(!pass_total, "total_n_lt_min", "stratum_group_lt_min")))
}

emit_tier3_primary_outputs <- function(meta_res, betas, targets, curr_anno, prefix, out_dir,
                                       group_con_in, group_test_in, clean_con, clean_test,
                                       max_points, pval_thresh, lfc_thresh,
                                       meta_method = "", i2_median = NA_real_,
                                       tier3_stats = NULL) {
  if (is.null(meta_res) || nrow(meta_res) == 0) return(NULL)
  if (!("CpG" %in% colnames(meta_res))) {
    meta_res$CpG <- rownames(meta_res)
  }
  res <- meta_res
  con_mask <- targets$primary_group == clean_con
  test_mask <- targets$primary_group == clean_test
  mean_beta_con <- rowMeans(betas[, con_mask, drop = FALSE], na.rm = TRUE)
  mean_beta_test <- rowMeans(betas[, test_mask, drop = FALSE], na.rm = TRUE)
  delta_beta <- mean_beta_test - mean_beta_con
  res$Mean_Beta_Con <- mean_beta_con[res$CpG]
  res$Mean_Beta_Test <- mean_beta_test[res$CpG]
  res$Delta_Beta <- delta_beta[res$CpG]
  if (!is.null(curr_anno)) {
    res <- merge(res, curr_anno, by.x = "CpG", by.y = "CpG", all.x = TRUE)
  }
  clean_gene_names <- function(g_str) {
    if (is.na(g_str) || g_str == "") return("")
    genes <- unlist(strsplit(as.character(g_str), ";"))
    genes <- unique(genes)
    if (length(genes) > 2) genes <- genes[1:2]
    final_str <- paste(genes, collapse = ";")
    if (nchar(final_str) > 15) final_str <- paste0(substr(final_str, 1, 12), "...")
    final_str
  }
  if ("Gene" %in% colnames(res)) {
    res$Gene <- sapply(res$Gene, clean_gene_names)
  }
  
  file_prefix <- paste0(prefix, "_Tier3_Primary")
  meta_method_disp <- ifelse(nzchar(meta_method), meta_method, "fixed")
  i2_disp <- ifelse(is.finite(i2_median), sprintf("%.3f", i2_median), "NA")
  note_lines <- c(
    "Tier3 confounding detected: primary inference uses stratified EWAS + meta-analysis.",
    sprintf("Tier3 meta-analysis method: %s (I2 median=%s).", meta_method_disp, i2_disp),
    "These Tier3_Primary outputs are the recommended results for interpretation.",
    "Standard (non-stratified) outputs are provided as sensitivity checks only."
  )
  min_total_warn <- suppressWarnings(as.numeric(if (!is.null(config_settings$tier3_meta$min_total_warn))
    config_settings$tier3_meta$min_total_warn else 20))
  if (!is.finite(min_total_warn)) min_total_warn <- 20
  min_stratum_warn <- suppressWarnings(as.numeric(if (!is.null(config_settings$tier3_meta$min_stratum_warn))
    config_settings$tier3_meta$min_stratum_warn else 6))
  if (!is.finite(min_stratum_warn)) min_stratum_warn <- 6
  if (!is.null(tier3_stats)) {
    if (is.finite(tier3_stats$total_n) && tier3_stats$total_n < min_total_warn) {
      note_lines <- c(note_lines,
                      sprintf("WARNING: Tier3 total N=%d below recommended minimum %d; interpret with caution.",
                              tier3_stats$total_n, min_total_warn))
    }
    if (is.finite(tier3_stats$min_stratum_n) && tier3_stats$min_stratum_n < min_stratum_warn) {
      note_lines <- c(note_lines,
                      sprintf("WARNING: Tier3 smallest stratum N=%d below recommended minimum %d; power may be limited.",
                              tier3_stats$min_stratum_n, min_stratum_warn))
    }
  }
  writeLines(note_lines, file.path(out_dir, paste0(file_prefix, "_Notes.txt")))
  
  p_vals <- res$P.Value
  p_vals <- p_vals[!is.na(p_vals)]
  n_p <- length(p_vals)
  lambda_val <- NA_real_
  if (n_p > 1) {
    lambda_val <- compute_genomic_lambda(p_vals)
    expected <- -log10(ppoints(n_p))
    observed <- -log10(pmax(sort(p_vals), .Machine$double.xmin))
    qq_df <- data.frame(Expected = expected, Observed = observed)
    if (nrow(qq_df) > max_points) {
      idx <- unique(c(1:1000, seq(1001, n_p, length.out = max_points)))
      qq_df <- qq_df[idx, ]
    }
    p_qq <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      geom_point(alpha = 0.5, size = 1) +
      theme_minimal() +
      labs(x = "Expected -log10(P)", y = "Observed -log10(P)") +
      ggtitle(paste(prefix, "Tier3 Primary Q-Q (Lambda =", round(lambda_val, 3), ")"))
    save_interactive_plot(p_qq, paste0(file_prefix, "_QQPlot.html"), out_dir)
    save_static_plot(p_qq, paste0(file_prefix, "_QQPlot.png"), out_dir, width = 5, height = 4)
  }
  
  plot_res <- res[order(res$P.Value), ]
  if (nrow(plot_res) > max_points) {
    plot_res <- plot_res[1:max_points, ]
  }
  plot_res$diffexpressed <- "NO"
  delta_pass <- delta_beta_pass(plot_res$Delta_Beta)
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC > lfc_thresh & delta_pass] <- "UP"
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC < -lfc_thresh & delta_pass] <- "DOWN"
  subtitle_str <- paste0("Control: ", group_con_in, " vs Test: ", group_test_in, " (Tier3 stratified)")
  p_vol <- ggplot(plot_res, aes(x = logFC, y = -log10(P.Value), color = diffexpressed,
                                text = paste("CpG:", CpG,
                                             "<br>Gene:", Gene,
                                             "<br>DeltaBeta:", round(Delta_Beta, 4),
                                             "<br>Region:", Region,
                                             "<br>Island:", Island_Context))) +
    geom_point(alpha = 0.5) + theme_minimal() +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
    labs(title = paste(prefix, "Tier3 Primary Volcano"), subtitle = subtitle_str, color = "Status")
  save_interactive_plot(p_vol, paste0(file_prefix, "_Volcano.html"), out_dir)
  save_static_plot(p_vol, paste0(file_prefix, "_Volcano.png"), out_dir, width = 6, height = 5)
  
  if ("chr" %in% colnames(plot_res) && "pos" %in% colnames(plot_res)) {
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
      chr_labels <- names(chr_vals)[match(axis_set$chr_num, chr_vals)]
      p_man <- ggplot(plot_res_man, aes(x = pos_cum, y = -log10(P.Value), color = as.factor(chr_num),
                                        text = paste("CpG:", CpG,
                                                     "<br>Gene:", Gene,
                                                     "<br>DeltaBeta:", round(Delta_Beta, 4),
                                                     "<br>Region:", Region,
                                                     "<br>Island:", Island_Context))) +
        geom_point(alpha = 0.7, size = 1) +
        scale_x_continuous(label = chr_labels, breaks = axis_set$center) +
        theme_minimal() +
        theme(legend.position = "none",
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)) +
        xlab("Chromosome") +
        ggtitle(paste(prefix, "Tier3 Primary Manhattan"))
      save_interactive_plot(p_man, paste0(file_prefix, "_Manhattan.html"), out_dir)
      save_static_plot(p_man, paste0(file_prefix, "_Manhattan.png"), out_dir, width = 8, height = 4)
    }
  }
  
  res <- res[order(res$P.Value, na.last = TRUE), ]
  write.csv(res, file.path(out_dir, paste0(file_prefix, "_DMPs.csv")), row.names = FALSE)
  save_datatable(head(res, min(10000, nrow(res))), paste0(file_prefix, "_DMPs.html"), out_dir)
  
  list(res = res, lambda = lambda_val)
}

run_pooled_batch_model <- function(betas, targets, batch_col, covariates, prefix, out_dir) {
  if (is.null(batch_col) || !(batch_col %in% colnames(targets))) return(NULL)
  if (length(unique(targets[[batch_col]])) < 2) return(NULL)
  covars <- intersect(covariates, colnames(targets))
  formula_str <- paste("~ 0 + primary_group +", batch_col)
  if (length(covars) > 0) formula_str <- paste(formula_str, "+", paste(covars, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = targets)
  colnames(design) <- gsub("primary_group", "", colnames(design))
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  group_cols <- make.names(levels(targets$primary_group))
  group_cols <- intersect(group_cols, colnames(design))
  if (length(group_cols) < 2) return(NULL)
  design_ld <- drop_linear_dependencies(design, group_cols = group_cols)
  design <- design_ld$mat
  if (ncol(design) < 2) return(NULL)
  m_vals <- logit_offset(betas)
  mval_out_path <- file.path(out_dir, paste0(prefix, "_MvalueMatrix.tsv.gz"))
  write_matrix_tsv_gz(m_vals, mval_out_path)
  fit <- lmFit(m_vals, design)
  if (!is.null(fit$df.residual) && all(is.finite(fit$df.residual)) && all(fit$df.residual <= 0)) {
    message("  Pooled batch model skipped: no residual degrees of freedom.")
    return(NULL)
  }
  contrast_str <- paste0(group_cols[2], "-", group_cols[1])
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- safe_ebayes(fit2, "pooled_batch_model")
  if (is.null(fit2)) return(NULL)
  res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  res$CpG <- rownames(res)
  write.csv(res, file.path(out_dir, paste0(prefix, "_Pooled_Batch_DMPs.csv")), row.names = FALSE)
  save_datatable(head(res, min(10000, nrow(res))), paste0(prefix, "_Pooled_Batch_DMPs.html"), out_dir)
  res
}

run_overlap_restricted_analysis <- function(betas, targets, batch_col, group_col, covariates, prefix, out_dir) {
  overlap_batches <- identify_overlap_batches(targets, batch_col, group_col)
  if (length(overlap_batches) == 0) return(NULL)
  idx <- targets[[batch_col]] %in% overlap_batches
  sub_targets <- targets[idx, , drop = FALSE]
  sub_betas <- betas[, idx, drop = FALSE]
  if (!(group_col %in% colnames(sub_targets))) return(NULL)
  sub_targets[[group_col]] <- as.factor(sub_targets[[group_col]])
  group_counts <- table(sub_targets[[group_col]])
  if (length(group_counts) < 2) {
    message("  Overlap restricted analysis skipped: <2 groups after overlap filter.")
    return(NULL)
  }
  if (any(group_counts < 2)) {
    message("  Overlap restricted analysis skipped: insufficient samples per group in overlap strata.")
    return(NULL)
  }
  covars <- intersect(covariates, colnames(sub_targets))
  formula_str <- "~ 0 + primary_group"
  if (length(covars) > 0) formula_str <- paste(formula_str, "+", paste(covars, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = sub_targets)
  colnames(design) <- gsub("primary_group", "", colnames(design))
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  group_cols <- make.names(levels(sub_targets$primary_group))
  group_cols <- intersect(group_cols, colnames(design))
  if (length(group_cols) < 2) return(NULL)
  design_ld <- drop_linear_dependencies(design, group_cols = group_cols)
  design <- design_ld$mat
  if (ncol(design) < 2) return(NULL)
  m_vals <- logit_offset(sub_betas)
  fit <- lmFit(m_vals, design)
  if (!is.null(fit$df.residual) && all(is.finite(fit$df.residual)) && all(fit$df.residual <= 0)) {
    message("  Overlap restricted analysis skipped: no residual degrees of freedom.")
    return(NULL)
  }
  contrast_str <- paste0(group_cols[2], "-", group_cols[1])
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- safe_ebayes(fit2, "overlap_restricted")
  if (is.null(fit2)) return(NULL)
  res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  res$CpG <- rownames(res)
  write.csv(res, file.path(out_dir, paste0(prefix, "_Overlap_Restricted_DMPs.csv")), row.names = FALSE)
  save_datatable(head(res, min(10000, nrow(res))), paste0(prefix, "_Overlap_Restricted_DMPs.html"), out_dir)
  res
}

run_limma_dmp <- function(betas, design, group_levels, m_vals = NULL) {
  if (is.null(design) || ncol(design) < 2) return(NULL)
  if (length(group_levels) < 2) return(NULL)
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  group_levels <- make.names(group_levels)
  group_levels <- intersect(group_levels, colnames(design))
  if (length(group_levels) < 2) return(NULL)
  contrast_str <- paste0(group_levels[2], "-", group_levels[1])
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  if (is.null(m_vals)) {
    m_vals <- logit_offset(betas)
  }
  fit <- lmFit(m_vals, design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- safe_ebayes(fit2, "limma_dmp")
  if (is.null(fit2)) return(NULL)
  res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  res$CpG <- rownames(res)
  res
}

run_lambda_guard <- function(betas, targets, group_col, prefix, out_dir, max_points, pval_thresh, lfc_thresh) {
  if (is.null(betas) || is.null(targets)) return(list(status = "failed", lambda = NA_real_))
  if (!(group_col %in% colnames(targets))) return(list(status = "failed", lambda = NA_real_))
  targets_use <- targets
  targets_use[[group_col]] <- as.factor(targets_use[[group_col]])
  if (length(levels(targets_use[[group_col]])) < 2) {
    return(list(status = "failed", lambda = NA_real_))
  }
  design <- model.matrix(as.formula(paste("~ 0 +", group_col)), data = targets_use)
  colnames(design) <- gsub(group_col, "", colnames(design), fixed = TRUE)
  design <- drop_linear_dependencies(design, group_cols = colnames(design))$mat
  if (ncol(design) < 2) return(list(status = "failed", lambda = NA_real_))
  group_cols <- colnames(design)
  res <- tryCatch(run_limma_dmp(betas, design, group_cols), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(list(status = "failed", lambda = NA_real_))
  p_vals <- res$P.Value
  p_vals <- p_vals[!is.na(p_vals)]
  if (length(p_vals) < 2) return(list(status = "failed", lambda = NA_real_))
  lambda_val <- compute_genomic_lambda(p_vals)
  if (!is.finite(lambda_val)) return(list(status = "failed", lambda = NA_real_))
  
  expected <- -log10(ppoints(length(p_vals)))
  observed <- -log10(pmax(sort(p_vals), .Machine$double.xmin))
  qq_df <- data.frame(Expected = expected, Observed = observed)
  if (nrow(qq_df) > max_points) {
    idx <- unique(c(1:1000, seq(1001, nrow(qq_df), length.out = max_points)))
    qq_df <- qq_df[idx, ]
  }
  p_qq <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_point(alpha = 0.5, size = 1) +
    theme_minimal() +
    labs(x = "Expected -log10(P)", y = "Observed -log10(P)") +
    ggtitle(paste(prefix, "Lambda Guard Q-Q (Lambda =", round(lambda_val, 3), ")"))
  save_interactive_plot(p_qq, paste0(prefix, "_LambdaGuard_QQPlot.html"), out_dir)
  save_static_plot(p_qq, paste0(prefix, "_LambdaGuard_QQPlot.png"), out_dir, width = 5, height = 4)
  
  res_out <- res[order(res$P.Value), ]
  write.csv(res_out, file.path(out_dir, paste0(prefix, "_LambdaGuard_DMPs.csv")), row.names = FALSE)
  save_datatable(head(res_out, min(10000, nrow(res_out))), paste0(prefix, "_LambdaGuard_DMPs.html"), out_dir)
  
  n_sig <- sum(res_out$adj.P.Val < pval_thresh & abs(res_out$logFC) > lfc_thresh, na.rm = TRUE)
  metrics <- data.frame(
    metric = c("lambda_guard_lambda", "lambda_guard_n_sig", "lambda_guard_n_cpgs", "lambda_guard_n_samples"),
    value = c(lambda_val, n_sig, nrow(res_out), nrow(targets_use))
  )
  write.csv(metrics, file.path(out_dir, paste0(prefix, "_LambdaGuard_Metrics.csv")), row.names = FALSE)
  
  list(status = "ok", lambda = lambda_val, n_sig = n_sig)
}

run_dmrff_with_inputs <- function(dmr_df, betas, targets, curr_anno, prefix, out_dir,
                                  dmr_p_cutoff, dmr_min_cpgs, con_label, test_label) {
  dmr_res <- data.frame()
  status <- "ok"
  reason <- ""
  if (is.null(dmr_df) || nrow(dmr_df) == 0) {
    return(list(res = dmr_res, status = "skipped_missing_estimates", reason = "missing_estimates"))
  }
  if (is.null(curr_anno) || !all(c("chr", "pos") %in% colnames(curr_anno))) {
    return(list(res = dmr_res, status = "skipped_missing_annotation", reason = "missing_annotation"))
  }
  if (nrow(targets) < 4) {
    return(list(res = dmr_res, status = "skipped_small_n", reason = "n_samples<4"))
  }
  if (!all(c("CpG", "logFC", "SE", "P.Value") %in% colnames(dmr_df))) {
    return(list(res = dmr_res, status = "skipped_missing_columns", reason = "missing_columns"))
  }
  dmr_df <- dmr_df[!is.na(dmr_df$CpG), , drop = FALSE]
  idx <- match(dmr_df$CpG, rownames(betas))
  keep <- !is.na(idx)
  if (sum(keep) < 100) {
    return(list(res = dmr_res, status = "skipped_too_few_cpgs", reason = "insufficient_matching_cpgs"))
  }
  dmr_df <- dmr_df[keep, , drop = FALSE]
  idx <- idx[keep]
  betas_dmr <- betas[idx, , drop = FALSE]
  dmr_anno <- curr_anno[match(dmr_df$CpG, rownames(curr_anno)), ]
  valid <- is.finite(dmr_df$logFC) & is.finite(dmr_df$SE) & is.finite(dmr_df$P.Value) &
    !is.na(dmr_anno$chr) & !is.na(dmr_anno$pos)
  if (sum(valid) < 100) {
    return(list(res = dmr_res, status = "skipped_too_few_cpgs", reason = "insufficient_valid_cpgs"))
  }
  dmr_df <- dmr_df[valid, , drop = FALSE]
  betas_dmr <- betas_dmr[valid, , drop = FALSE]
  dmr_anno <- dmr_anno[valid, , drop = FALSE]
  con_mask <- targets$primary_group == con_label
  test_mask <- targets$primary_group == test_label
  dmr_mean_con <- rowMeans(betas_dmr[, con_mask, drop = FALSE], na.rm = TRUE)
  dmr_mean_test <- rowMeans(betas_dmr[, test_mask, drop = FALSE], na.rm = TRUE)
  dmr_est <- dmr_df$logFC
  dmr_se <- dmr_df$SE
  dmr_p <- dmr_df$P.Value
  if (length(dmr_est) < 100) {
    return(list(res = dmr_res, status = "skipped_too_few_cpgs", reason = "fewer_than_100_cpgs"))
  }

	  annotate_dmr <- function(df, anno_tbl) {
	    if (nrow(df) == 0) return(df)
	    anno_tbl$Gene <- ifelse(is.na(anno_tbl$Gene), "", anno_tbl$Gene)
	    anno_tbl$Region <- ifelse(is.na(anno_tbl$Region), "", anno_tbl$Region)
	    collapse_unique_tokens <- function(vals, max_n = 5) {
	      vals <- as.character(vals)
	      vals[is.na(vals)] <- ""
	      vals <- vals[nzchar(vals)]
	      if (length(vals) == 0) return("")
	      toks <- unlist(strsplit(vals, "[;,]"))
	      toks <- trimws(toks)
	      toks <- toks[nzchar(toks)]
	      if (length(toks) == 0) return("")
	      toks <- unique(toks)
	      paste(head(toks, max_n), collapse = ";")
	    }
	    gene_region <- function(chr, start, end, field) {
	      idx <- which(anno_tbl$chr == chr & anno_tbl$pos >= start & anno_tbl$pos <= end)
	      if (length(idx) == 0) return("")
	      collapse_unique_tokens(anno_tbl[idx, field], max_n = 5)
	    }
	    df$Genes <- mapply(gene_region, df$chr, df$start, df$end, MoreArgs = list(field = "Gene"))
	    df$Regions <- mapply(gene_region, df$chr, df$start, df$end, MoreArgs = list(field = "Region"))
	    df
	  }

  compute_dmr_delta_beta <- function(df, anno_tbl, mean_con, mean_test) {
    if (nrow(df) == 0) return(df)
    anno_tbl$CpG <- rownames(anno_tbl)
    dmr_delta <- numeric(nrow(df))
    for (i in seq_len(nrow(df))) {
      chr <- df$chr[i]
      start <- df$start[i]
      end <- df$end[i]
      idx <- which(anno_tbl$chr == chr & anno_tbl$pos >= start & anno_tbl$pos <= end)
      if (length(idx) > 0) {
        cpgs <- anno_tbl$CpG[idx]
        cpgs <- cpgs[!is.na(cpgs)]
        if (length(cpgs) > 0) {
          dcon <- mean_con[cpgs]
          dtest <- mean_test[cpgs]
          dmr_delta[i] <- mean(dtest - dcon, na.rm = TRUE)
        } else {
          dmr_delta[i] <- NA_real_
        }
      } else {
        dmr_delta[i] <- NA_real_
      }
    }
    df$Delta_Beta <- dmr_delta
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
      p.cutoff = dmr_p_cutoff
    )

    if (nrow(dmr_res) > 0) {
      if (dmr_min_cpgs > 1) {
        dmr_before <- nrow(dmr_res)
        dmr_res <- dmr_res[dmr_res$n >= dmr_min_cpgs, , drop = FALSE]
        message(paste("    - Filtered DMRs by min CpGs >=", dmr_min_cpgs, ":", dmr_before, "->", nrow(dmr_res)))
      }
      if (nrow(dmr_res) > 0) {
        dmr_res <- annotate_dmr(dmr_res, curr_anno)
        dmr_res <- compute_dmr_delta_beta(dmr_res, dmr_anno, dmr_mean_con, dmr_mean_test)
        dmr_res <- dmr_res[order(dmr_res$p.value, dmr_res$p.adjust, na.last = TRUE), ]
        message(paste("    - Found", nrow(dmr_res), "DMRs."))

        dmr_filename <- paste0(prefix, "_DMRs_Table.html")
        save_datatable(head(dmr_res, 3000), dmr_filename, out_dir, sort_col = "p.value")
        write.csv(dmr_res, file.path(out_dir, paste0(prefix, "_DMRs.csv")), row.names = FALSE)

	        # Prefer Delta_Beta for volcano x-axis (interpretable, bounded). Fall back to dmrff estimate.
	        x_col <- if ("Delta_Beta" %in% colnames(dmr_res) && any(is.finite(dmr_res$Delta_Beta))) "Delta_Beta" else "estimate"
	        x_label <- if (identical(x_col, "Delta_Beta")) "Delta beta (mean test - control)" else "dmrff estimate"
	        dmr_res$x_effect <- dmr_res[[x_col]]
	        
	        dmr_res$signif_status <- "NO"
	        dmr_res$signif_status[dmr_res$p.value < dmr_p_cutoff & is.finite(dmr_res$x_effect) & dmr_res$x_effect > 0] <- "UP"
	        dmr_res$signif_status[dmr_res$p.value < dmr_p_cutoff & is.finite(dmr_res$x_effect) & dmr_res$x_effect < 0] <- "DOWN"
	        p_dmr_volcano <- ggplot(dmr_res, aes(x = x_effect, y = -log10(p.value), color = signif_status,
	                                             text = paste(
	                                               "Region:", paste0(chr, ":", start, "-", end),
	                                               "<br>Genes:", ifelse(nzchar(Genes), Genes, "NA"),
	                                               "<br>Regions:", ifelse(nzchar(Regions), Regions, "NA"),
	                                               "<br>nCpG:", n,
	                                               "<br>Estimate (dmrff):", signif(estimate, 4),
	                                               "<br>DeltaBeta:", ifelse(is.finite(Delta_Beta), signif(Delta_Beta, 4), "NA"),
	                                               "<br>P:", signif(p.value, 3),
	                                               "<br>FDR:", signif(p.adjust, 3)
	                                             ))) +
	          geom_point(alpha = 0.75, size = 2) +
	          scale_color_manual(values = c("DOWN" = "#3B82F6", "NO" = "#9CA3AF", "UP" = "#EF4444")) +
	          theme_minimal() +
	          labs(title = paste(prefix, "DMR Volcano Plot"),
	               x = x_label, y = "-log10(p)", color = "Status")

        save_interactive_plot(p_dmr_volcano, paste0(prefix, "_DMR_Volcano.html"), out_dir)
        save_static_plot(p_dmr_volcano, paste0(prefix, "_DMR_Volcano.png"), out_dir, width = 6, height = 4)

        if (nrow(dmr_res) > 0) {
          top_label_n <- min(DMR_LABEL_TOP_N, nrow(dmr_res))
          top_dmrs <- dmr_res[1:top_label_n, , drop = FALSE]
          p_dmr_scatter <- ggplot(top_dmrs, aes(x = estimate, y = -log10(p.value), label = Genes)) +
            geom_point(color = "#e67e22") +
            ggrepel::geom_text_repel(size = 3, max.overlaps = 10) +
            theme_minimal() +
            labs(title = paste(prefix, "Top DMRs"),
                 x = "dmrff estimate", y = "-log10(p)")
          save_interactive_plot(p_dmr_scatter, paste0(prefix, "_DMR_Top_Labelled.html"), out_dir)
          save_static_plot(p_dmr_scatter, paste0(prefix, "_DMR_Top_Labelled.png"), out_dir, width = 6, height = 4)
        }

	        if (nrow(dmr_res) > 0 && all(c("chr", "start", "end") %in% colnames(dmr_res))) {
	          chr_map <- c(as.character(1:22), "X", "Y", "M")
	          chr_vals <- 1:25
	          names(chr_vals) <- chr_map
	          clean_chr <- gsub("chr", "", dmr_res$chr)
	          dmr_res$chr_num <- chr_vals[clean_chr]
	          plot_dmr_man <- dmr_res[!is.na(dmr_res$chr_num) & !is.na(dmr_res$start) & !is.na(dmr_res$end), ]
	          plot_dmr_man <- plot_dmr_man[order(plot_dmr_man$chr_num, plot_dmr_man$start), ]
	          chr_len <- tapply(plot_dmr_man$end, plot_dmr_man$chr_num, max)
	          chr_len <- chr_len[!is.na(chr_len)]
	          if (length(chr_len) > 0) {
	            cum_lengths <- cumsum(as.numeric(chr_len))
	            offsets <- c(0, cum_lengths[-length(cum_lengths)])
	            names(offsets) <- names(chr_len)
	            plot_dmr_man$pos_cum <- plot_dmr_man$start + offsets[as.character(plot_dmr_man$chr_num)]
	            axis_set_dmr <- plot_dmr_man %>%
	              group_by(chr_num) %>%
	              summarize(center = (min(pos_cum) + max(pos_cum)) / 2)
	            chr_labels_dmr <- names(chr_vals)[match(axis_set_dmr$chr_num, chr_vals)]
	            plot_dmr_man$chr_parity <- as.factor(plot_dmr_man$chr_num %% 2)
	            p_dmr_manhattan <- ggplot(plot_dmr_man, aes(x = pos_cum, y = -log10(p.value),
	                                                       color = chr_parity,
	                                                       text = paste("Region:", paste0(chr, ":", start, "-", end),
	                                                                    "<br>Genes:", ifelse(nzchar(Genes), Genes, "NA"),
	                                                                    "<br>Regions:", ifelse(nzchar(Regions), Regions, "NA"),
	                                                                    "<br>nCpG:", n,
	                                                                    "<br>Estimate (dmrff):", signif(estimate, 4),
	                                                                    "<br>DeltaBeta:", ifelse(is.finite(Delta_Beta), signif(Delta_Beta, 4), "NA"),
	                                                                    "<br>P:", signif(p.value, 3),
	                                                                    "<br>FDR:", signif(p.adjust, 3)))) +
	              geom_point(alpha = 0.75, size = 1.2) +
	              scale_color_manual(values = c("0" = "#64748B", "1" = "#0EA5E9")) +
	              scale_x_continuous(label = chr_labels_dmr, breaks = axis_set_dmr$center) +
	            theme_minimal() +
	            theme(legend.position = "none",
	                  panel.grid.major.x = element_blank(),
	                  panel.grid.minor.x = element_blank(),
	                  axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)) +
	            labs(title = paste(prefix, "DMR Manhattan Plot"),
	                 x = "Chromosome", y = "-log10(p)")
	            save_interactive_plot(p_dmr_manhattan, paste0(prefix, "_DMR_Manhattan.html"), out_dir)
	            save_static_plot(p_dmr_manhattan, paste0(prefix, "_DMR_Manhattan.png"), out_dir, width = 8, height = 4)
	          }
	        }

        top_dmrs <- head(dmr_res, 50)
        if (nrow(top_dmrs) > 0) {
          dmr_means <- matrix(NA_real_, nrow = nrow(top_dmrs), ncol = ncol(betas_dmr))
          rownames(dmr_means) <- paste0(top_dmrs$chr, ":", top_dmrs$start, "-", top_dmrs$end)
          colnames(dmr_means) <- colnames(betas_dmr)
          for (i in seq_len(nrow(top_dmrs))) {
            idx <- which(dmr_anno$chr == top_dmrs$chr[i] &
                           dmr_anno$pos >= top_dmrs$start[i] &
                           dmr_anno$pos <= top_dmrs$end[i])
            if (length(idx) > 0) {
              cpgs_in_region <- dmr_df$CpG[idx]
              if (length(cpgs_in_region) > 0) {
                dmr_means[i, ] <- colMeans(betas_dmr[match(cpgs_in_region, rownames(betas_dmr)), , drop = FALSE], na.rm = TRUE)
              }
            }
          }
          keep_rows <- rowSums(is.finite(dmr_means)) > 0
          if (any(keep_rows)) {
            dmr_means <- dmr_means[keep_rows, , drop = FALSE]
            heatmap_df <- as.data.frame(dmr_means)
            heatmap_df$Region <- rownames(dmr_means)
            heatmap_dt <- data.table::as.data.table(heatmap_df)
	            heatmap_df_melt <- data.table::melt(heatmap_dt, id.vars = "Region",
	                                                variable.name = "Sample",
	                                                value.name = "Methylation")
	            heatmap_df_melt$Group <- NA_character_
	            if ("primary_group" %in% colnames(targets)) {
	              heatmap_df_melt$Group <- as.character(targets$primary_group[match(heatmap_df_melt$Sample, targets[[gsm_col]])])
	              if (any(is.na(heatmap_df_melt$Group))) {
	                  alt <- as.character(targets$primary_group[match(heatmap_df_melt$Sample, rownames(targets))])
	                  heatmap_df_melt$Group[is.na(heatmap_df_melt$Group)] <- alt[is.na(heatmap_df_melt$Group)]
	              }
	            }
	            
	            # Order samples by group (Control -> Test) when possible.
	            sample_ids <- colnames(betas_dmr)
	            grp_vec <- as.character(targets$primary_group[match(sample_ids, targets[[gsm_col]])])
	            if (any(is.na(grp_vec))) {
	              alt <- as.character(targets$primary_group[match(sample_ids, rownames(targets))])
	              grp_vec[is.na(grp_vec)] <- alt[is.na(grp_vec)]
	            }
	            grp_levels <- unique(as.character(targets$primary_group))
	            if (!is.null(con_label) && !is.null(test_label) &&
	                con_label %in% grp_levels && test_label %in% grp_levels) {
	              grp_levels <- c(con_label, test_label, setdiff(grp_levels, c(con_label, test_label)))
	            }
	            if (length(grp_levels) >= 2 && !all(is.na(grp_vec))) {
	              sample_ids <- c(sample_ids[grp_vec == grp_levels[1]],
	                              sample_ids[grp_vec == grp_levels[2]],
	                              sample_ids[!(grp_vec %in% grp_levels[1:2])])
	            }
	            heatmap_df_melt$Sample <- factor(heatmap_df_melt$Sample, levels = sample_ids)
	            
	            # Tooltip: show region and genes for easier interpretation.
	            region_key <- paste0(top_dmrs$chr, ":", top_dmrs$start, "-", top_dmrs$end)
	            genes_map <- top_dmrs$Genes; names(genes_map) <- region_key
	            heatmap_df_melt$Genes <- genes_map[as.character(heatmap_df_melt$Region)]
	            heatmap_df_melt$text <- paste0(
	              "Sample: ", as.character(heatmap_df_melt$Sample),
	              "<br>Group: ", ifelse(is.na(heatmap_df_melt$Group), "NA", heatmap_df_melt$Group),
	              "<br>Region: ", heatmap_df_melt$Region,
	              "<br>Genes: ", ifelse(is.na(heatmap_df_melt$Genes) | !nzchar(heatmap_df_melt$Genes), "NA", heatmap_df_melt$Genes),
	              "<br>Mean beta: ", ifelse(is.finite(heatmap_df_melt$Methylation), signif(heatmap_df_melt$Methylation, 4), "NA")
	            )
	            
	            p_heat_dmr <- ggplot(heatmap_df_melt, aes(x = Sample, y = Region, fill = Methylation, text = text)) +
	              geom_tile() +
	              scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444", midpoint = 0.5) +
	              theme_minimal() +
	              theme(axis.text.x = element_blank(),
	                    axis.ticks.x = element_blank(),
	                    axis.text.y = element_text(size = 6),
	                    panel.grid = element_blank()) +
	              labs(title = paste(prefix, "Top 50 DMRs Heatmap"), x = "Samples", y = "DMR Region")
	            
	            anno_df <- unique(heatmap_df_melt[, c("Sample", "Group")])
	            okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
	            anno_levels <- grp_levels
	            if (length(anno_levels) == 0) anno_levels <- sort(unique(as.character(anno_df$Group)))
	            anno_df$Group <- factor(anno_df$Group, levels = anno_levels)
	            group_palette <- setNames(okabe_ito[seq_len(min(length(okabe_ito), length(anno_levels)))], anno_levels)
	            if (length(anno_levels) > length(okabe_ito)) {
	              extra <- grDevices::hcl.colors(length(anno_levels) - length(okabe_ito), palette = "Dark 3")
	              group_palette <- c(group_palette, setNames(extra, anno_levels[(length(okabe_ito) + 1):length(anno_levels)]))
	            }
	            p_anno <- ggplot(anno_df, aes(x = Sample, y = 1, fill = Group, text = paste("Group:", Group))) +
	              geom_tile() +
	              scale_fill_manual(values = group_palette, drop = FALSE) +
	              theme_void() +
	              theme(legend.position = "bottom",
	                    axis.text = element_blank(),
	                    axis.title = element_blank(),
	                    plot.margin = margin(0, 0, 0, 0))
	            
	            p_heat_dmr_ly <- ggplotly(p_heat_dmr, tooltip = "text")
	            p_anno_ly <- ggplotly(p_anno, tooltip = "text")
	            final_heatmap <- subplot(p_anno_ly, p_heat_dmr_ly, nrows = 2, heights = c(0.06, 0.94),
	                                     shareX = TRUE, titleX = FALSE) %>%
	              layout(title = paste(prefix, "Top 50 DMRs Heatmap"), hovermode = "closest")
	            
	            saveWidget(final_heatmap, file = file.path(out_dir, paste0(prefix, "_DMR_Heatmap.html")), selfcontained = TRUE)
	            save_static_plot(p_heat_dmr, paste0(prefix, "_DMR_Heatmap.png"), out_dir, width = 8, height = 6)
	          }
	        }
      } else {
        status <<- "ok_empty"
        reason <<- "no_dmrs_after_filter"
      }
    } else {
      status <<- "ok_empty"
      reason <<- "no_dmrs"
    }
  }, error = function(e) {
    status <<- "failed"
    reason <<- e$message
    dmr_res <<- data.frame()
  })

  list(res = dmr_res, status = status, reason = reason)
}

run_method_sensitivity <- function(betas, targets_base, design_base, targets_sva, design_sva,
                                   batch_col, covariates, methods, prefix, out_dir,
                                   pval_thresh, lfc_thresh, max_probes = 0) {
  if (length(methods) == 0) return(NULL)
  betas <- subset_top_variable(betas, max_probes)
  if (!is.null(design_base)) {
    design_base <- drop_linear_dependencies(design_base, group_cols = make.names(levels(targets_base$primary_group)))$mat
  }
  if (!is.null(design_sva)) {
    design_sva <- drop_linear_dependencies(design_sva, group_cols = make.names(levels(targets_sva$primary_group)))$mat
  }
  res_list <- list()
  for (m in methods) {
    if (m == "sva") {
      res <- run_limma_dmp(betas, design_sva, make.names(levels(targets_sva$primary_group)))
      if (!is.null(res)) res_list[[m]] <- res
      next
    }
    design_use <- design_base
    targets_use <- targets_base
    betas_use <- betas
    mvals_use <- NULL
    if (m %in% c("combat", "limma")) {
      bc_res <- apply_batch_correction(betas, m, batch_col, covariates, targets_use, design_use)
      betas_use <- bc_res$betas
      if (!is.null(bc_res$mvals)) mvals_use <- bc_res$mvals
    }
    res <- run_limma_dmp(betas_use, design_use, make.names(levels(targets_use$primary_group)), m_vals = mvals_use)
    if (!is.null(res)) res_list[[m]] <- res
  }
  if (length(res_list) < 2) return(NULL)
  methods_used <- names(res_list)
  pair_metrics <- data.frame()
  for (i in seq_len(length(methods_used) - 1)) {
    for (j in (i + 1):length(methods_used)) {
      m1 <- methods_used[i]
      m2 <- methods_used[j]
      a <- res_list[[m1]][, c("CpG", "logFC", "adj.P.Val")]
      b <- res_list[[m2]][, c("CpG", "logFC", "adj.P.Val")]
      merged <- merge(a, b, by = "CpG", suffixes = c(paste0(".", m1), paste0(".", m2)))
      corr <- suppressWarnings(cor(merged[[paste0("logFC.", m1)]], merged[[paste0("logFC.", m2)]], use = "complete.obs"))
      sig1 <- merged[[paste0("adj.P.Val.", m1)]] < pval_thresh & abs(merged[[paste0("logFC.", m1)]]) > lfc_thresh
      sig2 <- merged[[paste0("adj.P.Val.", m2)]] < pval_thresh & abs(merged[[paste0("logFC.", m2)]]) > lfc_thresh
      overlap <- sum(sig1 & sig2, na.rm = TRUE)
      jaccard <- if (sum(sig1 | sig2, na.rm = TRUE) > 0) overlap / sum(sig1 | sig2, na.rm = TRUE) else NA_real_
      pair_metrics <- rbind(pair_metrics, data.frame(method_a = m1, method_b = m2,
                                                     logFC_correlation = corr, jaccard_overlap = jaccard,
                                                     n_overlap = overlap))
      flip_idx <- (sign(merged[[paste0("logFC.", m1)]]) != sign(merged[[paste0("logFC.", m2)]])) & (sig1 | sig2)
      if (any(flip_idx, na.rm = TRUE)) {
        flip_df <- merged[flip_idx, , drop = FALSE]
        out_name <- paste0(prefix, "_Method_Flips_", m1, "_vs_", m2, ".csv")
        write.csv(flip_df, file.path(out_dir, out_name), row.names = FALSE)
      }
    }
  }
  write.csv(pair_metrics, file.path(out_dir, paste0(prefix, "_Method_Sensitivity.csv")), row.names = FALSE)
  list(res_list = res_list, pair_metrics = pair_metrics)
}

run_batch_pseudo_voi <- function(betas, targets, batch_col, covariates, prefix, out_dir) {
  if (is.null(batch_col) || !(batch_col %in% colnames(targets))) return(NULL)
  if (length(unique(targets[[batch_col]])) < 2) return(NULL)
  covars <- intersect(covariates, colnames(targets))
  formula_str <- paste("~ 0 +", batch_col)
  if (length(covars) > 0) formula_str <- paste(formula_str, "+", paste(covars, collapse = " + "))
  if ("primary_group" %in% colnames(targets)) formula_str <- paste(formula_str, "+ primary_group")
  design <- model.matrix(as.formula(formula_str), data = targets)
  m_vals <- logit_offset(betas)
  fit <- lmFit(m_vals, design)
  batch_cols <- grep(paste0("^", batch_col), colnames(design))
  if (length(batch_cols) == 0) return(NULL)
  contrast_batch <- diag(ncol(design))[, batch_cols, drop = FALSE]
  fit_con <- contrasts.fit(fit, contrast_batch)
  fit_con <- safe_ebayes(fit_con, "batch_pseudo_voi")
  if (is.null(fit_con)) return(NULL)
  p_batch <- fit_con$F.p.value
  fdr_batch <- p.adjust(p_batch, "BH")
  sig_count <- sum(fdr_batch < 0.05, na.rm = TRUE)
  out <- data.frame(metric = c("batch_pseudo_voi_sig_count", "batch_pseudo_voi_min_p"),
                    value = c(sig_count, min(p_batch, na.rm = TRUE)))
  write.csv(out, file.path(out_dir, paste0(prefix, "_Batch_PseudoVOI.csv")), row.names = FALSE)
  out
}

profile_input_columns <- function(df) {
  cols <- colnames(df)
  out <- lapply(cols, function(col) {
    vals <- df[[col]]
    vals_chr <- if (is.character(vals)) vals else as.character(vals)
    vals_chr <- trimws(vals_chr)
    missing <- is.na(vals) | vals_chr == ""
    non_missing <- vals_chr[!missing]
    dtype <- if (is.numeric(vals) || is.integer(vals)) "numeric" else "categorical"
    data.frame(
      column = col,
      dtype_guess = dtype,
      missing_count = sum(missing),
      missing_frac = ifelse(length(vals) > 0, sum(missing) / length(vals), 0),
      unique_non_missing = length(unique(non_missing)),
      example_values = paste(head(non_missing, 3), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

write_input_profile <- function(targets, out_dir) {
  if (is.null(targets) || nrow(targets) == 0) return(invisible(NULL))
  prof <- profile_input_columns(targets)
  write.csv(prof, file.path(out_dir, "Input_Profile.csv"), row.names = FALSE)
  high_missing <- prof$column[prof$missing_frac >= 0.5]
  constant_cols <- prof$column[prof$unique_non_missing <= 1]
  if (length(high_missing) > 0) {
    message(paste("Warning: high missingness columns (>=50%):", paste(head(high_missing, 8), collapse = ", ")))
  }
  if (length(constant_cols) > 0) {
    message(paste("Warning: constant columns will be ignored:", paste(head(constant_cols, 8), collapse = ", ")))
  }
  grp_tbl <- table(targets$primary_group, useNA = "ifany")
  grp_df <- data.frame(group = names(grp_tbl), n = as.integer(grp_tbl), stringsAsFactors = FALSE)
  write.csv(grp_df, file.path(out_dir, "Input_Group_Distribution.csv"), row.names = FALSE)
  invisible(prof)
}

# --- 1. Load and Prepare Data ---

seed_value <- suppressWarnings(as.integer(Sys.getenv("ILLUMETA_SEED", "")))
if (!is.finite(seed_value) || seed_value <= 0) seed_value <- 12345
set.seed(seed_value) # Ensure reproducibility

message("Loading configuration...")
first_line <- tryCatch(readLines(config_file, n = 1, warn = FALSE), error = function(e) "")
if (!length(first_line) || !grepl("\t", first_line)) {
  stop("configure.tsv must be tab-delimited (TSV). CSV/other delimiters are not supported.")
}
targets <- tryCatch(
  read.delim(config_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8"),
  error = function(e) {
    tryCatch(
      read.delim(config_file, stringsAsFactors = FALSE),
      error = function(e2) {
        stop(paste("Cannot read configure.tsv:", conditionMessage(e2)))
      }
    )
  }
)

# Force Batch/ID columns to be factors (categorical)
if ("Sentrix_ID" %in% colnames(targets)) targets$Sentrix_ID <- as.factor(targets$Sentrix_ID)
if ("Sentrix_Position" %in% colnames(targets)) targets$Sentrix_Position <- as.factor(targets$Sentrix_Position)

if (!"primary_group" %in% colnames(targets)) stop("configure.tsv missing 'primary_group' column.")
targets$primary_group <- trimws(targets$primary_group) # remove whitespace
write_input_profile(targets, out_dir)

arg_tissue_norm <- normalize_tissue(opt$tissue)
if (!is.na(arg_tissue_norm)) tissue_use <- arg_tissue_norm
if (tolower(tissue_use) == "auto") {
  inferred <- infer_tissue_from_config(targets)
  if (!is.null(inferred)) {
    tissue_use <- inferred$tissue
    tissue_source <- paste0("config:", inferred$source)
    message(sprintf("Auto tissue detected from %s: %s (raw: %s)", inferred$source, tissue_use, inferred$raw))
  } else {
    tissue_source <- "auto"
    message("Auto tissue selection: no usable tissue column found; using reference-free mode.")
  }
}

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
  # Fallback: match basenames starting with sample_id followed by '_' or end of string.
  # This avoids GSM1234 matching GSM12345 (substring collision).
  pattern <- paste0("^", gsub("([.+?^${}()|\\[\\]\\\\])", "\\\\\\1", sample_id), "(_|$)")
  prefix_match <- basenames[grepl(pattern, basename(basenames))]
  if (length(prefix_match) > 0) return(prefix_match[1])
  return(NA_character_)
}

strip_idat_suffix <- function(path) {
  for (sfx in c("_Grn.idat.gz", "_Red.idat.gz", "_Grn.idat", "_Red.idat")) {
    if (endsWith(path, sfx)) return(substr(path, 1, nchar(path) - nchar(sfx)))
  }
  path
}

resolve_basename <- function(candidate, sample_id, basenames, project_dir, idat_dir) {
  cand <- candidate
  if (!is.na(cand) && nzchar(cand)) cand <- strip_idat_suffix(cand)
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

idat_has <- function(base, suffix) {
  file.exists(paste0(base, suffix)) || file.exists(paste0(base, suffix, ".gz"))
}
idat_pairs <- data.frame(
  SampleID = targets$SampleID,
  Basename = targets$Basename,
  Has_Grn = vapply(targets$Basename, idat_has, logical(1), suffix = "_Grn.idat"),
  Has_Red = vapply(targets$Basename, idat_has, logical(1), suffix = "_Red.idat"),
  Group = as.character(targets$primary_group),
  stringsAsFactors = FALSE
)
write.csv(idat_pairs, file.path(out_dir, "Preflight_IDAT_Pairs.csv"), row.names = FALSE)
missing_pair_idx <- which(!(idat_pairs$Has_Grn & idat_pairs$Has_Red))
if (length(missing_pair_idx) > 0) {
  bad_ids <- idat_pairs$SampleID[missing_pair_idx]
  stop(paste0("Missing IDAT pairs for: ", paste(bad_ids, collapse = ", "),
              ". Ensure both _Grn.idat and _Red.idat (or .gz) exist."))
}

dup_basename_count <- sum(duplicated(targets$Basename))
if (dup_basename_count > 0) {
  dup_ids <- targets$SampleID[duplicated(targets$Basename)]
  stop(paste0("Duplicate Basename entries detected (", dup_basename_count,
              " duplicates). Affected sample IDs: ",
              paste(head(dup_ids, 10), collapse = ", "),
              ". Each sample must map to a unique IDAT pair."))
}

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

n_targets_raw <- nrow(targets)

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

# CRF: sample size tiering
crf_cfg <- config_settings$crf
crf_enabled <- !is.null(crf_cfg) && isTRUE(crf_cfg$enabled)
crf_tier_info <- resolve_crf_tier(n_con, n_test, crf_cfg)
crf_tier <- crf_tier_info$tier
combat_par_prior <- crf_tier_info$combat_par_prior
combat_allowed <- crf_tier_info$combat_allowed
mmc_methods <- crf_tier_info$mmc_methods
if (crf_enabled && !isTRUE(crf_tier_info$sva_allowed)) {
  disable_sva <- TRUE
  message(sprintf("  CRF tier '%s' disables SVA for stability.", crf_tier))
}
log_decision("crf", "sample_tier", crf_tier,
             reason = "sample_size_adaptive",
             metrics = list(total_n = crf_tier_info$total_n,
                            min_per_group = crf_tier_info$min_per_group))
write.csv(data.frame(tier = crf_tier,
                     total_n = crf_tier_info$total_n,
                     min_per_group = crf_tier_info$min_per_group,
                     warnings = paste(crf_tier_info$warnings, collapse = "; ")),
          file.path(out_dir, "CRF_Sample_Tier.csv"), row.names = FALSE)

# Display tier-specific warnings prominently in console
if (length(crf_tier_info$warnings) > 0 && crf_tier %in% c("minimal", "small")) {
  message("\n", strrep("=", 70))
  message(sprintf("  CRF SAMPLE SIZE WARNING (Tier: %s, n=%d)", toupper(crf_tier), crf_tier_info$total_n))
  message(strrep("=", 70))
  for (w in crf_tier_info$warnings) {
    message("  * ", w)
  }
  message(strrep("=", 70), "\n")
}

preflight_summary <- data.frame(
  metric = c("Total_samples_config", "Samples_group_filtered", "Control_group_n", "Test_group_n",
             "Duplicate_basenames", "IDAT_pairs_missing", "IDAT_modal_size"),
  value = c(n_targets_raw, nrow(targets), n_con, n_test,
            dup_basename_count, length(missing_pair_idx),
            ifelse(is.na(mode_size), "NA", format(mode_size, scientific = FALSE))),
  stringsAsFactors = FALSE
)
write.csv(preflight_summary, file.path(out_dir, "Preflight_Summary.csv"), row.names = FALSE)

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
clean_names <- make.unique(c(clean_con, clean_test))
if (!identical(clean_names, c(clean_con, clean_test))) {
  message(sprintf("WARNING: Group labels collided after sanitization; using unique labels: %s vs %s",
                  clean_names[1], clean_names[2]))
  log_decision("group", "sanitized_collision",
               paste(clean_con, clean_test, sep = " vs "),
               reason = paste(clean_names[1], clean_names[2], sep = " vs "))
  clean_con <- clean_names[1]
  clean_test <- clean_names[2]
}

targets$primary_group[tolower(targets$primary_group) == con_lower] <- clean_con
targets$primary_group[tolower(targets$primary_group) == test_lower] <- clean_test

targets$primary_group <- factor(targets$primary_group, levels = c(clean_con, clean_test))


# --- 2. Minfi Analysis ---

message("--- Starting Minfi Analysis ---")
if (force_idat) {
  message("  Forcing IDAT read despite differing array sizes (force=TRUE).")
}
rgSet <- tryCatch(
  read.metharray.exp(targets = targets, force = force_idat),
  error = function(e) {
    stop(paste("Failed to read IDAT files via minfi:",
               conditionMessage(e),
               "\nCheck that IDAT files exist and are not corrupted."))
  }
)
array_type <- detect_array_type(rgSet)
if (!is.na(array_type)) {
  message(sprintf("Detected array type: %s", array_type))
}
array_label <- normalize_array_label(array_type)
if (!is.na(array_type) && array_label == "EPICv2") {
  epicv2_pkgs <- c(
    "IlluminaHumanMethylationEPICv2manifest",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
  )
  missing_epicv2 <- epicv2_pkgs[!vapply(epicv2_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_epicv2) > 0) {
    stop(paste0(
      "EPICv2 array detected but required packages are missing: ",
      paste(missing_epicv2, collapse = ", "),
      ". Install EPICv2 dependencies (e.g. rerun r_scripts/setup_env.R with ILLUMETA_REQUIRE_EPICV2=1) and retry."
    ))
  }
}
if (!is.na(array_type) && array_label == "EPICv2" && tissue_use != "Auto") {
  message(sprintf("EPICv2 detected; reference-based cell composition may not be available for tissue '%s'.", tissue_use))
}

cross_reactive_probes <- NULL
cross_reactive_source <- ""
cross_reactive_version <- ""
cross_reactive_active <- isTRUE(cross_reactive_enabled) && !isTRUE(unsafe_skip_cross_reactive)
cross_reactive_required <- cross_reactive_active && array_label %in% c("450K", "EPIC")
cross_reactive_best_effort <- cross_reactive_active && array_label == "EPICv2"
if (cross_reactive_active) {
  cross_loaded <- resolve_cross_reactive_probes(cross_reactive_path,
                                                use_maxprobes = cross_reactive_use_maxprobes,
                                                array_type = array_type,
                                                local_dir = cross_reactive_local_dir,
                                                allow_missing = !cross_reactive_required)
  if (!is.null(cross_loaded) && !is.null(cross_loaded$probes) && length(cross_loaded$probes) > 0) {
    cross_reactive_probes <- cross_loaded$probes
    cross_reactive_source <- cross_loaded$source
    cross_reactive_version <- ifelse(is.null(cross_loaded$version), "", cross_loaded$version)
    message(sprintf("Cross-reactive probe list loaded (%d probes; source: %s).",
                    length(cross_reactive_probes), cross_reactive_source))
    log_decision("qc", "cross_reactive_list", "loaded",
                 reason = cross_reactive_source,
                 metrics = list(n = length(cross_reactive_probes)))
  } else if (cross_reactive_required) {
    stop("Cross-reactive probe filtering is mandatory for 450K/EPIC arrays, but no list was found. ",
         "Install the maxprobes package or provide a local list via config cross_reactive.path or references/probe_blacklists. ",
         "To bypass (unsafe), rerun with --unsafe-skip-cross-reactive.")
  } else if (cross_reactive_best_effort) {
    message("WARNING: EPICv2 cross-reactive list not found (best-effort). Proceeding without filtering.")
    log_decision("qc", "cross_reactive_list", "missing", reason = "EPICv2_best_effort")
  } else {
    message("Cross-reactive probe filtering enabled but no list found; skipping.")
    log_decision("qc", "cross_reactive_list", "missing", reason = "not_available")
  }
} else {
  message("Cross-reactive probe filtering disabled (unsafe).")
  log_decision("qc", "cross_reactive_list", "skipped", reason = "unsafe_skip")
}

# Persist processed/normalized matrices for all samples (before sample/probe QC)
message("Saving pre-filter processed matrices (all samples, full probe set for GEO)...")
mSet_prefilter <- preprocessNoob(rgSet)
beta_minfi_prefilter <- getBeta(mSet_prefilter)
if (ncol(beta_minfi_prefilter) != nrow(targets)) {
  stop("Failed to align metadata to prefilter Minfi beta matrix: sample counts differ.")
}
# Preserve explicit sample IDs (no underscore trimming)
colnames(beta_minfi_prefilter) <- targets[[gsm_col]]
stopifnot("Minfi prefilter beta column order does not match targets" =
  identical(colnames(beta_minfi_prefilter), targets[[gsm_col]]))
prefilter_beta_path <- file.path(out_dir, "Minfi_BetaMatrix_PreFilter.tsv.gz")
write_matrix_tsv_gz(beta_minfi_prefilter, prefilter_beta_path)
prefilter_mval_path <- file.path(out_dir, "Minfi_MvalueMatrix_PreFilter.tsv.gz")
write_matrix_tsv_gz(logit_offset(beta_minfi_prefilter), prefilter_mval_path)
prefilter_meth <- getMeth(mSet_prefilter)
prefilter_unmeth <- getUnmeth(mSet_prefilter)
prefilter_detp <- detectionP(rgSet)
prefilter_probe_set <- rownames(beta_minfi_prefilter)
if (!all(prefilter_probe_set %in% rownames(prefilter_detp))) {
  missing_detp <- setdiff(prefilter_probe_set, rownames(prefilter_detp))
  warning(sprintf("Prefilter detection P missing %d probes; aligning to available probes.", length(missing_detp)))
}
prefilter_detp <- prefilter_detp[prefilter_probe_set, , drop = FALSE]
if (!all(rownames(prefilter_meth) == prefilter_probe_set)) {
  prefilter_meth <- prefilter_meth[prefilter_probe_set, , drop = FALSE]
}
if (!all(rownames(prefilter_unmeth) == prefilter_probe_set)) {
  prefilter_unmeth <- prefilter_unmeth[prefilter_probe_set, , drop = FALSE]
}
if (ncol(prefilter_meth) != nrow(targets) || ncol(prefilter_unmeth) != nrow(targets)) {
  stop("Failed to align metadata to prefilter Minfi signal matrices: sample counts differ.")
}
if (ncol(prefilter_detp) != nrow(targets)) {
  stop("Failed to align metadata to prefilter detection P matrix: sample counts differ.")
}
colnames(prefilter_meth) <- targets[[gsm_col]]
colnames(prefilter_unmeth) <- targets[[gsm_col]]
colnames(prefilter_detp) <- targets[[gsm_col]]
prefilter_meth_path <- file.path(out_dir, "Minfi_MethylatedSignal_PreFilter.tsv.gz")
prefilter_unmeth_path <- file.path(out_dir, "Minfi_UnmethylatedSignal_PreFilter.tsv.gz")
write_matrix_tsv_gz(prefilter_meth, prefilter_meth_path)
write_matrix_tsv_gz(prefilter_unmeth, prefilter_unmeth_path)
prefilter_detp_path <- file.path(out_dir, "Minfi_DetectionP_PreFilter.tsv.gz")
write_matrix_tsv_gz(prefilter_detp, prefilter_detp_path)
rm(mSet_prefilter, beta_minfi_prefilter, prefilter_meth, prefilter_unmeth, prefilter_detp, prefilter_probe_set)
invisible(gc())

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
        message(sprintf("WARNING: Removing %d samples due to low signal intensity (< %.1f log2):",
                        length(bad_samples), QC_MEDIAN_INTENSITY_THRESHOLD))
        message(paste(bad_samples, collapse = ", "))
        targets <- targets[-bad_samples_idx, ]
        rgSet <- rgSet[, -bad_samples_idx]
        detP <- detP[, -bad_samples_idx]
        if (nrow(targets) < 2) stop("Too few samples remaining after QC.")
    }
} else {
    message("Performing Sample QC: intensity filter disabled (--qc_intensity_threshold<=0); skipping intensity-based removal.")
}
samples_failed_qc <- sum(sample_fail) + length(bad_samples_idx)
samples_failed_sex <- 0
sex_mismatch_count <- 0
sex_check_column_used <- ""
sex_check_res <- NULL
if (isTRUE(sex_check_enabled)) {
  message("Running sex mismatch QC...")
  sex_check_res <- run_sex_mismatch_check(rgSet, targets, gsm_col, out_dir,
                                          action = sex_check_action, column = sex_check_column)
  if (!is.null(sex_check_res)) {
    if (!is.null(sex_check_res$mismatch_count)) sex_mismatch_count <- sex_check_res$mismatch_count
    if (!is.null(sex_check_res$column)) sex_check_column_used <- sex_check_res$column
  }
  log_decision("qc", "sex_mismatch_action", sex_check_action,
               reason = ifelse(sex_check_enabled, "enabled", "skipped"),
               metrics = list(mismatch_count = sex_mismatch_count))
  if (isTRUE(sex_check_action == "stop") && sex_mismatch_count > 0) {
    stop(sprintf("Sex mismatch detected in %d sample(s). Fix metadata or rerun with --sex-mismatch-action drop|ignore.",
                 sex_mismatch_count))
  }
  if (isTRUE(sex_check_action == "drop") &&
      !is.null(sex_check_res) &&
      !is.null(sex_check_res$mismatch_samples) &&
      length(sex_check_res$mismatch_samples) > 0) {
    drop_ids <- intersect(sex_check_res$mismatch_samples, colnames(rgSet))
    if (length(drop_ids) > 0) {
      message(sprintf("  - Dropping %d sample(s) due to sex mismatch.", length(drop_ids)))
      keep <- !(colnames(rgSet) %in% drop_ids)
      rgSet <- rgSet[, keep]
      detP <- detP[, keep, drop = FALSE]
      targets <- targets[keep, , drop = FALSE]
      samples_failed_sex <- length(drop_ids)
      samples_failed_qc <- samples_failed_qc + samples_failed_sex
    }
  }
  if (isTRUE(sex_check_action == "ignore") && sex_mismatch_count > 0) {
    message("WARNING: Sex mismatches were ignored; results may be unreliable.")
  }
}

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
cell_est <- estimate_cell_counts_safe(
  rgSet,
  tissue = tissue_use,
  out_dir = out_dir,
  cell_reference = cell_reference,
  cell_reference_platform = cell_reference_platform,
  array_type = array_type
)
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
    summary_row <- data.frame(
      Pipeline = "Minfi",
      Method = cell_est$reference,
      Tissue = tissue_use,
      K = ifelse(grepl("RefFree", cell_est$reference, ignore.case = TRUE), cell_ref_free_k, NA_integer_),
      Sample_Count = nrow(cell_df),
      Cell_Types = length(cell_cols),
      stringsAsFactors = FALSE
    )
    append_cell_deconv_summary(out_dir, summary_row)
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
        write.csv(assoc_results, file.path(out_dir, "Cell_vs_Group_Association.csv"), row.names = FALSE)
        
        # Check for high confounding (Eta-squared threshold is configurable)
        high_confound <- assoc_results[assoc_results$Eta_Squared > cell_confound_eta2_threshold, ]
        if (nrow(high_confound) > 0) {
            message("WARNING: Strong association detected between Cell Composition and Group. Risk of Over-correction!")
            for (i in 1:nrow(high_confound)) {
                message(sprintf("  - %s: Eta^2 = %.3f (P = %.3g)", high_confound$CellType[i], high_confound$Eta_Squared[i], high_confound$P_Value[i]))
            }
            message("  > Consider removing these covariates if they represent biological pathology rather than noise.")
            if (cell_confound_action == "drop") {
              drop_cells <- unique(high_confound$CellType)
              cell_covariates <- setdiff(cell_covariates, drop_cells)
              config_settings$do_not_adjust <- unique(c(config_settings$do_not_adjust, drop_cells))
              message(sprintf("  > Dropping confounded cell covariates (action=drop): %s", paste(drop_cells, collapse = ", ")))
            }
            if (cell_adjustment_action == "stop") {
              stop("Cell composition strongly associated with group (high Eta^2). Stop requested (cell_adjustment.on_high_eta2=stop).")
            }
        } else {
            message(sprintf("  - Cell composition vs Group association checked: No strong confounding detected (all Eta^2 <= %.2f).",
                            cell_confound_eta2_threshold))
        }
    }
  } else {
    message("  - Cell composition estimated but sample IDs did not match metadata; skipping merge.")
  }
}

targets_global <- targets

# B. Normalization
message("Preprocessing (Noob)...")
mSet <- preprocessNoob(rgSet)
gmSet <- mapToGenome(mSet)

# C. Probe-Level QC
message("Performing Probe QC...")
nrow_raw <- nrow(gmSet)
probe_filter_log <- data.frame(step = "raw", removed = 0, remaining = nrow_raw,
                               stringsAsFactors = FALSE)

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
probe_filter_log <- rbind(probe_filter_log,
                          data.frame(step = "detectionP", removed = probes_failed_detection, remaining = nrow(gmSet)))

n_cross_probes <- 0
cross_reactive_removed_ids <- character(0)
if (isTRUE(cross_reactive_active)) {
  before_cross <- nrow(gmSet)
  cr_res <- apply_cross_reactive_to_gmSet(gmSet, cross_reactive_probes, use_maxprobes = cross_reactive_use_maxprobes)
  gmSet <- cr_res$gmSet
  cross_reactive_removed_ids <- cr_res$removed
  n_cross_probes <- before_cross - nrow(gmSet)
  if (n_cross_probes > 0) {
    if (nzchar(cr_res$source)) cross_reactive_source <- cr_res$source
    source_disp <- ifelse(nzchar(cross_reactive_source), cross_reactive_source, "list")
    message(sprintf("  - Removed %d cross-reactive probes (source: %s).", n_cross_probes, source_disp))
    log_decision("qc", "cross_reactive_removed", as.character(n_cross_probes),
                 reason = source_disp)
  } else {
    message("  - Cross-reactive probe filtering applied but no probes were removed.")
  }
  if (length(cross_reactive_removed_ids) > 0) {
    if (is.null(cross_reactive_probes) || length(cross_reactive_probes) == 0) {
      cross_reactive_probes <- cross_reactive_removed_ids
    } else {
      cross_reactive_probes <- unique(c(cross_reactive_probes, cross_reactive_removed_ids))
    }
  }
} else if (isTRUE(cross_reactive_enabled)) {
  message("  - Cross-reactive probe filtering disabled (unsafe).")
}
probe_filter_log <- rbind(probe_filter_log,
                          data.frame(step = "cross_reactive", removed = n_cross_probes, remaining = nrow(gmSet)))

before_snps <- nrow(gmSet)
gmSet <- dropLociWithSnps(gmSet, snps=c("SBE","CpG"), maf=SNP_MAF_THRESHOLD)
n_snp_probes <- before_snps - nrow(gmSet)
message(paste("  - SNP filtering applied (MAF threshold:", SNP_MAF_THRESHOLD, "). Current probes:", nrow(gmSet)))
probe_filter_log <- rbind(probe_filter_log,
                          data.frame(step = "snp", removed = n_snp_probes, remaining = nrow(gmSet)))

ann <- getAnnotation(gmSet)
keep_sex <- !(ann$chr %in% c("chrX", "chrY"))
before_sex <- nrow(gmSet)
gmSet <- gmSet[keep_sex, ]
n_sex_probes <- before_sex - nrow(gmSet)
message(paste("  - Removed Sex Chromosome probes. Final probes:", nrow(gmSet)))
probe_filter_log <- rbind(probe_filter_log,
                          data.frame(step = "sex_chr", removed = n_sex_probes, remaining = nrow(gmSet)))

qc_report <- data.frame(
  metric = c("Total_samples_input", "Samples_failed_QC", "Samples_passed_QC",
             "Samples_failed_sex_mismatch", "Sex_mismatch_samples",
             "Total_probes_raw", "Probes_failed_detection", "Probes_cross_reactive",
             "Probes_with_SNPs", "Probes_sex_chromosomes", "Probes_final"),
  value = c(n_samples_input, samples_failed_qc, nrow(targets),
            samples_failed_sex, sex_mismatch_count,
            nrow_raw, probes_failed_detection, n_cross_probes,
            n_snp_probes, n_sex_probes, nrow(gmSet))
)
write.csv(qc_report, file.path(out_dir, "QC_Summary.csv"), row.names = FALSE)
write.csv(probe_filter_log, file.path(out_dir, "Probe_Filter_Summary.csv"), row.names = FALSE)
if (length(cross_reactive_removed_ids) == 0) {
  cr_out <- data.frame(
    CpG = character(0),
    Array = character(0),
    Source = character(0),
    Source_Version = character(0),
    Retrieved_On = character(0),
    stringsAsFactors = FALSE
  )
} else {
  cr_out <- data.frame(
    CpG = cross_reactive_removed_ids,
    Array = ifelse(is.na(array_type), "", array_type),
    Source = ifelse(nzchar(cross_reactive_source), cross_reactive_source,
                    ifelse(isTRUE(cross_reactive_active), "unknown", "skipped")),
    Source_Version = ifelse(nzchar(cross_reactive_version), cross_reactive_version, ""),
    Retrieved_On = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
}
write.csv(cr_out, file.path(out_dir, "CrossReactive_Removed_Probes.csv"), row.names = FALSE)

beta_minfi <- getBeta(gmSet)
if (ncol(beta_minfi) != nrow(targets)) {
  stop("Failed to align metadata to Minfi beta matrix: sample counts differ.")
}
# Preserve explicit sample IDs (no underscore trimming)
colnames(beta_minfi) <- targets[[gsm_col]]
stopifnot("Minfi beta column order does not match targets" =
  identical(colnames(beta_minfi), targets[[gsm_col]]))

anno_data <- getAnnotation(gmSet)
anno <- anno_data[, c("chr", "pos", "Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
anno_df <- as.data.frame(anno)
anno_df$CpG <- rownames(anno_df)

colnames(anno_df)[colnames(anno_df) == "UCSC_RefGene_Name"] <- "Gene"
colnames(anno_df)[colnames(anno_df) == "UCSC_RefGene_Group"] <- "Region"
colnames(anno_df)[colnames(anno_df) == "Relation_to_Island"] <- "Island_Context"


# --- 3. Sesame Analysis ---
message("--- Starting Sesame Analysis ---")
beta_sesame_strict <- NULL
beta_sesame_native <- NULL
sesame_strict_before <- 0
sesame_native_before <- 0
sesame_typeinorm_enabled <- isTRUE(opt$sesame_typeinorm)
sesame_typeinorm_disabled <- !sesame_typeinorm_enabled && !isTRUE(opt$skip_sesame)
sesame_dyebias_mode <- if (opt$skip_sesame) "skipped" else if (sesame_typeinorm_enabled) "dyeBiasCorrTypeINorm" else "dyeBiasL"
sesame_dyebias_thread_error <- FALSE
sesame_dyebias_missing_controls <- FALSE
sesame_dyebias_typeinorm_failed <- FALSE
sesame_dyebias_note <- ""
update_sesame_dyebias_mode <- function(new_mode) {
  # Higher = more degraded preprocessing mode (used for transparent reporting).
  rank <- c(skipped = 0, dyeBiasCorrTypeINorm = 1, dyeBiasL = 2, dyeBiasCorr = 2, noob_only = 3)
  if (!(new_mode %in% names(rank))) return(NULL)
  if (!(sesame_dyebias_mode %in% names(rank))) {
    sesame_dyebias_mode <<- new_mode
    return(NULL)
  }
  if (rank[[new_mode]] > rank[[sesame_dyebias_mode]]) {
    sesame_dyebias_mode <<- new_mode
  }
}
if (opt$skip_sesame) {
  message("Sesame analysis skipped (--skip-sesame).")
} else {
  if (isTRUE(sesame_typeinorm_disabled)) {
    message("  Sesame: TypeINorm disabled by default for stability (use --sesame_typeinorm to enable).")
  }
  sesame_thread_env <- c(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  sesame_thread_env_prev <- Sys.getenv(names(sesame_thread_env), unset = NA_character_)
  do.call(Sys.setenv, as.list(sesame_thread_env))
  message("  Sesame: forcing single-thread BLAS/OMP to avoid pthread errors.")
  tryCatch({
    ssets <- lapply(targets$Basename, function(x) readIDATpair(x))
    sesame_preprocess <- function(x) {
      sdf <- x
      # Align with SeSAMe openSesame "QCEPB" (stable under thread constraints):
      # qualityMask -> channel inference -> linear dye bias -> pOOBAH masking -> noob.
      sdf <- tryCatch(qualityMask(sdf, verbose = FALSE), error = function(e) {
        message("  Sesame qualityMask failed; continuing without design mask: ", e$message)
        sdf
      })
      sdf <- tryCatch(inferInfiniumIChannel(sdf), error = function(e) sdf)

      if (isTRUE(sesame_typeinorm_enabled)) {
        sdf <- tryCatch(
          {
            update_sesame_dyebias_mode("dyeBiasCorrTypeINorm")
            dyeBiasCorrTypeINorm(sdf)
          },
          error = function(e) {
            msg <- conditionMessage(e)
            if (grepl("pthread_create", msg, ignore.case = TRUE)) {
              sesame_dyebias_thread_error <<- TRUE
              message("  Sesame dyeBiasCorrTypeINorm disabled after thread error; using dyeBiasL.")
            } else {
              message("  Sesame dyeBiasCorrTypeINorm failed; using dyeBiasL: ", msg)
              sesame_dyebias_typeinorm_failed <<- TRUE
            }
            tryCatch({
              update_sesame_dyebias_mode("dyeBiasL")
              dyeBiasL(sdf)
            }, error = function(e2) {
              message("  Sesame dyeBiasL failed; proceeding without dye bias correction: ", e2$message)
              update_sesame_dyebias_mode("noob_only")
              sdf
            })
          }
        )
      } else {
        sdf <- tryCatch({
          update_sesame_dyebias_mode("dyeBiasL")
          dyeBiasL(sdf)
        }, error = function(e) {
          message("  Sesame dyeBiasL failed; proceeding without dye bias correction: ", e$message)
          update_sesame_dyebias_mode("noob_only")
          sdf
        })
      }

      sdf <- tryCatch(pOOBAH(sdf), error = function(e) {
        message("  Sesame pOOBAH masking failed; continuing without pOOBAH: ", e$message)
        sdf
      })
      sdf <- tryCatch(noob(sdf), error = function(e) {
        message("  Sesame noob failed; continuing with current signals: ", e$message)
        sdf
      })
      sdf
    }
    sdf_list <- lapply(ssets, sesame_preprocess)
    betas_sesame_list <- lapply(sdf_list, getBetas)
    common_probes_sesame <- Reduce(intersect, lapply(betas_sesame_list, names))
    beta_sesame_raw <- do.call(cbind, lapply(betas_sesame_list, function(x) x[common_probes_sesame]))
    colnames(beta_sesame_raw) <- targets[[gsm_col]]
    stopifnot("SeSAMe beta column order does not match targets" =
      identical(colnames(beta_sesame_raw), targets[[gsm_col]]))
    if (nrow(beta_sesame_raw) == 0) {
      stop("Sesame preprocessing produced no probes after alignment.")
    }

    # Persist pre-filter Sesame matrices (before NA filtering and Minfi-alignment)
    prefilter_beta_path <- file.path(out_dir, "Sesame_BetaMatrix_PreFilter.tsv.gz")
    write_matrix_tsv_gz(beta_sesame_raw, prefilter_beta_path)
    prefilter_mval_path <- file.path(out_dir, "Sesame_MvalueMatrix_PreFilter.tsv.gz")
    write_matrix_tsv_gz(logit_offset(beta_sesame_raw), prefilter_mval_path)

    if (isTRUE(cross_reactive_active) && !is.null(cross_reactive_probes) && length(cross_reactive_probes) > 0) {
      cr_res <- apply_cross_reactive_to_matrix(beta_sesame_raw, cross_reactive_probes)
      beta_sesame_raw <- cr_res$mat
      n_cross <- length(cr_res$removed)
      source_disp <- ifelse(nzchar(cross_reactive_source), cross_reactive_source, "list")
      message(sprintf("  - Sesame: removed %d cross-reactive probes (source: %s).", n_cross, source_disp))
    }

    sesame_cell_est <- estimate_sesame_cell_counts(sdf_list, targets[[gsm_col]], out_dir = out_dir)
    if (!is.null(sesame_cell_est)) {
      message("  - Sesame cell composition estimated using ", sesame_cell_est$reference)
    }

    # Native Sesame: keep pOOBAH masking and avoid Minfi-driven probe filtering
    sesame_native_before <- nrow(beta_sesame_raw)
    native_res <- prepare_sesame_native(beta_sesame_raw, SESAME_NATIVE_NA_MAX_FRAC)
    beta_sesame_native <- native_res$mat
    if (!is.null(beta_sesame_native) && nrow(beta_sesame_native) == 0) beta_sesame_native <- NULL
    if (!is.null(native_res$impute_method)) {
      sesame_native_impute_method <<- native_res$impute_method
    }
    if (is.null(sesame_cell_est) && !is.null(beta_sesame_native)) {
      sesame_cell_est <- estimate_cell_counts_epidish(beta_sesame_native, tissue_use, out_dir = out_dir)
      if (!is.null(sesame_cell_est)) {
        message("  - Sesame cell composition estimated using ", sesame_cell_est$reference)
      }
    }
    if (is.null(sesame_cell_est) && !is.null(beta_sesame_native) && tolower(tissue_use) == "placenta") {
      sesame_cell_est <- estimate_cell_counts_planet(beta_sesame_native, out_dir = out_dir)
      if (!is.null(sesame_cell_est)) {
        message("  - Sesame cell composition estimated using ", sesame_cell_est$reference)
      }
    }
    if (is.null(sesame_cell_est) && !is.null(beta_sesame_native)) {
      sesame_cell_est <- estimate_cell_counts_reffree(beta_sesame_native, out_dir = out_dir,
                                                      K_latent = cell_ref_free_k,
                                                      max_probes = cell_ref_free_max_probes)
      if (!is.null(sesame_cell_est)) {
        message("  - Sesame cell composition estimated using ", sesame_cell_est$reference)
      }
    }
    if (!is.null(sesame_cell_est)) {
      sesame_cell_reference <<- sesame_cell_est$reference
      if (!is.null(sesame_cell_est$counts)) {
        cell_cols_sesame <- grep("^Cell_", colnames(sesame_cell_est$counts), value = TRUE)
        summary_row <- data.frame(
          Pipeline = "Sesame",
          Method = sesame_cell_est$reference,
          Tissue = tissue_use,
          K = ifelse(grepl("RefFree", sesame_cell_est$reference, ignore.case = TRUE), cell_ref_free_k, NA_integer_),
          Sample_Count = nrow(sesame_cell_est$counts),
          Cell_Types = length(cell_cols_sesame),
          stringsAsFactors = FALSE
        )
        append_cell_deconv_summary(out_dir, summary_row)
      }
    }
    if (is.null(beta_sesame_native)) {
      message("Sesame native: no probes retained after NA filtering.")
    } else {
      message(sprintf(
        "Sesame native: %d probes -> %d (dropped all-NA: %d, dropped NA>%.2f: %d, imputed: %d; method: %s).",
        sesame_native_before, nrow(beta_sesame_native), native_res$dropped_all_na,
        SESAME_NATIVE_NA_MAX_FRAC, native_res$dropped_na, native_res$imputed,
        native_res$impute_method
      ))
    }

    # Strict Sesame: harmonize to Minfi footprint for conservative cross-validation
    beta_sesame_strict <- na.omit(beta_sesame_raw)
    sesame_strict_before <- nrow(beta_sesame_strict)
    qc_probes <- rownames(beta_minfi)
    shared_qc <- intersect(rownames(beta_sesame_strict), qc_probes)
    if (length(shared_qc) > 0) {
      beta_sesame_strict <- beta_sesame_strict[shared_qc, , drop = FALSE]
    }
    if (nrow(beta_sesame_strict) == 0) beta_sesame_strict <- NULL
    if (is.null(beta_sesame_strict)) {
      message("Sesame strict: no probes retained after QC alignment.")
    } else {
      message(sprintf(
        "Sesame strict: %d probes (after NA omit) -> %d after QC alignment.",
        sesame_strict_before, nrow(beta_sesame_strict)
      ))
    }

    rm(beta_sesame_raw, betas_sesame_list, sdf_list, ssets)
    invisible(gc())
  }, error = function(e) {
    message("Sesame analysis failed; continuing with Minfi only: ", e$message)
    beta_sesame_strict <<- NULL
    beta_sesame_native <<- NULL
  })
  for (nm in names(sesame_thread_env_prev)) {
    val <- sesame_thread_env_prev[[nm]]
    if (is.na(val) || val == "") {
      Sys.unsetenv(nm)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(val, nm)))
    }
  }
}

# --- 4. Intersection ---
if (is.null(beta_sesame_strict)) {
  message("Intersection (strict): skipped (Sesame not available).")
} else {
  common_cpgs <- intersect(rownames(beta_minfi), rownames(beta_sesame_strict))
  message(paste("Intersection (strict):", length(common_cpgs), "CpGs common to Minfi and Sesame."))
}
if (is.null(beta_sesame_native)) {
  message("Intersection (native): skipped (Sesame native not available).")
} else {
  common_cpgs <- intersect(rownames(beta_minfi), rownames(beta_sesame_native))
  message(paste("Intersection (native):", length(common_cpgs), "CpGs common to Minfi and Sesame."))
}

targets_sesame <- NULL
if (!is.null(sesame_cell_est)) {
  cell_df <- sesame_cell_est$counts
  cell_cols <- grep("^Cell_", colnames(cell_df), value = TRUE)
  if (length(cell_cols) > 0) {
    targets_sesame <- targets_global
    drop_cols <- grep("^Cell_", colnames(targets_sesame), value = TRUE)
    if (length(drop_cols) > 0) {
      targets_sesame <- targets_sesame[, setdiff(colnames(targets_sesame), drop_cols), drop = FALSE]
    }
    clean_id <- function(x) sub("\\.idat.*$", "", basename(as.character(x)))
    match_idx <- match(clean_id(targets_sesame[[gsm_col]]), clean_id(cell_df$SampleID))
    if (all(is.na(match_idx)) && "Basename" %in% colnames(targets_sesame)) {
      match_idx <- match(clean_id(targets_sesame$Basename), clean_id(cell_df$SampleID))
    }
    matched <- !is.na(match_idx)
    if (any(matched)) {
      for (cc in cell_cols) {
        targets_sesame[[cc]] <- NA_real_
        targets_sesame[[cc]][matched] <- cell_df[[cc]][match_idx[matched]]
      }
      write.csv(targets_sesame[, c(gsm_col, cell_cols), drop = FALSE],
                file.path(out_dir, "cell_counts_sesame_merged.csv"), row.names = FALSE)
      message(paste("  - Sesame cell covariates applied:", paste(cell_cols, collapse = ", ")))
    } else {
      targets_sesame <- NULL
      message("  - Sesame cell composition estimated but sample IDs did not match metadata; using Minfi covariates.")
    }
  }
}
if (is.null(targets_sesame) && !is.null(beta_sesame_native)) {
  message(sprintf("  - Sesame pipelines using %s cell covariates.", sesame_cell_reference))
}

# --- Pipeline Function ---

#' Main DMP/DMR Analysis Pipeline
#'
#' Core analysis function that executes the complete methylation differential
#' analysis workflow: preprocessing, batch correction, DMP detection, DMR
#' identification, and comprehensive reporting.
#'
#' @param betas Beta value matrix (probes x samples) with normalized values
#' @param prefix Pipeline identifier (e.g., "Minfi_Noob", "Sesame_pOOBAH")
#' @param annotation_df Probe annotation data frame with chr, pos, gene info
#' @param targets_override Optional metadata to use instead of global targets
#'
#' @return List containing:
#'   - dmp_result: Limma DMP results data frame
#'   - dmr_result: DMR regions if computed
#'   - batch_method: Applied batch correction method
#'   - crf_result: CRF assessment outcome
#'   - stats: Summary statistics for dashboard
#'
#' @details
#' Pipeline steps:
#' 1. Probe filtering (low range, zero variance, cross-reactive)
#' 2. PCA computation and visualization
#' 3. Batch effect detection and correction selection
#' 4. DMP analysis with limma (eBayes/robust)
#' 5. Lambda inflation guard check
#' 6. CRF robustness assessment
#' 7. DMR analysis (dmrff)
#' 8. Epigenetic clock estimation (if enabled)
#' 9. HTML/CSV output generation
#'
#' @seealso run_crf_assessment, run_intersection, select_batch_strategy
run_pipeline <- function(betas, prefix, annotation_df, targets_override = NULL) {
  message(paste("Running pipeline for:", prefix))
  # Work on a local copy of targets to avoid altering the shared metadata across pipelines
  targets <- if (is.null(targets_override)) targets_global else targets_override
  
  best_method <- "none"
  batch_evaluated <- FALSE
  batch_candidates_df <- NULL
  best_candidate_score <- NA_real_
  pvca_ci_before <- NA_real_
  pvca_ci_after <- NA_real_
  crf_report <- NULL
  raw_res <- NULL
  raw_delta_beta <- NULL
  lambda_guard_status <- "disabled"
  lambda_guard_action <- "none"
  lambda_guard_lambda <- NA_real_
  lambda_guard_threshold <- NA_real_
  tier3_mode <- FALSE
  tier3_primary <- FALSE
  tier3_primary_lambda <- NA_real_
  tier3_primary_res <- NULL
  tier3_meta_method <- ""
  tier3_meta_i2_median <- NA_real_
  tier3_low_power <- FALSE
  tier3_ineligible <- FALSE
  tier3_eligibility <- NULL
  dmr_status <- "not_run"
  dmr_reason <- ""
  dmr_strategy <- "standard"
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
	  
	  # Add group context to tooltips and order samples by group for readability.
	  grp_map <- setNames(as.character(targets$primary_group), targets[[gsm_col]])
	  dist_df$Group1 <- grp_map[as.character(dist_df$Sample1)]
	  dist_df$Group2 <- grp_map[as.character(dist_df$Sample2)]
	  sample_ids <- colnames(beta_clust)
	  grp_vec <- grp_map[sample_ids]
	  grp_levels <- unique(as.character(targets$primary_group))
	  if (length(grp_levels) == 2) {
	      sample_ids <- c(sample_ids[grp_vec == grp_levels[1]], sample_ids[grp_vec == grp_levels[2]])
	  }
	  dist_df$Sample1 <- factor(as.character(dist_df$Sample1), levels = sample_ids)
	  dist_df$Sample2 <- factor(as.character(dist_df$Sample2), levels = sample_ids)
	  
	  subtitle_dist <- ""
	  if (length(grp_levels) >= 2) {
	      subtitle_dist <- paste0("Samples ordered by group: ", paste(grp_levels, collapse = " -> "))
	  }
	  p_dist <- ggplot(dist_df, aes(x=Sample1, y=Sample2, fill=Distance,
	                                text=paste0("Sample1: ", Sample1,
	                                            "<br>Group1: ", ifelse(is.na(Group1), "NA", Group1),
	                                            "<br>Sample2: ", Sample2,
	                                            "<br>Group2: ", ifelse(is.na(Group2), "NA", Group2),
	                                            "<br>Distance: ", signif(Distance, 4)))) +
	      geom_tile() +
	      scale_fill_gradient(low="red", high="white") + 
	      theme_minimal() +
	      theme(axis.text.x = element_text(angle=90, hjust=1, size=6), axis.text.y = element_text(size=6)) +
	      labs(title = paste(prefix, "Sample Distance Matrix (Top 1000 Var Probes)"),
	           subtitle = subtitle_dist)
	      
	  save_interactive_plot(p_dist, paste0(prefix, "_Sample_Clustering_Distance.html"), out_dir)
	  save_static_plot(p_dist, paste0(prefix, "_Sample_Clustering_Distance.png"), out_dir, width = 7, height = 6)

  pca_res <- prcomp(t(betas), scale. = TRUE)
  pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = targets$primary_group)
  p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) + geom_point(size=3) + theme_minimal() + ggtitle(paste(prefix, "PCA (Raw)"))
  save_interactive_plot(p_pca, paste0(prefix, "_PCA_Before.html"), out_dir)
  save_static_plot(p_pca, paste0(prefix, "_PCA_Before.png"), out_dir, width = 5, height = 4)
  
  # Auto-detect covariates via PC association (top PCs)
  covariates <- character(0)
  forced_covariates <- character(0)
  covar_log_path <- file.path(out_dir, paste0(prefix, "_AutoCovariates.csv"))
  drop_log <- data.frame(Variable=character(0), Reason=character(0))
  covar_log_df <- data.frame(Variable = character(0), MinP = numeric(0), MinPC = integer(0), Type = character(0))

  if (include_clock_covariates) {
    clock_cov <- compute_epigenetic_clocks(betas_full, targets, id_col = gsm_col, prefix = prefix,
                                           out_dir = out_dir, tissue = tissue_use, return_covariates = TRUE)
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
  
  covar_info <- auto_detect_covariates(targets, gsm_col, pca_res$x,
                                       alpha = AUTO_COVARIATE_ALPHA,
                                       max_pcs = MAX_PCS_FOR_COVARIATE_DETECTION)
  covar_log_df <- covar_info$log
  covar_report <- build_covariate_candidate_report(targets, covar_log_df, group_col = "primary_group",
                                                   group_assoc_p_threshold = auto_cov_group_p)
  if (nrow(covar_report) > 0) {
    covar_log_df <- covar_report
  }
  if (nrow(covar_report) > 0) {
    covar_report$Selected <- FALSE
    covar_report$Decision_Reason <- ""
  }
  if (isTRUE(auto_cov_enabled)) {
    covariates <- covar_info$selected
    if (length(covariates) > 0) {
      filt <- filter_covariates(targets, covariates, group_col = "primary_group")
      covariates <- filt$keep
      if (nrow(filt$dropped) > 0) drop_log <- rbind(drop_log, filt$dropped)
    }
    if (auto_cov_exclude_group && nrow(covar_report) > 0) {
      group_flag <- covar_report$Variable[which(covar_report$Group_Assoc_Flag %in% TRUE)]
      drop_group <- intersect(covariates, group_flag)
      if (length(drop_group) > 0) {
        covariates <- setdiff(covariates, drop_group)
        drop_log <- rbind(drop_log, data.frame(Variable = drop_group, Reason = "auto_covariate_group_association"))
        for (cv in drop_group) {
          log_decision("auto_covariates", cv, "excluded", reason = "group_associated")
        }
      }
    }
    if (is.finite(auto_cov_max_cor) && nrow(covar_report) > 0) {
      high_cor <- covar_report$Variable[which(abs(covar_report$Max_Correlation) >= auto_cov_max_cor)]
      drop_cor <- intersect(covariates, high_cor)
      if (length(drop_cor) > 0) {
        covariates <- setdiff(covariates, drop_cor)
        drop_log <- rbind(drop_log, data.frame(Variable = drop_cor, Reason = "auto_covariate_collinearity"))
        for (cv in drop_cor) {
          log_decision("auto_covariates", cv, "excluded", reason = "collinearity")
        }
      }
    }
  } else {
    message("  Auto covariate selection disabled by flag/config (report still generated).")
  }
  if (nrow(covar_report) > 0) {
    covar_report$Selected <- covar_report$Variable %in% covariates
    covar_report$Decision_Reason <- ifelse(covar_report$Selected, "selected_auto", "")
    if (auto_cov_exclude_group) {
      idx <- covar_report$Variable %in% drop_log$Variable[drop_log$Reason == "auto_covariate_group_association"]
      covar_report$Decision_Reason[idx] <- "excluded_group_association"
    }
    if (is.finite(auto_cov_max_cor)) {
      idx <- covar_report$Variable %in% drop_log$Variable[drop_log$Reason == "auto_covariate_collinearity"]
      covar_report$Decision_Reason[idx] <- "excluded_collinearity"
    }
    write.csv(covar_report, covar_log_path, row.names = FALSE)
    for (i in seq_len(nrow(covar_report))) {
      cv <- covar_report$Variable[i]
      if (isTRUE(covar_report$Selected[i])) {
        log_decision("auto_covariates", cv, "selected", reason = "pc_association")
      } else if (nzchar(covar_report$Decision_Reason[i])) {
        log_decision("auto_covariates", cv, "excluded", reason = covar_report$Decision_Reason[i])
      }
    }
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

  gov <- apply_covariate_governance(targets, covariates, forced_covariates, config_settings, drop_log)
  targets <- gov$targets
  covariates <- gov$covariates
  forced_covariates <- gov$forced_covariates
  drop_log <- gov$drop_log
  
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
  if (length(covariates) > 0) {
    post_single <- drop_single_level_covariates(targets, covariates)
    covariates <- post_single$keep
    if (nrow(post_single$dropped) > 0) {
      drop_log <- rbind(drop_log, post_single$dropped)
      message(paste("  Dropped single-level covariates after NA filtering:",
                    paste(unique(post_single$dropped$Variable), collapse = ", ")))
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

  covariate_sets <- generate_covariate_sets(targets, covariates, forced_covariates, covar_log_df)
  
  # PVCA before model fitting to quantify variance explained by group/batch factors
  pvca_factors <- unique(c("primary_group", covariates))
  run_pvca_assessment(betas, targets, pvca_factors, prefix, out_dir, threshold = 0.6, max_probes = 5000, sample_ids = colnames(betas))
  if (!clock_computed) {
    compute_epigenetic_clocks(betas, targets, id_col = gsm_col, prefix = prefix, out_dir = out_dir, tissue = tissue_use)
  }

  # Batch-VOI confounding detection and candidate selection
  batch_candidates <- detect_batch_candidates(targets, "primary_group", config_settings)
  conf_df <- assess_batch_confounding(targets, batch_candidates, "primary_group", config_settings)
  if (nrow(conf_df) > 0) {
    write.csv(conf_df, file.path(out_dir, paste0(prefix, "_Batch_Confounding.csv")), row.names = FALSE)
    for (bc in conf_df$batch) {
      plot_confounding_heatmap(targets, bc, "primary_group", prefix, out_dir)
    }
    plot_confounding_summary(conf_df, prefix, out_dir)
  }
  tier3_batches <- conf_df$batch[conf_df$tier == 3]
  tier3_batch <- if (length(tier3_batches) > 0) tier3_batches[1] else NULL
  tier3_mode <- length(tier3_batches) > 0
  if (length(tier3_batches) > 0) {
    log_decision("confounding", "tier3_batch", paste(tier3_batches, collapse = ","),
                 reason = "non_identifiable_overlap")
  }
  if (tier3_mode && !is.null(tier3_batch)) {
    tier3_eligibility <- check_tier3_eligibility(targets, tier3_batch, "primary_group",
                                                 min_total_n = tier3_min_total_n,
                                                 min_per_group_per_stratum = tier3_min_per_group_per_stratum,
                                                 out_dir = out_dir)
    log_decision("tier3", "eligibility",
                 ifelse(isTRUE(tier3_eligibility$eligible), "pass", "fail"),
                 reason = tier3_eligibility$reason,
                 metrics = list(total_n = tier3_eligibility$total_n,
                                min_stratum_n = tier3_eligibility$min_stratum_n))
    if (!isTRUE(tier3_eligibility$eligible)) {
      tier3_ineligible <- TRUE
      msg <- sprintf("Tier3 confounding detected but eligibility failed (total_n=%s, min_stratum_n=%s).",
                     tier3_eligibility$total_n, tier3_eligibility$min_stratum_n)
      if (tier3_on_fail == "stop") {
        stop(paste0(msg, " Stop for safety. Use tier3.on_fail=skip to proceed without Tier3."))
      } else {
        message("WARNING: ", msg, " Proceeding without Tier3 stratified/meta-analysis.")
        tier3_mode <- FALSE
      }
    }
  }
  conf_df$effect_strength <- ifelse(conf_df$voi_type == "continuous", conf_df$r2, conf_df$cramer_v)
  eligible_df <- conf_df[conf_df$tier < 3, , drop = FALSE]
  if (nrow(eligible_df) > 0) {
    eligible_df <- eligible_df[order(eligible_df$tier, -eligible_df$effect_strength), ]
  }
  batch_col <- NULL
  if (nzchar(batch_col_override)) {
    if (!(batch_col_override %in% colnames(targets))) {
      stop(sprintf("Batch override column '%s' not found in metadata.", batch_col_override))
    }
    batch_col <- batch_col_override
    log_decision("confounding", "batch_candidate_override", batch_col, reason = "manual_override")
  } else if (nrow(eligible_df) > 0) {
    batch_col <- eligible_df$batch[1]
  }
  batch_tier <- if (!is.null(batch_col) && batch_col %in% conf_df$batch) conf_df$tier[match(batch_col, conf_df$batch)] else NA_integer_
  if (!is.null(batch_col)) {
    message(paste("  Selected batch candidate:", batch_col, "(tier", batch_tier, ")"))
    log_decision("confounding", "batch_candidate", batch_col,
                 reason = paste("tier", batch_tier))
  } else {
    message("  No eligible batch candidate found after confounding checks.")
  }

  # PVCA-based confounding index (batch/group) on the pre-correction matrix.
  # CI ~ 0: little batch; CI ~ 1: batch comparable to group; CI > 1: batch dominates (risk).
  if (!is.null(batch_col) && batch_col %in% colnames(targets)) {
    pvca_batch_df <- run_pvca_assessment(
      betas,
      targets,
      unique(c("primary_group", batch_col)),
      paste0(prefix, "_WithBatch"),
      out_dir,
      threshold = 0.6,
      max_probes = 5000,
      sample_ids = colnames(betas)
    )
    pvca_ci_before <- compute_pvca_confounding_index(pvca_batch_df, batch_col)
    if (is.finite(pvca_ci_before)) {
      message(sprintf("  PVCA confounding index (batch/group) before correction: %.3f", pvca_ci_before))
    }
  }

  if (!is.null(batch_col) && nrow(targets) >= 4) {
    message(paste("  Evaluating batch correction methods using batch factor:", batch_col))
    covariates_full <- covariates
    sel <- select_batch_strategy(
      betas = betas,
      targets = targets,
      batch_col = batch_col,
      batch_tier = batch_tier,
      covariate_sets = covariate_sets,
      group_col = "primary_group",
      config_settings = config_settings,
      scoring_preset = scoring_preset,
      perm_n = perm_n,
      vp_top = vp_top,
      prefix = prefix,
      out_dir = out_dir,
      method_pool = if (crf_enabled) mmc_methods else NULL
    )
    best_method <- sel$best_method
    if (nzchar(batch_method_override)) {
      best_method <- batch_method_override
      log_decision("confounding", "batch_method_override", best_method, reason = "manual_override")
      message(paste("  Batch method override applied:", best_method))
    }
    if (!is.null(sel$candidates)) {
      batch_evaluated <- TRUE
      batch_candidates_df <- sel$candidates
      if (nrow(batch_candidates_df) > 0) {
        best_candidate_score <- batch_candidates_df$total_score[1]
      }
      out_batch_path <- file.path(out_dir, paste0(prefix, "_BatchMethodComparison.csv"))
      write.csv(sel$candidates, out_batch_path, row.names = FALSE)
    }
    if (length(sel$best_covariates) > 0) {
      covariates <- sel$best_covariates
      drop_covs <- setdiff(covariates_full, covariates)
      if (length(drop_covs) > 0) {
        drop_log <- rbind(drop_log, data.frame(Variable = drop_covs, Reason = "covariate_set_selection"))
      }
      log_decision("covariates", "selected_set", paste(covariates, collapse = ","), reason = "batch_opt")
    }
    message(paste("  Selected batch method:", best_method))
    log_decision("batch_opt", "method", best_method, reason = "scored_selection")
  } else if (!is.null(batch_col)) {
    message("  Skipping batch method comparison (too few samples for stable evaluation).")
  } else {
    message("  No eligible batch factor found for batch method comparison.")
  }
  if (nzchar(batch_method_override) && is.null(batch_col)) {
    message("  Batch method override requested but no batch column available; ignoring.")
  } else if (nzchar(batch_method_override) && !is.null(batch_col) && best_method != batch_method_override) {
    best_method <- batch_method_override
    log_decision("confounding", "batch_method_override", best_method, reason = "manual_override_post_eval")
    message(paste("  Batch method override applied:", best_method))
  }
  if (crf_enabled && best_method == "combat" && !isTRUE(combat_allowed)) {
    message("  CRF guard: ComBat disabled for this sample size tier; falling back to none.")
    log_decision("crf", "combat_guard", "disabled", reason = crf_tier)
    best_method <- "none"
  }
  
  # Build Design Matrix
  formula_str <- "~ 0 + primary_group"
  if (length(covariates) > 0) formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = targets)
  colnames(design) <- gsub("primary_group", "", colnames(design))
  colnames(design) <- make.unique(make.names(colnames(design)))
  group_cols <- make.names(levels(targets$primary_group))
  design_base <- design
  targets_base <- targets
  
  svobj <- NULL
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
  sv_mat_all <- NULL
  if (length(sv_cols) > 0) {
    sv_mat_all <- as.data.frame(targets[, sv_cols, drop = FALSE])
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
  if (!is.null(sv_mat_all) && ncol(sv_mat_all) > 0 &&
      length(sv_cols) < UNDERCORRECTION_GUARD_MIN_SV && use_sva_in_model && !disable_sva) {
    needed_sv <- UNDERCORRECTION_GUARD_MIN_SV - length(sv_cols)
    max_sv_possible <- max_terms_total - length(covariates)
    if (max_sv_possible < needed_sv) {
      drop_needed <- needed_sv - max_sv_possible
      drop_candidates <- setdiff(rank_covariates(targets, covariates, covar_log_df), forced_covariates)
      if (length(drop_candidates) >= drop_needed) {
        drop_covs <- tail(drop_candidates, drop_needed)
        drop_log <- rbind(drop_log, data.frame(Variable = drop_covs, Reason = "undercorrection_guard_drop_covariates"))
        message(paste("  Undercorrection guard: dropping covariates to allow SVs:", paste(drop_covs, collapse = ", ")))
        covariates <- setdiff(covariates, drop_covs)
        design <- design[, setdiff(colnames(design), drop_covs), drop = FALSE]
        targets <- targets[, setdiff(colnames(targets), drop_covs), drop = FALSE]
      } else {
        message("  Undercorrection guard: unable to free space for SVs; keeping covariates.")
      }
      max_sv_possible <- max_terms_total - length(covariates)
    }
    available_sv <- setdiff(colnames(sv_mat_all), sv_cols)
    add_n <- min(needed_sv, max_sv_possible - length(sv_cols), length(available_sv))
    if (add_n > 0) {
      add_sv <- head(available_sv, add_n)
      design <- cbind(design, sv_mat_all[, add_sv, drop = FALSE])
      targets <- cbind(targets, sv_mat_all[, add_sv, drop = FALSE])
      sv_cols <- c(sv_cols, add_sv)
      message(paste("  Undercorrection guard: added SVs:", paste(add_sv, collapse = ", ")))
    } else if (needed_sv > 0) {
      message("  Undercorrection guard: unable to add SVs within cap.")
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
    
    clean_betas <- tryCatch(
      removeBatchEffect(betas, covariates=cov_mat, design=design_keep),
      error = function(e) {
        message(sprintf("  [!] removeBatchEffect for PCA correction failed: %s", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(clean_betas)) {
      message("  Skipping corrected PCA plot.")
    } else {
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
  }

  betas_pre_correction <- betas
  lambda_raw <- NA_real_
  if (!is.null(design_base)) {
    raw_res <- run_limma_dmp(betas_pre_correction, design_base, make.names(levels(targets_base$primary_group)))
    if (!is.null(raw_res) && nrow(raw_res) > 0) {
      pvals_raw <- raw_res$P.Value
      lambda_raw <- compute_genomic_lambda(pvals_raw)
    }
  }
  if (!is.null(betas_pre_correction) && nrow(betas_pre_correction) > 0) {
    con_mask_raw <- targets$primary_group == clean_con
    test_mask_raw <- targets$primary_group == clean_test
    if (any(con_mask_raw) && any(test_mask_raw)) {
      mean_beta_con_raw <- rowMeans(betas_pre_correction[, con_mask_raw, drop = FALSE], na.rm = TRUE)
      mean_beta_test_raw <- rowMeans(betas_pre_correction[, test_mask_raw, drop = FALSE], na.rm = TRUE)
      raw_delta_beta <- mean_beta_test_raw - mean_beta_con_raw
    }
  }

  # Apply the selected batch correction method to the modeling matrix
  betas_for_model <- betas
  mvals_for_model <- NULL
  applied_batch_method <- "none"
  if (!is.null(batch_col) && (batch_col %in% colnames(targets)) && length(unique(targets[[batch_col]])) > 1) {
    bc_res <- apply_batch_correction(betas, best_method, batch_col, all_covariates, targets, design)
    betas_for_model <- bc_res$betas
    if (!is.null(bc_res$mvals)) mvals_for_model <- bc_res$mvals
    applied_batch_method <- bc_res$method
    if (applied_batch_method != "none") {
      message(paste("  Batch correction applied for modeling using:", applied_batch_method))
    } else {
      message("  Batch correction not applied to modeling matrix (method=none or unsupported).")
    }
  }
  betas <- betas_for_model

  # PVCA-based confounding index (batch/group) after the applied batch correction (modeling matrix).
  if (!is.null(batch_col) && batch_col %in% colnames(targets) && length(unique(targets[[batch_col]])) > 1) {
    pvca_batch_after_df <- run_pvca_assessment(
      betas,
      targets,
      unique(c("primary_group", batch_col)),
      paste0(prefix, "_AfterCorrection_WithBatch"),
      out_dir,
      threshold = 0.6,
      max_probes = 5000,
      sample_ids = colnames(betas)
    )
    pvca_ci_after <- compute_pvca_confounding_index(pvca_batch_after_df, batch_col)
    if (is.finite(pvca_ci_after)) {
      message(sprintf("  PVCA confounding index (batch/group) after correction: %.3f", pvca_ci_after))
    }
  }

  # Precompute the probe subset used for permutation calibration (top-variable CpGs).
  # This is reused to compute an "observed lambda" on the same probe set so obs/null comparisons are meaningful.
  top_calib_idx <- integer(0)
  if (is.finite(vp_top) && vp_top > 0 && !is.null(betas) && nrow(betas) > 0) {
    top_calib_idx <- head(order(apply(betas, 1, var), decreasing = TRUE), min(vp_top, nrow(betas)))
  }

  # Persist processed matrices used for modeling (GEO processed files)
  beta_out_path <- file.path(out_dir, paste0(prefix, "_BetaMatrix.tsv.gz"))
  write_matrix_tsv_gz(betas, beta_out_path)

  # Effect sizes on the beta scale (paper-friendly interpretability)
  con_mask <- targets$primary_group == clean_con
  test_mask <- targets$primary_group == clean_test
  mean_beta_con <- rowMeans(betas[, con_mask, drop = FALSE], na.rm = TRUE)
  mean_beta_test <- rowMeans(betas[, test_mask, drop = FALSE], na.rm = TRUE)
  delta_beta <- mean_beta_test - mean_beta_con
  
  groups <- levels(targets$primary_group)
  if (length(groups) < 2) {
      message("Skipping stats: Not enough groups.")
      branch_summary <- list(
        pipeline = prefix,
        n_samples = nrow(targets),
        n_cpgs = nrow(betas),
        batch_candidate = ifelse(is.null(batch_col), "", batch_col),
        batch_tier = ifelse(is.na(batch_tier), "", batch_tier),
        best_method = best_method
      )
      return(list(res = NULL, n_con = n_con_local, n_test = n_test_local, n_samples = nrow(targets), summary = branch_summary))
  }
  contrast_str <- paste0(groups[2], "-", groups[1])
  message(paste("  - Contrast:", contrast_str))
  
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  
  # B. Stats
  m_vals <- if (!is.null(mvals_for_model)) mvals_for_model else logit_offset(betas)
  mval_out_path <- file.path(out_dir, paste0(prefix, "_MvalueMatrix.tsv.gz"))
  write_matrix_tsv_gz(m_vals, mval_out_path)
  fit <- tryCatch(lmFit(m_vals, design), error = function(e) {
    message(sprintf("  [!] lmFit failed: %s", conditionMessage(e)))
    NULL
  })
  if (is.null(fit)) {
    message("  Skipping stats: lmFit failed; no DMP/DMR results generated for this pipeline.")
    branch_summary <- list(
      pipeline = prefix,
      n_samples = nrow(targets),
      n_cpgs = nrow(betas),
      batch_candidate = ifelse(is.null(batch_col), "", batch_col),
      batch_tier = ifelse(is.na(batch_tier), "", batch_tier),
      best_method = best_method,
      lambda_guard_status = lambda_guard_status,
      lambda_guard_action = lambda_guard_action,
      lambda_guard_threshold = lambda_guard_threshold,
      lambda_guard_lambda = lambda_guard_lambda,
      tier3_mode = tier3_mode,
      tier3_ineligible = tier3_ineligible,
      tier3_low_power = tier3_low_power,
      tier3_batch = ifelse(is.null(tier3_batch), "", tier3_batch),
      primary_result_mode = "analysis_failed",
      tier3_primary_lambda = NA_real_,
      tier3_meta_method = "",
      tier3_meta_i2_median = NA_real_,
      dmr_status = "skipped",
      dmr_reason = "lmfit_failed",
      dmr_strategy = "standard",
      caf_score = NA_real_,
      best_candidate_score = best_candidate_score,
      covariates_used = used_covariates,
      sv_used = sv_cols,
      permutations = perm_n
    )
    return(branch_summary)
  }
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- safe_ebayes(fit2, paste(prefix, "primary_model"))
  if (is.null(fit2)) {
    message("  Skipping stats: eBayes failed; no DMP/DMR results generated for this pipeline.")
    branch_summary <- list(
      pipeline = prefix,
      n_samples = nrow(targets),
      n_cpgs = nrow(betas),
      batch_candidate = ifelse(is.null(batch_col), "", batch_col),
      batch_tier = ifelse(is.na(batch_tier), "", batch_tier),
      best_method = best_method,
      lambda_guard_status = lambda_guard_status,
      lambda_guard_action = lambda_guard_action,
      lambda_guard_threshold = lambda_guard_threshold,
      lambda_guard_lambda = lambda_guard_lambda,
      tier3_mode = tier3_mode,
      tier3_ineligible = tier3_ineligible,
      tier3_low_power = tier3_low_power,
      tier3_batch = ifelse(is.null(tier3_batch), "", tier3_batch),
      primary_result_mode = "analysis_failed",
      tier3_primary_lambda = NA_real_,
      tier3_meta_method = "",
      tier3_meta_i2_median = NA_real_,
      dmr_status = "skipped",
      dmr_reason = "ebayes_failed",
      dmr_strategy = "standard",
      caf_score = NA_real_,
      best_candidate_score = best_candidate_score,
      covariates_used = used_covariates,
      sv_used = sv_cols,
      permutations = perm_n,
      perm_ks_p_median = NA_real_,
      perm_lambda_median = NA_real_
    )
    return(list(res = NULL, n_con = n_con_local, n_test = n_test_local, n_samples = nrow(targets), summary = branch_summary))
  }
  res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  res$CpG <- rownames(res)

  res$Mean_Beta_Con <- mean_beta_con[res$CpG]
  res$Mean_Beta_Test <- mean_beta_test[res$CpG]
  res$Delta_Beta <- delta_beta[res$CpG]

  effect_metrics <- compute_effect_size_consistency(
    raw_res = raw_res,
    corr_res = res,
    raw_delta_beta = raw_delta_beta,
    corr_delta_beta = delta_beta,
    pval_thresh = pval_thresh,
    lfc_thresh = lfc_thresh,
    top_n = 1000
  )
  if (!is.null(effect_metrics)) {
    effect_df <- data.frame(metric = names(effect_metrics),
                            value = vapply(effect_metrics, function(x) ifelse(is.null(x), NA, x), numeric(1)),
                            stringsAsFactors = FALSE)
    write.csv(effect_df, file.path(out_dir, paste0(prefix, "_Signal_Preservation.csv")), row.names = FALSE)
  }
  marker_metrics <- evaluate_marker_retention(res, marker_list, pval_thresh, lfc_thresh, raw_res = raw_res)
  if (!is.null(marker_metrics)) {
    marker_df <- data.frame(metric = names(marker_metrics),
                            value = vapply(marker_metrics, function(x) ifelse(is.null(x), NA, x), numeric(1)),
                            stringsAsFactors = FALSE)
    write.csv(marker_df, file.path(out_dir, paste0(prefix, "_Known_Marker_Summary.csv")), row.names = FALSE)
  }

  sig_delta_pass <- delta_beta_pass(res$Delta_Beta)
  sig_up <- sum(res$adj.P.Val < pval_thresh & res$logFC > lfc_thresh & sig_delta_pass, na.rm = TRUE)
  sig_down <- sum(res$adj.P.Val < pval_thresh & res$logFC < -lfc_thresh & sig_delta_pass, na.rm = TRUE)
  no_signal <- (sig_up + sig_down) == 0
  no_signal_reason <- if (no_signal) "no_significant_dmps_after_thresholds" else ""
  
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
  
  lambda_val <- compute_genomic_lambda(p_vals)
  message(paste("  Genomic Inflation Factor (Lambda):", round(lambda_val, 3)))

  # Observed lambda computed on the same probe set used for permutation calibration (vp_top).
  lambda_vp_top <- NA_real_
  if (length(top_calib_idx) > 0) {
    top_cpgs <- rownames(betas)[top_calib_idx]
    p_top <- tryCatch(res[top_cpgs, "P.Value"], error = function(e) NA_real_)
    p_top <- p_top[!is.na(p_top)]
    if (length(p_top) > 0) {
      lambda_vp_top <- compute_genomic_lambda(p_top)
    }
  }

  lambda_guard_cfg <- config_settings$lambda_guard
  if (!is.null(lambda_guard_cfg) && isTRUE(lambda_guard_cfg$enabled)) {
    lambda_guard_action <- if (!is.null(lambda_guard_cfg$action)) as.character(lambda_guard_cfg$action) else "warn_simplify"
    valid_actions <- c("none", "warn", "warn_simplify", "simplify")
    if (!lambda_guard_action %in% valid_actions) lambda_guard_action <- "warn_simplify"
    lambda_guard_threshold <- suppressWarnings(as.numeric(lambda_guard_cfg$threshold))
    lambda_guard_min_samples <- suppressWarnings(as.integer(lambda_guard_cfg$min_samples))
    if (!is.finite(lambda_guard_min_samples)) lambda_guard_min_samples <- 0L
    if (!is.finite(lambda_val)) {
      lambda_guard_status <- "missing_lambda"
    } else if (!is.finite(lambda_guard_threshold)) {
      lambda_guard_status <- "skipped_no_threshold"
    } else if (nrow(targets) < lambda_guard_min_samples) {
      lambda_guard_status <- "skipped_small_n"
    } else if (lambda_val <= lambda_guard_threshold) {
      lambda_guard_status <- "ok"
    } else {
      lambda_guard_status <- "triggered"
      log_decision("lambda_guard", "triggered", sprintf("%.3f", lambda_val),
                   reason = "lambda_above_threshold",
                   metrics = list(threshold = lambda_guard_threshold,
                                  n_samples = nrow(targets),
                                  action = lambda_guard_action))
      if (lambda_guard_action %in% c("warn", "warn_simplify", "simplify")) {
        message(sprintf("  WARNING: Lambda guard triggered (lambda=%.3f, threshold=%.2f).",
                        lambda_val, lambda_guard_threshold))
      }
      if (lambda_guard_action %in% c("warn_simplify", "simplify")) {
        guard_res <- run_lambda_guard(betas_pre_correction, targets, "primary_group", prefix, out_dir,
                                      max_points, pval_thresh, lfc_thresh)
        if (!is.null(guard_res) && guard_res$status == "ok") {
          lambda_guard_status <- "simplified"
          lambda_guard_lambda <- guard_res$lambda
        } else {
          lambda_guard_status <- "simplify_failed"
        }
      }
    }
  }
  
  expected <- -log10(ppoints(n_p))
  observed <- -log10(pmax(sort(p_vals), .Machine$double.xmin))
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
  delta_pass <- delta_beta_pass(plot_res$Delta_Beta)
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC > lfc_thresh & delta_pass] <- "UP"
  plot_res$diffexpressed[plot_res$adj.P.Val < pval_thresh & plot_res$logFC < -lfc_thresh & delta_pass] <- "DOWN"
  
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
	          grp_levels <- unique(as.character(targets$primary_group))
	          if (length(grp_levels) == 0) {
	              grp_levels <- sort(unique(as.character(anno_df$Group)))
	          }
	          anno_df$Group <- factor(anno_df$Group, levels = grp_levels)
	          okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
	          group_palette <- setNames(okabe_ito[seq_len(min(length(okabe_ito), length(grp_levels)))], grp_levels)
	          if (length(grp_levels) > length(okabe_ito)) {
	              extra <- grDevices::hcl.colors(length(grp_levels) - length(okabe_ito), palette = "Dark 3")
	              group_palette <- c(group_palette, setNames(extra, grp_levels[(length(okabe_ito) + 1):length(grp_levels)]))
	          }
	          p_anno <- ggplot(anno_df, aes(x=Sample, y=1, fill=Group, text=paste("Group:", Group))) +
	              geom_tile() +
	              scale_fill_manual(values = group_palette, drop = FALSE) +
	              theme_void() +
	              theme(legend.position = "bottom",
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
  if (tier3_mode) {
    message("  Tier3 confounding detected: deferring DMR analysis to stratified/meta results.")
    dmr_res <- data.frame()
    dmr_status <- "deferred_tier3"
    dmr_reason <- "tier3_non_identifiable_batch"
    dmr_strategy <- "tier3_meta"
  } else if (nrow(targets) < 4) {
    message("  Skipping DMR analysis (dmrff) due to very small sample size (<4).")
    dmr_res <- data.frame()
    dmr_status <- "skipped_small_n"
    dmr_reason <- "n_samples<4"
  } else {
    message("  Running DMR analysis (dmrff)...")
    dmr_df <- data.frame(
      CpG = rownames(fit2$coefficients),
      logFC = fit2$coefficients[, 1],
      SE = fit2$stdev.unscaled[, 1] * fit2$sigma,
      P.Value = fit2$p.value[, 1],
      stringsAsFactors = FALSE
    )
    dmr_out <- run_dmrff_with_inputs(
      dmr_df = dmr_df,
      betas = betas,
      targets = targets,
      curr_anno = curr_anno,
      prefix = prefix,
      out_dir = out_dir,
      dmr_p_cutoff = dmr_p_cutoff,
      dmr_min_cpgs = DMR_MIN_CPGS,
      con_label = clean_con,
      test_label = clean_test
    )
    dmr_res <- dmr_out$res
    dmr_status <- dmr_out$status
    dmr_reason <- dmr_out$reason
  }

  evaluate_batch <- function(data_betas, meta, label, file_suffix, allowed_vars) {
     allowed <- unique(c("primary_group", allowed_vars))
     allowed <- intersect(allowed, colnames(meta))
     keep_cols <- vapply(allowed, function(nm) {
         if (nm == "primary_group") return(TRUE)
         if (grepl("^SV", nm)) return(length(unique(meta[[nm]])) > 1)
         vals <- normalize_meta_vals(meta[[nm]])
         vals <- vals[!is.na(vals)]
         uniq <- length(unique(vals))
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
         val <- normalize_meta_vals(test_meta[[m]])
         cc <- !is.na(val)
         if (sum(cc) < 3) next
         val_cc <- val[cc]
         if (length(unique(val_cc)) < 2) next
         for (k in 1:n_pcs) {
             if (is.numeric(val_cc)) {
                 fit <- tryCatch(lm(pcs[cc, k] ~ val_cc), error = function(e) NULL)
                 if (!is.null(fit)) {
                   coefs <- summary(fit)$coefficients
                   if (nrow(coefs) >= 2) {
                       pvals[m,k] <- coefs[2,4]
                   } else {
                       pvals[m,k] <- NA
                   }
                 }
             } else {
                 pvals[m,k] <- safe_kruskal_p(pcs[cc, k], val_cc)
             }
         }
     }
     
     csv_name <- paste0(prefix, "_Batch_Evaluation_", file_suffix, "_Table.csv")
     write.csv(pvals, file.path(out_dir, csv_name))
     
     pvals_melt <- as.data.frame(as.table(pvals))
     colnames(pvals_melt) <- c("Variable", "PC", "P_Value")
     pvals_melt$LogP <- -log10(pmax(pvals_melt$P_Value, .Machine$double.xmin))
     
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
  perm_n_probes <- NA_real_
  if (perm_n > 0) {
      message(paste("  Running", perm_n, "permutations for null DMP counts (top", vp_top, "CpGs)..."))
      top_perm <- top_calib_idx
      if (length(top_perm) == 0) {
        message("  Permutation test skipped: vp_top <= 0 or no probes available for calibration.")
      } else {
        perm_n_probes <- length(top_perm)
        m_perm <- logit_offset(betas[top_perm, , drop=FALSE])
      perm_results <- data.frame(run = integer(0), sig_count = integer(0), min_p = numeric(0),
                                 ks_p = numeric(0), lambda = numeric(0))
      perm_group_cols <- make.names(levels(targets$primary_group))
      perm_start <- Sys.time()
      
      for (i in seq_len(perm_n)) {
          perm_labels <- targets$primary_group
          if (!is.null(batch_col) && batch_col %in% colnames(targets)) {
              idx_list <- split(seq_len(nrow(targets)), targets[[batch_col]])
              for (idx in idx_list) {
                  perm_labels[idx] <- sample(perm_labels[idx])
              }
          } else {
              perm_labels <- sample(perm_labels)
          }
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
          # Find group columns by name rather than assuming fixed positions
          grp_idx <- match(perm_group_cols, colnames(perm_design))
          grp_idx <- grp_idx[!is.na(grp_idx)]
          if (length(grp_idx) < 2) next

          cont_vec <- rep(0, ncol(perm_design))
          cont_vec[grp_idx[2]] <- 1; cont_vec[grp_idx[1]] <- -1
          cm_perm <- matrix(cont_vec, ncol = 1)
          rownames(cm_perm) <- colnames(perm_design)
          
          fitp <- lmFit(m_perm, perm_design)
          fitp2 <- contrasts.fit(fitp, cm_perm)
          fitp2 <- safe_ebayes(fitp2, "perm_null")
          if (is.null(fitp2)) next
          resp <- topTable(fitp2, coef = 1, number = Inf, adjust.method = "BH")
          
          delta_pass_perm <- rep(TRUE, nrow(resp))
          if (delta_beta_thresh > 0) {
              perm_con_mask <- perm_labels == clean_con
              perm_test_mask <- perm_labels == clean_test
              if (sum(perm_con_mask) > 0 && sum(perm_test_mask) > 0) {
                  perm_delta <- rowMeans(betas[top_perm, perm_test_mask, drop=FALSE], na.rm = TRUE) -
                      rowMeans(betas[top_perm, perm_con_mask, drop=FALSE], na.rm = TRUE)
                  perm_delta <- perm_delta[rownames(resp)]
                  delta_pass_perm <- delta_beta_pass(perm_delta)
              } else {
                  delta_pass_perm <- rep(FALSE, nrow(resp))
              }
          }
	          sig <- sum(resp$adj.P.Val < pval_thresh & abs(resp$logFC) > lfc_thresh & delta_pass_perm, na.rm=TRUE)
		          pvals <- suppressWarnings(as.numeric(resp$P.Value))
		          pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
		          if (length(pvals) < 2) next
		          minp <- suppressWarnings(min(pvals, na.rm = TRUE))
		          ks_p <- tryCatch(ks.test(pvals, "punif")$p.value, error = function(e) NA_real_)
		          # NOTE: keep the main-model lambda_val intact; use a dedicated variable for permutation lambdas.
		          perm_lambda <- compute_genomic_lambda(pvals)
		          perm_results <- rbind(perm_results, data.frame(run = i, sig_count = sig, min_p = minp,
		                                                         ks_p = ks_p, lambda = perm_lambda))
          
          # Progress every 10 permutations (and at the end)
          if (i %% 10 == 0 || i == perm_n) {
              elapsed <- as.numeric(difftime(Sys.time(), perm_start, units = "mins"))
              eta <- (elapsed / i) * (perm_n - i)
	              message(sprintf("    [perm %d/%d] elapsed: %.1f min, ETA: %.1f min", i, perm_n, elapsed, eta))
	          }
	      }
	      perm_mean_sig <- mean(perm_results$sig_count)
      perm_summary <- data.frame(
        stat = c("mean_sig", "median_sig", "max_sig", "min_p_median", "ks_p_median", "lambda_median",
                 "n_probes", "mean_sig_rate"),
        value = c(perm_mean_sig, median(perm_results$sig_count),
                  max(perm_results$sig_count), median(perm_results$min_p),
                  median(perm_results$ks_p, na.rm = TRUE), median(perm_results$lambda, na.rm = TRUE),
                  perm_n_probes,
                  ifelse(is.finite(perm_mean_sig) && perm_n_probes > 0, perm_mean_sig / perm_n_probes, NA))
      )
	      write.csv(perm_results, file.path(out_dir, paste0(prefix, "_Permutation_Results.csv")), row.names = FALSE)
	      write.csv(perm_summary, file.path(out_dir, paste0(prefix, "_Permutation_Summary.csv")), row.names = FALSE)
      }
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
        vp_fit <- fit_variance_partition_safe(m_vp, form_vp, meta_vp,
                                              vp_cfg = config_settings$variance_partition,
                                              label = "variancePartition")
        vp_res <- vp_fit$result
        
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
  
  if (isTRUE(tier3_mode) && !is.null(tier3_batch)) {
    msg <- paste("Tier3 confounding detected for", tier3_batch, "- running stratified/meta analyses.")
    message(paste("  ", msg))
    log_decision("tier3", "fallback", tier3_batch, reason = "non_identifiable_batch")
    ident_lines <- c(
      paste("Tier3 confounding detected for batch variable:", tier3_batch),
      "Primary VOI effect is not globally identifiable under batch removal.",
      "Primary inference uses stratified EWAS + meta-analysis across overlap strata.",
      "Tier3_Primary outputs are the recommended results for interpretation.",
      "Pooled batch model and overlap-restricted analysis are reported as sensitivity checks."
    )
    writeLines(ident_lines, file.path(out_dir, paste0(prefix, "_Identifiability.txt")))
    meta_out <- run_stratified_meta_analysis(
      betas_pre_correction, targets, tier3_batch, "primary_group",
      used_covariates, forced_covariates, prefix, out_dir,
      pval_thresh, lfc_thresh, delta_beta_thresh
    )
    meta_res <- if (!is.null(meta_out)) meta_out$res else NULL
    if (!is.null(meta_out)) {
      if (!is.null(meta_out$method)) tier3_meta_method <- meta_out$method
      if (!is.null(meta_out$i2_median)) tier3_meta_i2_median <- meta_out$i2_median
    }
    min_total_warn <- suppressWarnings(as.numeric(if (!is.null(config_settings$tier3_meta$min_total_warn))
      config_settings$tier3_meta$min_total_warn else 20))
    if (!is.finite(min_total_warn)) min_total_warn <- 20
    min_stratum_warn <- suppressWarnings(as.numeric(if (!is.null(config_settings$tier3_meta$min_stratum_warn))
      config_settings$tier3_meta$min_stratum_warn else 6))
    if (!is.finite(min_stratum_warn)) min_stratum_warn <- 6
    if (!is.null(meta_out)) {
      if (is.finite(meta_out$total_n) && meta_out$total_n < min_total_warn) tier3_low_power <- TRUE
      if (is.finite(meta_out$min_stratum_n) && meta_out$min_stratum_n < min_stratum_warn) tier3_low_power <- TRUE
      if (tier3_low_power) {
        log_decision("tier3", "low_power", "true",
                     reason = "below_min_n",
                     metrics = list(total_n = meta_out$total_n,
                                    min_stratum_n = meta_out$min_stratum_n,
                                    min_total_warn = min_total_warn,
                                    min_stratum_warn = min_stratum_warn))
      }
    }
    tier3_out <- emit_tier3_primary_outputs(meta_res, betas_pre_correction, targets, curr_anno,
                                            prefix, out_dir, group_con_in, group_test_in,
                                            clean_con, clean_test, max_points, pval_thresh, lfc_thresh,
                                            meta_method = tier3_meta_method, i2_median = tier3_meta_i2_median,
                                            tier3_stats = meta_out)
    if (!is.null(tier3_out) && !is.null(tier3_out$res)) {
      tier3_primary_res <- tier3_out$res
      tier3_primary_lambda <- tier3_out$lambda
      tier3_primary <- TRUE
      res <- tier3_primary_res
      log_decision("tier3", "primary_result", "stratified_meta", reason = "tier3_primary")
    }
    pooled_res <- run_pooled_batch_model(betas_pre_correction, targets, tier3_batch, used_covariates, prefix, out_dir)
    overlap_res <- run_overlap_restricted_analysis(betas_pre_correction, targets, tier3_batch, "primary_group",
                                                   used_covariates, prefix, out_dir)
  } else if (!is.null(tier3_batch) && !isTRUE(tier3_mode)) {
    reason <- if (!is.null(tier3_eligibility$reason)) tier3_eligibility$reason else "ineligible"
    log_decision("tier3", "skipped", "ineligible", reason = reason)
    message("  Tier3 analyses skipped due to ineligible strata.")
  }

  if (tier3_mode) {
    if (!is.null(tier3_primary_res)) {
      tier3_cols_ok <- all(c("CpG", "logFC", "SE", "P.Value") %in% colnames(tier3_primary_res))
      if (!isTRUE(tier3_cols_ok)) {
        dmr_status <- "skipped_tier3_missing_columns"
        dmr_reason <- "missing_columns"
        dmr_strategy <- "tier3_meta"
      } else {
        tier3_prefix <- paste0(prefix, "_Tier3_Primary")
        tier3_df <- tier3_primary_res[, c("CpG", "logFC", "SE", "P.Value")]
        dmr_out <- run_dmrff_with_inputs(
          dmr_df = tier3_df,
          betas = betas_pre_correction,
          targets = targets,
          curr_anno = curr_anno,
          prefix = tier3_prefix,
          out_dir = out_dir,
          dmr_p_cutoff = dmr_p_cutoff,
          dmr_min_cpgs = DMR_MIN_CPGS,
          con_label = clean_con,
          test_label = clean_test
        )
        if (dmr_out$status == "ok") {
          dmr_status <- "ran_tier3_meta"
        } else if (dmr_out$status == "ok_empty") {
          dmr_status <- "ran_tier3_meta_empty"
        } else {
          dmr_status <- dmr_out$status
        }
        dmr_reason <- dmr_out$reason
        dmr_strategy <- "tier3_meta"
      }
    } else {
      dmr_status <- "skipped_tier3_missing_meta"
      dmr_reason <- "tier3_meta_unavailable"
      dmr_strategy <- "tier3_meta"
    }
  }

  mmc_res <- NULL
  if (crf_enabled) {
    mmc_methods_use <- mmc_methods
    if (is.null(batch_col) || !(batch_col %in% colnames(targets)) || length(unique(targets[[batch_col]])) < 2) {
      mmc_methods_use <- setdiff(mmc_methods_use, c("combat", "limma"))
    }
    if (length(mmc_methods_use) >= 2) {
      max_mmc_probes <- suppressWarnings(as.integer(config_settings$crf$max_probes_mmc))
      if (!is.finite(max_mmc_probes)) max_mmc_probes <- 0
      mmc_res <- run_method_sensitivity(betas_pre_correction, targets_base, design_base, targets, design,
                                        batch_col, used_covariates, mmc_methods_use,
                                        prefix, out_dir, pval_thresh, lfc_thresh,
                                        max_probes = max_mmc_probes)
    }
  } else if (!is.null(batch_candidates_df) && nrow(batch_candidates_df) > 0) {
    top_methods <- batch_candidates_df$method[order(-batch_candidates_df$total_score)]
    top_methods <- unique(top_methods)
    top_methods <- head(top_methods, 3)
    run_method_sensitivity(betas_pre_correction, targets_base, design_base, targets, design,
                           batch_col, used_covariates, top_methods,
                           prefix, out_dir, pval_thresh, lfc_thresh)
  }

  if (crf_enabled) {
    crf_report <- run_crf_assessment(
      betas_raw = betas_pre_correction,
      betas_corr = betas,
      targets = targets,
      covariates = used_covariates,
      tier_info = crf_tier_info,
      batch_col = batch_col,
      applied_batch_method = applied_batch_method,
      out_dir = out_dir,
      prefix = prefix,
      pval_thresh = pval_thresh,
      lfc_thresh = lfc_thresh,
      config_settings = config_settings,
      rgSet = rgSet,
      mmc_res = mmc_res,
      pvca_ci_before = pvca_ci_before,
      pvca_ci_after = pvca_ci_after
    )
    if (!is.null(crf_report)) {
      log_decision("crf", "report", "generated", reason = crf_tier)
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
  primary_result_mode <- ifelse(tier3_ineligible,
                                "tier3_ineligible",
                                ifelse(tier3_primary,
                                       ifelse(tier3_low_power, "tier3_low_power", "tier3_stratified_meta"),
                                       "standard_model"))
  tier3_batch_val <- ifelse(is.null(tier3_batch), "", tier3_batch)

  caf_res <- NULL
  caf_report_lines <- NULL
	  if (isTRUE(caf_enabled)) {
      caf_ncs_score <- NA_real_
      if (!is.null(crf_report) && !is.null(crf_report$ncs_best_score)) {
        caf_ncs_score <- crf_report$ncs_best_score
      }
	    caf_res <- compute_caf_scores(
	      perm_summary = perm_summary,
	      perm_n_probes = perm_n_probes,
	      effect_metrics = effect_metrics,
	      marker_metrics = marker_metrics,
	      batch_before = batch_before,
	      batch_after = batch_after,
	      weights = caf_weights,
	      target_fpr = caf_target_fpr,
	      observed_lambda_all = lambda_val,
	      observed_lambda_vp_top = lambda_vp_top,
	      ncs_score = caf_ncs_score
	    )
	    caf_report_lines <- build_caf_report_lines(prefix, applied_batch_method, n_con_local, n_test_local,
	                                               perm_n, perm_n_probes, caf_res, marker_metrics, caf_target_fpr,
	                                               observed_lambda_all = lambda_val,
	                                               observed_lambda_vp_top = lambda_vp_top)
	    if (!is.null(caf_res)) {
	      caf_df <- data.frame(
	        metric = c("calibration_score", "preservation_score", "batch_removal_score", "ncs_score", "cai",
	                   "null_fpr", "null_lambda", "null_ks_p",
	                   "lambda_observed_all", "lambda_observed_vp_top", "lambda_ratio",
	                   "marker_retention", "effect_concordance",
	                   "batch_reduction_ratio",
	                   "weight_calibration", "weight_preservation", "weight_batch", "weight_ncs"),
	        value = c(caf_res$calibration_score, caf_res$preservation_score, caf_res$batch_reduction, caf_res$ncs_score,
	                  caf_res$cai, caf_res$null_fpr, caf_res$null_lambda, caf_res$null_ks_p,
	                  caf_res$observed_lambda_all, caf_res$observed_lambda_vp_top, caf_res$lambda_ratio,
	                  caf_res$marker_retention, caf_res$effect_concordance, caf_res$batch_reduction,
	                  caf_res$weights[["calibration"]], caf_res$weights[["preservation"]], caf_res$weights[["batch"]],
	                  caf_res$weights[["ncs"]]),
	        stringsAsFactors = FALSE
	      )
      write.csv(caf_df, file.path(out_dir, paste0(prefix, "_CAF_Summary.csv")), row.names = FALSE)
    }
    if (!is.null(caf_report_lines)) {
      writeLines(caf_report_lines, file.path(out_dir, paste0(prefix, "_CAF_Report.txt")))
    }
  }
  
	  metrics <- data.frame(
	    metric = c("pipeline", "lambda", "lambda_guard_status", "lambda_guard_action",
	               "lambda_guard_threshold", "lambda_guard_lambda",
	               "primary_result_mode", "no_signal", "no_signal_reason",
	               "tier3_batch", "tier3_primary_lambda", "tier3_low_power",
	               "tier3_ineligible", "tier3_meta_method", "tier3_meta_i2_median",
	               "batch_method_applied", "batch_candidate", "batch_tier", "n_samples", "n_cpgs",
	               "n_covariates_used", "covariates_used", "n_sv_used", "sv_used",
	               "batch_sig_p_lt_0.05_before", "batch_min_p_before",
		               "batch_sig_p_lt_0.05_after", "batch_min_p_after",
		               "group_min_p_before", "group_min_p_after",
		               "pvca_ci_before", "pvca_ci_after",
		               "dropped_covariates",
		               "perm_mean_sig", "perm_max_sig", "perm_ks_p_median", "perm_lambda_median", "vp_primary_group_mean",
		               "lambda_vp_top",
		               "caf_score", "caf_calibration_score", "caf_preservation_score", "caf_batch_score", "caf_ncs_score",
	               "lambda_ratio"),
	    value = c(prefix, lambda_val, lambda_guard_status, lambda_guard_action,
	              lambda_guard_threshold, lambda_guard_lambda,
	              primary_result_mode, no_signal, no_signal_reason,
              tier3_batch_val, tier3_primary_lambda, tier3_low_power,
              tier3_ineligible, tier3_meta_method, tier3_meta_i2_median,
              applied_batch_method,
              ifelse(is.null(batch_col), "", batch_col),
              ifelse(is.na(batch_tier), "", batch_tier),
              nrow(targets), nrow(betas), length(used_covariates),
              paste(used_covariates, collapse=";"),
              length(sv_cols), paste(sv_cols, collapse=";"),
	              batch_before$sig_count, batch_before$min_p,
	              batch_after$sig_count, batch_after$min_p,
	              grp_min_before, grp_min_after,
	              pvca_ci_before, pvca_ci_after,
	              paste(unique(drop_log$Variable), collapse=";"),
	              if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat=="mean_sig"] else NA,
	              if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat=="max_sig"] else NA,
		              if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat=="ks_p_median"] else NA,
	              if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat=="lambda_median"] else NA,
	              vp_primary,
	              lambda_vp_top,
	              if (!is.null(caf_res)) caf_res$cai else NA,
	              if (!is.null(caf_res)) caf_res$calibration_score else NA,
	              if (!is.null(caf_res)) caf_res$preservation_score else NA,
	              if (!is.null(caf_res)) caf_res$batch_reduction else NA,
	              if (!is.null(caf_res)) caf_res$ncs_score else NA,
	              if (!is.null(caf_res)) caf_res$lambda_ratio else NA)
	  )
  if (!is.null(effect_metrics)) {
    eff_df <- data.frame(metric = names(effect_metrics),
                         value = vapply(effect_metrics, function(x) ifelse(is.null(x), NA, x), numeric(1)),
                         stringsAsFactors = FALSE)
    metrics <- rbind(metrics, eff_df)
  }
  if (!is.null(marker_metrics)) {
    mk_df <- data.frame(metric = names(marker_metrics),
                        value = vapply(marker_metrics, function(x) ifelse(is.null(x), NA, x), numeric(1)),
                        stringsAsFactors = FALSE)
    metrics <- rbind(metrics, mk_df)
  }
  metrics_path <- file.path(out_dir, paste0(prefix, "_Metrics.csv"))
  write.csv(metrics, metrics_path, row.names = FALSE)

  ablation_df <- data.frame(
    metric = c("lambda", "batch_sig_p_lt_0.05", "batch_min_p"),
    raw = c(lambda_raw, batch_before$sig_count, batch_before$min_p),
    corrected = c(lambda_val, batch_after$sig_count, batch_after$min_p)
  )
  write.csv(ablation_df, file.path(out_dir, paste0(prefix, "_Ablation_Summary.csv")), row.names = FALSE)

  if (!is.null(batch_col)) {
    run_batch_pseudo_voi(betas, targets, batch_col, used_covariates, prefix, out_dir)
  }
  
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

  branch_summary <- list(
    pipeline = prefix,
    n_samples = nrow(targets),
    n_cpgs = nrow(betas),
    batch_candidate = ifelse(is.null(batch_col), "", batch_col),
    batch_tier = ifelse(is.na(batch_tier), "", batch_tier),
    best_method = best_method,
    lambda_guard_status = lambda_guard_status,
    lambda_guard_action = lambda_guard_action,
    lambda_guard_threshold = lambda_guard_threshold,
    lambda_guard_lambda = lambda_guard_lambda,
    tier3_mode = tier3_mode,
    tier3_ineligible = tier3_ineligible,
    tier3_low_power = tier3_low_power,
    tier3_batch = tier3_batch_val,
    primary_result_mode = primary_result_mode,
    no_signal = no_signal,
    no_signal_reason = no_signal_reason,
    tier3_primary_lambda = tier3_primary_lambda,
    tier3_meta_method = tier3_meta_method,
    tier3_meta_i2_median = tier3_meta_i2_median,
    dmr_status = dmr_status,
    dmr_reason = dmr_reason,
    dmr_strategy = dmr_strategy,
    caf_score = if (!is.null(caf_res)) caf_res$cai else NA_real_,
    best_candidate_score = best_candidate_score,
    covariates_used = used_covariates,
    sv_used = sv_cols,
    permutations = perm_n,
    perm_ks_p_median = if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat == "ks_p_median"] else NA,
    perm_lambda_median = if (!is.null(perm_summary)) perm_summary$value[perm_summary$stat == "lambda_median"] else NA,
    lambda_ratio = if (!is.null(caf_res)) caf_res$lambda_ratio else NA
  )
  return(list(res = res, n_con = n_con_local, n_test = n_test_local, n_samples = nrow(targets), summary = branch_summary))
}

#' Dual-Pipeline Consensus Intersection Analysis
#'
#' Identifies high-confidence differentially methylated positions by requiring
#' significance in both Minfi and Sesame pipelines. Applies Fisher's method for
#' combined p-values and generates concordance visualizations.
#'
#' @param res_minfi Data frame of Minfi DMP results (CpG, logFC, P.Value, adj.P.Val)
#' @param res_sesame Data frame of Sesame DMP results with matching columns
#' @param prefix Output file prefix (e.g., "Intersect_Noob_pOOBAH")
#' @param sesame_label Label for Sesame pipeline variant in plots
#'
#' @return List containing:
#'   - up: Count of consensus hypermethylated CpGs
#'   - down: Count of consensus hypomethylated CpGs
#'
#' @details
#' Consensus criteria:
#' - adj.P.Val < threshold in BOTH pipelines
#' - Concordant logFC direction (both positive or both negative)
#' - Delta beta above minimum threshold in both pipelines
#'
#' Generates: Consensus_DMPs.csv, LogFC_Concordance plot, Overlap_Summary plot
#'
#' @seealso run_pipeline for individual pipeline analysis
run_intersection <- function(res_minfi, res_sesame, prefix, sesame_label) {
  consensus_counts <- list(up = 0, down = 0)
  if (is.null(res_minfi) || is.null(res_sesame)) {
    message(sprintf("Warning: Consensus (%s) skipped because one or more pipelines produced no results.", prefix))
    return(consensus_counts)
  }
  tryCatch({
    keep_cols <- intersect(
      c("CpG", "Gene", "chr", "pos", "Region", "Island_Context", "logFC", "Delta_Beta", "P.Value", "adj.P.Val"),
      colnames(res_minfi)
    )
    a <- res_minfi[, keep_cols, drop = FALSE]
    b <- res_sesame[, intersect(c("CpG", "logFC", "Delta_Beta", "P.Value", "adj.P.Val"), colnames(res_sesame)), drop = FALSE]
    concord <- merge(a, b, by = "CpG", suffixes = c(".Minfi", ".Sesame"))
    concord <- concord[is.finite(concord$logFC.Minfi) & is.finite(concord$logFC.Sesame), , drop = FALSE]
    delta_pass_minfi <- delta_beta_pass(concord$Delta_Beta.Minfi)
    delta_pass_sesame <- delta_beta_pass(concord$Delta_Beta.Sesame)
    p1 <- pmin(pmax(concord$P.Value.Minfi, .Machine$double.xmin), 1)
    p2 <- pmin(pmax(concord$P.Value.Sesame, .Machine$double.xmin), 1)
    concord$P.Fisher <- pchisq(-2 * (log(p1) + log(p2)), df = 4, lower.tail = FALSE)
    concord$adj.P.Fisher <- p.adjust(concord$P.Fisher, "BH")

    is_up <- concord$adj.P.Val.Minfi < pval_thresh &
      concord$adj.P.Val.Sesame < pval_thresh &
      concord$logFC.Minfi > lfc_thresh &
      concord$logFC.Sesame > lfc_thresh &
      delta_pass_minfi &
      delta_pass_sesame

    is_down <- concord$adj.P.Val.Minfi < pval_thresh &
      concord$adj.P.Val.Sesame < pval_thresh &
      concord$logFC.Minfi < -lfc_thresh &
      concord$logFC.Sesame < -lfc_thresh &
      delta_pass_minfi &
      delta_pass_sesame

    consensus_counts$up <- sum(is_up, na.rm = TRUE)
    consensus_counts$down <- sum(is_down, na.rm = TRUE)

    consensus_df <- concord[is_up | is_down, , drop = FALSE]
    consensus_df$logFC_mean <- rowMeans(consensus_df[, c("logFC.Minfi", "logFC.Sesame")], na.rm = TRUE)
    # Use Fisher combined p-value as the consensus ranking statistic.
    # Selection itself is still based on per-pipeline significance + direction
    # (is_up / is_down defined above).
    consensus_df$P.Value <- consensus_df$P.Fisher
    # Use genome-wide BH-corrected Fisher p-values (computed on the full
    # concordance table, not just the consensus subset) for proper FDR control.
    consensus_df$adj.P.Val <- consensus_df$adj.P.Fisher
    # Retain selection-rule p-values explicitly for transparency/back-compat.
    consensus_df$P.Value.max <- pmax(consensus_df$P.Value.Minfi, consensus_df$P.Value.Sesame, na.rm = TRUE)
    consensus_df$adj.P.Val.max <- pmax(consensus_df$adj.P.Val.Minfi, consensus_df$adj.P.Val.Sesame, na.rm = TRUE)
    consensus_df$P.Value.selection <- consensus_df$P.Value.max
    consensus_df$adj.P.Val.selection <- consensus_df$adj.P.Val.max
    consensus_df <- consensus_df[order(consensus_df$P.Value, consensus_df$adj.P.Val), , drop = FALSE]

    out_cons_csv <- file.path(out_dir, paste0(prefix, "_Consensus_DMPs.csv"))
    write.csv(consensus_df, out_cons_csv, row.names = FALSE)
    save_datatable(head(consensus_df, min(10000, nrow(consensus_df))), paste0(prefix, "_Consensus_DMPs.html"), out_dir)
    message(paste("Consensus DMPs saved to", out_cons_csv, "(n=", nrow(consensus_df), ")."))

    # Concordance plot (logFC Minfi vs Sesame)
    concord$Consensus <- (is_up | is_down)
    plot_df <- concord
    if (nrow(plot_df) > max_points) {
      set.seed(seed_value)
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
        title = sprintf("Minfi vs %s logFC concordance", sesame_label),
        subtitle = sprintf("Pearson r = %.3f (points subsampled to max_plots=%d for rendering)", r_val, max_points),
        x = "log2FC (Minfi)",
        y = "log2FC (Sesame)",
        color = "Consensus"
      )
    save_interactive_plot(p_conc, paste0(prefix, "_LogFC_Concordance.html"), out_dir)
    save_static_plot(p_conc, paste0(prefix, "_LogFC_Concordance.png"), out_dir, width = 6, height = 6)

    # Overlap summary plot (counts)
    sig_minfi <- concord$adj.P.Val.Minfi < pval_thresh & abs(concord$logFC.Minfi) > lfc_thresh & delta_pass_minfi
    sig_sesame <- concord$adj.P.Val.Sesame < pval_thresh & abs(concord$logFC.Sesame) > lfc_thresh & delta_pass_sesame
    n_both <- sum(is_up | is_down, na.rm = TRUE)
    overlap_df <- data.frame(
      Category = c("Minfi only", paste0(sesame_label, " only"), "Both (consensus)"),
      Count = c(sum(sig_minfi, na.rm = TRUE) - n_both, sum(sig_sesame, na.rm = TRUE) - n_both, n_both)
    )
    fill_vals <- c("#3498db", "#2ecc71", "#e74c3c")
    names(fill_vals) <- c("Minfi only", paste0(sesame_label, " only"), "Both (consensus)")
    p_ov <- ggplot(overlap_df, aes(x = Category, y = Count, fill = Category, text = Count)) +
      geom_col() +
      scale_fill_manual(values = fill_vals) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = sprintf("Significant DMP overlap (Minfi vs %s)", sesame_label), y = "Count", x = NULL)
    save_interactive_plot(p_ov, paste0(prefix, "_Significant_Overlap.html"), out_dir)
    save_static_plot(p_ov, paste0(prefix, "_Significant_Overlap.png"), out_dir, width = 7, height = 4)

    discordant <- concord[
      (sig_minfi & !sig_sesame) | (!sig_minfi & sig_sesame) |
        (sign(concord$logFC.Minfi) != sign(concord$logFC.Sesame)),
      , drop = FALSE
    ]
    if (nrow(discordant) > 0) {
      write.csv(discordant, file.path(out_dir, paste0(prefix, "_Discordant_Probes.csv")), row.names = FALSE)
      enrich_fields <- c("Region.Minfi", "Island_Context.Minfi")
      for (fld in enrich_fields) {
        if (fld %in% colnames(discordant)) {
          vals <- as.character(discordant[[fld]])
          tbl <- sort(table(vals), decreasing = TRUE)
          enrich_df <- data.frame(Category = names(tbl), Count = as.integer(tbl),
                                  Fraction = as.numeric(tbl) / nrow(discordant))
          out_name <- paste0(prefix, "_Discordant_Enrichment_", gsub("[^A-Za-z0-9]", "_", fld), ".csv")
          write.csv(enrich_df, file.path(out_dir, out_name), row.names = FALSE)
        }
      }
    }
    logfc_corr <- suppressWarnings(cor(concord$logFC.Minfi, concord$logFC.Sesame, use = "complete.obs"))
    jaccard <- if (sum(sig_minfi | sig_sesame) > 0) n_both / sum(sig_minfi | sig_sesame) else NA_real_
    comp_metrics <- data.frame(
      metric = c("logFC_correlation", "jaccard_overlap", "n_both"),
      value = c(logfc_corr, jaccard, n_both)
    )
    write.csv(comp_metrics, file.path(out_dir, paste0(prefix, "_Comparison_Metrics.csv")), row.names = FALSE)

  }, error = function(e) {
    message("Warning: Consensus (Intersection) outputs failed: ", e$message)
  })
  consensus_counts
}

minfi_out <- run_pipeline(beta_minfi, "Minfi", anno_df)
res_minfi <- if (!is.null(minfi_out)) minfi_out$res else NULL
minfi_summary <- if (!is.null(minfi_out)) minfi_out$summary else NULL
if (!is.null(minfi_summary)) {
  write_yaml_safe(minfi_summary, file.path(results_dirs$minfi, "branch_summary.yaml"))
}
rm(minfi_out, beta_minfi)
invisible(gc())

sesame_out <- NULL
res_sesame <- NULL
if (!is.null(beta_sesame_strict)) {
  sesame_out <- run_pipeline(beta_sesame_strict, "Sesame", anno_df, targets_override = targets_sesame)
  res_sesame <- if (!is.null(sesame_out)) sesame_out$res else NULL
  sesame_summary <- if (!is.null(sesame_out)) sesame_out$summary else NULL
  if (!is.null(sesame_summary)) {
    write_yaml_safe(sesame_summary, file.path(results_dirs$sesame, "branch_summary.yaml"))
  }
  rm(sesame_out, beta_sesame_strict)
  invisible(gc())
} else {
  message("Sesame strict pipeline skipped (no Sesame betas available).")
}

sesame_native_out <- NULL
res_sesame_native <- NULL
if (!is.null(beta_sesame_native)) {
  sesame_native_out <- run_pipeline(beta_sesame_native, "Sesame_Native", anno_df, targets_override = targets_sesame)
  res_sesame_native <- if (!is.null(sesame_native_out)) sesame_native_out$res else NULL
  sesame_native_summary <- if (!is.null(sesame_native_out)) sesame_native_out$summary else NULL
  if (!is.null(sesame_native_summary)) {
    write_yaml_safe(sesame_native_summary, file.path(results_dirs$sesame_native, "branch_summary.yaml"))
  }
  rm(sesame_native_out, beta_sesame_native)
  invisible(gc())
} else {
  message("Sesame native pipeline skipped (no Sesame native betas available).")
}

sesame_dyebias_notes <- character(0)
if (isTRUE(sesame_typeinorm_disabled)) sesame_dyebias_notes <- c(sesame_dyebias_notes, "typeinorm_disabled")
if (isTRUE(sesame_dyebias_thread_error)) sesame_dyebias_notes <- c(sesame_dyebias_notes, "pthread_error")
if (isTRUE(sesame_dyebias_typeinorm_failed)) sesame_dyebias_notes <- c(sesame_dyebias_notes, "typeinorm_failed")
if (isTRUE(sesame_dyebias_missing_controls)) sesame_dyebias_notes <- c(sesame_dyebias_notes, "missing_controls")
if (exists("sesame_cache_refresh_attempted", inherits = FALSE) && isTRUE(sesame_cache_refresh_attempted)) {
  sesame_dyebias_notes <- c(sesame_dyebias_notes, "cache_refresh")
}
sesame_dyebias_note <- if (length(sesame_dyebias_notes) > 0) paste(unique(sesame_dyebias_notes), collapse = ";") else ""

# --- Consensus (Intersection) between Minfi and Sesame ---
consensus_counts <- run_intersection(res_minfi, res_sesame, "Intersection", "Sesame (Strict)")
consensus_native_counts <- run_intersection(res_minfi, res_sesame_native, "Intersection_Native", "Sesame (Native)")

get_sig_counts <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(list(up = 0, down = 0))
  delta_pass <- delta_beta_pass(res_df$Delta_Beta)
  up <- sum(res_df$adj.P.Val < pval_thresh & res_df$logFC > lfc_thresh & delta_pass, na.rm=TRUE)
  down <- sum(res_df$adj.P.Val < pval_thresh & res_df$logFC < -lfc_thresh & delta_pass, na.rm=TRUE)
  return(list(up=up, down=down))
}

minfi_counts <- get_sig_counts(res_minfi)
sesame_counts <- get_sig_counts(res_sesame)
sesame_native_counts <- get_sig_counts(res_sesame_native)
intersect_counts <- consensus_counts
intersect_native_counts <- consensus_native_counts
summary_n_con <- n_con
summary_n_test <- n_test

primary_branch <- "Minfi"
primary_summary <- NULL
primary_result_mode <- ""
primary_tier3_batch <- ""
primary_tier3_meta_method <- ""
primary_tier3_meta_i2_median <- NA_real_
primary_tier3_ineligible <- FALSE
primary_tier3_low_power <- FALSE
primary_caf_score <- NA_real_
primary_dmr_status <- ""
primary_dmr_reason <- ""
primary_no_signal <- FALSE
primary_no_signal_reason <- ""
primary_lambda_guard_status <- ""
primary_lambda_guard_action <- ""
primary_lambda_guard_lambda <- NA_real_
primary_lambda_guard_threshold <- NA_real_
primary_branch_override <- ""
primary_branch_reason <- ""

tryCatch({
  branch_df <- data.frame(
    branch = c("Minfi", "Sesame", "Sesame_Native"),
    best_candidate_score = c(
      if (!is.null(minfi_summary)) minfi_summary$best_candidate_score else NA,
      if (!is.null(sesame_summary)) sesame_summary$best_candidate_score else NA,
      if (!is.null(sesame_native_summary)) sesame_native_summary$best_candidate_score else NA
    ),
    stringsAsFactors = FALSE
  )
  branch_df$excluded <- FALSE
  branch_df$excluded_reason <- ""
  sesame_issue <- isTRUE(sesame_dyebias_thread_error) || isTRUE(sesame_dyebias_missing_controls) ||
    sesame_dyebias_mode %in% c("noob_only")
  if (isTRUE(sesame_issue) && !is.null(minfi_summary)) {
    branch_df$best_candidate_score[branch_df$branch %in% c("Sesame", "Sesame_Native")] <- NA
    branch_df$excluded[branch_df$branch %in% c("Sesame", "Sesame_Native")] <- TRUE
    primary_branch_override <- "Minfi"
    sesame_issue_notes <- character(0)
    if (isTRUE(sesame_dyebias_thread_error)) sesame_issue_notes <- c(sesame_issue_notes, "pthread_error")
    if (isTRUE(sesame_dyebias_missing_controls)) sesame_issue_notes <- c(sesame_issue_notes, "missing_controls")
    if (sesame_dyebias_mode %in% c("noob_only")) sesame_issue_notes <- c(sesame_issue_notes, "noob_only")
    primary_branch_reason <- paste0("sesame_fallback:", paste(unique(sesame_issue_notes), collapse = ";"))
    branch_df$excluded_reason[branch_df$branch %in% c("Sesame", "Sesame_Native")] <- primary_branch_reason
    log_decision("consensus", "primary_branch_override", "Minfi", reason = primary_branch_reason)
  }
  missing_mask <- is.na(branch_df$best_candidate_score) & !nzchar(branch_df$excluded_reason)
  branch_df$excluded[missing_mask] <- TRUE
  branch_df$excluded_reason[missing_mask] <- "missing_score"
  branch_df_use <- branch_df[!is.na(branch_df$best_candidate_score), , drop = FALSE]
  primary_branch <- if (nrow(branch_df_use) > 0) branch_df_use$branch[which.max(branch_df_use$best_candidate_score)] else "Minfi"
  branch_df$primary <- branch_df$branch == primary_branch
  write.csv(branch_df, file.path(results_dirs$consensus, "comparison_metrics.csv"), row.names = FALSE)
  log_decision("consensus", "primary_branch", primary_branch, reason = "best_candidate_score")
  writeLines(primary_branch, file.path(results_dirs$consensus, "primary_branch.txt"))

  primary_summary <- switch(primary_branch,
                            Minfi = minfi_summary,
                            Sesame = sesame_summary,
                            Sesame_Native = sesame_native_summary,
                            NULL)
  if (!is.null(primary_summary)) {
    if (!is.null(primary_summary$primary_result_mode)) primary_result_mode <- primary_summary$primary_result_mode
    if (!is.null(primary_summary$tier3_batch)) primary_tier3_batch <- primary_summary$tier3_batch
    if (!is.null(primary_summary$tier3_meta_method)) primary_tier3_meta_method <- primary_summary$tier3_meta_method
    if (!is.null(primary_summary$tier3_meta_i2_median)) primary_tier3_meta_i2_median <- primary_summary$tier3_meta_i2_median
    if (!is.null(primary_summary$tier3_ineligible)) primary_tier3_ineligible <- primary_summary$tier3_ineligible
    if (!is.null(primary_summary$tier3_low_power)) primary_tier3_low_power <- primary_summary$tier3_low_power
    if (!is.null(primary_summary$caf_score)) primary_caf_score <- primary_summary$caf_score
    if (!is.null(primary_summary$dmr_status)) primary_dmr_status <- primary_summary$dmr_status
    if (!is.null(primary_summary$dmr_reason)) primary_dmr_reason <- primary_summary$dmr_reason
    if (!is.null(primary_summary$no_signal)) primary_no_signal <- primary_summary$no_signal
    if (!is.null(primary_summary$no_signal_reason)) primary_no_signal_reason <- primary_summary$no_signal_reason
    if (!is.null(primary_summary$lambda_guard_status)) primary_lambda_guard_status <- primary_summary$lambda_guard_status
    if (!is.null(primary_summary$lambda_guard_action)) primary_lambda_guard_action <- primary_summary$lambda_guard_action
    if (!is.null(primary_summary$lambda_guard_lambda)) primary_lambda_guard_lambda <- primary_summary$lambda_guard_lambda
    if (!is.null(primary_summary$lambda_guard_threshold)) primary_lambda_guard_threshold <- primary_summary$lambda_guard_threshold
  }

  summary_payload <- list(
    n_con = summary_n_con,
    n_test = summary_n_test,
    beginner_safe = beginner_safe,
    minfi_up = minfi_counts$up,
    minfi_down = minfi_counts$down,
    sesame_up = sesame_counts$up,
    sesame_down = sesame_counts$down,
    sesame_native_up = sesame_native_counts$up,
    sesame_native_down = sesame_native_counts$down,
    intersect_up = intersect_counts$up,
    intersect_down = intersect_counts$down,
    intersect_native_up = intersect_native_counts$up,
    intersect_native_down = intersect_native_counts$down,
    primary_branch = primary_branch,
    primary_result_mode = primary_result_mode,
    primary_tier3_batch = primary_tier3_batch,
    primary_tier3_meta_method = primary_tier3_meta_method,
    primary_tier3_meta_i2_median = primary_tier3_meta_i2_median,
    primary_tier3_ineligible = primary_tier3_ineligible,
    primary_tier3_low_power = primary_tier3_low_power,
    primary_no_signal = primary_no_signal,
    primary_no_signal_reason = primary_no_signal_reason,
    primary_caf_score = primary_caf_score,
    primary_dmr_status = primary_dmr_status,
    primary_dmr_reason = primary_dmr_reason,
    primary_lambda_guard_status = primary_lambda_guard_status,
    primary_lambda_guard_action = primary_lambda_guard_action,
    primary_lambda_guard_lambda = primary_lambda_guard_lambda,
    primary_lambda_guard_threshold = primary_lambda_guard_threshold,
    primary_branch_override = primary_branch_override,
    primary_branch_reason = primary_branch_reason,
    sesame_dyebias_mode = sesame_dyebias_mode,
    sesame_dyebias_note = sesame_dyebias_note,
    crf_enabled = crf_enabled,
    crf_sample_tier = crf_tier
  )
  
  write_json(summary_payload, file.path(out_dir, "summary.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  message("Summary statistics saved to summary.json")

  tryCatch({
    ledger_path <- file.path(out_dir, "decision_ledger.tsv")
    write.table(decision_log, ledger_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Decision ledger saved to decision_ledger.tsv")
  }, error = function(e) {
    message("Warning: failed to write decision_ledger.tsv: ", e$message)
  })

  consensus_files <- c(
    "Intersection_Consensus_DMPs.csv",
    "Intersection_Native_Consensus_DMPs.csv",
    "Intersection_Discordant_Probes.csv",
    "Intersection_Native_Discordant_Probes.csv",
    "Intersection_Comparison_Metrics.csv",
    "Intersection_Native_Comparison_Metrics.csv"
  )
  for (cf in consensus_files) {
    src <- file.path(out_dir, cf)
    if (file.exists(src)) {
      file.copy(src, results_dirs$consensus, overwrite = TRUE)
    }
  }

  caf_report_src <- file.path(out_dir, paste0(primary_branch, "_CAF_Report.txt"))
  if (file.exists(caf_report_src)) {
    caf_report_dest <- file.path(out_dir, "Correction_Adequacy_Report.txt")
    file.copy(caf_report_src, caf_report_dest, overwrite = TRUE)
    file.copy(caf_report_dest, results_dirs$reports, overwrite = TRUE)
  }
  caf_summary_src <- file.path(out_dir, paste0(primary_branch, "_CAF_Summary.csv"))
  if (file.exists(caf_summary_src)) {
    caf_summary_dest <- file.path(out_dir, "Correction_Adequacy_Summary.csv")
    file.copy(caf_summary_src, caf_summary_dest, overwrite = TRUE)
    file.copy(caf_summary_dest, results_dirs$reports, overwrite = TRUE)
  }

  crf_report_src <- file.path(out_dir, paste0(primary_branch, "_CRF_Report.txt"))
  if (file.exists(crf_report_src)) {
    crf_report_dest <- file.path(out_dir, "Correction_Robustness_Report.txt")
    file.copy(crf_report_src, crf_report_dest, overwrite = TRUE)
    file.copy(crf_report_dest, results_dirs$reports, overwrite = TRUE)
  }
  for (crf_name in c("CRF_MMC_Summary.csv", "CRF_NCS_Summary.csv", "CRF_RSS_Summary.csv")) {
    src <- file.path(out_dir, paste0(primary_branch, "_", crf_name))
    if (file.exists(src)) {
      dest <- file.path(out_dir, crf_name)
      file.copy(src, dest, overwrite = TRUE)
      file.copy(dest, results_dirs$reports, overwrite = TRUE)
    }
  }
  
  sva_rule <- ifelse(disable_sva, "disabled", "best_method_or_no_batch")
  config_yaml_resolved <- ifelse(file.exists(config_yaml_path), config_yaml_path, "")
  params_payload <- list(
    pval_threshold = pval_thresh,
    lfc_threshold = lfc_thresh,
    delta_beta_threshold = delta_beta_thresh,
    beginner_safe_delta_beta = beginner_safe_delta_beta,
    beginner_safe = beginner_safe,
    snp_maf = SNP_MAF_THRESHOLD,
    qc_intensity_threshold = QC_MEDIAN_INTENSITY_THRESHOLD,
    detection_p_threshold = QC_DETECTION_P_THRESHOLD,
    qc_sample_fail_frac = QC_SAMPLE_DET_FAIL_FRAC,
    cross_reactive_enabled = cross_reactive_enabled,
    unsafe_skip_cross_reactive = unsafe_skip_cross_reactive,
    cross_reactive_active = cross_reactive_active,
    cross_reactive_source = ifelse(nzchar(cross_reactive_source), cross_reactive_source, ""),
    cross_reactive_version = ifelse(nzchar(cross_reactive_version), cross_reactive_version, ""),
    cross_reactive_path = ifelse(nzchar(cross_reactive_path), cross_reactive_path, ""),
    cross_reactive_local_dir = ifelse(nzchar(cross_reactive_local_dir), cross_reactive_local_dir, ""),
    cross_reactive_count = ifelse(is.null(cross_reactive_probes), 0, length(cross_reactive_probes)),
    sex_check_enabled = sex_check_enabled,
    sex_check_action = sex_check_action,
    sex_check_column = ifelse(nzchar(sex_check_column_used), sex_check_column_used, ""),
    sex_mismatch_count = sex_mismatch_count,
    sex_mismatch_dropped = samples_failed_sex,
    auto_covariates_enabled = auto_cov_enabled,
    auto_covariates_exclude_group_associated = auto_cov_exclude_group,
    auto_covariates_group_assoc_p_threshold = auto_cov_group_p,
    auto_covariates_max_cor = auto_cov_max_cor,
    batch_override_column = ifelse(nzchar(batch_col_override), batch_col_override, ""),
    batch_override_method = ifelse(nzchar(batch_method_override), batch_method_override, ""),
    tissue = tissue_use,
    tissue_source = tissue_source,
    array_type = ifelse(is.na(array_type), "", array_type),
    cell_reference = cell_reference,
    cell_reference_platform = cell_reference_platform,
    cell_ref_free_k = cell_ref_free_k,
    cell_ref_free_max_probes = cell_ref_free_max_probes,
    cell_confound_eta2_threshold = cell_confound_eta2_threshold,
    cell_confound_action = cell_confound_action,
    cell_adjustment_on_high_eta2 = cell_adjustment_action,
    auto_covariate_alpha = AUTO_COVARIATE_ALPHA,
    auto_covariate_max_pcs = MAX_PCS_FOR_COVARIATE_DETECTION,
    sesame_native_na_max_frac = SESAME_NATIVE_NA_MAX_FRAC,
    sesame_native_impute_method = sesame_native_impute_method,
    sesame_native_knn_k = SESAME_NATIVE_KNN_K,
    sesame_native_knn_max_rows = SESAME_NATIVE_KNN_MAX_ROWS,
    sesame_native_knn_ref_rows = SESAME_NATIVE_KNN_REF_ROWS,
    sesame_cell_reference = sesame_cell_reference,
    sesame_typeinorm_enabled = sesame_typeinorm_enabled,
    overcorrection_guard_ratio = OVERCORRECTION_GUARD_RATIO,
    undercorrection_guard_min_sv = UNDERCORRECTION_GUARD_MIN_SV,
    sva_enabled = !disable_sva,
    sva_inclusion_rule = sva_rule,
    clock_covariates_enabled = include_clock_covariates,
    sv_group_p_threshold = SV_GROUP_P_THRESHOLD,
    sv_group_eta2_threshold = SV_GROUP_ETA2_THRESHOLD,
    pvca_min_samples = PVCA_MIN_SAMPLES,
    permutations = perm_n,
    crf_enabled = crf_enabled,
    crf_sample_tier = crf_tier,
    crf_mmc_methods = if (!is.null(mmc_methods)) paste(mmc_methods, collapse = ",") else "",
    combat_par_prior = combat_par_prior,
    combat_allowed = combat_allowed,
    crf_ncs_types = if (!is.null(crf_tier_info$ncs_types)) paste(crf_tier_info$ncs_types, collapse = ",") else "",
    crf_rss_mode = if (!is.null(crf_tier_info$rss_mode)) crf_tier_info$rss_mode else "",
    crf_rss_iterations = if (!is.null(crf_tier_info$rss_iterations)) crf_tier_info$rss_iterations else NA,
    crf_rss_top_k = if (!is.null(crf_tier_info$rss_top_k)) paste(crf_tier_info$rss_top_k, collapse = ",") else "",
    crf_max_probes_mmc = if (!is.null(config_settings$crf$max_probes_mmc)) config_settings$crf$max_probes_mmc else NA,
    crf_max_probes_rss = if (!is.null(config_settings$crf$max_probes_rss)) config_settings$crf$max_probes_rss else NA,
    crf_housekeeping_path = if (!is.null(config_settings$crf$housekeeping_path)) config_settings$crf$housekeeping_path else "",
    vp_top = vp_top,
    dmr_maxgap = DMR_MAXGAP,
    dmr_p_cutoff = dmr_p_cutoff,
    dmr_min_cpgs = DMR_MIN_CPGS,
    logit_offset = LOGIT_OFFSET,
    seed = seed_value,
    preset = preset_name,
    config_yaml = config_yaml_resolved,
    min_overlap_per_cell = config_settings$min_overlap_per_cell,
    confounding_tier1_v = config_settings$confounding$tier1_v,
    confounding_tier2_v = config_settings$confounding$tier2_v,
    confounding_tier1_r2 = config_settings$confounding$tier1_r2,
    confounding_tier2_r2 = config_settings$confounding$tier2_r2,
    calibration_ks_p = config_settings$calibration$ks_p,
    caf_enabled = caf_enabled,
    caf_target_fpr = caf_target_fpr,
    caf_weights = caf_weights,
    lambda_guard_enabled = ifelse(!is.null(config_settings$lambda_guard$enabled),
                                  config_settings$lambda_guard$enabled, NA),
    lambda_guard_threshold = ifelse(!is.null(config_settings$lambda_guard$threshold),
                                    config_settings$lambda_guard$threshold, NA),
    lambda_guard_min_samples = ifelse(!is.null(config_settings$lambda_guard$min_samples),
                                      config_settings$lambda_guard$min_samples, NA),
    lambda_guard_action = ifelse(!is.null(config_settings$lambda_guard$action),
                                 config_settings$lambda_guard$action, ""),
    variance_partition_autoscale_numeric = ifelse(!is.null(config_settings$variance_partition$autoscale_numeric),
                                                  config_settings$variance_partition$autoscale_numeric, NA),
    variance_partition_autoscale_on_fail = ifelse(!is.null(config_settings$variance_partition$autoscale_on_fail),
                                                  config_settings$variance_partition$autoscale_on_fail, NA),
    scoring_weights = scoring_preset$weights,
    scoring_guards = scoring_preset$guards,
    scoring_top_k = scoring_preset$top_k,
    tier3_meta_method = ifelse(!is.null(config_settings$tier3_meta$method),
                               config_settings$tier3_meta$method, "auto"),
    tier3_meta_i2_threshold = ifelse(!is.null(config_settings$tier3_meta$i2_threshold),
                                     config_settings$tier3_meta$i2_threshold, NA),
    tier3_meta_min_total_n = ifelse(!is.null(config_settings$tier3_meta$min_total_n),
                                    config_settings$tier3_meta$min_total_n, NA),
    tier3_meta_min_per_group_per_stratum = ifelse(!is.null(config_settings$tier3_meta$min_per_group_per_stratum),
                                                  config_settings$tier3_meta$min_per_group_per_stratum, NA),
    tier3_meta_on_fail = ifelse(!is.null(config_settings$tier3_meta$on_fail),
                                config_settings$tier3_meta$on_fail, ""),
    tier3_meta_min_total_warn = ifelse(!is.null(config_settings$tier3_meta$min_total_warn),
                                       config_settings$tier3_meta$min_total_warn, NA),
    tier3_meta_min_stratum_warn = ifelse(!is.null(config_settings$tier3_meta$min_stratum_warn),
                                         config_settings$tier3_meta$min_stratum_warn, NA)
  )
  write_json(params_payload, file.path(out_dir, "analysis_parameters.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
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

config_used <- list(
  preset = preset_name,
  config_yaml = ifelse(file.exists(config_yaml_path), config_yaml_path, ""),
  beginner_safe = beginner_safe,
  min_overlap_per_cell = config_settings$min_overlap_per_cell,
  logit_offset = LOGIT_OFFSET,
  qc = config_settings$qc,
  dmr = config_settings$dmr,
  confounding = config_settings$confounding,
  batch_candidate_patterns = config_settings$batch_candidate_patterns,
  forced_adjust = config_settings$forced_adjust,
  allowed_adjust = config_settings$allowed_adjust,
  do_not_adjust = config_settings$do_not_adjust,
  batch_override = config_settings$batch_override,
  calibration = config_settings$calibration,
  unsafe = config_settings$unsafe,
  caf = config_settings$caf,
  auto_covariates = config_settings$auto_covariates,
  lambda_guard = config_settings$lambda_guard,
  variance_partition = config_settings$variance_partition,
  tier3_meta = config_settings$tier3_meta,
  cell_adjustment = config_settings$cell_adjustment,
  cell_deconv = config_settings$cell_deconv,
  cross_reactive = config_settings$cross_reactive,
  sex_check = config_settings$sex_check,
  scoring_preset = scoring_preset,
  permutations = perm_n
)
write_yaml_safe(config_used, file.path(out_dir, "config_used.yaml"))
if (nrow(decision_log) > 0) {
  write.table(decision_log, file.path(out_dir, "decision_ledger.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
}
file.copy(file.path(out_dir, "analysis_parameters.json"), results_dirs$logs, overwrite = TRUE)
file.copy(file.path(out_dir, "sessionInfo.txt"), results_dirs$logs, overwrite = TRUE)
file.copy(file.path(out_dir, "config_used.yaml"), results_dirs$logs, overwrite = TRUE)
if (file.exists(file.path(out_dir, "decision_ledger.tsv"))) {
  file.copy(file.path(out_dir, "decision_ledger.tsv"), results_dirs$logs, overwrite = TRUE)
}

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
  array_type_str <- ifelse(is.na(array_type) || !nzchar(array_type), "NA", array_type)
  array_label_str <- ifelse(nzchar(array_label), array_label, array_type_str)
  genome_build <- if (array_label_str == "EPICv2") {
    "hg38"
  } else if (array_label_str %in% c("EPIC", "450K")) {
    "hg19"
  } else {
    "minfi annotation"
  }
  cell_ref_str <- if (nzchar(cell_reference)) {
    if (nzchar(cell_reference_platform)) {
      sprintf("custom (%s; platform=%s)", cell_reference, cell_reference_platform)
    } else {
      sprintf("custom (%s; platform=auto)", cell_reference)
    }
  } else {
    "default references (fallback to RefFreeEWAS when unavailable)"
  }
  sig_line <- sprintf("- Significance thresholds: BH FDR < %.3f and |log2FC| > %.2f", pval_thresh, lfc_thresh)
  if (delta_beta_thresh > 0) {
    sig_line <- paste0(sig_line, sprintf(" and |DeltaBeta| >= %.3f", delta_beta_thresh))
  }
  primary_result_mode_disp <- ifelse(nzchar(primary_result_mode), primary_result_mode, "unknown")
  primary_tier3_batch_disp <- ifelse(nzchar(primary_tier3_batch), primary_tier3_batch, "none")
  primary_tier3_meta_method_disp <- ifelse(nzchar(primary_tier3_meta_method), primary_tier3_meta_method, "fixed")
  primary_tier3_meta_i2_disp <- ifelse(is.finite(primary_tier3_meta_i2_median),
                                       sprintf("%.3f", primary_tier3_meta_i2_median), "NA")
  tier3_low_power_line <- if (primary_result_mode_disp == "tier3_low_power") {
    "WARNING: Tier3 stratified meta-analysis flagged low power; interpret results with extreme caution."
  } else {
    NULL
  }
  tier3_ineligible_line <- if (isTRUE(primary_tier3_ineligible) || primary_result_mode_disp == "tier3_ineligible") {
    "WARNING: Tier3 confounding detected but eligibility failed; stratified/meta-analysis was not run."
  } else {
    NULL
  }
  no_signal_line <- if (isTRUE(primary_no_signal)) {
    "WARNING: No significant DMPs detected at the configured thresholds; results likely underpowered or overly stringent."
  } else {
    NULL
  }
  primary_dmr_status_disp <- ifelse(nzchar(primary_dmr_status), primary_dmr_status, "unknown")
  primary_dmr_reason_disp <- ifelse(nzchar(primary_dmr_reason), paste0(" (reason: ", primary_dmr_reason, ")"), "")
  lambda_guard_status_disp <- ifelse(nzchar(primary_lambda_guard_status), primary_lambda_guard_status, "unknown")
  lambda_guard_action_disp <- ifelse(nzchar(primary_lambda_guard_action), primary_lambda_guard_action, "unknown")
  lambda_guard_threshold_disp <- ifelse(is.finite(primary_lambda_guard_threshold), sprintf("%.3f", primary_lambda_guard_threshold), "NA")
  lambda_guard_lambda_disp <- ifelse(is.finite(primary_lambda_guard_lambda), sprintf("%.3f", primary_lambda_guard_lambda), "NA")
  primary_branch_override_disp <- ifelse(nzchar(primary_branch_override), primary_branch_override, "none")
  primary_branch_reason_disp <- ifelse(nzchar(primary_branch_reason), paste0(" (", primary_branch_reason, ")"), "")
  sesame_dyebias_note_disp <- ifelse(nzchar(sesame_dyebias_note), paste0(" (", sesame_dyebias_note, ")"), "")
  cross_reactive_line <- if (isTRUE(cross_reactive_active) && !is.null(cross_reactive_probes) &&
                             length(cross_reactive_probes) > 0) {
    source_disp <- ifelse(nzchar(cross_reactive_source), cross_reactive_source, "list")
    sprintf("- Cross-reactive probes: removed using %s (n=%d). References: Pidsley 2016; McCartney 2016; Chen 2013.", source_disp, n_cross_probes)
  } else if (isTRUE(cross_reactive_active)) {
    "- Cross-reactive probes: filtering enabled but list unavailable (EPICv2 best-effort or missing local list)."
  } else {
    "- WARNING: Cross-reactive probe filtering was skipped (unsafe)."
  }
  sex_check_line <- if (isTRUE(sex_check_enabled)) {
    col_disp <- ifelse(nzchar(sex_check_column_used), sex_check_column_used, "auto")
    sprintf("- Sex mismatch QC: predicted sex via minfi::getSex() vs metadata column '%s' (action=%s; mismatches=%d).",
            col_disp, sex_check_action, sex_mismatch_count)
  } else {
    "- Sex mismatch QC: disabled or not available."
  }
  beginner_line <- if (beginner_safe) {
    sprintf("- Beginner-safe mode: enabled (min total n >= %d; |DeltaBeta| >= %.3f).", MIN_TOTAL_SIZE_STOP, delta_beta_thresh)
  } else {
    NULL
  }
  methods_lines <- c(
    "# IlluMeta analysis methods (auto-generated)",
    "",
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("- Config: `%s`", config_file),
    sprintf("- Output: `%s`", out_dir),
    sprintf("- Groups: control=`%s` (n=%s), test=`%s` (n=%s)", group_con_in, summary_n_con, group_test_in, summary_n_test),
    beginner_line,
    sprintf("- Tissue: `%s` (source: %s)", tissue_use, tissue_source),
    sprintf("- Array type detected: %s", array_label_str),
    sprintf("- Genomic coordinates: standardized to minfi annotation (%s) for cross-pipeline comparison.", genome_build),
    sprintf("- Cell reference: %s", cell_ref_str),
    sprintf("- Sesame cell composition: %s.", sesame_cell_reference),
    sprintf("- Cell adjustment guard: action=%s when Eta^2 > %.2f.", cell_adjustment_action, cell_confound_eta2_threshold),
    sig_line,
    "",
    "## Overview",
    "IlluMeta performs an end-to-end DNA methylation analysis from raw Illumina IDAT files, running minfi and sesame independently and reporting Sesame in both strict (Minfi-aligned) and native (pOOBAH-preserving) views, plus consensus (intersection) call sets.",
    "",
    "## Data input",
    "- Raw IDATs are read from the project `idat/` directory.",
    sprintf("- Samples are defined by `%s` and `primary_group` in `configure.tsv`.", gsm_col),
    "- Auto-grouping (if used) is heuristic; users must verify group labels/counts before interpretation.",
    "",
    "## Sample-level QC",
    sprintf("- Detection P-value: samples with > %.0f%% probes failing (P > %.3g) are excluded.", QC_SAMPLE_DET_FAIL_FRAC * 100, QC_DETECTION_P_THRESHOLD),
    sprintf("- Signal intensity QC: %s.", qc_intensity_str),
    "- Mixed-array safeguard: samples whose IDAT array size deviates from the modal array size are excluded (unless `--force_idat`).",
    sex_check_line,
    "",
    "## Probe-level QC (minfi-derived probe set)",
    sprintf("- Detection P-value filter: probes are retained only if P < %.3g in all retained samples.", QC_DETECTION_P_THRESHOLD),
    cross_reactive_line,
    sprintf("- SNP filtering: probes overlapping common SNPs (SBE/CpG; MAF >= %.2f) are removed using `minfi::dropLociWithSnps()`.", SNP_MAF_THRESHOLD),
    "- Sex chromosome probes (chrX/chrY) are removed.",
    "",
    "## Normalization (two pipelines)",
    "- **minfi**: `preprocessNoob()` followed by `mapToGenome()`; beta values are extracted with `getBeta()`.",
    "- **sesame**: `qualityMask()` + `inferInfiniumIChannel()` + dye bias correction (default: `dyeBiasL()`; optional `dyeBiasCorrTypeINorm()` when `--sesame_typeinorm` is enabled with fallback to `dyeBiasL()`), followed by `pOOBAH()` masking and `noob()` background correction; beta values are extracted with `getBetas()`.",
    "- **sesame strict view**: probes with any masked values are removed and then intersected with the minfi QC probe set for conservative cross-validation.",
    sprintf("- **sesame native view**: probes with <= %.2f missingness are retained; remaining NAs are imputed via KNN (impute::impute.knn, k=%d) when available, with row-mean fallback if KNN is unavailable.", SESAME_NATIVE_NA_MAX_FRAC, SESAME_NATIVE_KNN_K),
    "",
    "## Covariates and batch control",
    sprintf("- Automatic covariate discovery: metadata variables associated with the top PCs (base alpha=%.3f) are considered; per-variable significance is Bonferroni-corrected by the number of finite PCs tested for that variable (effective alpha = alpha / n_tested_PCs), then filtered for stability/confounding.", AUTO_COVARIATE_ALPHA),
    sprintf("- Auto covariate guard: variables associated with the group (p < %.1g) are excluded by default.", auto_cov_group_p),
    "- Auto-selected covariates may include mediators; verify biological plausibility before interpretation.",
    "- Small-n safeguard: if covariate count would eliminate residual degrees of freedom, covariates are capped using PC-association/variance ranking to preserve model stability.",
    sprintf("- Overcorrection guard: total covariates + SVs are capped at %.0f%% of sample size (excess terms are dropped to preserve power).", OVERCORRECTION_GUARD_RATIO * 100),
    sprintf("- Undercorrection guard: if SVA detects hidden structure, at least %d SV(s) are retained when possible (dropping non-forced covariates first).", UNDERCORRECTION_GUARD_MIN_SV),
    sprintf("- Cell composition: when tissue is `Auto`, IlluMeta attempts to infer tissue from metadata/heuristics; if unresolved, reference-free deconvolution is performed via RefFreeEWAS (K=%d latent components). For `Placenta`, IlluMeta uses planet::plCellCpGsThird with minfi::projectCellType (Houseman) when available. Custom references via `--cell_reference` override defaults when provided.", cell_ref_free_k),
    "- Reference-based deconvolution can be biased if disease alters methylation states; interpret cell-adjusted models cautiously and consider RefFreeEWAS outputs.",
    "- Sesame pipelines attempt Sesame-native cell composition; if unavailable, they try planet (Placenta), EpiDISH (Blood), or RefFreeEWAS on Sesame betas, and fall back to Minfi-derived covariates if needed.",
    sprintf("- Surrogate variable analysis (SVA): enabled unless `--disable_sva`; SVs are estimated on top-variable probes and included in the model only when selected as the best batch strategy or when no batch factor is evaluated (to avoid double correction). SVs strongly associated with the group (P < %.1e or Eta^2 > %.2f) are excluded to avoid over-correction.", SV_GROUP_P_THRESHOLD, SV_GROUP_ETA2_THRESHOLD),
    "- Epigenetic clock covariates: when `--include_clock_covariates` is enabled, clock outputs are merged into metadata and considered by auto covariate selection (clocks with missing/constant values are excluded).",
    "- Batch method comparison: batch candidates are screened for VOI confounding (Cramer's V / R2 and overlap tiers). Tier 3 triggers stratified/meta-analysis fallbacks when eligible; Tier 0-2 proceed to scoring.",
    sprintf("- Tier3 eligibility constraints: min_total_n=%d, min_per_group_per_stratum=%d (on_fail=%s).",
            tier3_min_total_n, tier3_min_per_group_per_stratum, tier3_on_fail),
    sprintf("- Optimization preset: `%s` (weights on batch removal, biology preservation, calibration, stability).", preset_name),
    sprintf("- CRF sample tier: %s (total_n=%d; min_per_group=%d).",
            crf_tier, crf_tier_info$total_n, crf_tier_info$min_per_group),
    "- Calibration: permutation tests shuffle labels within batch strata to assess p-value uniformity (KS) and inflation.",
    "- Lambda guard is a heuristic inflation check; interpret with EWAS correlation structure in mind.",
    sprintf("- Correction Adequacy Framework (CAF): combines calibration (target FPR=%.2f), signal preservation, and batch removal into a single CAI; see `Correction_Adequacy_Report.txt`.", caf_target_fpr),
    "- Decision ledger: automated decisions (covariates, batch choice, method) are logged with reasons.",
    "",
    "## Differential methylation (DMP)",
    sprintf("- Beta values are transformed to M-values with a logit offset of %.6f.", LOGIT_OFFSET),
    "- Differential methylation is tested using limma (`lmFit`/`eBayes(robust=TRUE)`) with a group contrast (test - control).",
    "- Multiple testing is controlled by Benjamini-Hochberg FDR.",
    "",
    "## Differentially methylated regions (DMR)",
    sprintf("- DMRs are called with `dmrff` (maxgap=%d; p.cutoff=%.3f; min_cpgs=%d).", DMR_MAXGAP, dmr_p_cutoff, DMR_MIN_CPGS),
    "",
    "## Consensus (intersection) call set",
    "- Consensus DMPs are defined as CpGs significant in **both** minfi and sesame with the **same direction** under the same thresholds.",
    "- Intersection is intended as a high-confidence subset; pipeline-specific results may capture additional true positives and are reported as sensitivity/discovery sets.",
    "- Consensus table ranking uses Fisher's combined probability test (chi-squared, df=4) with genome-wide Benjamini-Hochberg FDR (`P.Value`/`adj.P.Val`). Selection-rule p-values are retained as `P.Value.selection` and `adj.P.Val.selection` (also mirrored in `*.max` legacy columns).",
    "- Consensus is computed for both the strict (Minfi-aligned) and native Sesame views.",
    "- Consensus outputs: `Intersection_Consensus_DMPs.*` and `Intersection_Native_Consensus_DMPs.*`, plus concordance/overlap plots.",
    "- Primary branch is selected by the optimization score (see `results/consensus/primary_branch.txt`); the other branch is reported as sensitivity.",
    sprintf("- Primary branch override: %s%s.", primary_branch_override_disp, primary_branch_reason_disp),
    sprintf("- Primary inference mode: %s (tier3_batch=%s).", primary_result_mode_disp, primary_tier3_batch_disp),
    tier3_low_power_line,
    tier3_ineligible_line,
    no_signal_line,
    sprintf("- Tier3 meta-analysis (primary): method=%s, I2_median=%s.", primary_tier3_meta_method_disp, primary_tier3_meta_i2_disp),
    sprintf("- Lambda guard (primary): status=%s, action=%s, threshold=%s, lambda_guard_lambda=%s.",
            lambda_guard_status_disp, lambda_guard_action_disp, lambda_guard_threshold_disp, lambda_guard_lambda_disp),
    sprintf("- DMR status (primary): %s%s.", primary_dmr_status_disp, primary_dmr_reason_disp),
    sprintf("- Sesame dye bias normalization: %s%s.", sesame_dyebias_mode, sesame_dyebias_note_disp),
    "",
    "## Reproducibility artifacts",
    "- `analysis_parameters.json`: run parameters and thresholds.",
    "- `sessionInfo.txt`: full R session/package versions.",
    "- `code_version.txt`: git commit hash (when available).",
    "- `config_used.yaml`: resolved config and preset details.",
    "- `decision_ledger.tsv`: auditable record of automated decisions.",
    "",
    "## QC summary (from QC_Summary.csv)",
    sprintf("- Total_samples_input: %s", qc_val("Total_samples_input", NA)),
    sprintf("- Samples_failed_QC: %s", qc_val("Samples_failed_QC", NA)),
    sprintf("- Samples_passed_QC: %s", qc_val("Samples_passed_QC", NA)),
    sprintf("- Samples_failed_sex_mismatch: %s", qc_val("Samples_failed_sex_mismatch", NA)),
    sprintf("- Sex_mismatch_samples: %s", qc_val("Sex_mismatch_samples", NA)),
    sprintf("- Total_probes_raw: %s", qc_val("Total_probes_raw", NA)),
    sprintf("- Probes_cross_reactive: %s", qc_val("Probes_cross_reactive", NA)),
    sprintf("- Probes_final: %s", qc_val("Probes_final", NA))
  )
  writeLines(methods_lines, file.path(out_dir, "methods.md"))
  message("Methods summary saved to methods.md")
  file.copy(file.path(out_dir, "methods.md"), results_dirs$reports, overwrite = TRUE)
}, error = function(e) {
  message("Warning: failed to write methods.md: ", e$message)
})

message("Analysis Complete.")
