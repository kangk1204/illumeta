# Prefer a project-local library unless the caller overrides R_LIBS_USER
lib_dir <- Sys.getenv("R_LIBS_USER")
if (lib_dir == "") {
    lib_dir <- file.path(getwd(), ".r-lib")
}
if (!dir.exists(lib_dir)) {
    dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
    message(paste("Created project library directory:", lib_dir))
}
.libPaths(c(lib_dir, .libPaths()))
message(paste("Using R library:", .libPaths()[1]))

# Warn early if a conda toolchain is forced via ~/.Renviron or CC/CXX envs
toolchain <- c(Sys.getenv("CC"), Sys.getenv("CXX"), Sys.getenv("CXX11"), Sys.getenv("CXX14"), Sys.getenv("CXX17"), Sys.getenv("CXX20"))
if (any(grepl("conda", toolchain, ignore.case = TRUE))) {
    stop("CC/CXX point to a conda toolchain. Disable those overrides (e.g., `export R_ENVIRON_USER=/dev/null && unset CC CXX CXX11 CXX14 CXX17 CXX20 FC F77`) and rerun.")
}

ensure_writable_lib <- function() {
    writable <- file.access(.libPaths()[1], mode = 2) == 0
    if (writable) return(invisible(TRUE))
    
    user_lib <- Sys.getenv("R_LIBS_USER")
    if (user_lib == "") {
        user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
    }
    if (!dir.exists(user_lib)) {
        dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
        message(paste("Created user library directory:", user_lib))
    }
    .libPaths(c(user_lib, .libPaths()))
    message(paste("System library not writable. Using user library:", user_lib))
}

cleanup_stale_locks <- function(libdir = .libPaths()[1]) {
    locks <- Sys.glob(file.path(libdir, "00LOCK*"))
    if (length(locks) == 0) return(invisible())
    message("Removing stale package locks: ", paste(basename(locks), collapse = ", "))
    for (lk in locks) {
        tryCatch(unlink(lk, recursive = TRUE, force = TRUE), error = function(e) {
            message("  Warning: could not remove ", lk, " (", conditionMessage(e), ")")
        })
    }
}

ensure_cmake_available <- function() {
    cmake <- Sys.which("cmake")
    if (cmake != "") return(invisible(cmake))
    message("ERROR: 'cmake' not found in PATH. Packages nloptr/lme4/variancePartition need it.")
    message("Install cmake then rerun setup_env.R. Examples:")
    message("  - Ubuntu/WSL: sudo apt-get install -y cmake")
    message("  - macOS (Homebrew): brew install cmake")
    message("  - No sudo: python3 -m pip install cmake  # ensure ~/.local/bin or your venv bin is on PATH")
    quit(status = 1)
}

ensure_writable_lib()
cleanup_stale_locks()
ensure_cmake_available()

# Prefer conda/system libxml2/icu/curl if available (for xml2/XML builds)
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = Sys.getenv("CONDA_DEFAULT_ENV", unset = ""))
if (nzchar(conda_prefix) && basename(conda_prefix) == Sys.getenv("CONDA_DEFAULT_ENV")) {
    conda_prefix <- Sys.getenv("CONDA_PREFIX")
}
conda_lib <- file.path(conda_prefix, "lib")
conda_pkgconfig <- file.path(conda_lib, "pkgconfig")
add_path <- function(var, path) {
    cur <- Sys.getenv(var, unset = "")
    parts <- if (nzchar(cur)) strsplit(cur, ":", fixed = TRUE)[[1]] else character(0)
    if (!(path %in% parts) && dir.exists(path)) {
        Sys.setenv(`var` = paste(c(path, parts), collapse = ":"))
    }
}
if (dir.exists(conda_lib)) {
    if (!grepl(conda_lib, Sys.getenv("LD_LIBRARY_PATH", ""))) {
        Sys.setenv(LD_LIBRARY_PATH = paste(conda_lib, Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
        message(paste("Added to LD_LIBRARY_PATH:", conda_lib))
    }
}
if (dir.exists(conda_pkgconfig)) {
    if (!grepl(conda_pkgconfig, Sys.getenv("PKG_CONFIG_PATH", ""))) {
        Sys.setenv(PKG_CONFIG_PATH = paste(conda_pkgconfig, Sys.getenv("PKG_CONFIG_PATH"), sep = ":"))
        message(paste("Added to PKG_CONFIG_PATH:", conda_pkgconfig))
    }
}

suppressPackageStartupMessages({
  if (!require("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
    library(optparse)
  }
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
    library(remotes)
  }
})

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

# Core BioC packages
bioc_pkgs <- c(
    "minfi", 
    "sesame", 
    "limma", 
    "GEOquery", 
    "Biobase",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICv2manifest", # EPIC v2 manifest required
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "sesameData", # Required for Sesame annotation caching
    "sva", # Added for Surrogate Variable Analysis
    "variancePartition",
    "pvca",
    "FlowSorted.Blood.EPIC",   # Cell type reference (EPIC/EPICv2)
    "FlowSorted.Blood.450k"    # Cell type reference (450k)
)

# EPIC v2 annotation names (accept any)
epicv2_annos <- c(
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38", # Bioc name
    "IlluminaHumanMethylationEPICv2anno.ilm10b5.hg38", # alt Bioc name seen before
    "IlluminaHumanMethylationEPICanno.ilm10b5.hg38" # GitHub legacy name
)

# CRAN packages
cran_pkgs <- c(
    "optparse", 
    "dplyr", 
    "stringr", 
    "ggplot2", 
    "plotly", 
    "DT", 
    "htmlwidgets", 
    "data.table",
    "qqman",
    "ggrepel",
    "remotes", # Added for install_github
    "reformulas",
    "lme4"
)

install_bioc_safe <- function(pkgs) {
    miss <- setdiff(pkgs, rownames(installed.packages()))
    if (length(miss) == 0) return(invisible())
    tryCatch({
        BiocManager::install(miss, update = FALSE, ask = FALSE, lib = .libPaths()[1])
    }, error = function(e) {
        message("Warning: BiocManager::install failed for: ", paste(miss, collapse = ", "))
        message("  Detail: ", conditionMessage(e))
    })
}

install_cran_safe <- function(pkgs) {
    miss <- setdiff(pkgs, rownames(installed.packages()))
    if (length(miss) == 0) return(invisible())
    tryCatch({
        install.packages(miss, repos = "http://cran.us.r-project.org", lib = .libPaths()[1])
    }, error = function(e) {
        message("Warning: install.packages failed for: ", paste(miss, collapse = ", "))
        message("  Detail: ", conditionMessage(e))
    })
}

install_bioc_safe(bioc_pkgs)
install_cran_safe(cran_pkgs)

ensure_min_version <- function(pkg, min_ver, installer_fn, lib = .libPaths()[1]) {
    cur_ver <- tryCatch(packageVersion(pkg), error = function(e) NULL)
    if (is.null(cur_ver) || cur_ver < as.package_version(min_ver)) {
        tryCatch(installer_fn(pkg, lib = lib), error = function(e) {
            message(paste("Warning: version check install failed for", pkg, "-", e$message))
        })
    }
}

# Refresh packages implicated in variancePartition/reformulas changes
ensure_min_version("reformulas", "0.3.0", function(p, lib) install.packages(p, repos = "http://cran.us.r-project.org", lib = lib))
ensure_min_version("lme4", "1.1-35", function(p, lib) install.packages(p, repos = "http://cran.us.r-project.org", lib = lib))
ensure_min_version("variancePartition", "1.30.0", function(p, lib) BiocManager::install(p, update = TRUE, ask = FALSE, lib = lib))

# Install dmrff from GitHub
if (!requireNamespace("dmrff", quietly = TRUE)) {
    message("Installing dmrff from GitHub...")
    tryCatch({
        remotes::install_github("perishky/dmrff", upgrade = "never", lib = .libPaths()[1])
    }, error = function(e) {
        message("Warning: Failed to install dmrff from GitHub.")
        message(paste("Error details:", e$message))
    })
}

install_epicv2_manifest <- function() {
  if (requireNamespace("IlluminaHumanMethylationEPICv2manifest", quietly = TRUE)) return(invisible(TRUE))
  message("Attempting EPIC v2 manifest from Bioconductor...")
  install_bioc_safe("IlluminaHumanMethylationEPICv2manifest")
  if (requireNamespace("IlluminaHumanMethylationEPICv2manifest", quietly = TRUE)) return(invisible(TRUE))
  # Try GitHub variants
  manifest_repos <- c("achilleasNP/IlluminaHumanMethylationEPICv2manifest",
                      "achilleasNP/IlluminaHumanMethylationEPICmanifest")
  for (repo in manifest_repos) {
    message(paste("Attempting EPIC v2 manifest from GitHub:", repo))
    tryCatch({
      remotes::install_github(repo, upgrade = "never", lib = .libPaths()[1])
      if (requireNamespace("IlluminaHumanMethylationEPICv2manifest", quietly = TRUE)) return(invisible(TRUE))
    }, error = function(e) {
      message("Warning: Failed to install EPIC v2 manifest from GitHub.")
      message(paste("Error details:", e$message))
    })
  }
}

install_epicv2_annos <- function() {
  # Try Bioc names
  tryCatch({
    install_bioc_safe("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
    install_bioc_safe("IlluminaHumanMethylationEPICv2anno.ilm10b5.hg38")
  }, error = function(e) {})
  # Try GitHub legacy name
  if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b5.hg38", quietly = TRUE)) {
    message("Attempting EPIC v2 annotation (achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38) from GitHub...")
    tryCatch({
      remotes::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38", upgrade = "never", lib = .libPaths()[1])
    }, error = function(e) {
      message("Warning: Failed to install EPIC v2 annotation from GitHub (ilm10b5).")
      message(paste("Error details:", e$message))
    })
  }
}

install_epicv2_manifest()
install_epicv2_annos()

# --- Configure Local Cache (Project Root/cache) ---
# To save space on internal drive, we use a local 'cache' folder in the project directory.
# We assume this script is run from the project root (where 'cache' should be).
local_cache_dir <- file.path(getwd(), "cache")

if (!dir.exists(local_cache_dir)) {
    dir.create(local_cache_dir, recursive = TRUE)
    message(paste("Created local cache directory:", local_cache_dir))
}

# Force BiocFileCache and other tools to use this directory via Environment Variables
Sys.setenv(R_USER_CACHE_DIR = local_cache_dir)
Sys.setenv(XDG_CACHE_HOME = local_cache_dir)
message(paste("Set R_USER_CACHE_DIR and XDG_CACHE_HOME to:", local_cache_dir))

# Set Options for Cache *BEFORE* loading packages that might use it
options(ExperimentHub.cache = local_cache_dir)
options(AnnotationHub.cache = local_cache_dir)
message(paste("Setting ExperimentHub/AnnotationHub cache to:", local_cache_dir))

# --- Pre-download Sesame Annotation Data ---
# This ensures data is available locally to prevent runtime download errors
message("Pre-downloading/Caching Sesame annotation data (HM450, EPIC, EPICv2)...")

tryCatch({
    library(sesame)
    library(sesameData)
    
    # Cache specific platforms
    # Recent versions of sesameData use `sesameDataCache()` to cache recommended data.
    # Specific keys like "HM450" might be deprecated or require specific filenames.
    # Calling sesameDataCache() without arguments caches the default/essential datasets.
    
    message("  - Running sesameDataCache() for essential data...")
    sesameDataCache()
    
    # Attempt to ensure specific platforms are present if not covered by default
    # We use tryCatch for each to avoid stopping the script if one key is wrong/already present
    
    platforms <- c("HM450", "EPIC", "EPICv2")
    
    for (plt in platforms) {
        tryCatch({
            message(paste("  - Attempting explicit cache for:", plt))
            sesameDataCache(plt)
        }, error = function(e) {
            message(paste("    (Note: specific cache for", plt, "skipped or failed. It might be included in defaults or key is different.)"))
        })
    }
    
    message("Sesame data caching complete.")
    
}, error = function(e) {
    message("Warning: Failed to run sesameDataCache. Analysis may try to download files at runtime.")
    message(paste("Error detail:", e$message))
})

# --- Verify required packages load ---
required_pkgs <- c(
  "xml2", "XML",
  "minfi", "sesame", "limma", "dmrff", "GEOquery",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "sesameData", "sva", "variancePartition"
)

failed <- character(0)
errors <- list()
for (pkg in required_pkgs) {
  ok <- tryCatch({
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    errors[[pkg]] <<- conditionMessage(e)
    FALSE
  })
  if (!ok) failed <- c(failed, pkg)
}

if (length(failed) > 0) {
  message("ERROR: Some required R packages failed to install or load:")
  for (pkg in failed) {
    msg <- errors[[pkg]]
    message(paste0("  - ", pkg, ": ", ifelse(is.null(msg), "unknown error", msg)))
  }
  message("\nCommon fixes on Ubuntu:")
  message("  sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev libicu-dev")
  message("  export R_LIBS_USER=\"$HOME/R/library\" && mkdir -p \"$R_LIBS_USER\"")
  message("  conda deactivate  # if a conda R is interfering with system libs")
  message("Then rerun: Rscript r_scripts/setup_env.R")
  quit(status = 1)
} else {
  message("All required dependencies installed and loadable.")
}

# EPIC v2 annotation check (must have at least one)
epicv2_anno_ok <- any(vapply(epicv2_annos, function(pkg) {
  tryCatch(requireNamespace(pkg, quietly = TRUE), error = function(e) FALSE)
}, logical(1)))

if (!epicv2_anno_ok) {
  message("ERROR: EPIC v2 annotation package not installed. At least one of the following is required:")
  for (pkg in epicv2_annos) message(paste0("  - ", pkg))
  message("Tried Bioconductor (IlluminaHumanMethylationEPICv2anno.20a1.hg38 / ...ilm10b5...) and GitHub (achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38).")
  message("If unavailable for your R/Bioc version, try upgrading Bioc or installing from source when it becomes available.")
  quit(status = 1)
}

message("EPIC v2 manifest/annotation available. Setup complete.")
