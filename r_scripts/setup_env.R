# Prefer a project-local library by default for reproducibility.
# To respect an externally-set R_LIBS_USER, set ILLUMETA_RESPECT_R_LIBS_USER=1.
rp <- function(x) {
    tryCatch(normalizePath(x, winslash = "/", mustWork = FALSE), error = function(e) x)
}

r_version <- getRversion()
r_major_minor <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))

respect_existing <- Sys.getenv("ILLUMETA_RESPECT_R_LIBS_USER", unset = "") == "1"
allow_external_lib <- Sys.getenv("ILLUMETA_ALLOW_EXTERNAL_LIB", unset = "") == "1"
existing_lib <- Sys.getenv("R_LIBS_USER")
lib_root <- file.path(getwd(), ".r-lib")
if (!respect_existing) {
    is_macos <- Sys.info()[["sysname"]] == "Darwin"
    proj_path <- rp(getwd())
    project_on_external <- is_macos && grepl("^/Volumes/", proj_path)
    if (project_on_external && !allow_external_lib) {
        lib_root <- file.path(Sys.getenv("HOME"), ".illumeta", "r-lib")
        message(paste("Project is on external volume; using local R library:", lib_root,
                      "(set ILLUMETA_ALLOW_EXTERNAL_LIB=1 to use .r-lib)"))
    }
}
if (respect_existing && nzchar(existing_lib)) {
    lib_dir <- existing_lib
} else {
    lib_dir <- file.path(lib_root, paste0("R-", r_major_minor))
    if (nzchar(existing_lib) && existing_lib != lib_dir && !respect_existing) {
        message(paste("Overriding R_LIBS_USER:", existing_lib, "->", lib_dir,
                      "(set ILLUMETA_RESPECT_R_LIBS_USER=1 to keep)"))
    }
    Sys.setenv(R_LIBS_USER = lib_dir)
}
if (!dir.exists(lib_dir)) {
    dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
    message(paste("Created R library directory:", lib_dir))
}
base_libs <- unique(c(.Library, .Library.site))
.libPaths(unique(c(lib_dir, base_libs)))
message(paste("Using R library:", .libPaths()[1]))
tryCatch({
    writeLines(lib_dir, con = file.path(getwd(), ".illumeta_r_lib_path"))
}, error = function(e) {
    message("Warning: could not write .illumeta_r_lib_path (", conditionMessage(e), ")")
})

# Prefer HTTPS CRAN (some environments block plain HTTP)
cran_repo <- Sys.getenv("ILLUMETA_CRAN_REPO", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))
if (capabilities("libcurl")) {
    options(download.file.method = "libcurl")
}
options(timeout = max(600, getOption("timeout", 60)))
require_epicv2 <- Sys.getenv("ILLUMETA_REQUIRE_EPICV2", unset = "") == "1"
install_devtools <- Sys.getenv("ILLUMETA_INSTALL_DEVTOOLS", unset = "") == "1"
install_clocks <- Sys.getenv("ILLUMETA_INSTALL_CLOCKS", unset = "") == "1"
download_retries <- suppressWarnings(as.integer(Sys.getenv("ILLUMETA_DOWNLOAD_RETRIES", unset = "2")))
if (is.na(download_retries) || download_retries < 0) download_retries <- 2

# Warn early if a conda toolchain is forced while R is NOT from conda.
# This mismatch is a common cause of "C compiler cannot create executables" and link errors.
toolchain <- c(Sys.getenv("CC"), Sys.getenv("CXX"), Sys.getenv("CXX11"), Sys.getenv("CXX14"), Sys.getenv("CXX17"), Sys.getenv("CXX20"))
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
r_home <- rp(R.home())
conda_prefix_norm <- rp(conda_prefix)
conda_toolchain <- any(grepl("conda", toolchain, ignore.case = TRUE))
r_in_conda <- nzchar(conda_prefix_norm) && startsWith(r_home, conda_prefix_norm)
if (conda_toolchain && !r_in_conda) {
    stop(
        "Detected a conda compiler toolchain (CC/CXX contains 'conda') but R is not from conda.\n",
        "This mixed setup often breaks package compilation.\n",
        "Fix options:\n",
        "  (1) Use the conda environment that provides R (recommended): conda env create -f environment.yml\n",
        "  (2) Or deactivate conda / unset compiler overrides:\n",
        "      conda deactivate\n",
        "      export R_ENVIRON_USER=/dev/null && unset CC CXX CXX11 CXX14 CXX17 CXX20 FC F77\n"
    )
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

ensure_fortran_available <- function() {
    fc <- Sys.which("gfortran")
    if (fc != "") return(invisible(fc))
    sysname <- Sys.info()[["sysname"]]
    message("Warning: 'gfortran' not found in PATH. Some packages (e.g., quadprog, lme4) may fail to compile.")
    if (identical(sysname, "Darwin")) {
        message("  - macOS (Homebrew): brew install gcc")
        message("  - Conda: conda install -c conda-forge gfortran")
    } else if (identical(sysname, "Linux")) {
        message("  - Ubuntu/WSL: sudo apt-get install -y gfortran")
    }
    return(invisible(FALSE))
}

ensure_writable_lib()
cleanup_stale_locks()
ensure_cmake_available()
ensure_fortran_available()

ensure_openmp_flags <- function() {
    if (Sys.info()[["sysname"]] != "Darwin") return(invisible(FALSE))
    libomp_inc <- "/opt/homebrew/opt/libomp/include"
    libomp_lib <- "/opt/homebrew/opt/libomp/lib"
    if (!dir.exists(libomp_inc) || !dir.exists(libomp_lib)) {
        message("Warning: libomp not found. Some Bioconductor packages (e.g., SparseArray) may fail to compile.")
        message("  - macOS (Homebrew): brew install libomp")
        return(invisible(FALSE))
    }
    cflags <- Sys.getenv("PKG_CFLAGS", unset = "")
    cppflags <- Sys.getenv("CPPFLAGS", unset = "")
    cflags_global <- Sys.getenv("CFLAGS", unset = "")
    libs <- Sys.getenv("PKG_LIBS", unset = "")
    ldflags <- Sys.getenv("LDFLAGS", unset = "")
    openmp_cflags <- paste("-Xclang -fopenmp", paste0("-I", libomp_inc))
    openmp_libs <- paste(paste0("-L", libomp_lib), "-lomp")
    shlib_openmp <- paste(openmp_cflags, openmp_libs)
    if (!grepl("-fopenmp", cflags, fixed = TRUE)) {
        cflags <- paste(cflags, openmp_cflags)
    }
    if (!grepl("-fopenmp", cppflags, fixed = TRUE)) {
        cppflags <- paste(cppflags, openmp_cflags)
    }
    if (!grepl("-fopenmp", cflags_global, fixed = TRUE)) {
        cflags_global <- paste(cflags_global, openmp_cflags)
    }
    if (!grepl("-lomp", libs, fixed = TRUE)) {
        libs <- paste(libs, openmp_libs)
    }
    if (!grepl("-lomp", ldflags, fixed = TRUE)) {
        ldflags <- paste(ldflags, openmp_libs)
    }
    Sys.setenv(PKG_CFLAGS = trimws(cflags))
    Sys.setenv(CPPFLAGS = trimws(cppflags))
    Sys.setenv(CFLAGS = trimws(cflags_global))
    Sys.setenv(PKG_LIBS = trimws(libs))
    Sys.setenv(LDFLAGS = trimws(ldflags))
    Sys.setenv(SHLIB_OPENMP_CFLAGS = trimws(shlib_openmp))
    Sys.setenv(SHLIB_OPENMP_CXXFLAGS = trimws(shlib_openmp))
    Sys.setenv(SHLIB_OPENMP_FCFLAGS = trimws(shlib_openmp))
    Sys.setenv(SHLIB_OPENMP_LDFLAGS = trimws(openmp_libs))
    makevars_path <- file.path(tempdir(), "Makevars.illumeta")
    writeLines(c(
        paste0("SHLIB_OPENMP_CFLAGS=", trimws(shlib_openmp)),
        paste0("SHLIB_OPENMP_CXXFLAGS=", trimws(shlib_openmp)),
        paste0("SHLIB_OPENMP_FCFLAGS=", trimws(shlib_openmp)),
        paste0("SHLIB_OPENMP_LDFLAGS=", trimws(openmp_libs))
    ), con = makevars_path)
    Sys.setenv(R_MAKEVARS_USER = makevars_path)
    message("Using custom Makevars for OpenMP: ", makevars_path)
    message("Configured OpenMP flags for macOS (libomp).")
    return(invisible(TRUE))
}

ensure_openmp_flags()

# Prefer conda/system libxml2/icu/curl if available (for xml2/XML builds)
conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (!nzchar(conda_prefix) || !dir.exists(conda_prefix)) {
    conda_prefix <- ""
}
use_conda_libs <- Sys.getenv("ILLUMETA_USE_CONDA_LIBS", unset = "") == "1"
conda_lib <- file.path(conda_prefix, "lib")
conda_pkgconfig <- file.path(conda_lib, "pkgconfig")
add_path <- function(var, path) {
    if (!dir.exists(path)) return(FALSE)
    cur <- Sys.getenv(var, unset = "")
    parts <- if (nzchar(cur)) strsplit(cur, ":", fixed = TRUE)[[1]] else character(0)
    if (path %in% parts) return(FALSE)
    new_val <- paste(c(path, parts), collapse = ":")
    do.call(Sys.setenv, setNames(list(new_val), var))
    return(TRUE)
}
if (nzchar(conda_prefix) && (use_conda_libs || r_in_conda)) {
    if (add_path("LD_LIBRARY_PATH", conda_lib)) {
        message(paste("Added to LD_LIBRARY_PATH:", conda_lib))
    }
    if (add_path("PKG_CONFIG_PATH", conda_pkgconfig)) {
        message(paste("Added to PKG_CONFIG_PATH:", conda_pkgconfig))
    }
} else if (nzchar(conda_prefix) && !r_in_conda) {
    message("Note: CONDA_PREFIX is set but R is not from conda; skipping conda libs.",
            " Set ILLUMETA_USE_CONDA_LIBS=1 to override.")
}

sys_which_any <- function(candidates) {
    for (cand in candidates) {
        path <- Sys.which(cand)
        if (nzchar(path)) return(path)
    }
    return("")
}

r_cmd_config <- function(key) {
    cmd <- file.path(R.home("bin"), "R")
    out <- tryCatch(system2(cmd, c("CMD", "config", key), stdout = TRUE, stderr = TRUE), error = function(e) character(0))
    if (length(out) == 0) return("")
    trimws(out[1])
}

append_makevars_lines <- function(lines) {
    if (length(lines) == 0) return(invisible(FALSE))
    makevars_path <- Sys.getenv("R_MAKEVARS_USER", unset = "")
    if (!nzchar(makevars_path)) {
        makevars_path <- file.path(tempdir(), "Makevars.illumeta")
        Sys.setenv(R_MAKEVARS_USER = makevars_path)
    }
    existing <- character(0)
    if (file.exists(makevars_path)) {
        existing <- tryCatch(readLines(makevars_path, warn = FALSE), error = function(e) character(0))
    }
    new_lines <- lines[!lines %in% existing]
    if (length(new_lines) > 0) {
        cat(paste0(new_lines, "\n"), file = makevars_path, append = TRUE)
        message("Updated Makevars for compiler settings: ", makevars_path)
    }
    return(invisible(TRUE))
}

ensure_c17_compiler <- function() {
    if (nzchar(Sys.getenv("CC17", unset = ""))) return(invisible(TRUE))
    cc <- Sys.getenv("CC", unset = "")
    if (!nzchar(cc)) cc <- r_cmd_config("CC")
    if (!nzchar(cc)) cc <- sys_which_any(c("clang", "gcc"))
    if (!nzchar(cc)) {
        message("Warning: CC17 not set and no C compiler found. C17 packages may fail to compile.")
        return(invisible(FALSE))
    }
    Sys.setenv(CC17 = cc)

    c17flags <- Sys.getenv("C17FLAGS", unset = "")
    if (!nzchar(c17flags)) c17flags <- r_cmd_config("C17FLAGS")
    if (!nzchar(c17flags)) c17flags <- Sys.getenv("CFLAGS", unset = "")
    if (!grepl("-std=", c17flags, fixed = FALSE)) {
        c17flags <- paste(c17flags, "-std=gnu17")
    }
    Sys.setenv(C17FLAGS = trimws(c17flags))
    append_makevars_lines(c(
        paste0("CC17=", cc),
        paste0("C17FLAGS=", trimws(c17flags))
    ))
    message("Configured CC17/C17FLAGS for C17 packages.")
    return(invisible(TRUE))
}

ensure_cxx17_compiler <- function() {
    if (nzchar(Sys.getenv("CXX17", unset = ""))) return(invisible(TRUE))
    cxx <- Sys.getenv("CXX", unset = "")
    if (!nzchar(cxx)) cxx <- r_cmd_config("CXX")
    if (!nzchar(cxx)) cxx <- sys_which_any(c("clang++", "g++"))
    if (!nzchar(cxx)) {
        message("Warning: CXX17 not set and no C++ compiler found. C++17 packages may fail to compile.")
        return(invisible(FALSE))
    }
    Sys.setenv(CXX17 = cxx)

    cxx17flags <- Sys.getenv("CXX17FLAGS", unset = "")
    if (!nzchar(cxx17flags)) cxx17flags <- r_cmd_config("CXX17FLAGS")
    if (!nzchar(cxx17flags)) cxx17flags <- Sys.getenv("CXXFLAGS", unset = "")
    if (!grepl("-std=", cxx17flags, fixed = FALSE)) {
        cxx17flags <- paste(cxx17flags, "-std=gnu++17")
    }
    Sys.setenv(CXX17FLAGS = trimws(cxx17flags))
    append_makevars_lines(c(
        paste0("CXX17=", cxx),
        paste0("CXX17FLAGS=", trimws(cxx17flags))
    ))
    message("Configured CXX17/CXX17FLAGS for C++17 packages.")
    return(invisible(TRUE))
}

ensure_c17_compiler()
ensure_cxx17_compiler()

clean_mismatched <- Sys.getenv("ILLUMETA_CLEAN_MISMATCHED_RLIB", unset = "") == "1"
pkg_built_map <- NULL

built_major_minor <- function(built) {
    if (is.na(built) || !nzchar(built)) return(NA_character_)
    parts <- strsplit(built, "\\.", fixed = FALSE)[[1]]
    if (length(parts) < 2) return(NA_character_)
    paste(parts[1], parts[2], sep = ".")
}

get_pkg_built_map <- function() {
    if (!is.null(pkg_built_map)) return(pkg_built_map)
    ip <- tryCatch(installed.packages(lib.loc = .libPaths()[1]), error = function(e) NULL)
    if (is.null(ip) || !"Built" %in% colnames(ip)) {
        pkg_built_map <<- character(0)
        return(pkg_built_map)
    }
    built <- ip[, "Built"]
    names(built) <- rownames(ip)
    pkg_built_map <<- built
    pkg_built_map
}

is_pkg_mismatch <- function(pkg) {
    built_map <- get_pkg_built_map()
    if (!pkg %in% names(built_map)) return(FALSE)
    built_mm <- built_major_minor(built_map[[pkg]])
    !is.na(built_mm) && built_mm != r_major_minor
}

remove_pkg_dir <- function(pkg, lib = .libPaths()[1]) {
    pkg_path <- file.path(lib, pkg)
    if (!dir.exists(pkg_path)) return(invisible(FALSE))
    message("Removing package:", pkg)
    tryCatch(unlink(pkg_path, recursive = TRUE, force = TRUE), error = function(e) {
        message("  Warning: could not remove ", pkg, " (", conditionMessage(e), ")")
    })
    return(invisible(TRUE))
}

mismatch_pkgs <- function() {
    built_map <- get_pkg_built_map()
    if (length(built_map) == 0) return(character(0))
    built_mm <- vapply(built_map, built_major_minor, character(1))
    names(built_mm)[!is.na(built_mm) & built_mm != r_major_minor]
}

mismatched <- mismatch_pkgs()
if (length(mismatched) > 0) {
    preview <- paste(head(mismatched, 8), collapse = ", ")
    suffix <- if (length(mismatched) > 8) " ..." else ""
    message("Detected packages built for a different R version in ", .libPaths()[1], ": ", preview, suffix)
    if (clean_mismatched) {
        message("Cleaning mismatched packages (ILLUMETA_CLEAN_MISMATCHED_RLIB=1).")
        for (pkg in mismatched) remove_pkg_dir(pkg)
        pkg_built_map <- NULL
    } else {
        message("Set ILLUMETA_CLEAN_MISMATCHED_RLIB=1 to purge these and reinstall cleanly.")
    }
}

suppressPackageStartupMessages({
  if (is_pkg_mismatch("optparse")) remove_pkg_dir("optparse")
  if (!require("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = cran_repo, lib = .libPaths()[1])
    library(optparse)
  }
  if (is_pkg_mismatch("remotes")) remove_pkg_dir("remotes")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = cran_repo, lib = .libPaths()[1])
    library(remotes)
  }
})

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo, lib = .libPaths()[1])
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("ERROR: Failed to install BiocManager. Cannot install Bioconductor packages.")
    quit(status = 1)
}

infer_bioc_version <- function() {
    if (r_version >= "4.5") return("3.22")
    if (r_version >= "4.4") return("3.20")
    if (r_version >= "4.3") return("3.18")
    return("")
}

ensure_bioc_version <- function() {
    target <- Sys.getenv("ILLUMETA_BIOC_VERSION", unset = "")
    if (!nzchar(target)) {
        target <- infer_bioc_version()
    }
    if (!nzchar(target)) return(invisible(""))

    current <- tryCatch(as.character(BiocManager::version()), error = function(e) "")
    if (!nzchar(current) || current != target) {
        message(paste("Setting Bioconductor version to", target, "for R", as.character(r_version)))
        tryCatch({
            need_force <- FALSE
            if (nzchar(current)) {
                need_force <- tryCatch(package_version(current) > package_version(target), error = function(e) FALSE)
            }
            install_args <- list(version = target, ask = FALSE, update = FALSE)
            if (need_force && "force" %in% names(formals(BiocManager::install))) {
                install_args$force <- TRUE
            }
            do.call(BiocManager::install, install_args)
        }, error = function(e) {
            message("Warning: Failed to set Bioconductor version automatically.")
            message("  Detail: ", conditionMessage(e))
        })
    }
    options(BiocManager.version = target)
    current_after <- tryCatch(as.character(BiocManager::version()), error = function(e) "")
    if (nzchar(current_after) && current_after != target) {
        message("Warning: Bioconductor version is still ", current_after, "; expected ", target, ".")
    }
    return(invisible(target))
}

configure_bioc_repos <- function(target = "") {
    repos <- tryCatch(BiocManager::repositories(), error = function(e) NULL)
    if ((is.null(repos) || length(repos) == 0) && nzchar(target)) {
        repos <- tryCatch(BiocManager::repositories(version = target), error = function(e) NULL)
    }
    if (is.null(repos) || length(repos) == 0) {
        message("Warning: Failed to configure Bioconductor repositories; using CRAN only.")
        if (nzchar(target)) {
            message("  Hint: set ILLUMETA_BIOC_VERSION (e.g., 3.22 for R 4.5) and rerun.")
        }
        return(invisible(FALSE))
    }
    repos["CRAN"] <- cran_repo
    options(repos = repos)
    message(paste("Using Bioconductor repositories (Bioc", as.character(BiocManager::version()), ") + CRAN:", cran_repo))
    return(invisible(TRUE))
}
bioc_target <- ensure_bioc_version()
configure_bioc_repos(bioc_target)

# Core BioC packages
bioc_pkgs <- c(
    "minfi", 
    "sesame", 
    "limma", 
    "GEOquery", 
    "Biobase",
    "illuminaio",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "sesameData", # Required for Sesame annotation caching
    "sva", # Added for Surrogate Variable Analysis
    "variancePartition",
    "pvca",
    "FlowSorted.Blood.EPIC",   # Cell type reference (EPIC/EPICv2)
    "FlowSorted.Blood.450k"   # Cell type reference (450k)
)

clock_pkgs <- c(
    "wateRmelon",        # Epigenetic clocks (Horvath/Hannum/PhenoAge)
    "planet",            # Placental DNA methylation clocks
    "methylclock",       # General epigenetic clocks
    "methylclockData"    # Data for methylclock
)

# CRAN packages
cran_pkgs <- c(
    "optparse", 
    "xml2",
    "XML",
    "rlang",
    "stringi",
    "glue",
    "digest",
    "matrixStats",
    "fs",
    "quadprog",
    "dplyr", 
    "stringr", 
    "ggplot2", 
    "plotly", 
    "DT", 
    "htmlwidgets", 
    "data.table",
    "qqman",
    "ggrepel",
    "remotes",
    "lme4",
    # Add other core tidyverse components explicitly if needed by methylclock
    # For now, rely on archived tidyverse for simplicity
    "tibble",
    "purrr",
    "readr",
    "tidyr",
    "forcats"
)

bioc_available_cache <- NULL
get_bioc_available <- function() {
    if (!is.null(bioc_available_cache)) return(bioc_available_cache)
    bioc_available_cache <<- tryCatch(BiocManager::available(), error = function(e) {
        message("Warning: Could not fetch Bioconductor package index: ", conditionMessage(e))
        NULL
    })
    bioc_available_cache
}

is_bioc_available <- function(pkg) {
    avail <- get_bioc_available()
    if (is.null(avail)) return(NA)
    pkg %in% avail
}

bioc_version <- tryCatch(as.character(BiocManager::version()), error = function(e) "unknown")

clean_downloaded_packages <- function() {
    dl_dir <- file.path(tempdir(), "downloaded_packages")
    if (!dir.exists(dl_dir)) return(invisible(FALSE))
    tryCatch({
        unlink(dl_dir, recursive = TRUE, force = TRUE)
        message("Cleared temporary downloaded_packages: ", dl_dir)
    }, error = function(e) {
        message("  Warning: could not clear downloaded_packages (", conditionMessage(e), ")")
    })
    return(invisible(TRUE))
}

is_loadable <- function(pkg) {
    if (is_pkg_mismatch(pkg)) return(FALSE)
    ok <- tryCatch(requireNamespace(pkg, quietly = TRUE), error = function(e) FALSE)
    isTRUE(ok)
}

remove_broken_pkg <- function(pkg, lib = .libPaths()[1]) {
    if (remove_pkg_dir(pkg, lib = lib)) {
        pkg_built_map <<- NULL
    }
}

install_bioc_safe <- function(pkgs) {
    miss <- pkgs[!vapply(pkgs, is_loadable, logical(1))]
    if (length(miss) == 0) return(invisible())
    for (pkg in miss) remove_broken_pkg(pkg)
    message(paste("Installing Bioconductor packages:", paste(miss, collapse = ", ")))
    tryCatch({
        BiocManager::install(miss, update = FALSE, ask = FALSE, lib = .libPaths()[1])
    }, error = function(e) {
        message("Warning: BiocManager::install failed for: ", paste(miss, collapse = ", "))
        message("  Detail: ", conditionMessage(e))
    })
    miss_after <- miss[!vapply(miss, is_loadable, logical(1))]
    if (length(miss_after) > 0 && download_retries > 0) {
        for (i in seq_len(download_retries)) {
            clean_downloaded_packages()
            message("Retrying Bioconductor packages (attempt ", i, "/", download_retries, "): ",
                    paste(miss_after, collapse = ", "))
            tryCatch({
                BiocManager::install(miss_after, update = FALSE, ask = FALSE, lib = .libPaths()[1])
            }, error = function(e) {
                message("Warning: BiocManager::install retry failed for: ", paste(miss_after, collapse = ", "))
                message("  Detail: ", conditionMessage(e))
            })
            miss_after <- miss_after[!vapply(miss_after, is_loadable, logical(1))]
            if (length(miss_after) == 0) break
        }
    }
}

install_cran_safe <- function(pkgs) {
    miss <- pkgs[!vapply(pkgs, is_loadable, logical(1))]
    if (length(miss) == 0) return(invisible())
    for (pkg in miss) remove_broken_pkg(pkg)
    message(paste("Installing CRAN packages:", paste(miss, collapse = ", ")))
    tryCatch({
        install.packages(miss, repos = cran_repo, lib = .libPaths()[1])
    }, error = function(e) {
        message("Warning: install.packages failed for: ", paste(miss, collapse = ", "))
        message("  Detail: ", conditionMessage(e))
    })
    miss_after <- miss[!vapply(miss, is_loadable, logical(1))]
    if (length(miss_after) > 0 && download_retries > 0) {
        for (i in seq_len(download_retries)) {
            clean_downloaded_packages()
            message("Retrying CRAN packages (attempt ", i, "/", download_retries, "): ",
                    paste(miss_after, collapse = ", "))
            tryCatch({
                install.packages(miss_after, repos = cran_repo, lib = .libPaths()[1])
            }, error = function(e) {
                message("Warning: install.packages retry failed for: ", paste(miss_after, collapse = ", "))
                message("  Detail: ", conditionMessage(e))
            })
            miss_after <- miss_after[!vapply(miss_after, is_loadable, logical(1))]
            if (length(miss_after) == 0) break
        }
    }
}

install_archived_safe <- function(pkg, version) {
    if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
    message(paste0("Installing archived ", pkg, " (", version, ") to avoid pkgdown/ragg issues..."))
    tryCatch({
        remotes::install_version(pkg, version = version, repos = cran_repo, dependencies = NA, upgrade = "never", lib = .libPaths()[1])
    }, error = function(e) {
        message("Warning: Failed to install archived ", pkg, ".")
        message("  Detail: ", conditionMessage(e))
    })
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Warning: Falling back to latest CRAN version of ", pkg, ".")
        install_cran_safe(pkg)
    }
}

install_cran_safe(cran_pkgs)
if (install_devtools) {
    install_archived_safe("devtools", "2.4.3")
    install_archived_safe("tidyverse", "1.3.2")
} else {
    message("Skipping devtools/tidyverse install (set ILLUMETA_INSTALL_DEVTOOLS=1 to enable).")
}

install_bioc_safe(bioc_pkgs)
if (install_clocks) {
    message("Installing optional clock packages (set ILLUMETA_INSTALL_CLOCKS=0 to skip)...")
    install_bioc_safe(clock_pkgs)
} else {
    message("Skipping optional clock packages (set ILLUMETA_INSTALL_CLOCKS=1 to install).")
}

ensure_min_version <- function(pkg, min_ver, installer_fn, lib = .libPaths()[1]) {
    cur_ver <- tryCatch(packageVersion(pkg), error = function(e) NULL)
    if (is.null(cur_ver) || cur_ver < as.package_version(min_ver)) {
        tryCatch(installer_fn(pkg, lib = lib), error = function(e) {
            message(paste("Warning: version check install failed for", pkg, "-", e$message))
        })
    }
}

# Refresh packages implicated in variancePartition/reformulas changes
ensure_min_version("reformulas", "0.3.0", function(p, lib) install.packages(p, repos = cran_repo, lib = lib))
ensure_min_version("lme4", "1.1-35", function(p, lib) install.packages(p, repos = cran_repo, lib = lib))
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

# Install RefFreeEWAS from CRAN Archive (version 2.2)
if (!requireNamespace("RefFreeEWAS", quietly = TRUE)) {
    message("Installing RefFreeEWAS (v2.2) from CRAN Archive...")
    tryCatch({
        remotes::install_version("RefFreeEWAS", version = "2.2", repos = "https://cloud.r-project.org", lib = .libPaths()[1], upgrade = "never")
    }, error = function(e) {
        message("Warning: Failed to install RefFreeEWAS from CRAN Archive.")
        message(paste("Error details:", e$message))
    })
}

install_epicv2_manifest <- function() {
  if (requireNamespace("IlluminaHumanMethylationEPICv2manifest", quietly = TRUE)) return(invisible(TRUE))
  avail <- is_bioc_available("IlluminaHumanMethylationEPICv2manifest")
  if (identical(avail, FALSE)) {
      message(paste("Note: EPIC v2 manifest is not available for Bioconductor", bioc_version, "(R", r_version, ")."))
      message("To enable EPIC v2 support, use R 4.4+ with Bioconductor 3.19+ (see README).")
      if (require_epicv2) {
          message("ERROR: EPIC v2 manifest required but not available for this R/Bioconductor version.")
          quit(status = 1)
      }
      return(invisible(FALSE))
  }
  message("Installing EPIC v2 manifest (IlluminaHumanMethylationEPICv2manifest)...")
  install_bioc_safe("IlluminaHumanMethylationEPICv2manifest")
  if (requireNamespace("IlluminaHumanMethylationEPICv2manifest", quietly = TRUE)) return(invisible(TRUE))
  if (require_epicv2) {
      message("ERROR: EPIC v2 manifest failed to install. Check Bioconductor configuration.")
      quit(status = 1)
  }
  message("Warning: EPIC v2 manifest not installed. EPIC v2 arrays will be unsupported.")
  return(invisible(FALSE))
}

install_epicv2_annos <- function() {
  # Install the official Bioconductor annotation package for EPIC v2
  if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38", quietly = TRUE)) {
      avail <- is_bioc_available("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
      if (identical(avail, FALSE)) {
          message(paste("Note: EPIC v2 annotation is not available for Bioconductor", bioc_version, "(R", r_version, ")."))
          message("To enable EPIC v2 support, use R 4.4+ with Bioconductor 3.19+ (see README).")
          if (require_epicv2) {
              message("ERROR: EPIC v2 annotation required but not available for this R/Bioconductor version.")
              quit(status = 1)
          }
          return(invisible(FALSE))
      }
      message("Installing EPIC v2 annotation (IlluminaHumanMethylationEPICv2anno.20a1.hg38)...")
      install_bioc_safe("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
  }
  if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38", quietly = TRUE)) {
      if (require_epicv2) {
          message("ERROR: EPIC v2 annotation failed to install. Check Bioconductor configuration.")
          quit(status = 1)
      }
      message("Warning: EPIC v2 annotation not installed. EPIC v2 arrays will be unsupported.")
      return(invisible(FALSE))
  }
  return(invisible(TRUE))
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
epicv2_pkgs <- c(
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
)
required_pkgs <- c(
  "xml2", "XML",
  "lme4", "reformulas", "illuminaio",
  "minfi", "sesame", "limma", "dmrff", "GEOquery",
  "Biobase",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "sesameData", "sva", "variancePartition", "pvca",
  "RefFreeEWAS"
)
required_pkgs <- c(required_pkgs, epicv2_pkgs)
if (!require_epicv2) {
  required_pkgs <- setdiff(required_pkgs, epicv2_pkgs)
}
if (install_clocks) {
  required_pkgs <- c(required_pkgs, clock_pkgs)
}

failed <- character(0)
errors <- new.env(parent = emptyenv())
installed_pkgs <- rownames(installed.packages())
for (pkg in required_pkgs) {
  ok <- tryCatch({
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    errors[[pkg]] <- conditionMessage(e)
    FALSE
  })
  if (!ok) {
    if (is.null(errors[[pkg]])) {
      errors[[pkg]] <- if (pkg %in% installed_pkgs) "failed to load" else "not installed"
    }
    failed <- c(failed, pkg)
  }
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

epicv2_installed <- all(sapply(epicv2_pkgs, function(pkg) requireNamespace(pkg, quietly = TRUE)))
if (epicv2_installed) {
  message("EPIC v2 manifest/annotation available.")
} else if (require_epicv2) {
  message("ERROR: EPIC v2 packages are required but not installed.")
  message("Please ensure Bioconductor is configured correctly and use R 4.4+ with Bioconductor 3.19+.")
  quit(status = 1)
} else {
  message("Note: EPIC v2 support is not installed; EPIC v2 arrays will be skipped.")
  message("To enable, use R 4.4+ with Bioconductor 3.19+ and set ILLUMETA_REQUIRE_EPICV2=1.")
}

message("Setup complete.")
