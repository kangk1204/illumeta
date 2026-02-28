#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

options(timeout = max(3600, getOption("timeout", 60)))
options(download.file.method = "libcurl")

option_list <- list(
  make_option(c("-g", "--gse"), type="character", default=NULL, 
              help="GEO Series ID (e.g., GSE86831)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="Output directory", metavar="character"),
  make_option(c("--platform"), type="character", default="",
              help="Optional platform ID (e.g., GPL21145) to force selection when multiple platforms exist",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Timestamped messages for easier runtime tracking
ts_message <- function(..., domain = NULL, appendLF = TRUE) {
  base::message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " ")),
                domain = domain, appendLF = appendLF)
}
message <- ts_message

if (is.null(opt$gse)){
  print_help(opt_parser)
  stop("GSE ID must be supplied", call.=FALSE)
}

gse_id <- opt$gse
out_dir <- opt$out
idat_dir <- file.path(out_dir, "idat")

if (!dir.exists(idat_dir)) {
  dir.create(idat_dir, recursive = TRUE)
}

message(paste("Fetching metadata for:", gse_id))

retry_attempts <- suppressWarnings(as.integer(Sys.getenv("ILLUMETA_DOWNLOAD_RETRIES", "3")))
if (!is.finite(retry_attempts) || retry_attempts < 1) retry_attempts <- 3
retry_wait <- suppressWarnings(as.numeric(Sys.getenv("ILLUMETA_DOWNLOAD_WAIT", "3")))
if (!is.finite(retry_wait) || retry_wait < 0) retry_wait <- 3
retry_backoff <- suppressWarnings(as.numeric(Sys.getenv("ILLUMETA_DOWNLOAD_BACKOFF", "1")))
if (!is.finite(retry_backoff) || retry_backoff < 1) retry_backoff <- 1

retry_run <- function(fn, label = "request", attempts = retry_attempts, wait = retry_wait, backoff = retry_backoff) {
  last_err <- NULL
  last_msg <- ""
  for (i in seq_len(attempts)) {
    res <- tryCatch(fn(), error = function(e) { last_err <<- e; NULL })
    if (!is.null(res)) {
      return(res)
    }
    if (!is.null(last_err) && !is.null(last_err$message) && nzchar(last_err$message)) {
      last_msg <- last_err$message
    } else if (last_msg == "") {
      last_msg <- "returned NULL"
    }
    message(sprintf("[%s] attempt %d/%d failed: %s", label, i, attempts, last_msg))
    if (i < attempts) Sys.sleep(wait * (backoff ^ (i - 1)))
  }
  stop(sprintf("%s failed after %d attempts: %s", label, attempts, last_msg))
}

# Fetch GEO series (possibly multi-platform)
gse_list <- retry_run(function() getGEO(gse_id, GSEMatrix = TRUE), label = "getGEO")
platform_override <- toupper(trimws(opt$platform))
if (length(gse_list) > 1) {
  message(sprintf("Multiple platforms detected (%d). Selecting the one with IDAT evidence...", length(gse_list)))
  chosen_idx <- NULL
  if (nzchar(platform_override)) {
    match_idx <- which(vapply(gse_list, function(eset) {
      ann <- tryCatch(annotation(eset), error = function(e) "")
      toupper(trimws(as.character(ann))) == platform_override
    }, logical(1)))
    if (length(match_idx) > 0) {
      chosen_idx <- match_idx[1]
      message(sprintf("  Platform override requested: %s (selected index %d).", platform_override, chosen_idx))
    } else {
      message(sprintf("WARNING: platform override %s not found; falling back to IDAT evidence.", platform_override))
    }
  }
  if (is.null(chosen_idx)) {
    idat_counts <- vapply(gse_list, function(eset) {
      pdat <- pData(eset)
      supp_cols <- grep("supplementary_file", colnames(pdat), value = TRUE)
      if (length(supp_cols) == 0) return(0L)
      vals <- as.character(unlist(pdat[, supp_cols, drop = FALSE]))
      sum(grepl("\\.idat(\\.gz)?$", vals, ignore.case = TRUE))
    }, integer(1))
    best_idx <- if (all(idat_counts == 0)) 1 else which.max(idat_counts)
    if (all(idat_counts == 0)) {
      message("  No platform showed explicit IDAT links; defaulting to the first.")
    } else {
      message(sprintf("  Choosing platform %d (IDAT hits: %d)", best_idx, idat_counts[best_idx]))
    }
    chosen_idx <- best_idx
  }
  gse <- gse_list[[chosen_idx]]
} else {
  gse <- gse_list[[1]]
  if (nzchar(platform_override)) {
    ann <- tryCatch(annotation(gse), error = function(e) "")
    if (!identical(toupper(trimws(as.character(ann))), platform_override)) {
      message(sprintf("WARNING: platform override %s does not match single-platform annotation (%s).",
                      platform_override, as.character(ann)))
    }
  }
}

meta <- pData(gse)

# Check for IDAT files in supplementary info
# Usually typically listed in 'supplementary_file' columns
supp_cols <- grep("supplementary_file", colnames(meta), value = TRUE)
has_idat <- FALSE

if (length(supp_cols) > 0) {
    # Pattern to detect IDAT files in GEO supplementary metadata
    # Accept both compressed and uncompressed IDATs
    has_idat <- any(apply(meta[, supp_cols, drop=FALSE], 1, function(x) any(grepl("\\.idat(\\.gz)?$", x, ignore.case=TRUE))))
}

if (!has_idat) {
  message("WARNING: No IDAT files found in metadata (supplementary_file).")
  message("This pipeline requires raw IDAT files.")
}

message("Checking for existing IDAT files...")
existing_idats <- list.files(idat_dir, pattern = "idat(\\.gz)?$", ignore.case = TRUE)

if (length(existing_idats) > 0) {
    message(paste("Found", length(existing_idats), "IDAT files already in", idat_dir))
    message("Skipping download step.")
} else {
    message("Attempting to download raw files (IDATs)...")

    geo_group_from_gse <- function(gse_id) {
        num <- suppressWarnings(as.integer(gsub("[^0-9]", "", gse_id)))
        if (is.na(num)) return(NULL)
        grp_num <- floor(num / 1000)
        unique(c(
            paste0("GSE", grp_num, "nnn"),
            sprintf("GSE%03dnnn", grp_num)
        ))
    }

    manual_download_suppl <- function(gse_id, out_dir, pattern = "[Rr][Aa][Ww]\\.tar(\\.gz)?$") {
        grps <- geo_group_from_gse(gse_id)
        if (is.null(grps) || length(grps) == 0) return(NULL)
        dl_dir <- file.path(out_dir, gse_id)
        if (!dir.exists(dl_dir)) dir.create(dl_dir, recursive = TRUE)
        for (grp in grps) {
            base_url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl", grp, gse_id)
            filelist_url <- paste0(base_url, "/filelist.txt")
            tmp <- tempfile(fileext = ".txt")
            files <- character(0)
            ok <- tryCatch({
                utils::download.file(filelist_url, tmp, quiet = TRUE, mode = "wb")
                TRUE
            }, error = function(e) FALSE)
            if (ok) {
                lines <- readLines(tmp, warn = FALSE)
                files <- sub("\\s+.*$", "", lines)
                files <- files[grepl(pattern, files)]
            }
            if (length(files) == 0) {
                files <- c(paste0(gse_id, "_RAW.tar"), paste0(gse_id, "_RAW.tar.gz"))
            }
            downloaded <- character(0)
            for (fname in unique(files)) {
                url <- paste0(base_url, "/", fname)
                dest <- file.path(dl_dir, fname)
                ok <- tryCatch({
                    utils::download.file(url, dest, quiet = TRUE, mode = "wb")
                    file.exists(dest) && file.info(dest)$size > 0
                }, error = function(e) FALSE)
                if (ok) {
                    downloaded <- c(downloaded, dest)
                } else if (file.exists(dest)) {
                    unlink(dest)
                }
            }
            if (length(downloaded) > 0) {
                return(downloaded)
            }
        }
        NULL
    }
    
    tryCatch({
        # Strategy 1: Try to download individual IDAT files
        files_info <- tryCatch(
            retry_run(
                function() getGEOSuppFiles(gse_id, baseDir = out_dir, makeDirectory = TRUE, filter_regex = "idat(\\\\.gz)?$"),
                label = "getGEOSuppFiles(idat[.gz])"
            ),
            error = function(e) {
                message("IDAT download attempt failed: ", e$message)
                NULL
            }
        )
        
        downloaded_dir <- file.path(out_dir, gse_id)
        
        # Check if we got any IDATs directly
        files <- list.files(downloaded_dir, pattern = "idat(\\.gz)?$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
        
        if (length(files) == 0) {
            message("No individual IDAT files found. Checking for RAW.tar bundle...")
            
            # Strategy 2: Download RAW.tar
            # Note: filter_regex might need to be broad or specific. usually "RAW.tar" matches "GSEXXXX_RAW.tar"
            files_info_tar <- tryCatch(
                retry_run(
                    function() getGEOSuppFiles(gse_id, baseDir = out_dir, makeDirectory = TRUE, filter_regex = "[Rr][Aa][Ww]\\\\.tar(\\\\.gz)?$"),
                    label = "getGEOSuppFiles(RAW.tar)"
                ),
                error = function(e) {
                    message("RAW tar download attempt failed: ", e$message)
                    NULL
                }
            )
            
            tar_files <- list.files(downloaded_dir, pattern = "RAW\\.tar(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
            if (length(tar_files) == 0) {
                manual_tar <- manual_download_suppl(gse_id, out_dir)
                if (!is.null(manual_tar)) {
                    tar_files <- manual_tar
                    message("Downloaded RAW tar bundle via direct HTTPS fallback.")
                }
            }
            
            if (length(tar_files) > 0) {
                message(paste("Found RAW tar bundle:", tar_files[1]))
                message("Extracting TAR file...")
                # Untar to downloaded_dir
                untar(tar_files[1], exdir = downloaded_dir)
                
                # Now check for IDATs again (they might be compressed or not)
                files <- list.files(downloaded_dir, pattern = "idat", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
            }
        }
        
        if (length(files) > 0) {
            # Copy files to idat_dir (safer than rename across filesystems)
            # We strictly look for .idat or .idat.gz
            valid_idats <- files[grepl("idat(\\.gz)?$", files, ignore.case=TRUE)]
            
            if (length(valid_idats) > 0) {
                target_paths <- file.path(idat_dir, basename(valid_idats))
                copied <- file.copy(valid_idats, target_paths, overwrite = TRUE)
                if (!all(copied)) {
                    bad <- basename(valid_idats)[!copied]
                    stop("Failed to copy IDAT files: ", paste(bad, collapse = ", "))
                }
                message(paste("Downloaded and copied", length(valid_idats), "IDAT files to", idat_dir))
                
                # Clean up GSE folder and tar files if desired
                unlink(downloaded_dir, recursive = TRUE)
            } else {
                 stop("Files extracted but no .idat files found.")
            }
        } else {
             message("No IDAT files were downloaded via getGEOSuppFiles (neither individual nor TAR).")
             stop("No IDAT files found.")
        }
        
    }, error = function(e) {
        message("Error during download or no IDATs available: ", e$message)
        quit(status=1)
    })
}

# Prepare configure.tsv
message("Creating configure.tsv...")

# Clean up metadata
keep_cols <- colnames(meta)[!grepl("contact|data_row_count|channel_count|status|submission_date|last_update_date", colnames(meta))]
simple_meta <- meta[, keep_cols, drop=FALSE]

# Automatic cleaning of "Key: Value" format in columns
# Many GEO characteristic columns are like "disease: preeclampsia"
for (col in colnames(simple_meta)) {
    if (is.character(simple_meta[[col]]) || is.factor(simple_meta[[col]])) {
        # Check if the column values look like "Key: Value"
        # We check the first few non-NA values
        vals <- na.omit(as.character(simple_meta[[col]]))
        if (length(vals) > 0) {
             # Simple heuristic: if all/most contain ": ", we strip the prefix
             if (all(grepl(": ", vals))) {
                 # Remove prefix (everything up to first ": ")
                 simple_meta[[col]] <- sub("^[^:]+: ", "", as.character(simple_meta[[col]]))
             }
        }
    }
}

# Extract Sentrix_ID and Sentrix_Position from supplementary_file paths
# We assume 'supplementary_file' column contains the IDAT path
supp_col <- grep("supplementary_file", colnames(simple_meta), value = TRUE)[1]

if (!is.na(supp_col)) {
    # Typical pattern: ..._SentrixID_SentrixPosition_Grn.idat.gz
    # or ..._SentrixID_SentrixPosition.idat.gz
    # Regex to capture (SentrixID)_(RxxCxx)
    # We look for 10+ digits followed by _RxxCxx
    
    paths <- simple_meta[[supp_col]]
    
    # Extract Position (RxxCxx)
    # Pattern: _(R[0-9]{2}C[0-9]{2})_?
    pos_matches <- str_extract(paths, "_R[0-9]{2}C[0-9]{2}")
    simple_meta$Sentrix_Position <- sub("^_", "", pos_matches)
    
    # Extract Sentrix ID (digits before Position)
    # Pattern: ([0-9]{10,})_R[0-9]{2}C[0-9]{2}
    # We capture the digits group
    id_matches <- str_match(paths, "([0-9]{10,})_R[0-9]{2}C[0-9]{2}")
    if (ncol(id_matches) >= 2) {
        simple_meta$Sentrix_ID <- id_matches[, 2]
    } else {
        simple_meta$Sentrix_ID <- NA
    }
}

# Add primary_group
simple_meta <- data.frame(primary_group = "", simple_meta, check.names = FALSE)

# Drop duplicate columns with identical values (keep the first)
encode_col <- function(x) paste0(ifelse(is.na(x), "<NA>", as.character(x)), collapse = "|")
col_enc <- vapply(simple_meta, encode_col, character(1))
dup_cols <- names(simple_meta)[duplicated(col_enc)]
if (length(dup_cols) > 0) {
    keep_cols <- setdiff(names(simple_meta), dup_cols)
    simple_meta <- simple_meta[, keep_cols, drop = FALSE]
    message("Dropped duplicate columns (identical values): ", paste(dup_cols, collapse = ", "))
}

# Save original (full) metadata snapshot
config_orig_path <- file.path(out_dir, "configure_original.tsv")
write.table(simple_meta, config_orig_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Original metadata saved to:", config_orig_path))

# Retain all columns except degenerate ones; always keep geo_accession for ID mapping
geo_cols <- grep("geo_accession", colnames(simple_meta), value = TRUE, ignore.case = TRUE)
dropped <- data.frame(Variable=character(0), Reason=character(0))
keep_cols <- unique(c("primary_group", geo_cols))
for (nm in colnames(simple_meta)) {
    if (nm %in% keep_cols) next
    uniq <- length(unique(simple_meta[[nm]]))
    if (uniq < 2) {
        dropped <- rbind(dropped, data.frame(Variable=nm, Reason="single_level"))
        next
    }
    if (uniq >= nrow(simple_meta)) {
        dropped <- rbind(dropped, data.frame(Variable=nm, Reason="all_unique"))
        next
    }
    keep_cols <- c(keep_cols, nm)
}
filtered_meta <- simple_meta[, keep_cols, drop=FALSE]

if (nrow(dropped) > 0) {
    message("Dropped columns (degenerate): ", paste(dropped$Variable, collapse=", "))
}

config_path <- file.path(out_dir, "configure.tsv")
write.table(filtered_meta, config_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Metadata saved.")
message("IMPORTANT: Fill in 'primary_group' in configure.tsv before analysis (or use illumeta.py --auto-group to populate it).")
