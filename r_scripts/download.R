#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

option_list <- list(
  make_option(c("-g", "--gse"), type="character", default=NULL, 
              help="GEO Series ID (e.g., GSE86831)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

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

# Fetch GEO series
gse <- getGEO(gse_id, GSEMatrix = TRUE)
if (length(gse) > 1) {
  message("Warning: Multiple platforms found. Using the first one.")
}
gse <- gse[[1]]

meta <- pData(gse)

# Check for IDAT files in supplementary info
# Usually typically listed in 'supplementary_file' columns
supp_cols <- grep("supplementary_file", colnames(meta), value = TRUE)
has_idat <- FALSE

if (length(supp_cols) > 0) {
    # specific check for .idat.gz
    # FIX: We use quadruple backslash in the tool input so it writes double backslash to the file.
    # The file content on disk must look like: grepl("\\.idat\\.gz$", ...)
    has_idat <- any(apply(meta[, supp_cols, drop=FALSE], 1, function(x) any(grepl("\\.idat\\.gz$", x, ignore.case=TRUE))))
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
    
    tryCatch({
        # Strategy 1: Try to download individual IDAT files
        files_info <- getGEOSuppFiles(gse_id, baseDir = out_dir, makeDirectory = TRUE, filter_regex = "idat.gz")
        
        downloaded_dir <- file.path(out_dir, gse_id)
        
        # Check if we got any IDATs directly
        files <- list.files(downloaded_dir, pattern = "idat.gz", full.names = TRUE, recursive = TRUE)
        
        if (length(files) == 0) {
            message("No individual IDAT files found. Checking for RAW.tar bundle...")
            
            # Strategy 2: Download RAW.tar
            # Note: filter_regex might need to be broad or specific. usually "RAW.tar" matches "GSEXXXX_RAW.tar"
            files_info_tar <- getGEOSuppFiles(gse_id, baseDir = out_dir, makeDirectory = TRUE, filter_regex = "RAW.tar")
            
            tar_files <- list.files(downloaded_dir, pattern = "RAW.tar", full.names = TRUE)
            
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
            # Move files to idat_dir
            # We strictly look for .idat or .idat.gz
            valid_idats <- files[grepl("idat(\\.gz)?$", files, ignore.case=TRUE)]
            
            if (length(valid_idats) > 0) {
                file.rename(valid_idats, file.path(idat_dir, basename(valid_idats)))
                message(paste("Downloaded and moved", length(valid_idats), "IDAT files to", idat_dir))
                
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
message("IMPORTANT: You must edit 'configure.tsv' and fill in the 'primary_group' column before running analysis.")
