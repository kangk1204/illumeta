suppressPackageStartupMessages(library(GEOquery))
options(timeout=3600, download.file.method='libcurl')
out <- 'projects/GSE306227/metadata'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
gse_list <- getGEO('GSE306227', GSEMatrix=TRUE)
cat('platforms\t', length(gse_list), '\n', sep='')
for (i in seq_along(gse_list)) {
  gse <- gse_list[[i]]
  meta <- pData(gse)
  ann <- annotation(gse)
  fn <- file.path(out, sprintf('GSE306227_pData_%s.tsv', ann))
  write.table(meta, fn, sep='\t', quote=FALSE, row.names=FALSE)
  cat('platform\t', i, '\t', ann, '\t', nrow(meta), '\t', ncol(meta), '\t', fn, '\n', sep='')
  print(colnames(meta))
}
