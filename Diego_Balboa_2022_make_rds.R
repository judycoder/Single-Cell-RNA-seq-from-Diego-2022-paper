#Sc-RNA seq test with one sample file from Diego Balboa 2022

library(Seurat)

data_dir <- "/Users/apple/Desktop/GSE167880_RAW"
files <- list.files(data_dir, full.names = TRUE)

mtx_files  <- files[grepl("matrix\\.mtx(\\.gz)?$", files)]
bar_files  <- files[grepl("barcodes\\.tsv(\\.gz)?$", files)]
feature_file <- "/Users/apple/Desktop/GSE167880_RAW/GSE167880_features.tsv"

prefix <- function(x) sub("_(matrix\\.mtx|barcodes\\.tsv)(\\.gz)?$", "", basename(x))
samples <- sort(unique(prefix(mtx_files)))

out_dir <- file.path(data_dir, "rds_output")
dir.create(out_dir, showWarnings = FALSE)

for (s in samples) {
  mtx  <- mtx_files[prefix(mtx_files) == s]
  bar  <- bar_files[prefix(bar_files) == s]
  feat <- feature_file
  
  if (length(mtx) != 1 || length(bar) != 1) {
    message("Skipping ", s, " (missing/duplicate files)")
    next
  }
  
  mat <- ReadMtx(mtx = mtx, cells = bar, features = feat)
  obj <- CreateSeuratObject(counts = mat, project = s)
  saveRDS(obj, file.path(out_dir, paste0(s, ".rds")))
}




