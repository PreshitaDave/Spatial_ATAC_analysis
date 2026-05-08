suppressPackageStartupMessages({
  library(data.table)
})

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
variant_root <- file.path(project_root, "Data/variant_calling")
out_dir <- file.path(project_root, "Data/tissue_barcodes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_barcode_order_simple <- function(somatic_dir) {
  # Read barcodes from ALL chromosome files and aggregate (not just chr1)
  cell_files <- list.files(somatic_dir, pattern = "\\.cell_snv\\.cellID\\.filter\\.csv$", full.names = TRUE)
  if (!length(cell_files)) {
    stop(sprintf("No *.cell_snv.cellID.filter.csv files found in: %s", somatic_dir), call. = FALSE)
  }
  all_barcodes <- character(0)
  for (f in cell_files) {
    x <- fread(f)
    all_barcodes <- c(all_barcodes, as.character(x[[ncol(x)]]))
  }
  unique(gsub("-1$", "", all_barcodes))
}

build_map_8_to_16 <- function(dataset) {
  f16 <- file.path(variant_root, sprintf("%s_cell_data.csv", dataset))
  f8 <- file.path(variant_root, sprintf("%s_cell_data_8bp.csv", dataset))

  if (file.exists(f16) && file.exists(f8)) {
    d16 <- fread(f16)
    d8 <- fread(f8)
    m <- data.table(bc8 = as.character(d8$cell), bc16 = as.character(d16$cell))
    m <- unique(m)
    return(setNames(m$bc16, m$bc8))
  }

  if (file.exists(f16)) {
    d16 <- fread(f16)
    bc16 <- as.character(d16$cell)
    bc8 <- substr(bc16, 1, 8)
    m <- data.table(bc8 = bc8, bc16 = bc16)
    m <- m[, .SD[1], by = bc8]
    if (!nrow(m)) {
      return(setNames(character(0), character(0)))
    }
    return(setNames(m$bc16, m$bc8))
  }

  setNames(character(0), character(0))
}

map_to_16 <- function(barcodes, map8to16) {
  bc <- as.character(barcodes)
  is8 <- nchar(bc) == 8
  mapped <- bc
  mapped[is8] <- ifelse(
    bc[is8] %in% names(map8to16),
    unname(map8to16[bc[is8]]),
    bc[is8]
  )
  gsub("-1$", "", mapped)
}

spatial <- fread(file.path(project_root, "Data/tissue_positions_list.csv"), header = FALSE)
setnames(spatial, c("barcode", "in_tissue", "array_row", "array_col", "x_spatial", "y_spatial"))
spatial <- spatial[in_tissue == 1]
spatial[, tissue := ifelse(y_spatial > 4000, "488B", "489")]

for (dataset in c("deepseq", "lowseq")) {
  cat(sprintf("[%s] Building tissue barcode lists for %s\n", format(Sys.time(), "%F %T"), dataset))

  somatic_dir <- file.path(variant_root, dataset, "somatic")
  somatic_barcodes <- read_barcode_order_simple(somatic_dir)

  map8to16 <- build_map_8_to_16(dataset)
  spatial_dt <- copy(spatial)
  spatial_dt[, barcode16 := map_to_16(barcode, map8to16)]

  for (tissue_name in c("488B", "489")) {
    tissue_barcodes <- unique(spatial_dt[tissue == tissue_name, barcode16])
    tissue_barcodes <- tissue_barcodes[tissue_barcodes %in% somatic_barcodes]
    tissue_barcodes <- tissue_barcodes[nzchar(tissue_barcodes)]

    out_file <- file.path(out_dir, sprintf("%s_%s.barcodes.tsv", dataset, tissue_name))
    writeLines(tissue_barcodes, out_file)

    cat(sprintf("[%s] %s %s: wrote %d barcodes to %s\n",
      format(Sys.time(), "%F %T"), dataset, tissue_name, length(tissue_barcodes), out_file
    ))
  }
}

cat(sprintf("[%s] Tissue barcode list generation complete\n", format(Sys.time(), "%F %T")))
