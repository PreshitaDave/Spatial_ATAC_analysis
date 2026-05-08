suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

project_root <- Sys.getenv("PROJECT_ROOT", "/projectnb/paxlab/presh/projects/spatial_atac")
barcode_dir <- file.path(project_root, "Data", "tissue_barcodes")
out_root <- file.path(project_root, "Data", "archr_tissue_no_edge")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

threads <- as.integer(Sys.getenv("NSLOTS", "4"))
addArchRThreads(threads = threads)

datasets <- c("deepseq", "lowseq")
tissues <- c("488B", "489")

archr_sources <- list(
  deepseq = file.path(project_root, "Data", "Save-ArchR-Project.rds"),
  lowseq = file.path(project_root, "Data", "lowseq_saveArchR", "Save-ArchR-Project.rds")
)

read_bc <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  x <- fread(path, header = FALSE, data.table = FALSE)
  bc <- as.character(x[[1]])
  bc <- bc[nzchar(bc)]
  unique(sub("-1$", "", bc))
}

norm_archr_cell <- function(cell_name) {
  sub("-1$", "", sub("#.*$", "", as.character(cell_name)))
}

summary_rows <- list()

for (dataset in datasets) {
  src <- archr_sources[[dataset]]
  if (!file.exists(src)) {
    stop(sprintf("Missing ArchR source project: %s", src), call. = FALSE)
  }

  log_msg("step", sprintf("Loading %s ArchR project from %s", dataset, src))
  proj <- readRDS(src)
  all_cells <- getCellNames(proj)
  all_cells_norm <- norm_archr_cell(all_cells)
  log_msg("step", sprintf("%s project cells: %d", dataset, length(all_cells)))

  for (tissue in tissues) {
    all_file <- file.path(barcode_dir, sprintf("%s_%s.barcodes.tsv", dataset, tissue))
    edge_file <- file.path(barcode_dir, sprintf("%s_%s.edge_effect.barcodes.tsv", dataset, tissue))
    no_edge_file <- file.path(barcode_dir, sprintf("%s_%s.no_edge_effect.barcodes.tsv", dataset, tissue))

    all_bc <- read_bc(all_file)
    edge_bc <- read_bc(edge_file)
    no_edge_bc <- read_bc(no_edge_file)

    if (!length(no_edge_bc) && length(all_bc)) {
      no_edge_bc <- setdiff(all_bc, edge_bc)
      writeLines(no_edge_bc, no_edge_file)
      log_msg("step", sprintf("Wrote derived no-edge barcodes: %s (%d)", no_edge_file, length(no_edge_bc)))
    }

    matched <- all_cells[all_cells_norm %in% no_edge_bc]
    out_dir <- file.path(out_root, sprintf("%s_%s_no_edge", dataset, tissue))

    if (length(matched) == 0L) {
      log_msg("warn", sprintf("No matching ArchR cells for %s %s (no-edge barcodes=%d)", dataset, tissue, length(no_edge_bc)))
      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        dataset = dataset,
        tissue = tissue,
        archr_total_cells = length(all_cells),
        no_edge_barcodes = length(no_edge_bc),
        matched_cells = 0L,
        saved_project = FALSE,
        output_dir = out_dir,
        output_rds = file.path(out_dir, "Save-ArchR-Project.rds")
      )
      next
    }

    log_msg("step", sprintf("Saving %s %s no-edge ArchR project (%d cells)", dataset, tissue, length(matched)))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    suppressMessages({
      subsetArchRProject(
        ArchRProj = proj,
        cells = matched,
        outputDirectory = out_dir,
        force = TRUE,
        threads = threads
      )
    })

    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      dataset = dataset,
      tissue = tissue,
      archr_total_cells = length(all_cells),
      no_edge_barcodes = length(no_edge_bc),
      matched_cells = length(matched),
      saved_project = file.exists(file.path(out_dir, "Save-ArchR-Project.rds")),
      output_dir = out_dir,
      output_rds = file.path(out_dir, "Save-ArchR-Project.rds")
    )
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
summary_file <- file.path(out_root, "archr_tissue_no_edge_summary.tsv")
fwrite(summary_dt, summary_file, sep = "\t")
log_msg("done", sprintf("Wrote summary: %s", summary_file))
