# ==============================================================================
# Export ArchR GeneScoreMatrix (subset to genes shared with the Xenium panel)
# for the gene-activity vs Xenium-expression loss evaluation notebook.
# ==============================================================================

suppressMessages({
  library(ArchR)
  library(Matrix)
  library(data.table)
})

# --- Parameters ---
archr_project_path <- "/projectnb/paxlab/presh/projects/spatial_atac/analysis/binsize_comparison/archr_projects/deepseq_488B_5000bp_binarizeFALSE_v2"
xenium_genes_path <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/Xenium_488B/giotto_output/xenium_genes.csv"
output_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/mosaicfield_outputs/gene_loss_inputs_v2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ==============================================================================
# 1. Load ArchR project, extract GeneScoreMatrix
# ==============================================================================

cat("Loading ArchR project...\n")
proj <- loadArchRProject(archr_project_path, showLogo = FALSE)
cat(sprintf("  %d cells\n", ncol(proj)))
cat(sprintf("  Available matrices: %s\n", paste(getAvailableMatrices(proj), collapse = ", ")))

cat("\nExtracting GeneScoreMatrix...\n")
gsm <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
archr_genes <- rowData(gsm)$name
cat(sprintf("  GeneScoreMatrix: %d genes x %d cells\n", nrow(gsm), ncol(gsm)))

# ==============================================================================
# 2. Determine shared gene space (ArchR genes ∩ Xenium panel genes)
# ==============================================================================

xenium_genes <- fread(xenium_genes_path)$gene
shared_genes <- intersect(archr_genes, xenium_genes)

cat(sprintf("\nGene overlap:\n"))
cat(sprintf("  ArchR GeneScoreMatrix genes: %d\n", length(archr_genes)))
cat(sprintf("  Xenium panel genes:          %d\n", length(xenium_genes)))
cat(sprintf("  Shared genes (used for loss):%d\n", length(shared_genes)))

if (length(shared_genes) < 10) {
  stop("Fewer than 10 shared genes found - check gene symbol formatting.")
}

# ==============================================================================
# 3. Subset GeneScoreMatrix to shared genes, export
# ==============================================================================

gene_idx <- match(shared_genes, archr_genes)
gs_mat <- assay(gsm)[gene_idx, , drop = FALSE]
rownames(gs_mat) <- shared_genes
cell_ids <- colnames(gsm)
barcodes <- sub(".*#", "", cell_ids)

if (!inherits(gs_mat, "sparseMatrix")) {
  gs_mat <- Matrix(gs_mat, sparse = TRUE)
}

cat("\nExporting files...\n")

# genes x cells, MTX (raw ArchR gene scores; log1p applied downstream in Python)
writeMM(gs_mat, file.path(output_dir, "archr_genescore_shared.mtx"))

fwrite(data.table(gene = shared_genes), file.path(output_dir, "archr_gene_names.csv"))
fwrite(data.table(cell_id = cell_ids), file.path(output_dir, "archr_cell_names.csv"))
fwrite(data.table(barcode = barcodes), file.path(output_dir, "archr_barcodes.csv"))
fwrite(data.table(gene = shared_genes), file.path(output_dir, "shared_genes.csv"))

cat(sprintf("\n%s\nExport complete: %s\n%s\n", strrep("=", 80), output_dir, strrep("=", 80)))
for (f in list.files(output_dir)) {
  cat(sprintf("  %s (%.0f KB)\n", f, file.size(file.path(output_dir, f)) / 1024))
}
