#!/usr/bin/env Rscript

# to run this script, first load R module and then execute the command below, replacing the paths and parameters as needed:
# module load R && Rscript /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/numbat/postprocess_numbat_results.R --out_dir /projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/results/lowseq/atac_only --plot_dir /projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/results/lowseq/atac_only/plots --iteration 2 --gtf /projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/reference/var220kb.rds --dataset lowseq


suppressPackageStartupMessages({
  library(optparse)
  library(numbat)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--out_dir", type = "character"),
  make_option("--plot_dir", type = "character", default = NULL),
  make_option("--iteration", type = "integer", default = 2),
  make_option("--gtf", type = "character", default = NULL),
  make_option("--dataset", type = "character", default = "sample")
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$out_dir)) {
  stop("--out_dir is required")
}

plot_dir <- if (is.null(args$plot_dir)) file.path(args$out_dir, "plots") else args$plot_dir
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

gtf_obj <- numbat::gtf_hg38
if (!is.null(args$gtf) && nzchar(args$gtf)) {
  gr_bins <- readRDS(args$gtf)
  if (inherits(gr_bins, "GRanges")) {
    gr_bins$gene <- with(as.data.frame(gr_bins), paste0(seqnames, ":", start, "-", end))
    gr_bins$gene_length <- IRanges::width(gr_bins)
    gr_bins$gene_start <- IRanges::start(gr_bins)
    gr_bins$gene_end <- IRanges::end(gr_bins)
    gr_bins$CHROM <- gsub("chr", "", as.character(GenomicRanges::seqnames(gr_bins)))
    gtf_obj <- as.data.frame(gr_bins) %>%
      filter(CHROM != "X") %>%
      select(-strand) %>%
      mutate(across(c(gene_start, gene_end, gene_length), as.integer))
    rownames(gtf_obj) <- gtf_obj$gene
    gtf_obj <- gtf_obj[, colnames(numbat::gtf_hg38)]
  }
}

nb <- Numbat$new(out_dir = args$out_dir, i = args$iteration, gtf = gtf_obj)

save_plot <- function(filename, expr, width = 10, height = 5) {
  plot_obj <- tryCatch(expr(), error = function(e) {
    message(sprintf("Skipping %s: %s", filename, e$message))
    NULL
  })
  if (!is.null(plot_obj)) {
    ggsave(file.path(plot_dir, filename), plot_obj, width = width, height = height, dpi = 300)
  }
}

clone_levels <- sort(unique(nb$clone_post$clone_opt))
pal_clone <- setNames(grDevices::hcl.colors(length(clone_levels), palette = "Dynamic"), clone_levels)

save_plot(
  "phylogeny_heatmap.png",
  function() nb$plot_phylo_heatmap(clone_bar = TRUE, p_min = 0.9, pal_clone = pal_clone),
  width = 11,
  height = 5
)

save_plot(
  "consensus_segments.png",
  function() nb$plot_consensus(),
  width = 14,
  height = 2.5
)

save_plot(
  "single_cell_phylogeny.png",
  function() nb$plot_sc_tree(label_size = 3, branch_width = 0.5, tip_length = 0.5, pal_clone = pal_clone, tip = TRUE),
  width = 8,
  height = 4
)

save_plot(
  "mutation_history.png",
  function() nb$plot_mut_history(pal = pal_clone),
  width = 8,
  height = 4
)

save_plot(
  "clone_profiles.png",
  function() nb$plot_clone_profile(),
  width = 10,
  height = 5
)

save_plot(
  "expression_roll.png",
  function() nb$plot_exp_roll(k = min(3, max(2, length(clone_levels)))),
  width = 10,
  height = 5
)

if (!is.null(nb$bulk_clones) && nrow(nb$bulk_clones) > 0) {
  save_plot(
    "bulk_clone_profiles.png",
    function() {
      bulk_use <- nb$bulk_clones
      if ("n_cells" %in% names(bulk_use)) {
        bulk_use <- bulk_use %>% filter(n_cells > 0)
      }
      plot_bulks(bulk_use, min_LLR = 10, legend = TRUE)
    },
    width = 13,
    height = 6
  )
}

cut_values <- seq_len(min(4, max(1, length(clone_levels) - 1)))
cut_panels <- lapply(cut_values, function(k) {
  nb$cutree(n_cut = k)
  nb$plot_phylo_heatmap(clone_bar = TRUE, p_min = 0.9, pal_clone = pal_clone) +
    ggtitle(sprintf("%s n_cut=%d", args$dataset, k))
})

if (length(cut_panels) > 0) {
  suppressPackageStartupMessages(library(patchwork))
  ggsave(
    file.path(plot_dir, "phylogeny_cutree_grid.png"),
    patchwork::wrap_plots(cut_panels),
    width = 12,
    height = 8,
    dpi = 300
  )
}

writeLines(
  c(
    sprintf("dataset=%s", args$dataset),
    sprintf("out_dir=%s", normalizePath(args$out_dir)),
    sprintf("plot_dir=%s", normalizePath(plot_dir)),
    sprintf("iteration=%s", args$iteration),
    sprintf("plots_generated=%s", length(list.files(plot_dir)))
  ),
  con = file.path(plot_dir, "plot_manifest.txt")
)

invisible(TRUE)