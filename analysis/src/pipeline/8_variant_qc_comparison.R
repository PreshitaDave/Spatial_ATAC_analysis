#!/usr/bin/env Rscript
# 8_variant_qc_comparison.R
# Per-chromosome QC: deepseq vs lowseq comparison using Beagle-phased VCFs
# Pipeline: samtools mpileup (.gl.vcf.gz) -> Beagle genotyping (.gp.vcf.gz) -> Beagle phasing (.phased.vcf.gz)
# This script uses .phased.vcf.gz (properly called + phased genotypes from 1000G panel)
# and .gp.vcf.gz (for GP-based quality metrics on chr19-22)
# Produces 9 PDFs in the comparison/phased/ folder

cat("=== Variant QC Comparison Script (Phased VCFs) ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

n_cores <- 6
cat("Using", n_cores, "cores for parallel loading.\n\n")

library(vcfR)
library(parallel)

# ---- Paths ----
base_dir  <- '/projectnb/paxlab/presh/projects/spatial_atac/Data'
deep_germ <- file.path(base_dir, 'variant_calling/deepseq/germline')
low_germ  <- file.path(base_dir, 'variant_calling/lowseq/germline')
ref_fai   <- file.path(base_dir, 'hg38_resources/Homo_sapiens_assembly38.fasta.fai')
comp_dir  <- file.path(base_dir, '../analysis/plots/variant_qc/comparison/phased')

chromosomes <- paste0('chr', 1:22)
# GP files (Beagle genotype posteriors) now available for all chromosomes
gp_chrs <- paste0('chr', 1:22)

dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Chromosome lengths from FASTA index ----
cat("[1] Reading chromosome lengths from FASTA index...\n")
fai <- read.table(ref_fai, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
chr_lengths <- setNames(fai$V2, fai$V1)

# ---- Helper functions ----

# Ts/Tv from REF/ALT columns (SNP sites only)
count_ts_tv <- function(vcf) {
  if (nrow(vcf@fix) == 0) return(c(ts = 0, tv = 0))
  ref <- vcf@fix[, "REF"]
  alt <- vcf@fix[, "ALT"]
  snp <- (nchar(ref) == 1 & nchar(alt) == 1)
  pair <- paste0(ref[snp], alt[snp])
  n_ts <- sum(pair %in% c("AG", "GA", "CT", "TC"))
  n_tv <- sum(!(pair %in% c("AG", "GA", "CT", "TC")))
  c(ts = n_ts, tv = n_tv)
}

# Count het / hom-alt / hom-ref from phased GT vector
count_gt_types <- function(gt) {
  het    <- sum(gt %in% c("0|1", "1|0"), na.rm = TRUE)
  hom    <- sum(gt == "1|1", na.rm = TRUE)
  homref <- sum(gt == "0|0", na.rm = TRUE)
  c(het = het, hom = hom, homref = homref,
    total = sum(!is.na(gt)))
}

# Which GT values carry an ALT allele
is_alt_gt <- function(gt) gt %in% c("0|1", "1|0", "1|1")

# ---- Load all chromosomes in parallel ----
cat("[2] Loading phased VCF data for all chromosomes in parallel (", n_cores, "cores)...\n")

load_one_chr <- function(chr) {
  deep_path <- file.path(deep_germ, paste0(chr, '.phased.vcf.gz'))
  low_path  <- file.path(low_germ,  paste0(chr, '.phased.vcf.gz'))
  if (!file.exists(deep_path) || !file.exists(low_path)) return(NULL)

  deep_vcf <- read.vcfR(deep_path, verbose = FALSE)
  low_vcf  <- read.vcfR(low_path,  verbose = FALSE)

  deep_gt <- extract.gt(deep_vcf, element = "GT")[, 1]
  low_gt  <- extract.gt(low_vcf,  element = "GT")[, 1]

  deep_gt_counts <- count_gt_types(deep_gt)
  low_gt_counts  <- count_gt_types(low_gt)

  deep_alt_idx <- is_alt_gt(deep_gt)
  low_alt_idx  <- is_alt_gt(low_gt)

  # Ts/Tv on ALT-carrying sites only
  deep_tstv <- count_ts_tv(deep_vcf[deep_alt_idx, ])
  low_tstv  <- count_ts_tv(low_vcf[low_alt_idx, ])

  # Overlap: ALT positions shared between datasets
  deep_pos_alt <- paste0(deep_vcf@fix[deep_alt_idx, "CHROM"], ":",
                         deep_vcf@fix[deep_alt_idx, "POS"])
  low_pos_alt  <- paste0(low_vcf@fix[low_alt_idx,  "CHROM"], ":",
                         low_vcf@fix[low_alt_idx,  "POS"])

  n_shared    <- length(intersect(deep_pos_alt, low_pos_alt))
  n_deep_only <- sum(!deep_pos_alt %in% low_pos_alt)
  n_low_only  <- sum(!low_pos_alt  %in% deep_pos_alt)

  list(
    deep_vcf = deep_vcf, low_vcf = low_vcf,
    deep_gt = deep_gt, low_gt = low_gt,
    deep_gt_counts = deep_gt_counts, low_gt_counts = low_gt_counts,
    deep_alt_idx = deep_alt_idx, low_alt_idx = low_alt_idx,
    deep_tstv = deep_tstv, low_tstv = low_tstv,
    n_shared = n_shared, n_deep_only = n_deep_only, n_low_only = n_low_only,
    chr_len = chr_lengths[chr]
  )
}

raw_results <- mclapply(chromosomes, load_one_chr, mc.cores = n_cores)
names(raw_results) <- chromosomes
results <- Filter(Negate(is.null), raw_results)
cat(sprintf("  Loaded %d chromosomes.\n\n", length(results)))

# ---- Aggregate per-chromosome stats ----
all_stats <- do.call(rbind, lapply(names(results), function(chr) {
  r   <- results[[chr]]
  dgc <- r$deep_gt_counts
  lgc <- r$low_gt_counts
  data.frame(
    chr         = chr,
    deep_total  = dgc["total"],  low_total  = lgc["total"],
    deep_alt    = dgc["het"] + dgc["hom"],
    low_alt     = lgc["het"] + lgc["hom"],
    deep_het    = dgc["het"],    low_het    = lgc["het"],
    deep_hom    = dgc["hom"],    low_hom    = lgc["hom"],
    deep_homref = dgc["homref"], low_homref = lgc["homref"],
    deep_ts = r$deep_tstv["ts"], deep_tv = r$deep_tstv["tv"],
    low_ts  = r$low_tstv["ts"],  low_tv  = r$low_tstv["tv"],
    n_shared    = r$n_shared,
    n_deep_only = r$n_deep_only,
    n_low_only  = r$n_low_only,
    row.names = NULL
  )
}))

all_stats$deep_tstv    <- all_stats$deep_ts / pmax(all_stats$deep_tv, 1)
all_stats$low_tstv     <- all_stats$low_ts  / pmax(all_stats$low_tv,  1)
all_stats$deep_het_hom <- all_stats$deep_het / pmax(all_stats$deep_hom, 1)
all_stats$low_het_hom  <- all_stats$low_het  / pmax(all_stats$low_hom,  1)

chr_nums <- sub("chr", "", all_stats$chr)

# ================================================================
# PLOTS
# ================================================================

# ---- PLOT 1: ALT variant counts + Ts/Tv per chromosome ----
cat("[PLOT 1/9] tstv_variant_count.pdf\n")
pdf(file.path(comp_dir, 'tstv_variant_count.pdf'), width = 16, height = 8)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))

bm <- rbind(all_stats$deep_alt, all_stats$low_alt)
colnames(bm) <- chr_nums
barplot(bm, beside = TRUE, col = c("steelblue", "darkorange"),
        main = "ALT-Carrying Variant Count per Chromosome (Beagle-genotyped)",
        xlab = "Chromosome", ylab = "Number of Variants", las = 1,
        legend.text = c("Deepseq", "Lowseq"),
        args.legend = list(x = "topright", bty = "n"))

tm <- rbind(all_stats$deep_tstv, all_stats$low_tstv)
colnames(tm) <- chr_nums
ylim_ts <- c(0, max(c(2.5, max(tm, na.rm = TRUE) * 1.1)))
barplot(tm, beside = TRUE, col = c("steelblue", "darkorange"),
        main = "Ts/Tv Ratio per Chromosome (~2.0 expected for WGS SNPs)",
        xlab = "Chromosome", ylab = "Ts/Tv", las = 1, ylim = ylim_ts,
        legend.text = c("Deepseq", "Lowseq"),
        args.legend = list(x = "topright", bty = "n"))
abline(h = 2.0, col = "red", lty = 2, lwd = 1.5)
text(1, 2.08, "Expected ~2.0", col = "red", adj = 0, cex = 0.8)
dev.off()

# ---- PLOT 2: Het vs HomAlt counts per chromosome ----
cat("[PLOT 2/9] het_hom_comparison.pdf\n")
pdf(file.path(comp_dir, 'het_hom_comparison.pdf'), width = 16, height = 8)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))

hm_d <- rbind(all_stats$deep_het, all_stats$deep_hom)
colnames(hm_d) <- chr_nums
barplot(hm_d, beside = TRUE, col = c("steelblue", "lightblue"),
        main = "Deepseq: Heterozygous vs Homozygous ALT per Chromosome",
        xlab = "Chromosome", ylab = "Count", las = 1,
        legend.text = c("Het (0|1)", "HomAlt (1|1)"),
        args.legend = list(x = "topright", bty = "n"))

hm_l <- rbind(all_stats$low_het, all_stats$low_hom)
colnames(hm_l) <- chr_nums
barplot(hm_l, beside = TRUE, col = c("darkorange", "moccasin"),
        main = "Lowseq: Heterozygous vs Homozygous ALT per Chromosome",
        xlab = "Chromosome", ylab = "Count", las = 1,
        legend.text = c("Het (0|1)", "HomAlt (1|1)"),
        args.legend = list(x = "topright", bty = "n"))
dev.off()

# ---- PLOT 3: Het:Hom ratio per chromosome ----
cat("[PLOT 3/9] het_hom_ratio.pdf\n")
pdf(file.path(comp_dir, 'het_hom_ratio.pdf'), width = 14, height = 6)
par(mar = c(5, 5, 4, 2))
rr <- rbind(all_stats$deep_het_hom, all_stats$low_het_hom)
colnames(rr) <- chr_nums
ylim_r <- c(0, max(c(3, max(rr, na.rm = TRUE) * 1.1)))
barplot(rr, beside = TRUE, col = c("steelblue", "darkorange"),
        main = "Het:Hom Ratio per Chromosome (~2.0 expected for diploid germline)",
        xlab = "Chromosome", ylab = "Het:Hom Ratio", las = 1, ylim = ylim_r,
        legend.text = c("Deepseq", "Lowseq"),
        args.legend = list(x = "topright", bty = "n"))
abline(h = 2.0, col = "red", lty = 2, lwd = 1.5)
text(1, 2.08, "Expected ~2.0", col = "red", adj = 0, cex = 0.8)
dev.off()

# ---- PLOT 4: ALT site overlap between deepseq and lowseq ----
cat("[PLOT 4/9] overlap_analysis.pdf\n")
pdf(file.path(comp_dir, 'overlap_analysis.pdf'), width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

overlap_mat <- rbind(all_stats$n_deep_only, all_stats$n_shared, all_stats$n_low_only)
colnames(overlap_mat) <- chr_nums
barplot(overlap_mat, beside = FALSE,
        col = c("steelblue", "mediumseagreen", "darkorange"),
        main = "ALT Site Overlap: Deepseq vs Lowseq per Chromosome",
        xlab = "Chromosome", ylab = "Number of Sites", las = 1,
        legend.text = c("Deep-only", "Shared", "Low-only"),
        args.legend = list(x = "top", bty = "n"))

frac_shared <- all_stats$n_shared /
  pmax(all_stats$n_shared + all_stats$n_deep_only + all_stats$n_low_only, 1)
barplot(frac_shared, names.arg = chr_nums, col = "mediumseagreen",
        main = "Fraction of ALT Sites Shared Between Deepseq and Lowseq",
        xlab = "Chromosome", ylab = "Fraction Shared", las = 1, ylim = c(0, 1))
abline(h = 0.8, col = "red", lty = 2, lwd = 1.5)
dev.off()

# ---- PLOT 5: Fraction of genotyped sites called HomRef ----
cat("[PLOT 5/9] homref_fraction.pdf\n")
pdf(file.path(comp_dir, 'homref_fraction.pdf'), width = 14, height = 6)
par(mar = c(5, 5, 4, 2))
frac_d <- all_stats$deep_homref / pmax(all_stats$deep_total, 1)
frac_l <- all_stats$low_homref  / pmax(all_stats$low_total,  1)
ff <- rbind(frac_d, frac_l)
colnames(ff) <- chr_nums
barplot(ff, beside = TRUE, col = c("steelblue", "darkorange"),
        main = "Fraction of 1000G Panel Sites Called HomRef (0|0) per Chromosome",
        xlab = "Chromosome", ylab = "Fraction HomRef", las = 1, ylim = c(0, 1),
        legend.text = c("Deepseq", "Lowseq"),
        args.legend = list(x = "topright", bty = "n"))
dev.off()

# ---- PLOT 6: Genotype distribution per chromosome (per page) ----
cat("[PLOT 6/9] genotype_distribution.pdf\n")
pdf(file.path(comp_dir, 'genotype_distribution.pdf'), width = 14, height = 6)
for (chr in names(results)) {
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  for (seq_type in c("deep", "low")) {
    gc <- if (seq_type == "deep") results[[chr]]$deep_gt_counts else results[[chr]]$low_gt_counts
    counts <- c(HomRef = gc["homref"], Het = gc["het"], HomAlt = gc["hom"])
    col    <- if (seq_type == "deep") c("lightblue", "steelblue", "navy") else
                                      c("moccasin",  "darkorange", "darkorange4")
    label  <- if (seq_type == "deep") "Deepseq" else "Lowseq"
    bp <- barplot(counts, col = col,
                  main = paste0(label, " - ", chr),
                  ylab = "Count", names.arg = c("Hom\nRef", "Het", "Hom\nAlt"))
    text(bp, counts / 2, labels = format(as.integer(counts), big.mark = ","), cex = 0.85)
  }
}
dev.off()

# ---- PLOT 7: ALT variant density along chromosome ----
cat("[PLOT 7/9] het_density_comparison.pdf\n")
pdf(file.path(comp_dir, 'het_density_comparison.pdf'), width = 14, height = 6)
for (chr in names(results)) {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  pos_d <- as.numeric(results[[chr]]$deep_vcf@fix[results[[chr]]$deep_alt_idx, "POS"])
  pos_l <- as.numeric(results[[chr]]$low_vcf@fix[results[[chr]]$low_alt_idx,  "POS"])
  if (length(pos_d) == 0 || length(pos_l) == 0) next
  chr_end    <- max(results[[chr]]$chr_len, max(pos_d), max(pos_l), na.rm = TRUE)
  breaks_seq <- seq(0, chr_end + 1e6, by = 1e6)
  h_d <- hist(pos_d, breaks = breaks_seq, plot = FALSE)
  h_l <- hist(pos_l, breaks = breaks_seq, plot = FALSE)
  ylim_max <- max(c(h_d$counts, h_l$counts), na.rm = TRUE) * 1.1

  plot(h_d$mids / 1e6, h_d$counts, type = "l", col = "steelblue", lwd = 1.5,
       main = paste0("Deepseq - ", chr),
       xlab = "Position (Mb)", ylab = "ALT Variants per 1 Mb", ylim = c(0, ylim_max))
  plot(h_l$mids / 1e6, h_l$counts, type = "l", col = "darkorange", lwd = 1.5,
       main = paste0("Lowseq - ", chr),
       xlab = "Position (Mb)", ylab = "ALT Variants per 1 Mb", ylim = c(0, ylim_max))
}
dev.off()

# ---- PLOTS 8–9: GP-based quality metrics (chr19-22 whole-chr GP files available) ----

extract_af <- function(vcf) {
  info <- vcf@fix[, "INFO"]
  af <- suppressWarnings(as.numeric(sub(".*\\bAF=([^;]+).*", "\\1", info)))
  af[!is.na(af)]
}

gp_to_gq <- function(gp_str) {
  # GQ = -10 * log10(1 - max_GP), capped at 99
  parts <- strsplit(gp_str, ",", fixed = TRUE)
  vapply(parts, function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (all(is.na(v))) return(NA_real_)
    best_p <- max(v, na.rm = TRUE)
    if (best_p >= 1.0) return(99.0)
    round(-10 * log10(1.0 - best_p), 1)
  }, numeric(1))
}

# Load GP data for all chromosomes in parallel
cat("[GP] Loading GP VCFs for all chromosomes in parallel...\n")
load_one_gp <- function(chr) {
  deep_gp_path <- file.path(deep_germ, paste0(chr, '.gp.vcf.gz'))
  low_gp_path  <- file.path(low_germ,  paste0(chr, '.gp.vcf.gz'))
  if (!file.exists(deep_gp_path) || !file.exists(low_gp_path)) return(NULL)
  deep_gp <- read.vcfR(deep_gp_path, verbose = FALSE)
  low_gp  <- read.vcfR(low_gp_path,  verbose = FALSE)
  list(deep_gp = deep_gp, low_gp = low_gp)
}
gp_results <- mclapply(gp_chrs, load_one_gp, mc.cores = n_cores)
names(gp_results) <- gp_chrs
gp_results <- Filter(Negate(is.null), gp_results)
cat(sprintf("  Loaded GP for %d chromosomes.\n\n", length(gp_results)))

cat("[PLOT 8/9] allele_freq_comparison.pdf\n")
pdf(file.path(comp_dir, 'allele_freq_comparison.pdf'), width = 14, height = 6)
for (chr in names(gp_results)) {
  cat("  AF:", chr, "\n")
  af_d <- extract_af(gp_results[[chr]]$deep_gp)
  af_l <- extract_af(gp_results[[chr]]$low_gp)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  if (length(af_d) > 0)
    hist(af_d, breaks = 50, col = "steelblue", border = "white",
         main = paste0("Deepseq - ", chr, " (n=", format(length(af_d), big.mark = ","), ")"),
         xlab = "Allele Frequency (from Beagle 1000G panel)", ylab = "Count")
  if (length(af_l) > 0)
    hist(af_l, breaks = 50, col = "darkorange", border = "white",
         main = paste0("Lowseq - ", chr, " (n=", format(length(af_l), big.mark = ","), ")"),
         xlab = "Allele Frequency (from Beagle 1000G panel)", ylab = "Count")
}
dev.off()

cat("[PLOT 9/9] gp_quality_comparison.pdf\n")
pdf(file.path(comp_dir, 'gp_quality_comparison.pdf'), width = 14, height = 6)
for (chr in names(gp_results)) {
  cat("  GQ:", chr, "\n")

  gp_d <- extract.gt(gp_results[[chr]]$deep_gp, element = "GP")[, 1]
  gp_l <- extract.gt(gp_results[[chr]]$low_gp,  element = "GP")[, 1]
  gq_d <- gp_to_gq(gp_d[!is.na(gp_d)]); gq_d <- gq_d[!is.na(gq_d)]
  gq_l <- gp_to_gq(gp_l[!is.na(gp_l)]); gq_l <- gq_l[!is.na(gq_l)]

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  if (length(gq_d) > 0) {
    hist(gq_d, breaks = 50, col = "steelblue", border = "white",
         main = paste0("Deepseq - ", chr, " (n=", format(length(gq_d), big.mark = ","), ")"),
         xlab = "Genotype Quality (GQ = -10*log10(1-max GP))", ylab = "Count")
    abline(v = median(gq_d), col = "red", lwd = 2, lty = 2)
    legend("topleft", paste0("Median GQ=", round(median(gq_d), 1)),
           col = "red", lty = 2, lwd = 2, bty = "n")
  }
  if (length(gq_l) > 0) {
    hist(gq_l, breaks = 50, col = "darkorange", border = "white",
         main = paste0("Lowseq - ", chr, " (n=", format(length(gq_l), big.mark = ","), ")"),
         xlab = "Genotype Quality (GQ = -10*log10(1-max GP))", ylab = "Count")
    abline(v = median(gq_l), col = "red", lwd = 2, lty = 2)
    legend("topleft", paste0("Median GQ=", round(median(gq_l), 1)),
           col = "red", lty = 2, lwd = 2, bty = "n")
  }
}
dev.off()

# ---- Summary ----
cat("\n[3] Summary tables:\n\n")
cat("=== Per-chromosome summary (phased genotypes) ===\n")
print(all_stats[, c("chr", "deep_alt", "low_alt", "deep_het", "low_het",
                    "deep_tstv", "low_tstv", "deep_het_hom", "low_het_hom",
                    "n_shared", "n_deep_only", "n_low_only")])

write.csv(all_stats, file.path(comp_dir, "summary_stats_phased.csv"), row.names = FALSE)

cat("\n=== ALL DONE ===\n")
cat("9 PDFs + 1 CSV saved to:", comp_dir, "\n")
cat("End time:", format(Sys.time()), "\n")
