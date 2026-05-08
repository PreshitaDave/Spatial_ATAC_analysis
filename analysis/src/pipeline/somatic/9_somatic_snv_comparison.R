#!/usr/bin/env Rscript
# 9_somatic_snv_comparison.R
# Compare deepseq vs lowseq putative somatic SNVs per chromosome
# Following Monopogen filtering workflow

cat("=== Somatic SNV Comparison: Deepseq vs Lowseq ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

library(parallel)
n_cores <- as.integer(Sys.getenv("NSLOTS", "6"))
cat("Using", n_cores, "cores\n")

# --- Paths ---
base_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling"
deep_somatic <- file.path(base_dir, "deepseq", "somatic")
low_somatic  <- file.path(base_dir, "lowseq", "somatic")
out_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/analysis/comparison/somatic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chromosomes <- paste0("chr", 1:22)

# --- Step 1: Load all putativeSNVs.csv ---
cat("\n--- Loading putativeSNVs for all chromosomes ---\n")

load_somatic <- function(dataset, chr) {
  f <- file.path(base_dir, dataset, "somatic", paste0(chr, ".putativeSNVs.csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$chr_name <- chr
  d$dataset <- dataset
  d$snv_id <- paste0(d$chr, ":", d$pos, ":", d$Ref_allele, ":", d$Alt_allele)
  d
}

deep_list <- mclapply(chromosomes, function(chr) load_somatic("deepseq", chr), mc.cores = n_cores)
low_list  <- mclapply(chromosomes, function(chr) load_somatic("lowseq", chr), mc.cores = n_cores)
names(deep_list) <- chromosomes
names(low_list)  <- chromosomes

cat("Loaded all chromosomes.\n")

# --- Step 2: Apply Monopogen filtering ---
# Following the Monopogen example:
# 1. Depth_ref > 5 AND Depth_alt > 5
# 2. BAF_alt < 0.5 (remove likely germline)
# 3. LDrefine_merged_score > 0.25 AND not NA

filter_somatic <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df <- df[df$Depth_ref > 5 & df$Depth_alt > 5, ]
  df <- df[df$BAF_alt < 0.5, ]
  df <- df[!is.na(df$LDrefine_merged_score) & df$LDrefine_merged_score > 0.25, ]
  df
}

deep_filt <- lapply(deep_list, filter_somatic)
low_filt  <- lapply(low_list, filter_somatic)

cat("Filtering complete.\n\n")

# --- Step 3: Per-chromosome summary ---
cat("--- Per-chromosome Summary ---\n")

chr_summary <- data.frame(
  chr = chromosomes,
  deep_raw = sapply(deep_list, nrow),
  low_raw = sapply(low_list, nrow),
  deep_filtered = sapply(deep_filt, nrow),
  low_filtered = sapply(low_filt, nrow),
  stringsAsFactors = FALSE
)

# Compute overlaps
chr_summary$overlap_raw <- sapply(chromosomes, function(chr) {
  length(intersect(deep_list[[chr]]$snv_id, low_list[[chr]]$snv_id))
})

chr_summary$overlap_filtered <- sapply(chromosomes, function(chr) {
  length(intersect(deep_filt[[chr]]$snv_id, low_filt[[chr]]$snv_id))
})

chr_summary$deep_only_filt <- chr_summary$deep_filtered - chr_summary$overlap_filtered
chr_summary$low_only_filt  <- chr_summary$low_filtered - chr_summary$overlap_filtered
chr_summary$jaccard_raw <- with(chr_summary, overlap_raw / (deep_raw + low_raw - overlap_raw))
chr_summary$jaccard_filt <- with(chr_summary, overlap_filtered / (deep_filtered + low_filtered - overlap_filtered))

cat(sprintf("%-6s %8s %8s %8s %8s %8s %8s %8s\n",
            "CHR", "D_raw", "L_raw", "D_filt", "L_filt", "Ovlp_r", "Ovlp_f", "Jacc_f"))
for (i in 1:nrow(chr_summary)) {
  cat(sprintf("%-6s %8d %8d %8d %8d %8d %8d %8.3f\n",
              chr_summary$chr[i], chr_summary$deep_raw[i], chr_summary$low_raw[i],
              chr_summary$deep_filtered[i], chr_summary$low_filtered[i],
              chr_summary$overlap_raw[i], chr_summary$overlap_filtered[i],
              chr_summary$jaccard_filt[i]))
}

# Save summary
write.csv(chr_summary, file.path(out_dir, "chr_summary.csv"), row.names = FALSE)

# --- Step 4: Combined data for filtered SNVs ---
deep_all <- do.call(rbind, deep_filt)
low_all  <- do.call(rbind, low_filt)

cat("\n--- Overall Summary ---\n")
cat("Deepseq filtered somatic SNVs:", nrow(deep_all), "\n")
cat("Lowseq  filtered somatic SNVs:", nrow(low_all), "\n")
overlap_all <- intersect(deep_all$snv_id, low_all$snv_id)
cat("Overlap (shared positions):", length(overlap_all), "\n")
cat("Deepseq-only:", nrow(deep_all) - length(overlap_all), "\n")
cat("Lowseq-only:", nrow(low_all) - length(overlap_all), "\n")
cat("Jaccard index:", round(length(overlap_all) / (nrow(deep_all) + nrow(low_all) - length(overlap_all)), 4), "\n")

# --- Step 5: Generate plots ---
cat("\n--- Generating plots ---\n")

pdf(file.path(out_dir, "somatic_comparison_plots.pdf"), width = 14, height = 10)

# =============================================
# PLOT 1: Per-chromosome SNV counts (raw vs filtered)
# =============================================
cat("  Plot 1: Per-chromosome SNV counts\n")
par(mfrow = c(1, 2), mar = c(7, 5, 3, 1))

# Raw counts
barplot_data_raw <- rbind(chr_summary$deep_raw, chr_summary$low_raw)
colnames(barplot_data_raw) <- gsub("chr", "", chromosomes)
bp <- barplot(barplot_data_raw, beside = TRUE, col = c("steelblue", "coral"),
              main = "Raw Putative Somatic SNVs per Chromosome",
              ylab = "Number of SNVs", xlab = "", las = 2, cex.names = 0.8)
legend("topright", c("Deepseq", "Lowseq"), fill = c("steelblue", "coral"), cex = 0.8)

# Filtered counts
barplot_data_filt <- rbind(chr_summary$deep_filtered, chr_summary$low_filtered)
colnames(barplot_data_filt) <- gsub("chr", "", chromosomes)
bp <- barplot(barplot_data_filt, beside = TRUE, col = c("steelblue", "coral"),
              main = "Filtered Somatic SNVs per Chromosome",
              ylab = "Number of SNVs", xlab = "", las = 2, cex.names = 0.8)
legend("topright", c("Deepseq", "Lowseq"), fill = c("steelblue", "coral"), cex = 0.8)

# =============================================
# PLOT 2: Overlap Venn-style bar chart per chromosome
# =============================================
cat("  Plot 2: Overlap breakdown per chromosome\n")
par(mfrow = c(1, 1), mar = c(7, 5, 3, 1))

overlap_data <- rbind(chr_summary$deep_only_filt, chr_summary$overlap_filtered, chr_summary$low_only_filt)
colnames(overlap_data) <- gsub("chr", "", chromosomes)
barplot(overlap_data, col = c("steelblue", "purple", "coral"),
        main = "Filtered Somatic SNV Overlap per Chromosome",
        ylab = "Number of SNVs", xlab = "", las = 2, cex.names = 0.8,
        legend.text = c("Deepseq-only", "Shared", "Lowseq-only"),
        args.legend = list(x = "topright", cex = 0.8))

# =============================================
# PLOT 3: Jaccard similarity per chromosome
# =============================================
cat("  Plot 3: Jaccard similarity\n")
par(mfrow = c(1, 2), mar = c(7, 5, 3, 1))

barplot(chr_summary$jaccard_raw, names.arg = gsub("chr", "", chromosomes),
        col = "gray60", main = "Jaccard Similarity (Raw)",
        ylab = "Jaccard Index", las = 2, cex.names = 0.8, ylim = c(0, max(chr_summary$jaccard_raw, na.rm = TRUE) * 1.2))

barplot(chr_summary$jaccard_filt, names.arg = gsub("chr", "", chromosomes),
        col = "darkgreen", main = "Jaccard Similarity (Filtered)",
        ylab = "Jaccard Index", las = 2, cex.names = 0.8, ylim = c(0, max(chr_summary$jaccard_filt, na.rm = TRUE) * 1.2))

# =============================================
# PLOT 4: SVM score distributions
# =============================================
cat("  Plot 4: SVM score distributions\n")
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

if (nrow(deep_all) > 0) {
  hist(deep_all$SVM_pos_score, breaks = 50, col = rgb(0.27, 0.51, 0.71, 0.6),
       main = "SVM Positive Score\n(Filtered Somatic SNVs)",
       xlab = "SVM Positive Score", ylab = "Frequency")
  if (nrow(low_all) > 0) {
    hist(low_all$SVM_pos_score, breaks = 50, col = rgb(1, 0.5, 0.31, 0.5), add = TRUE)
    legend("topleft", c("Deepseq", "Lowseq"), fill = c(rgb(0.27, 0.51, 0.71, 0.6), rgb(1, 0.5, 0.31, 0.5)))
  }
}

# =============================================
# PLOT 5: LD refinement score distributions
# =============================================
cat("  Plot 5: LD refinement scores\n")

if (nrow(deep_all) > 0) {
  hist(deep_all$LDrefine_merged_score, breaks = 50, col = rgb(0.27, 0.51, 0.71, 0.6),
       main = "LD Refinement Merged Score\n(Filtered Somatic SNVs)",
       xlab = "LDrefine Merged Score", ylab = "Frequency")
  if (nrow(low_all) > 0) {
    hist(low_all$LDrefine_merged_score, breaks = 50, col = rgb(1, 0.5, 0.31, 0.5), add = TRUE)
    legend("topleft", c("Deepseq", "Lowseq"), fill = c(rgb(0.27, 0.51, 0.71, 0.6), rgb(1, 0.5, 0.31, 0.5)))
  }
}

# =============================================
# PLOT 6: BAF distributions
# =============================================
cat("  Plot 6: BAF distributions\n")
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

if (nrow(deep_all) > 0) {
  hist(deep_all$BAF_alt, breaks = 50, col = "steelblue",
       main = "BAF Alt - Deepseq\n(Filtered)", xlab = "BAF Alt", ylab = "Frequency")
}

if (nrow(low_all) > 0) {
  hist(low_all$BAF_alt, breaks = 50, col = "coral",
       main = "BAF Alt - Lowseq\n(Filtered)", xlab = "BAF Alt", ylab = "Frequency")
}

# =============================================
# PLOT 7: Depth distributions
# =============================================
cat("  Plot 7: Depth distributions\n")
par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

if (nrow(deep_all) > 0) {
  hist(log10(deep_all$Depth_total + 1), breaks = 50, col = "steelblue",
       main = "Total Depth - Deepseq", xlab = "log10(Depth + 1)")
}
if (nrow(low_all) > 0) {
  hist(log10(low_all$Depth_total + 1), breaks = 50, col = "coral",
       main = "Total Depth - Lowseq", xlab = "log10(Depth + 1)")
}
if (nrow(deep_all) > 0) {
  hist(deep_all$Depth_alt, breaks = 50, col = "steelblue",
       main = "Alt Depth - Deepseq", xlab = "Alt Allele Depth")
}
if (nrow(low_all) > 0) {
  hist(low_all$Depth_alt, breaks = 50, col = "coral",
       main = "Alt Depth - Lowseq", xlab = "Alt Allele Depth")
}

# =============================================
# PLOT 8: Scatter - Deep vs Low filtered count per chr
# =============================================
cat("  Plot 8: Scatter deep vs low per chr\n")
par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))

max_val <- max(c(chr_summary$deep_filtered, chr_summary$low_filtered))
plot(chr_summary$deep_filtered, chr_summary$low_filtered,
     pch = 19, col = "darkblue", cex = 1.5,
     xlim = c(0, max_val * 1.1), ylim = c(0, max_val * 1.1),
     main = "Filtered Somatic SNVs: Deepseq vs Lowseq per Chr",
     xlab = "Deepseq (filtered count)", ylab = "Lowseq (filtered count)")
abline(0, 1, lty = 2, col = "gray50")
text(chr_summary$deep_filtered, chr_summary$low_filtered,
     labels = gsub("chr", "", chromosomes), pos = 3, cex = 0.7)
cor_val <- cor(chr_summary$deep_filtered, chr_summary$low_filtered)
legend("topleft", paste0("r = ", round(cor_val, 3)), bty = "n", cex = 1.2)

# =============================================
# PLOT 9: Mutation spectrum (Ti/Tv) comparison
# =============================================
cat("  Plot 9: Mutation spectrum\n")
par(mfrow = c(1, 2), mar = c(7, 5, 3, 1))

get_mutation_type <- function(ref, alt) {
  transitions <- c("AG", "GA", "CT", "TC")
  mut <- paste0(ref, alt)
  ifelse(mut %in% transitions, "Transition", "Transversion")
}

get_mut_class <- function(ref, alt) {
  paste0(ref, ">", alt)
}

if (nrow(deep_all) > 0 && nrow(low_all) > 0) {
  deep_all$mut_type <- get_mutation_type(deep_all$Ref_allele, deep_all$Alt_allele)
  low_all$mut_type  <- get_mutation_type(low_all$Ref_allele, low_all$Alt_allele)
  deep_all$mut_class <- get_mut_class(deep_all$Ref_allele, deep_all$Alt_allele)
  low_all$mut_class  <- get_mut_class(low_all$Ref_allele, low_all$Alt_allele)

  # Ti/Tv per chromosome
  tstv_deep <- sapply(chromosomes, function(chr) {
    d <- deep_filt[[chr]]
    if (is.null(d) || nrow(d) == 0) return(NA)
    mt <- get_mutation_type(d$Ref_allele, d$Alt_allele)
    ti <- sum(mt == "Transition")
    tv <- sum(mt == "Transversion")
    if (tv == 0) return(NA)
    ti / tv
  })

  tstv_low <- sapply(chromosomes, function(chr) {
    d <- low_filt[[chr]]
    if (is.null(d) || nrow(d) == 0) return(NA)
    mt <- get_mutation_type(d$Ref_allele, d$Alt_allele)
    ti <- sum(mt == "Transition")
    tv <- sum(mt == "Transversion")
    if (tv == 0) return(NA)
    ti / tv
  })

  tstv_data <- rbind(tstv_deep, tstv_low)
  colnames(tstv_data) <- gsub("chr", "", chromosomes)
  barplot(tstv_data, beside = TRUE, col = c("steelblue", "coral"),
          main = "Ti/Tv Ratio per Chromosome\n(Filtered Somatic)",
          ylab = "Ti/Tv Ratio", las = 2, cex.names = 0.8)
  abline(h = 2, lty = 2, col = "gray50")
  legend("topright", c("Deepseq", "Lowseq"), fill = c("steelblue", "coral"), cex = 0.8)

  # Mutation spectrum bar chart
  all_classes <- sort(unique(c(deep_all$mut_class, low_all$mut_class)))
  deep_spec <- table(factor(deep_all$mut_class, levels = all_classes))
  low_spec  <- table(factor(low_all$mut_class, levels = all_classes))
  deep_spec_pct <- deep_spec / sum(deep_spec) * 100
  low_spec_pct  <- low_spec / sum(low_spec) * 100
  spec_data <- rbind(deep_spec_pct, low_spec_pct)

  barplot(spec_data, beside = TRUE, col = c("steelblue", "coral"),
          main = "Mutation Spectrum (Filtered Somatic)",
          ylab = "Percentage (%)", las = 2, cex.names = 0.7)
  legend("topright", c("Deepseq", "Lowseq"), fill = c("steelblue", "coral"), cex = 0.8)
}

# =============================================
# PLOT 10: Per-chr overlap proportions
# =============================================
cat("  Plot 10: Overlap proportions\n")
par(mfrow = c(1, 2), mar = c(7, 5, 3, 1))

# What fraction of deep is shared with low?
deep_shared_pct <- ifelse(chr_summary$deep_filtered > 0,
                          chr_summary$overlap_filtered / chr_summary$deep_filtered * 100, 0)
low_shared_pct <- ifelse(chr_summary$low_filtered > 0,
                         chr_summary$overlap_filtered / chr_summary$low_filtered * 100, 0)

barplot(deep_shared_pct, names.arg = gsub("chr", "", chromosomes),
        col = "steelblue", main = "% of Deepseq SNVs Also in Lowseq",
        ylab = "Percentage (%)", las = 2, cex.names = 0.8, ylim = c(0, 100))

barplot(low_shared_pct, names.arg = gsub("chr", "", chromosomes),
        col = "coral", main = "% of Lowseq SNVs Also in Deepseq",
        ylab = "Percentage (%)", las = 2, cex.names = 0.8, ylim = c(0, 100))

# =============================================
# PLOT 11: Filtering effect (raw to filtered)
# =============================================
cat("  Plot 11: Filtering effect\n")
par(mfrow = c(1, 2), mar = c(7, 5, 3, 1))

deep_filt_pct <- chr_summary$deep_filtered / chr_summary$deep_raw * 100
low_filt_pct  <- chr_summary$low_filtered / chr_summary$low_raw * 100

filt_data <- rbind(deep_filt_pct, low_filt_pct)
colnames(filt_data) <- gsub("chr", "", chromosomes)
barplot(filt_data, beside = TRUE, col = c("steelblue", "coral"),
        main = "% SNVs Passing Filter per Chromosome",
        ylab = "% Passing", las = 2, cex.names = 0.8)
legend("topright", c("Deepseq", "Lowseq"), fill = c("steelblue", "coral"), cex = 0.8)

# Absolute numbers lost
deep_lost <- chr_summary$deep_raw - chr_summary$deep_filtered
low_lost  <- chr_summary$low_raw - chr_summary$low_filtered
lost_data <- rbind(deep_lost, low_lost)
colnames(lost_data) <- gsub("chr", "", chromosomes)
barplot(lost_data, beside = TRUE, col = c("steelblue", "coral"),
        main = "SNVs Removed by Filtering",
        ylab = "Count Removed", las = 2, cex.names = 0.8)
legend("topright", c("Deepseq", "Lowseq"), fill = c("steelblue", "coral"), cex = 0.8)

# =============================================
# PLOT 12: Shared SNV characteristics
# =============================================
cat("  Plot 12: Shared vs unique SNV characteristics\n")
par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

if (length(overlap_all) > 0) {
  deep_shared <- deep_all[deep_all$snv_id %in% overlap_all, ]
  deep_unique <- deep_all[!deep_all$snv_id %in% overlap_all, ]
  low_shared  <- low_all[low_all$snv_id %in% overlap_all, ]
  low_unique  <- low_all[!low_all$snv_id %in% overlap_all, ]

  # SVM scores shared vs unique
  boxplot(list("Deep\nShared" = deep_shared$SVM_pos_score,
               "Deep\nUnique" = deep_unique$SVM_pos_score,
               "Low\nShared" = low_shared$SVM_pos_score,
               "Low\nUnique" = low_unique$SVM_pos_score),
          col = c("steelblue", "lightblue", "coral", "lightyellow"),
          main = "SVM Score: Shared vs Unique", ylab = "SVM Positive Score")

  # LD scores shared vs unique
  boxplot(list("Deep\nShared" = deep_shared$LDrefine_merged_score,
               "Deep\nUnique" = deep_unique$LDrefine_merged_score,
               "Low\nShared" = low_shared$LDrefine_merged_score,
               "Low\nUnique" = low_unique$LDrefine_merged_score),
          col = c("steelblue", "lightblue", "coral", "lightyellow"),
          main = "LD Refinement Score: Shared vs Unique", ylab = "LDrefine Merged Score")

  # BAF shared vs unique
  boxplot(list("Deep\nShared" = deep_shared$BAF_alt,
               "Deep\nUnique" = deep_unique$BAF_alt,
               "Low\nShared" = low_shared$BAF_alt,
               "Low\nUnique" = low_unique$BAF_alt),
          col = c("steelblue", "lightblue", "coral", "lightyellow"),
          main = "BAF Alt: Shared vs Unique", ylab = "BAF Alt")

  # Depth shared vs unique
  boxplot(list("Deep\nShared" = log10(deep_shared$Depth_total + 1),
               "Deep\nUnique" = log10(deep_unique$Depth_total + 1),
               "Low\nShared" = log10(low_shared$Depth_total + 1),
               "Low\nUnique" = log10(low_unique$Depth_total + 1)),
          col = c("steelblue", "lightblue", "coral", "lightyellow"),
          main = "Total Depth: Shared vs Unique", ylab = "log10(Depth + 1)")
}

# =============================================
# PLOT 13: BAF correlation for shared SNVs
# =============================================
cat("  Plot 13: BAF correlation for shared SNVs\n")
par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))

if (length(overlap_all) > 0) {
  deep_shared_ord <- deep_all[deep_all$snv_id %in% overlap_all, ]
  low_shared_ord  <- low_all[low_all$snv_id %in% overlap_all, ]

  rownames(deep_shared_ord) <- deep_shared_ord$snv_id
  rownames(low_shared_ord)  <- low_shared_ord$snv_id

  common_ids <- intersect(deep_shared_ord$snv_id, low_shared_ord$snv_id)
  deep_shared_ord <- deep_shared_ord[common_ids, ]
  low_shared_ord  <- low_shared_ord[common_ids, ]

  plot(deep_shared_ord$BAF_alt, low_shared_ord$BAF_alt,
       pch = 19, cex = 0.5, col = rgb(0, 0, 0, 0.3),
       main = "BAF Alt Correlation: Shared Somatic SNVs",
       xlab = "Deepseq BAF Alt", ylab = "Lowseq BAF Alt",
       xlim = c(0, 0.5), ylim = c(0, 0.5))
  abline(0, 1, lty = 2, col = "red")
  baf_cor <- cor(deep_shared_ord$BAF_alt, low_shared_ord$BAF_alt, use = "complete.obs")
  legend("topleft", paste0("r = ", round(baf_cor, 3), "\nn = ", length(common_ids)),
         bty = "n", cex = 1.2)
}

dev.off()
cat("\nPlots saved to:", file.path(out_dir, "somatic_comparison_plots.pdf"), "\n")

# --- Save summary data ---
write.csv(deep_all, file.path(out_dir, "deepseq_filtered_somatic_all.csv"), row.names = FALSE)
write.csv(low_all, file.path(out_dir, "lowseq_filtered_somatic_all.csv"), row.names = FALSE)

cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Done ===\n")
