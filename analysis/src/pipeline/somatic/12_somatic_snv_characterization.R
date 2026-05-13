##############################################################################
# 12_somatic_snv_characterization.R
#
# Deep characterization of somatic SNVs across deepseq / lowseq and tissues
# (488B vs 489).  Builds on the variant sets collected by script 11.
#
# Analysis objectives
#   1. Reconstruct per-tissue variant sets from script-11 outputs
#   2. Characterise "only-deepseq", "only-lowseq", and "shared" variants:
#        depth, VAF, quality scores, Ti/Tv, trinucleotide context
#   3. Compare variants by tissue identity (488B vs 489)
#   4. Annotate every variant with the nearest gene (TxDb hg38)
#   5. Flag variants in breast-cancer / TNBC driver genes
#   6. Produce summary tables and plots
#
# Env vars (all optional)
#   COMP_DIR     : dir containing script-11 TSV outputs
#                  [default Data/05_results/variant_calling/tissue_variants/tables]
#   OUT_DIR      : table output directory
#                  [default Data/05_results/variant_calling/somatic_characterization/tables]
#   PLOT_DIR     : plot output directory
#                  [default analysis/plots/comparison/somatic_char]
#   NOTE_DIR     : interpretation notes output directory
#                  [default Data/05_results/variant_calling/somatic_characterization/notes]
#   DATASETS     : comma-separated  [default deepseq,lowseq]
#   CHR_START    : first chromosome number [default 1]
#   CHR_END      : last  chromosome number [default 22]
#   N_WORKERS    : parallel workers        [default NSLOTS or detectCores]
#   UPSTREAM_BP  : bp upstream of gene start to include [default 5000]
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# ── helpers ──────────────────────────────────────────────────────────────────

log_info <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(fmt, ...)))
  flush.console()
}

get_env_int <- function(name, default, min_value = NULL, allow_na = FALSE) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) {
    value <- default
  } else if (allow_na && toupper(raw) %in% c("NA", "NONE", "NULL")) {
    value <- NA_integer_
  } else {
    value <- suppressWarnings(as.integer(raw))
    if (is.na(value)) stop(sprintf("%s must be integer, got '%s'", name, raw), call. = FALSE)
  }
  if (!is.na(value) && !is.null(min_value) && value < min_value)
    stop(sprintf("%s must be >= %d, got %d", name, min_value, value), call. = FALSE)
  value
}

get_env_csv <- function(name, default) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) return(default)
  trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
}

get_env_num <- function(name, default, min_value = NULL) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  value <- if (!nzchar(raw)) default else suppressWarnings(as.numeric(raw))
  if (is.na(value)) stop(sprintf("%s must be numeric, got '%s'", name, raw), call. = FALSE)
  if (!is.null(min_value) && value < min_value)
    stop(sprintf("%s must be >= %s, got %s", name, format(min_value), format(value)), call. = FALSE)
  value
}

# ── configuration ─────────────────────────────────────────────────────────────

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
variant_root <- file.path(project_root, "Data/variant_calling")

comp_dir <- Sys.getenv("COMP_DIR",
  unset = file.path(project_root, "Data/05_results/variant_calling/tissue_variants/tables"))
out_dir  <- Sys.getenv("OUT_DIR",
  unset = file.path(project_root, "Data/05_results/variant_calling/somatic_characterization/tables"))
plot_dir <- Sys.getenv("PLOT_DIR",
  unset = file.path(project_root, "analysis/plots/comparison/somatic_char"))
note_dir <- Sys.getenv("NOTE_DIR",
  unset = file.path(project_root, "Data/05_results/variant_calling/somatic_characterization/notes"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(note_dir, recursive = TRUE, showWarnings = FALSE)

datasets     <- get_env_csv("DATASETS", c("deepseq", "lowseq"))
chr_start    <- get_env_int("CHR_START",   1L, min_value = 1L)
chr_end      <- get_env_int("CHR_END",    22L, min_value = 1L)
n_workers    <- get_env_int("N_WORKERS",
  {ns <- suppressWarnings(as.integer(Sys.getenv("NSLOTS", "")));
   if (!is.na(ns) && ns > 0L) ns else max(1L, parallel::detectCores(logical = FALSE), na.rm = TRUE)},
  min_value = 1L)
upstream_bp  <- get_env_int("UPSTREAM_BP", 5000L, min_value = 0L)

# Matched practical cutoffs applied identically to deepseq and lowseq.
min_dep       <- get_env_int("MIN_DEP", 20L, min_value = 1L)
min_alt_count <- get_env_int("MIN_ALT_COUNT", 5L, min_value = 1L)
min_alt_cells <- get_env_int("MIN_ALT_CELLS", 3L, min_value = 1L)
min_vaf       <- get_env_num("MIN_VAF", 0.03, min_value = 0)
min_svm       <- get_env_num("MIN_SVM", 0.90, min_value = 0)
min_qs        <- get_env_num("MIN_QS", 0.20, min_value = 0)

chr_numbers  <- seq.int(chr_start, chr_end)

log_info("somatic SNV characterization | out_dir=%s | plot_dir=%s | note_dir=%s | chrs=chr%d-chr%d | workers=%d",
         out_dir, plot_dir, note_dir, chr_start, chr_end, n_workers)
log_info("Matched filters | Dep>=%d alt_count>=%d alt_cells>=%d VAF>=%.3f SVM>=%.2f QS>=%.2f",
         min_dep, min_alt_count, min_alt_cells, min_vaf, min_svm, min_qs)

##############################################################################
# SECTION 1 – Reconstruct variant sets from script-11 outputs
##############################################################################

load_s11_sets <- function(dataset) {
  # Script 11 writes:
  #   {dataset}_same_488B_489.tsv  – variants present in BOTH tissues
  #   {dataset}_only_488B.tsv      – variants in 488B but NOT 489
  #   {dataset}_only_489.tsv       – variants in 489 but NOT 488B
  read_ids <- function(fname) {
    f <- file.path(comp_dir, fname)
    if (!file.exists(f)) {
      stop(sprintf("Script-11 output missing: %s\n  Run 11_comparing_tissue_variants.R first.", f),
           call. = FALSE)
    }
    dt <- fread(f)
    as.character(dt$variant_id)
  }

  same_ids    <- read_ids(sprintf("%s_same_488B_489.tsv", dataset))
  only488_ids <- read_ids(sprintf("%s_only_488B.tsv",     dataset))
  only489_ids <- read_ids(sprintf("%s_only_489.tsv",      dataset))

  list(
    v488  = union(same_ids, only488_ids),
    v489  = union(same_ids, only489_ids),
    same  = same_ids,
    only488 = only488_ids,
    only489 = only489_ids
  )
}

build_sets <- function(deep, low) {
  list(
    # Within 488B
    shared_488B         = intersect(deep$v488, low$v488),
    deep_only_488B      = setdiff(deep$v488,   low$v488),
    low_only_488B       = setdiff(low$v488,    deep$v488),

    # Within 489
    shared_489          = intersect(deep$v489, low$v489),
    deep_only_489       = setdiff(deep$v489,   low$v489),
    low_only_489        = setdiff(low$v489,    deep$v489),

    # Tissue identity (across both datasets)
    tissue_488B_only    = setdiff(union(deep$v488, low$v488),
                                  union(deep$v489, low$v489)),
    tissue_489_only     = setdiff(union(deep$v489, low$v489),
                                  union(deep$v488, low$v488)),
    tissue_both         = intersect(union(deep$v488, low$v488),
                                    union(deep$v489, low$v489)),

    # Within-dataset tissue-unique
    deep_488B_tissue_unique = deep$only488,
    deep_489_tissue_unique  = deep$only489,
    low_488B_tissue_unique  = low$only488,
    low_489_tissue_unique   = low$only489
  )
}

log_info("Loading script-11 variant sets …")
s11 <- setNames(lapply(datasets, load_s11_sets), datasets)

deep_raw <- s11[["deepseq"]]
low_raw  <- s11[["lowseq"]]
sets_raw <- build_sets(deep_raw, low_raw)

set_sizes <- data.table(
  set_name = names(sets_raw),
  n_variants = vapply(sets_raw, length, integer(1L))
)
log_info("Unfiltered variant set sizes (from script-11):\n%s", paste(capture.output(print(set_sizes)), collapse = "\n"))
fwrite(set_sizes, file.path(out_dir, "variant_set_sizes_unfiltered.tsv"), sep = "\t")

##############################################################################
# SECTION 2 – Build variant attribute table from allSNVs.csv + putativeSNVs.csv
##############################################################################

# All unique variant IDs we care about (union of unfiltered sets)
all_ids <- unique(unlist(sets_raw, use.names = FALSE))
log_info("Total unique variants to characterise: %s", format(length(all_ids), big.mark = ","))

# Parse variant IDs → data.table
parse_ids <- function(ids) {
  parts <- strsplit(ids, ":", fixed = TRUE)
  data.table(
    variant_id = ids,
    chr = vapply(parts, `[[`, character(1L), 1L),
    pos = as.integer(vapply(parts, `[[`, character(1L), 2L)),
    ref = vapply(parts, `[[`, character(1L), 3L),
    alt = vapply(parts, `[[`, character(1L), 4L)
  )
}
id_dt <- parse_ids(all_ids)
setkeyv(id_dt, c("chr", "pos", "ref", "alt"))

# Load allSNVs attributes for one dataset (deepseq preferred; fallback lowseq)
load_allsnvs_chr <- function(dataset, chr_num) {
  chr <- sprintf("chr%d", chr_num)
  f   <- file.path(variant_root, dataset, "somatic", sprintf("%s.allSNVs.csv", chr))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, chr := chr]
  dt
}

load_putative_chr <- function(dataset, chr_num) {
  chr <- sprintf("chr%d", chr_num)
  f   <- file.path(variant_root, dataset, "somatic", sprintf("%s.putativeSNVs.csv", chr))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  setnames(dt, c("Ref_allele", "Alt_allele"), c("ref", "alt"), skip_absent = TRUE)
  dt[, chr := chr]
  dt
}

load_attrs_for_dataset <- function(dataset) {
  log_info("Loading allSNVs attributes for %s …", dataset)
  all_list <- mclapply(chr_numbers, function(n) load_allsnvs_chr(dataset, n), mc.cores = n_workers)
  all_dt   <- rbindlist(Filter(Negate(is.null), all_list), fill = TRUE)
  if (!"ref" %in% names(all_dt)) setnames(all_dt, "Ref_allele", "ref", skip_absent = TRUE)
  if (!"alt" %in% names(all_dt)) setnames(all_dt, "Alt_allele", "alt", skip_absent = TRUE)
  setkeyv(all_dt, c("chr", "pos", "ref", "alt"))

  log_info("Loading putativeSNVs for %s …", dataset)
  put_list <- mclapply(chr_numbers, function(n) load_putative_chr(dataset, n), mc.cores = n_workers)
  put_dt   <- rbindlist(Filter(Negate(is.null), put_list), fill = TRUE)
  setkeyv(put_dt, c("chr", "pos", "ref", "alt"))

  # Merge: putative fields onto allSNVs
  keep_put <- setdiff(names(put_dt), c("chr", "pos", "ref", "alt", "Depth_total",
                                        "Depth_ref", "Depth_alt"))
  if (length(keep_put)) {
    all_dt <- merge(all_dt, put_dt[, c("chr", "pos", "ref", "alt", keep_put), with = FALSE],
                    by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
  }
  all_dt[, dataset := dataset]
  all_dt
}

apply_matched_filters <- function(dt, dataset_label) {
  if (!nrow(dt)) return(dt)

  if (!("VAF" %in% names(dt)) && all(c("dep_alt_new", "Dep") %in% names(dt))) {
    dt[, VAF := dep_alt_new / pmax(Dep, 1L)]
  }

  keep <- rep(TRUE, nrow(dt))

  if ("Dep" %in% names(dt)) keep <- keep & !is.na(dt$Dep) & dt$Dep >= min_dep
  if ("dep_alt_new" %in% names(dt)) keep <- keep & !is.na(dt$dep_alt_new) & dt$dep_alt_new >= min_alt_count
  if ("cell_alt" %in% names(dt)) keep <- keep & !is.na(dt$cell_alt) & dt$cell_alt >= min_alt_cells
  if ("VAF" %in% names(dt)) keep <- keep & !is.na(dt$VAF) & dt$VAF >= min_vaf
  if ("SVM_pos_score" %in% names(dt)) keep <- keep & !is.na(dt$SVM_pos_score) & dt$SVM_pos_score >= min_svm
  if ("QS" %in% names(dt)) keep <- keep & !is.na(dt$QS) & dt$QS >= min_qs

  out <- dt[keep]
  log_info("Applied matched filters to %s: %s -> %s rows (%.1f%% kept)",
           dataset_label,
           format(nrow(dt), big.mark = ","),
           format(nrow(out), big.mark = ","),
           100 * nrow(out) / max(1, nrow(dt)))
  out
}

# Load attributes for each dataset, then merge with our id_dt
attr_list <- lapply(datasets, function(ds) {
  a <- load_attrs_for_dataset(ds)
  # Keep only variants we care about
  a <- a[id_dt[, .(chr, pos, ref, alt, variant_id)], on = c("chr", "pos", "ref", "alt"), nomatch = 0L]
  a <- apply_matched_filters(a, ds)
  a
})
names(attr_list) <- datasets

# Rebuild per-tissue sets using variants that pass the same filters in each dataset.
deep_pass <- unique(attr_list[["deepseq"]]$variant_id)
low_pass  <- unique(attr_list[["lowseq"]]$variant_id)

deep <- list(
  v488    = intersect(deep_raw$v488, deep_pass),
  v489    = intersect(deep_raw$v489, deep_pass),
  same    = intersect(deep_raw$same, deep_pass),
  only488 = intersect(deep_raw$only488, deep_pass),
  only489 = intersect(deep_raw$only489, deep_pass)
)
low <- list(
  v488    = intersect(low_raw$v488, low_pass),
  v489    = intersect(low_raw$v489, low_pass),
  same    = intersect(low_raw$same, low_pass),
  only488 = intersect(low_raw$only488, low_pass),
  only489 = intersect(low_raw$only489, low_pass)
)
sets <- build_sets(deep, low)

set_sizes <- data.table(
  set_name = names(sets),
  n_variants = vapply(sets, length, integer(1L))
)
log_info("Filtered variant set sizes:\n%s", paste(capture.output(print(set_sizes)), collapse = "\n"))
fwrite(set_sizes, file.path(out_dir, "variant_set_sizes.tsv"), sep = "\t")

# Combine, keeping one row per variant × dataset (averaged if duplicates exist)
attr_all <- rbindlist(attr_list, fill = TRUE)

# Compute derived metrics per variant×dataset ─────────────────────────────
# Allele frequency in the bulk
if ("dep_alt_new" %in% names(attr_all) && "Dep" %in% names(attr_all)) {
  attr_all[, VAF := dep_alt_new / pmax(Dep, 1L)]
}

# Transition / Transversion
is_transition <- function(ref, alt) {
  pairs <- paste(toupper(ref), toupper(alt), sep = ">")
  pairs %in% c("A>G", "G>A", "C>T", "T>C")
}
attr_all[, ti_tv := {
  ok <- !is.na(ref) & !is.na(alt) & nchar(as.character(ref)) == 1 & nchar(as.character(alt)) == 1
  out <- rep(NA_character_, .N)
  if (any(ok)) {
    ref_ok <- ref[ok]
    alt_ok <- alt[ok]
    tri_idx <- is_transition(ref_ok, alt_ok)
    out[which(ok)[tri_idx]] <- "Ti"
    out[which(ok)[!tri_idx]] <- "Tv"
  }
  out
}]

# SBS trinucleotide context (C/T reference convention) ───────────────────
# Skippable: set env var SKIP_TRINUC=1 to skip if BSgenome step is too slow
SKIP_TRINUC <- Sys.getenv("SKIP_TRINUC", "0") %in% c("1", "TRUE", "true", "T", "yes")
if (SKIP_TRINUC) {
  log_info("Skipping trinucleotide context (SKIP_TRINUC=1)")
  attr_all[, tri_context := NA_character_]
} else {
log_info("Extracting trinucleotide context from BSgenome …")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# Work on the unique chr:pos:ref:alt coordinates
context_dt <- unique(attr_all[, .(chr, pos, ref, alt)])
context_dt[, tri_context := {
  # For each variant, get the -1/+1 flanking bases
  # pyrimidine convention: if ref is purine, take complement strand
  ref_upper <- toupper(ref)
  alt_upper <- toupper(alt)

  seqs <- tryCatch(
    as.character(getSeq(bsgenome, names = chr,
                        start = pos - 1L, end = pos + 1L)),
    error = function(e) rep(NA_character_, .N)
  )
  # Force to pyrimidine reference (C or T)
  comp <- c(A="T", T="A", C="G", G="C")
  rc_triplet <- function(t) {
    if (is.na(t)) return(NA_character_)
    bases <- strsplit(t, "")[[1]]
    paste(rev(comp[bases]), collapse = "")
  }
  out <- character(length(seqs))
  for (i in seq_along(seqs)) {
    tri <- seqs[i]
    if (is.na(tri) || nchar(tri) != 3L) { out[i] <- NA_character_; next }
    mid <- toupper(substr(tri, 2L, 2L))
    if (mid %in% c("A", "G")) tri <- rc_triplet(tri)
    out[i] <- tri
  }
  out
}, by = .(chr, pos, ref, alt)]

attr_all <- merge(attr_all, context_dt,
                  by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
} # end if (!SKIP_TRINUC)

##############################################################################
# SECTION 3 – Gene annotation (TxDb hg38 + org.Hs.eg.db)
##############################################################################

log_info("Building gene model from TxDb.Hsapiens.UCSC.hg38.knownGene …")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

# Gene body ranges + upstream promoter window
gene_gr <- suppressMessages(genes(txdb))
# Add upstream window for promoter coverage
promoter_gr <- GenomicRanges::promoters(gene_gr,
                                        upstream = upstream_bp,
                                        downstream = 0L)
combined_gr <- c(gene_gr, promoter_gr)

# Entrez → gene symbol map
entrez_ids    <- names(gene_gr)
sym_map_raw   <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys    = entrez_ids,
                        columns = c("SYMBOL", "GENENAME"),
                        keytype = "ENTREZID"))
sym_map <- as.data.table(sym_map_raw)
setnames(sym_map, c("ENTREZID", "SYMBOL", "GENENAME"), c("entrezid", "symbol", "genename"))
sym_map <- sym_map[!duplicated(entrezid)]

# Overlap variants with gene ranges ─────────────────────────────────────────
log_info("Annotating %s variant coordinates …", format(nrow(id_dt), big.mark = ","))

var_gr <- makeGRangesFromDataFrame(
  data.frame(seqnames = id_dt$chr,
             start    = id_dt$pos,
             end      = id_dt$pos,
             variant_id = id_dt$variant_id),
  keep.extra.columns = TRUE
)

hits <- findOverlaps(var_gr, combined_gr, select = "all")

# Collapse to one gene per variant (prefer gene body over promoter)
# combined_gr first half = gene_gr, second half = promoter_gr
n_genes <- length(gene_gr)
gene_hits_dt <- data.table(
  var_idx    = queryHits(hits),
  gene_idx   = subjectHits(hits) %% n_genes,   # maps promoter back to gene
  is_promoter = subjectHits(hits) > n_genes
)
gene_hits_dt[gene_idx == 0L, gene_idx := n_genes]   # modulo fix for exact n_genes

gene_hits_dt[, entrezid := names(gene_gr)[gene_idx]]
gene_hits_dt <- merge(gene_hits_dt, sym_map, by = "entrezid", all.x = TRUE)
gene_hits_dt[, variant_id := id_dt$variant_id[var_idx]]

# For each variant pick the gene with largest overlap (prefer body over promoter)
gene_hits_dt <- gene_hits_dt[order(variant_id, is_promoter)]
gene_annot   <- gene_hits_dt[!duplicated(variant_id),
                              .(variant_id, entrezid, symbol, genename, is_promoter)]

# Variants with no gene hit get NA
all_ids_dt   <- data.table(variant_id = all_ids)
gene_annot   <- merge(all_ids_dt, gene_annot, by = "variant_id", all.x = TRUE)

log_info("Variants with gene annotation: %d / %d",
         sum(!is.na(gene_annot$symbol)), nrow(gene_annot))

##############################################################################
# SECTION 4 – TNBC / breast-cancer driver gene list
##############################################################################

tnbc_genes <- list(
  # Classic tumour suppressors & drivers
  core_brca_tsg = c("TP53", "BRCA1", "BRCA2", "PTEN", "RB1"),

  # PI3K pathway
  pi3k_pathway  = c("PIK3CA", "PIK3R1", "AKT1", "AKT2", "AKT3", "MTOR",
                     "TSC1", "TSC2"),

  # Cell cycle / amplicons
  cell_cycle     = c("MYC", "CCND1", "CCNE1", "CDK4", "CDK6", "MDM2",
                     "CDKN2A", "CDKN1A", "RB1"),

  # TNBC-enriched / basal-like
  tnbc_basal     = c("KIT", "EGFR", "FGFR1", "FGFR2", "FGFR4",
                     "AR", "BRD4", "MET", "NOTCH1", "NOTCH2",
                     "JAK2", "STAT3", "BCL2", "BCL2L1"),

  # DNA damage response
  ddr            = c("ATM", "ATR", "CHEK1", "CHEK2", "BRIP1",
                     "PALB2", "RAD51C", "RAD51D", "NBN",
                     "RAD50", "MRE11"),

  # Epigenetic regulators
  epigenetic     = c("KMT2C", "KMT2D", "ARID1A", "ARID1B", "SMAD4",
                     "SETD2", "CREBBP", "EP300", "KDM6A"),

  # Hormone receptors & growth factor signalling
  hr_gf          = c("ESR1", "ERBB2", "ERBB3", "FGFR1", "IGF1R",
                     "VEGFA", "VEGFB", "VEGFC", "KDR"),

  # EMT / metastasis
  emt            = c("CDH1", "CDH2", "VIM", "FN1", "SNAI1", "SNAI2",
                     "ZEB1", "ZEB2", "TWIST1", "MMP2", "MMP9"),

  # Immune checkpoint
  immune         = c("CD274", "PDCD1LG2", "JAK1", "JAK2",
                     "STAT1", "B2M", "HLA-A", "HLA-B", "HLA-C"),

  # Other frequently mutated in breast cancer
  other_brca     = c("GATA3", "FOXA1", "TBX3", "RUNX1",
                     "MAP3K1", "MAP2K4", "NF1", "SF3B1",
                     "CBFB", "TP53BP1")
)
all_tnbc_genes <- unique(unlist(tnbc_genes))

gene_annot[, is_cancer_gene   := symbol %in% all_tnbc_genes]
gene_annot[, cancer_categories := {
  cats <- vapply(symbol, function(s) {
    if (is.na(s)) return("")
    hits <- names(tnbc_genes)[vapply(tnbc_genes, function(g) s %in% g, logical(1L))]
    paste(hits, collapse = ";")
  }, character(1L))
  cats
}]

##############################################################################
# SECTION 5 – Annotated variant master table
##############################################################################

log_info("Building master annotated variant table …")

# Assign each variant to the sets it belongs to
all_ids_filtered <- unique(unlist(sets, use.names = FALSE))
variant_membership <- data.table(variant_id = all_ids_filtered)
for (sname in names(sets)) {
  variant_membership[, (sname) := variant_id %in% sets[[sname]]]
}

# Summary: which named "cross-dataset" categories does each variant belong to?
core_cats <- c("shared_488B", "deep_only_488B", "low_only_488B",
                "shared_489",  "deep_only_489",  "low_only_489",
                "tissue_488B_only", "tissue_489_only", "tissue_both")

variant_membership[, cross_category := {
  cats <- apply(.SD, 1L, function(row) paste(names(row)[as.logical(row)], collapse = ";"))
  cats
}, .SDcols = core_cats]

master <- merge(id_dt[, .(variant_id, chr, pos, ref, alt)],
                gene_annot,  by = "variant_id", all.x = TRUE)
master <- merge(master, variant_membership, by = "variant_id", all = FALSE)

# Pull in deepseq attrs (preferred) then lowseq for missing
attrs_deep <- attr_list[["deepseq"]][, .(
  variant_id, Dep_deep = Dep,
  VAF_deep = if (exists("VAF", inherits = FALSE)) VAF else dep_alt_new / pmax(Dep, 1L),
  dep_alt_deep = dep_alt_new, dep_ref_deep = dep_ref_new,
  cell_alt_deep = cell_alt, cell_ref_deep = cell_ref,
  QS_deep = if ("QS" %in% names(.SD)) QS else NA_real_,
  ti_tv = if ("ti_tv" %in% names(.SD)) ti_tv else NA_character_,
  tri_context = if ("tri_context" %in% names(.SD)) tri_context else NA_character_,
  SVM_pos_score_deep   = if ("SVM_pos_score"   %in% names(.SD)) SVM_pos_score   else NA_real_,
  LDrefine_score_deep  = if ("LDrefine_merged_score" %in% names(.SD)) LDrefine_merged_score else NA_real_,
  BAF_alt_deep         = if ("BAF_alt" %in% names(.SD)) BAF_alt else NA_real_
)]

attrs_low <- attr_list[["lowseq"]][, .(
  variant_id, Dep_low = Dep,
  VAF_low = if (exists("VAF", inherits = FALSE)) VAF else dep_alt_new / pmax(Dep, 1L),
  dep_alt_low = dep_alt_new, dep_ref_low = dep_ref_new,
  cell_alt_low = cell_alt, cell_ref_low = cell_ref,
  QS_low = if ("QS" %in% names(.SD)) QS else NA_real_,
  SVM_pos_score_low   = if ("SVM_pos_score"   %in% names(.SD)) SVM_pos_score   else NA_real_,
  LDrefine_score_low  = if ("LDrefine_merged_score" %in% names(.SD)) LDrefine_merged_score else NA_real_,
  BAF_alt_low         = if ("BAF_alt" %in% names(.SD)) BAF_alt else NA_real_
)]

# De-dup: take one row per variant from each dataset (first occurrence)
attrs_deep <- attrs_deep[!duplicated(variant_id)]
attrs_low  <- attrs_low[!duplicated(variant_id)]

master <- merge(master, attrs_deep, by = "variant_id", all.x = TRUE)
master <- merge(master, attrs_low,  by = "variant_id", all.x = TRUE)

fwrite(master, file.path(out_dir, "master_annotated_variants.tsv"), sep = "\t")
log_info("Wrote master table: %d variants", nrow(master))

##############################################################################
# SECTION 6 – Per-set characteristic summary
##############################################################################

compute_set_summary <- function(ids, name) {
  sub <- master[variant_id %in% ids]
  if (!nrow(sub)) return(data.table(set_name = name, n = 0L))

  data.table(
    set_name              = name,
    n                     = nrow(sub),
    n_cancer_gene         = sum(sub$is_cancer_gene, na.rm = TRUE),
    pct_cancer_gene       = round(mean(sub$is_cancer_gene, na.rm = TRUE) * 100, 2),
    median_Dep_deep       = median(sub$Dep_deep,    na.rm = TRUE),
    median_Dep_low        = median(sub$Dep_low,     na.rm = TRUE),
    mean_VAF_deep         = round(mean(sub$VAF_deep, na.rm = TRUE), 4),
    mean_VAF_low          = round(mean(sub$VAF_low,  na.rm = TRUE), 4),
    mean_cell_alt_deep    = round(mean(sub$cell_alt_deep, na.rm = TRUE), 2),
    mean_cell_alt_low     = round(mean(sub$cell_alt_low,  na.rm = TRUE), 2),
    mean_QS_deep          = round(mean(sub$QS_deep,  na.rm = TRUE), 4),
    mean_QS_low           = round(mean(sub$QS_low,   na.rm = TRUE), 4),
    pct_Ti                = round(mean(sub$ti_tv == "Ti", na.rm = TRUE) * 100, 1),
    pct_Tv                = round(mean(sub$ti_tv == "Tv", na.rm = TRUE) * 100, 1),
    mean_SVM_deep         = round(mean(sub$SVM_pos_score_deep,  na.rm = TRUE), 4),
    mean_SVM_low          = round(mean(sub$SVM_pos_score_low,   na.rm = TRUE), 4),
    mean_BAF_deep         = round(mean(sub$BAF_alt_deep, na.rm = TRUE), 4),
    mean_BAF_low          = round(mean(sub$BAF_alt_low,  na.rm = TRUE), 4)
  )
}

set_summary <- rbindlist(lapply(names(sets), function(sn)
  compute_set_summary(sets[[sn]], sn)), fill = TRUE)

fwrite(set_summary, file.path(out_dir, "per_set_summary.tsv"), sep = "\t")
log_info("Wrote per-set summary")
print(set_summary[, .(set_name, n, pct_cancer_gene, mean_VAF_deep, mean_VAF_low,
                       pct_Ti, median_Dep_deep, median_Dep_low)])

##############################################################################
# SECTION 7 – Cancer gene hit tables per set
##############################################################################

log_info("Writing cancer gene hit tables …")

cancer_hits_all <- master[is_cancer_gene == TRUE,
                           .(variant_id, chr, pos, ref, alt, symbol, genename,
                             cancer_categories, cross_category,
                             Dep_deep, VAF_deep, Dep_low, VAF_low,
                             QS_deep, QS_low, ti_tv, tri_context)]

fwrite(cancer_hits_all[order(symbol, chr, pos)],
       file.path(out_dir, "cancer_gene_variants_all.tsv"), sep = "\t")

# Per-set cancer gene hits
for (sname in names(sets)) {
  hits_s <- cancer_hits_all[variant_id %in% sets[[sname]]]
  if (nrow(hits_s)) {
    fwrite(hits_s, file.path(out_dir, sprintf("cancer_gene_%s.tsv", sname)), sep = "\t")
  }
}
log_info("Cancer gene variants in key sets:")
for (sname in c("deep_only_488B", "low_only_488B", "shared_488B",
                 "deep_only_489",  "low_only_489",  "shared_489",
                 "tissue_488B_only", "tissue_489_only")) {
  n <- sum(master$variant_id %in% sets[[sname]] & master$is_cancer_gene, na.rm = TRUE)
  log_info("  %s : %d cancer-gene variants", sname, n)
}

# TNBC category breakdown for tissue-specific variants
tnbc_tissue_488 <- master[variant_id %in% sets$tissue_488B_only & is_cancer_gene == TRUE,
                            .(symbol, cancer_categories, cross_category,
                              Dep_deep, VAF_deep, Dep_low, VAF_low)]
tnbc_tissue_489 <- master[variant_id %in% sets$tissue_489_only & is_cancer_gene == TRUE,
                            .(symbol, cancer_categories, cross_category,
                              Dep_deep, VAF_deep, Dep_low, VAF_low)]
fwrite(tnbc_tissue_488, file.path(out_dir, "tnbc_tissue_488B_specific.tsv"), sep = "\t")
fwrite(tnbc_tissue_489, file.path(out_dir, "tnbc_tissue_489_specific.tsv"),  sep = "\t")

##############################################################################
# SECTION 8 – Per-gene hit counts table
##############################################################################

gene_counts <- master[!is.na(symbol), .N, by = .(symbol, genename, is_cancer_gene)]
setorder(gene_counts, -N)
fwrite(gene_counts, file.path(out_dir, "per_gene_variant_counts.tsv"), sep = "\t")

# Per gene × set
gene_set_counts <- lapply(names(sets), function(sname) {
  sub <- master[variant_id %in% sets[[sname]] & !is.na(symbol)]
  if (!nrow(sub)) return(NULL)
  ct <- sub[, .N, by = .(symbol)]
  ct[, set_name := sname]
  ct
})
gene_set_dt <- rbindlist(Filter(Negate(is.null), gene_set_counts))
fwrite(gene_set_dt, file.path(out_dir, "per_gene_per_set_counts.tsv"), sep = "\t")

##############################################################################
# SECTION 9 – Plots
##############################################################################

log_info("Generating plots …")

# Use a single stable category per variant for plotting (cross_category can be multi-label).
pick_plot_category <- function(x, targets) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  hit <- intersect(targets, parts)
  if (length(hit)) hit[1] else NA_character_
}

# Theme
theme_pub <- function() {
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
}

# Palette for cross-category
cat_pal <- c(
  shared_488B    = "#1b7837",
  deep_only_488B = "#762a83",
  low_only_488B  = "#c2a5cf",
  shared_489     = "#d6604d",
  deep_only_489  = "#4393c3",
  low_only_489   = "#92c5de",
  tissue_488B_only = "#2c7bb6",
  tissue_489_only  = "#d7191c",
  tissue_both      = "#636363"
)

target_cats <- c("shared_488B", "deep_only_488B", "low_only_488B",
                 "shared_489",  "deep_only_489",  "low_only_489")

save_empty_plot <- function(path, title, message, width, height) {
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.6, label = title,
             size = 6, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.45, label = message,
             size = 4.2, colour = "grey30") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
  ggsave(path, p_empty, width = width, height = height, dpi = 150, bg = "white")
}

draw_two_set_venn <- function(tissue, deep_only, low_only, shared, out_file) {
  total_deep <- deep_only + shared
  total_low  <- low_only + shared
  theta <- seq(0, 2 * pi, length.out = 240L)
  radius <- 1.8

  circle_dt <- rbindlist(list(
    data.table(set = "deepseq", x = -1 + radius * cos(theta), y = radius * sin(theta)),
    data.table(set = "lowseq",  x =  1 + radius * cos(theta), y = radius * sin(theta))
  ))

  lbl <- data.table(
    x = c(-1.9, 0, 1.9),
    y = c(0, 0, 0),
    txt = c(sprintf("Deep only\n%d", deep_only),
            sprintf("Shared\n%d", shared),
            sprintf("Low only\n%d", low_only))
  )

  p <- ggplot() +
    geom_polygon(data = circle_dt,
                 aes(x = x, y = y, group = set, fill = set),
                 alpha = 0.35, colour = "grey30", linewidth = 0.6) +
    geom_text(data = lbl, aes(x = x, y = y, label = txt), size = 4.5, fontface = "bold") +
    annotate("text", x = -1, y = 2.25, label = sprintf("Deep total: %d", total_deep), size = 3.8) +
    annotate("text", x = 1,  y = 2.25, label = sprintf("Low total: %d", total_low),  size = 3.8) +
    scale_fill_manual(values = c(deepseq = "#4393c3", lowseq = "#d6604d"), name = NULL) +
    coord_equal(xlim = c(-3.3, 3.3), ylim = c(-2.4, 2.7), clip = "off") +
    labs(title = sprintf("Deepseq vs Lowseq Variant Overlap (%s)", tissue)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    )

  ggsave(out_file, p, width = 8, height = 6, dpi = 150, bg = "white")
}

## 9.0  Venn diagrams (deep-only vs shared vs low-only), by tissue ───────────
set_sizes_file <- file.path(out_dir, "variant_set_sizes.tsv")
set_sizes_for_venn <- if (file.exists(set_sizes_file)) {
  fread(set_sizes_file)
} else {
  data.table(set_name = names(sets), n_variants = vapply(sets, length, integer(1L)))
}

get_set_size <- function(nm) {
  v <- set_sizes_for_venn[set_name == nm, n_variants]
  if (length(v) && !is.na(v[1])) as.integer(v[1]) else as.integer(length(sets[[nm]]))
}

draw_two_set_venn(
  tissue = "488B",
  deep_only = get_set_size("deep_only_488B"),
  low_only  = get_set_size("low_only_488B"),
  shared    = get_set_size("shared_488B"),
  out_file  = file.path(plot_dir, "venn_deep_low_488B.png")
)

draw_two_set_venn(
  tissue = "489",
  deep_only = get_set_size("deep_only_489"),
  low_only  = get_set_size("low_only_489"),
  shared    = get_set_size("shared_489"),
  out_file  = file.path(plot_dir, "venn_deep_low_489.png")
)

## 9.1  VAF distribution: deepseq vs lowseq, by cross-category ──────────────
plot_vaf_dt <- melt(
  master[cross_category %in% names(cat_pal)],
  id.vars = c("variant_id", "cross_category"),
  measure.vars = c("VAF_deep", "VAF_low"),
  variable.name = "seq_type", value.name = "VAF"
)[!is.na(VAF)]

plot_vaf_dt[, seq_type := fifelse(seq_type == "VAF_deep", "deepseq", "lowseq")]

# Focus on cross-dataset categories for comparison
plot_vaf_dt[, plot_category := vapply(cross_category, pick_plot_category,
                                      character(1L), targets = target_cats)]
plot_vaf_sub <- plot_vaf_dt[cross_category %in%
  c("shared_488B", "deep_only_488B", "low_only_488B",
    "shared_489",  "deep_only_489",  "low_only_489")]

# If exact membership filtering yields empty data (multi-label cross_category),
# fall back to the derived single plotting category.
if (!nrow(plot_vaf_sub)) {
  plot_vaf_sub <- plot_vaf_dt[!is.na(plot_category)]
}

if (nrow(plot_vaf_sub)) {
  plot_vaf_sub[, facet_category := if ("plot_category" %in% names(plot_vaf_sub)) plot_category else cross_category]
  plot_vaf_sub[, VAF_plot := pmin(pmax(VAF, 0), 1)]
  p_vaf <- ggplot(plot_vaf_sub, aes(x = VAF_plot, fill = facet_category)) +
    geom_histogram(bins = 60L, alpha = 0.75, position = "identity") +
    facet_grid(facet_category ~ seq_type, scales = "free_y") +
    scale_fill_manual(values = cat_pal, breaks = target_cats, drop = FALSE, name = "Category") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = "VAF distribution by dataset-category",
         x = "Variant Allele Frequency", y = "Count") +
    theme_pub()

  ggsave(file.path(plot_dir, "vaf_by_category.png"),
         p_vaf, width = 10, height = 12, dpi = 150, bg = "white")
} else {
  log_info("Skipping vaf_by_category plot: no rows after category filtering")
  save_empty_plot(
    file.path(plot_dir, "vaf_by_category.png"),
    "VAF distribution by dataset-category",
    "No VAF values available after filtering",
    width = 10,
    height = 12
  )
}

## 9.2  Depth distribution ────────────────────────────────────────────────────
dep_dt <- melt(
  master[cross_category %in% names(cat_pal) & !is.na(cross_category)],
  id.vars = c("variant_id", "cross_category"),
  measure.vars = c("Dep_deep", "Dep_low"),
  variable.name = "seq_type", value.name = "Depth"
)[!is.na(Depth)]
dep_dt[, seq_type := fifelse(seq_type == "Dep_deep", "deepseq", "lowseq")]
dep_dt[, plot_category := vapply(cross_category, pick_plot_category,
                                 character(1L), targets = target_cats)]
dep_sub <- dep_dt[cross_category %in%
  c("shared_488B", "deep_only_488B", "low_only_488B",
    "shared_489",  "deep_only_489",  "low_only_489")]

if (!nrow(dep_sub)) {
  dep_sub <- dep_dt[!is.na(plot_category)]
}

if (nrow(dep_sub)) {
  dep_sub[, facet_category := if ("plot_category" %in% names(dep_sub)) plot_category else cross_category]
  p_dep <- ggplot(dep_sub, aes(x = pmin(Depth, 500L), fill = facet_category)) +
    geom_histogram(bins = 50L, alpha = 0.75, position = "identity") +
    facet_grid(facet_category ~ seq_type, scales = "free_y") +
    scale_fill_manual(values = cat_pal, breaks = target_cats, drop = FALSE, name = "Category") +
    labs(title = "Coverage depth (capped 500×) by dataset-category",
         x = "Depth", y = "Count") +
    theme_pub()

  ggsave(file.path(plot_dir, "depth_by_category.png"),
         p_dep, width = 10, height = 12, dpi = 150, bg = "white")
} else {
  log_info("Skipping depth_by_category plot: no rows after category filtering")
  save_empty_plot(
    file.path(plot_dir, "depth_by_category.png"),
    "Coverage depth by dataset-category",
    "No depth values available after filtering",
    width = 10,
    height = 12
  )
}

## 9.3  Ti/Tv bar chart ───────────────────────────────────────────────────────
titv_dt <- master[!is.na(ti_tv) & cross_category %in% names(cat_pal),
                   .N, by = .(cross_category, ti_tv)]
if (!nrow(titv_dt)) {
  titv_dt <- master[!is.na(ti_tv), .(ti_tv, cross_category)]
  titv_dt[, plot_category := vapply(cross_category, pick_plot_category,
                                    character(1L), targets = names(cat_pal))]
  titv_dt <- titv_dt[!is.na(plot_category), .N, by = .(cross_category = plot_category, ti_tv)]
}
titv_dt[, pct := N / sum(N) * 100, by = cross_category]

if (nrow(titv_dt)) {
  p_titv <- ggplot(titv_dt, aes(x = cross_category, y = pct, fill = ti_tv)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c(Ti = "#2166ac", Tv = "#d6604d"), name = NULL) +
    coord_flip() +
    labs(title = "Ti/Tv ratio per variant set",
         x = NULL, y = "Percentage (%)") +
    theme_pub()

  ggsave(file.path(plot_dir, "titv_by_category.png"),
         p_titv, width = 9, height = 6, dpi = 150)
} else {
  log_info("Skipping titv_by_category plot: ti_tv unavailable")
}

## 9.4  Quality score comparison (deepseq vs lowseq) ──────────────────────────
qs_dt <- melt(
  master[cross_category %in%
    c("shared_488B", "deep_only_488B", "low_only_488B",
      "shared_489",  "deep_only_489",  "low_only_489")],
  id.vars = c("variant_id", "cross_category"),
  measure.vars = c("QS_deep", "QS_low"),
  variable.name = "seq_type", value.name = "QS"
)[!is.na(QS)]
qs_dt[, seq_type := fifelse(seq_type == "QS_deep", "deepseq", "lowseq")]
qs_dt[, plot_category := vapply(cross_category, pick_plot_category,
                                character(1L), targets = target_cats)]
if (!nrow(qs_dt)) {
  qs_dt <- melt(
    master[!is.na(cross_category)],
    id.vars = c("variant_id", "cross_category"),
    measure.vars = c("QS_deep", "QS_low"),
    variable.name = "seq_type", value.name = "QS"
  ) [!is.na(QS)]
  qs_dt[, seq_type := fifelse(seq_type == "QS_deep", "deepseq", "lowseq")]
  qs_dt[, plot_category := vapply(cross_category, pick_plot_category,
                                  character(1L), targets = target_cats)]
  qs_dt <- qs_dt[!is.na(plot_category)]
}

if (nrow(qs_dt)) {
  qs_dt[, facet_category := if ("plot_category" %in% names(qs_dt)) plot_category else cross_category]
  qs_dt[, QS_plot := pmin(pmax(QS, 0), 1)]
  p_qs <- ggplot(qs_dt, aes(x = QS, fill = facet_category)) +
    geom_histogram(data = qs_dt, aes(x = QS_plot, fill = facet_category),
                   bins = 50L, alpha = 0.75, position = "identity") +
    facet_grid(facet_category ~ seq_type, scales = "free_y") +
    scale_fill_manual(values = cat_pal, breaks = target_cats, drop = FALSE, name = "Category") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = "Quality score (QS) distribution by category",
         x = "QS", y = "Count") +
    theme_pub()

  ggsave(file.path(plot_dir, "qs_by_category.png"),
         p_qs, width = 10, height = 12, dpi = 150, bg = "white")
} else {
  log_info("Skipping qs_by_category plot: no QS values available")
  save_empty_plot(
    file.path(plot_dir, "qs_by_category.png"),
    "Quality score (QS) distribution",
    "No QS values available after filtering",
    width = 10,
    height = 12
  )
}

## 9.5  Top cancer genes per tissue ────────────────────────────────────────────
plot_top_genes <- function(set_name, label, top_n = 20L) {
  sub <- master[variant_id %in% sets[[set_name]] &
                  is_cancer_gene == TRUE & !is.na(symbol)]
  if (!nrow(sub)) {
    log_info("No cancer gene variants in %s — skipping plot", set_name)
    return(invisible(NULL))
  }
  ct <- sub[, .N, by = symbol][order(-N)][seq_len(min(.N, top_n))]
  ct[, symbol := factor(symbol, levels = rev(symbol))]
  p <- ggplot(ct, aes(x = symbol, y = N, fill = symbol)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
    labs(title = sprintf("Top cancer-gene variants — %s", label),
         x = NULL, y = "Variants") +
    theme_pub()
  fname <- file.path(plot_dir, sprintf("top_cancer_genes_%s.png", set_name))
  ggsave(fname, p, width = 7, height = 6, dpi = 150)
  invisible(p)
}

plot_top_genes("tissue_488B_only",  "Tissue 488B only")
plot_top_genes("tissue_489_only",   "Tissue 489 only")
plot_top_genes("deep_only_488B",    "Deepseq-only 488B")
plot_top_genes("low_only_488B",     "Lowseq-only 488B")
plot_top_genes("deep_only_489",     "Deepseq-only 489")
plot_top_genes("low_only_489",      "Lowseq-only 489")
plot_top_genes("shared_488B",       "Shared 488B")
plot_top_genes("shared_489",        "Shared 489")

## 9.6  Scatter: VAF deepseq vs VAF lowseq for shared variants ─────────────────
for (tissue in c("488B", "489")) {
  sname <- sprintf("shared_%s", tissue)
  sub <- master[variant_id %in% sets[[sname]] &
                  !is.na(VAF_deep) & !is.na(VAF_low)]
  if (!nrow(sub)) next
  p_sc <- ggplot(sub, aes(x = VAF_deep, y = VAF_low,
                            colour = is_cancer_gene)) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#d62728"),
                        labels = c("other", "cancer gene"), name = NULL) +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = sprintf("VAF deepseq vs lowseq — shared %s variants", tissue),
         x = "VAF (deepseq)", y = "VAF (lowseq)") +
    theme_pub()
  ggsave(file.path(plot_dir, sprintf("vaf_scatter_shared_%s.png", tissue)),
         p_sc, width = 7, height = 6, dpi = 150)
}

## 9.7  Tissue comparison: shared vs 488B-only vs 489-only (cancer genes) ──────
tissue_cancer <- master[is_cancer_gene == TRUE & !is.na(symbol)]
tissue_cancer[, tissue_cat := fifelse(
  variant_id %in% sets$tissue_488B_only, "488B only",
  fifelse(variant_id %in% sets$tissue_489_only, "489 only", "both tissues")
)]

gene_tissue_ct <- tissue_cancer[, .N, by = .(symbol, tissue_cat)]
top_genes_tissue <- gene_tissue_ct[, .(total = sum(N)), by = symbol][order(-total)][seq_len(min(.N, 30L))]
gene_tissue_ct2  <- gene_tissue_ct[symbol %in% top_genes_tissue$symbol]
gene_tissue_ct2[, symbol := factor(symbol, levels = rev(top_genes_tissue$symbol))]

p_tissue_cg <- ggplot(gene_tissue_ct2, aes(x = symbol, y = N, fill = tissue_cat)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("488B only" = "#2c7bb6",
                                "489 only"  = "#d7191c",
                                "both tissues" = "#636363"),
                    name = "Tissue presence") +
  labs(title = "Cancer-gene variants by tissue presence (top 30 genes)",
       x = NULL, y = "Variant count") +
  theme_pub()

ggsave(file.path(plot_dir, "cancer_genes_by_tissue.png"),
  p_tissue_cg, width = 9, height = 9, dpi = 150, bg = "white")

## 9.8  SVM score comparison (deep-only vs shared vs low-only per tissue) ─────
svm_dt <- melt(
  master[cross_category %in%
    c("shared_488B", "deep_only_488B", "low_only_488B",
      "shared_489",  "deep_only_489",  "low_only_489")],
  id.vars = c("variant_id", "cross_category"),
  measure.vars = c("SVM_pos_score_deep", "SVM_pos_score_low"),
  variable.name = "seq_type", value.name = "SVM"
)[!is.na(SVM)]
svm_dt[, seq_type := fifelse(grepl("deep", seq_type), "deepseq", "lowseq")]
svm_dt[, plot_category := vapply(cross_category, pick_plot_category,
                                 character(1L), targets = target_cats)]
svm_dt <- svm_dt[!is.na(plot_category)]

if (nrow(svm_dt)) {
  p_svm <- ggplot(svm_dt, aes(x = plot_category, y = SVM, fill = plot_category)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
    facet_wrap(~ seq_type) +
    scale_fill_manual(values = cat_pal, breaks = target_cats, drop = FALSE, name = "Category") +
    coord_flip() +
    labs(title = "SVM score by variant category",
         x = NULL, y = "SVM positive score") +
    theme_pub()

  ggsave(file.path(plot_dir, "svm_score_by_category.png"),
         p_svm, width = 11, height = 7, dpi = 150, bg = "white")
} else {
  log_info("Skipping svm_score_by_category plot: no SVM values available")
  save_empty_plot(
    file.path(plot_dir, "svm_score_by_category.png"),
    "SVM score by variant category",
    "No SVM values available after filtering",
    width = 11,
    height = 7
  )
}

##############################################################################
# SECTION 10 – Final summary and interpretation notes
##############################################################################

interp_summary <- function(set_summary_dt, sets_list, master_dt) {
  lines <- character(0)
  lines <- c(lines, "=== SOMATIC SNV CHARACTERIZATION SUMMARY ===\n")

  # Coverage-based interpretation for "only" sets
  for (tissue in c("488B", "489")) {
    deep_nm <- sprintf("deep_only_%s", tissue)
    low_nm  <- sprintf("low_only_%s", tissue)
    shr_nm  <- sprintf("shared_%s",   tissue)

    d_row <- set_summary_dt[set_name == deep_nm]
    l_row <- set_summary_dt[set_name == low_nm]
    s_row <- set_summary_dt[set_name == shr_nm]

    lines <- c(lines, sprintf("\n--- Tissue %s ---", tissue))
    lines <- c(lines, sprintf(
      "  Shared           : %d variants | deep VAF %.3f | low VAF %.3f | deep depth %.0f | low depth %.0f",
      s_row$n, s_row$mean_VAF_deep, s_row$mean_VAF_low,
      s_row$median_Dep_deep, s_row$median_Dep_low))
    lines <- c(lines, sprintf(
      "  Deepseq-only     : %d variants | deep VAF %.3f | (low VAF %.3f) | deep depth %.0f",
      d_row$n, d_row$mean_VAF_deep, d_row$mean_VAF_low, d_row$median_Dep_deep))
    lines <- c(lines, sprintf(
      "  Lowseq-only      : %d variants | (deep VAF %.3f) | low VAF %.3f | low depth %.0f",
      l_row$n, l_row$mean_VAF_deep, l_row$mean_VAF_low, l_row$median_Dep_low))

    # Interpretation heuristic
    if (!is.na(d_row$mean_VAF_deep) && !is.na(s_row$mean_VAF_deep)) {
      if (d_row$mean_VAF_deep < s_row$mean_VAF_deep - 0.05) {
        lines <- c(lines, "  -> Deepseq-only variants tend to have LOWER VAF than shared variants.")
        lines <- c(lines, "     Likely explanation: higher depth in deepseq allows detection of")
        lines <- c(lines, "     low-VAF / subclonal variants invisible in lowseq.")
      }
    }
    if (!is.na(l_row$mean_VAF_low) && !is.na(s_row$mean_VAF_low)) {
      if (l_row$mean_VAF_low > s_row$mean_VAF_low + 0.05) {
        lines <- c(lines, "  -> Lowseq-only variants tend to have HIGHER VAF in lowseq.")
        lines <- c(lines, "     Possible explanation: low-coverage stochastic inflation of VAF,")
        lines <- c(lines, "     or barcode-level noise not filtered by the somatic caller.")
      }
    }
  }

  # Tissue identity
  n488 <- sum(master_dt$variant_id %in% sets_list$tissue_488B_only)
  n489 <- sum(master_dt$variant_id %in% sets_list$tissue_489_only)
  nboth <- sum(master_dt$variant_id %in% sets_list$tissue_both)
  lines <- c(lines, sprintf(
    "\n--- Tissue specificity ---\n  488B-only: %d | 489-only: %d | both: %d",
    n488, n489, nboth))
  if (n488 > 50L || n489 > 50L) {
    lines <- c(lines, "  -> Substantial tissue-specific variant pools detected.")
    lines <- c(lines, "     These may represent different clonal compositions in each tissue slice.")
  }

  lines <- c(lines, "\n--- TNBC / breast-cancer gene highlights ---")
  for (sname in c("shared_488B", "deep_only_488B", "low_only_488B",
                  "shared_489", "deep_only_489", "low_only_489",
                  "tissue_488B_only", "tissue_489_only")) {
    row <- set_summary_dt[set_name == sname]
    if (nrow(row)) {
      lines <- c(lines, sprintf(
        "  %s: %d / %d cancer-gene variants (%.2f%%)",
        sname,
        row$n_cancer_gene,
        row$n,
        row$pct_cancer_gene
      ))
    }
  }

  tnbc488 <- sum(master_dt$variant_id %in% sets_list$tissue_488B_only & master_dt$is_cancer_gene, na.rm = TRUE)
  tnbc489 <- sum(master_dt$variant_id %in% sets_list$tissue_489_only & master_dt$is_cancer_gene, na.rm = TRUE)
  lines <- c(lines, sprintf("  Tissue-specific TNBC genes: 488B=%d | 489=%d", tnbc488, tnbc489))

  paste(lines, collapse = "\n")
}

interp <- tryCatch(
  interp_summary(set_summary, sets, master),
  error = function(e) {
    msg <- sprintf("Interpretation summary failed: %s", conditionMessage(e))
    log_info("%s", msg)
    paste(
      "=== SOMATIC SNV CHARACTERIZATION SUMMARY ===",
      "",
      msg,
      "Check per_set_summary.tsv and master_annotated_variants.tsv for details.",
      sep = "\n"
    )
  }
)
cat(interp, "\n")
writeLines(interp, file.path(note_dir, "interpretation_notes.txt"))
log_info("Wrote interpretation notes: %s",
         file.path(note_dir, "interpretation_notes.txt"))

log_info("Done. All outputs in: %s", out_dir)
