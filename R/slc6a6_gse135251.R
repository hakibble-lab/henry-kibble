# SLC6A6 (taurine transporter, TauT) expression in GSE135251
#
# GSE135251 (Govaere et al. 2020, Sci Transl Med): bulk liver RNA-seq from
# 206 biopsy-staged NAFLD patients + 10 controls. Samples are grouped as
# control, NAFL, NASH_F0-F1, NASH_F2, NASH_F3, NASH_F4, with fibrosis stage
# and NAFLD activity score (NAS) per sample.
#
# What this does:
#   1. Downloads the series metadata and per-sample raw counts from GEO
#   2. Builds a counts matrix and runs DESeq2 (~ group)
#   3. Reports SLC6A6 log2 fold changes for each disease group vs control
#   4. Tests SLC6A6 against fibrosis stage and NAS (Spearman, NAFLD only)
#   5. Writes tables and plots to output/GSE135251/
#
# One-time setup:
#   install.packages(c("BiocManager", "ggplot2"))
#   BiocManager::install(c("GEOquery", "DESeq2"))
#
# Run from the project root: source("R/slc6a6_gse135251.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
})

gse_id       <- "GSE135251"
gene_symbol  <- "SLC6A6"
gene_ensembl <- "ENSG00000131389"

data_dir <- file.path("data", gse_id)
out_dir  <- file.path("output", gse_id)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The raw counts tarball is large enough to exceed R's default 60 s timeout
options(timeout = max(3600, getOption("timeout")))

# ---- Sample metadata --------------------------------------------------------

gse <- getGEO(gse_id, destdir = data_dir, GSEMatrix = TRUE)[[1]]
pd  <- pData(gse)

# Find a characteristics column by pattern; fail loudly if GEO naming differs
pick_col <- function(patterns, required = TRUE) {
  for (p in patterns) {
    hit <- grep(p, colnames(pd), ignore.case = TRUE, value = TRUE)
    hit <- hit[grepl(":ch1$", hit)]
    if (length(hit) > 0) return(hit[1])
  }
  if (required) {
    stop("No column matching ", paste(patterns, collapse = " / "),
         ". Available: ", paste(grep(":ch1$", colnames(pd), value = TRUE),
                                collapse = ", "))
  }
  NA_character_
}

group_col    <- pick_col(c("group in paper", "^disease", "stage"))
fibrosis_col <- pick_col("fibrosis", required = FALSE)
nas_col      <- pick_col(c("\\bnas\\b", "activity"), required = FALSE)

meta <- data.frame(
  sample   = pd$geo_accession,
  group    = trimws(pd[[group_col]]),
  fibrosis = if (!is.na(fibrosis_col)) suppressWarnings(as.numeric(pd[[fibrosis_col]])) else NA,
  nas      = if (!is.na(nas_col)) suppressWarnings(as.numeric(pd[[nas_col]])) else NA,
  stringsAsFactors = FALSE
)

# Order groups by disease severity, control first
paper_levels <- c("control", "NAFL", "NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4")
present      <- unique(meta$group)
ordered      <- paper_levels[tolower(paper_levels) %in% tolower(present)]
ordered      <- present[match(tolower(ordered), tolower(present))]
levels_all   <- c(ordered, sort(setdiff(present, ordered)))
ref_level    <- levels_all[grepl("control", levels_all, ignore.case = TRUE)][1]
if (is.na(ref_level)) stop("Could not find a control group in: ", paste(present, collapse = ", "))
levels_all   <- c(ref_level, setdiff(levels_all, ref_level))
meta$group   <- factor(meta$group, levels = levels_all)
meta$nafld   <- meta$group != ref_level

message("Groups (", group_col, "):")
print(table(meta$group))

# ---- Raw counts -------------------------------------------------------------

supp <- getGEOSuppFiles(gse_id, baseDir = "data", fetch_files = FALSE)
raw_tar_name <- grep("_RAW\\.tar$", supp$fname, value = TRUE)
if (length(raw_tar_name) == 0) {
  stop("No _RAW.tar in supplementary files: ", paste(supp$fname, collapse = ", "))
}
raw_tar <- file.path(data_dir, raw_tar_name)
if (!file.exists(raw_tar)) {
  getGEOSuppFiles(gse_id, baseDir = "data", filter_regex = "_RAW\\.tar$")
}

raw_dir <- file.path(data_dir, "raw")
if (!dir.exists(raw_dir)) untar(raw_tar, exdir = raw_dir)

count_files <- list.files(raw_dir, pattern = "^GSM\\d+", full.names = TRUE)
if (length(count_files) == 0) stop("No GSM count files found in ", raw_dir)

read_counts <- function(path) {
  x <- read.delim(path, header = FALSE, stringsAsFactors = FALSE,
                  colClasses = c("character", "character"))
  # Drop a header row if present, and HTSeq summary rows (__no_feature etc.)
  x <- x[!is.na(suppressWarnings(as.numeric(x[[2]]))), ]
  x <- x[!startsWith(x[[1]], "__"), ]
  setNames(as.numeric(x[[2]]), sub("\\.\\d+$", "", x[[1]]))
}

counts_list <- lapply(count_files, read_counts)
names(counts_list) <- regmatches(basename(count_files),
                                 regexpr("GSM\\d+", basename(count_files)))
genes  <- Reduce(intersect, lapply(counts_list, names))
counts <- do.call(cbind, lapply(counts_list, `[`, genes))
rownames(counts) <- genes

shared <- intersect(meta$sample, colnames(counts))
if (length(shared) < nrow(meta)) {
  warning(nrow(meta) - length(shared), " samples in metadata have no counts file")
}
meta   <- meta[match(shared, meta$sample), ]
counts <- round(counts[, shared])
rownames(meta) <- meta$sample

gene_row <- if (gene_ensembl %in% rownames(counts)) gene_ensembl else gene_symbol
if (!gene_row %in% rownames(counts)) {
  stop(gene_symbol, " (", gene_ensembl, ") not found in counts row names, e.g. ",
       paste(head(rownames(counts)), collapse = ", "))
}

# ---- DESeq2 -----------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
dds <- DESeq(dds)

de <- do.call(rbind, lapply(setdiff(levels(meta$group), ref_level), function(lvl) {
  r <- results(dds, contrast = c("group", lvl, ref_level))
  data.frame(gene = gene_symbol, comparison = paste(lvl, "vs", ref_level),
             baseMean = r[gene_row, "baseMean"],
             log2FoldChange = r[gene_row, "log2FoldChange"],
             lfcSE = r[gene_row, "lfcSE"],
             pvalue = r[gene_row, "pvalue"],
             padj_genomewide = r[gene_row, "padj"])
}))
write.csv(de, file.path(out_dir, "slc6a6_deseq2_vs_control.csv"), row.names = FALSE)
message("\nDESeq2: ", gene_symbol, " vs ", ref_level)
print(de, digits = 3)

# ---- Per-sample expression --------------------------------------------------

vsd <- vst(dds, blind = FALSE)
meta$vst       <- assay(vsd)[gene_row, meta$sample]
meta$norm_count <- counts(dds, normalized = TRUE)[gene_row, meta$sample]
write.csv(meta, file.path(out_dir, "slc6a6_per_sample.csv"), row.names = FALSE)

summary_tbl <- do.call(rbind, lapply(split(meta, meta$group), function(d) {
  data.frame(group = d$group[1], n = nrow(d),
             median_norm_count = median(d$norm_count),
             mean_vst = mean(d$vst), sd_vst = sd(d$vst))
}))
write.csv(summary_tbl, file.path(out_dir, "slc6a6_group_summary.csv"), row.names = FALSE)
print(summary_tbl, digits = 3)

# ---- Association with histology (NAFLD patients only) -----------------------

nafld <- meta[meta$nafld, ]
assoc <- do.call(rbind, lapply(c("fibrosis", "nas"), function(v) {
  ok <- !is.na(nafld[[v]])
  if (sum(ok) < 10) return(NULL)
  ct <- cor.test(nafld$vst[ok], nafld[[v]][ok], method = "spearman", exact = FALSE)
  data.frame(gene = gene_symbol, variable = v, n = sum(ok),
             spearman_rho = unname(ct$estimate), pvalue = ct$p.value)
}))
kw <- kruskal.test(vst ~ group, data = meta)
assoc <- rbind(assoc, data.frame(gene = gene_symbol, variable = "group (Kruskal-Wallis)",
                                 n = nrow(meta), spearman_rho = NA, pvalue = kw$p.value))
write.csv(assoc, file.path(out_dir, "slc6a6_histology_association.csv"), row.names = FALSE)
message("\nAssociation with histology:")
print(assoc, digits = 3)

# ---- Plots ------------------------------------------------------------------

ylab <- paste(gene_symbol, "expression (VST)")

p_group <- ggplot(meta, aes(group, vst)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(x = NULL, y = ylab, title = paste(gene_symbol, "in", gse_id),
       subtitle = sprintf("Kruskal-Wallis p = %.2g", kw$p.value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "slc6a6_by_group.png"), p_group, width = 6, height = 4.5, dpi = 300)

if (any(!is.na(nafld$fibrosis))) {
  p_fib <- ggplot(nafld[!is.na(nafld$fibrosis), ], aes(factor(fibrosis), vst)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    labs(x = "Fibrosis stage", y = ylab, title = paste(gene_symbol, "vs fibrosis (NAFLD)")) +
    theme_bw()
  ggsave(file.path(out_dir, "slc6a6_by_fibrosis.png"), p_fib, width = 5, height = 4.5, dpi = 300)
}

if (any(!is.na(nafld$nas))) {
  p_nas <- ggplot(nafld[!is.na(nafld$nas), ], aes(nas, vst)) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black") +
    labs(x = "NAFLD activity score", y = ylab, title = paste(gene_symbol, "vs NAS (NAFLD)")) +
    theme_bw()
  ggsave(file.path(out_dir, "slc6a6_vs_nas.png"), p_nas, width = 5, height = 4.5, dpi = 300)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("\nDone. Results in ", out_dir)
