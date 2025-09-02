# ---- Setup ---------------------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Bats")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("sangerseqR","DECIPHER","Biostrings")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = TRUE)
}

library(sangerseqR)
library(DECIPHER)
library(Biostrings)

# Sanity check: ensure AlignSeqs is available
stopifnot("AlignSeqs" %in% getNamespaceExports("DECIPHER"))

# ---- Parameters ----------------------------------------------------------
folder <- "data/cytb"         # Folder with .ab1 files
trim_len <- 170               # initial trim on raw reads
align_trim_start <- 23        # trim off first 22 columns (primer/poly-A)
final_trim_len <- 148         # final retained length

# trim_len <- 550          # keep more from each raw read (or even nchar(raw_seq))
# align_trim_start <- 23   # same primer offset as before
# final_trim_len <- 300    # use a 300-bp window (or 297 if you prefer codon-clean)

output_raw <- "trimmed_170.fasta"
output_final_aligned <- "cytb_trimmed_for_tree.fasta"
known_seq_file <- "cytb_known_batseqs.fasta"
final_combined_output <- "cytb_known_and_unknown_trimmed.fasta"

# ---- Step 1: Read & hard-trim .ab1 files --------------------------------
ab1_files <- list.files(folder, pattern = "\\.ab1$", full.names = TRUE)
if (length(ab1_files) == 0) stop("No .ab1 files found in ", normalizePath(folder))

seqs_out <- vector("list", length(ab1_files))
names(seqs_out) <- tools::file_path_sans_ext(basename(ab1_files))

for (i in seq_along(ab1_files)) {
  abif <- read.abif(ab1_files[i])
  sanger <- sangerseq(abif)
  raw_seq <- as.character(primarySeq(sanger))
  trimmed <- substr(raw_seq, 1, trim_len)
  seqs_out[[i]] <- DNAString(trimmed)
}

seqs_DNA <- DNAStringSet(seqs_out)
writeXStringSet(seqs_DNA, filepath = output_raw)

# ---- Step 2: Align unknown sequences ------------------------------------
aligned <- DECIPHER::AlignSeqs(seqs_DNA)  # explicit namespace
# ---- Step 3: Window-trim to primerless region ---------------------------
aligned_trimmed <- subseq(aligned, start = align_trim_start, width = final_trim_len)
writeXStringSet(aligned_trimmed, filepath = output_final_aligned)

# ---- Step 4: Read known sequences & combine -----------------------------
if (!file.exists(known_seq_file)) stop("Known sequence file not found: ", known_seq_file)
knowns <- readDNAStringSet(known_seq_file)

# Remove gaps before re-aligning everything together (use DECIPHER's helper)
unknowns <- aligned_trimmed
unknowns_clean <- DECIPHER::RemoveGaps(unknowns)
knowns_clean   <- DECIPHER::RemoveGaps(knowns)

combined <- c(unknowns_clean, knowns_clean)

# ---- Step 5: Final MSA & consistent trimming ----------------------------
combined_aligned <- DECIPHER::AlignSeqs(combined)
final_trimmed <- subseq(combined_aligned, start = align_trim_start, width = final_trim_len)
writeXStringSet(final_trimmed, filepath = final_combined_output)

# ---- Step 6: Visual check -----------------------------------------------
DECIPHER::BrowseSeqs(combined_aligned, highlight = 1)


# ───────────────────────────── Setup ─────────────────────────────
setwd("~/Library/CloudStorage/Box-Box/Bats")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("Biostrings","DECIPHER")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=TRUE)
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")

library(Biostrings)
library(DECIPHER)
library(ape)

unknowns_only_fa <- "cytb_trimmed_for_tree.fasta"      # unknowns, aligned & 148 bp trimmed
known_fa         <- "cytb_known_batseqs.fasta"         # your reference panel (raw sequences)
combined_out_fa  <- "cytb_known_and_unknown_aligned.fasta"

# ─────────────── 0) Read inputs & basic checks ───────────────
stopifnot(file.exists(unknowns_only_fa), file.exists(known_fa))

unknown_aln <- readDNAStringSet(unknowns_only_fa)
stopifnot(length(unknown_aln) > 0)

known_raw <- readDNAStringSet(known_fa)
stopifnot(length(known_raw) > 0)

cat("Unknowns:", length(unknown_aln), "  Knowns:", length(known_raw), "\n")
cat("Unknown alignment widths (should be ~148):", paste(unique(width(unknown_aln)), collapse=", "), "\n")

# Clean any stray gaps from known headers/content and align knowns alone
known_raw <- DECIPHER::RemoveGaps(known_raw)
known_aln <- DECIPHER::AlignSeqs(known_raw)  # aligns & orients knowns

# ─────────────── 1) Combine alignments ───────────────
# Prefer AlignProfiles to preserve the unknowns' existing 148-bp window
if ("AlignProfiles" %in% getNamespaceExports("DECIPHER")) {
  combined_aln <- DECIPHER::AlignProfiles(unknown_aln, known_aln)
} else {
  # Fallback: align everything together (window may shift, but downstream p-distance is gap-robust)
  combined_aln <- DECIPHER::AlignSeqs(c(unknown_aln, known_raw))
}

# Save the combined alignment
writeXStringSet(combined_aln, combined_out_fa)
cat("Wrote:", combined_out_fa, "\n")

# Label unknown vs known from the unknowns file names
unk_names <- names(unknown_aln)
is_unknown <- names(combined_aln) %in% unk_names
is_known <- !is_unknown
cat("In combined alignment -> Unknowns:", sum(is_unknown), "  Knowns:", sum(is_known), "\n")

# Also write a label map to sanity-check
write.csv(
  data.frame(name = names(combined_aln),
             label = ifelse(is_unknown, "unknown", "known")),
  "cytb_name_labels.csv",
  row.names = FALSE
)
cat("Wrote: cytb_name_labels.csv (verify labels)\n")

# ─────────────────────── 2) Sanity checks ───────────────────────
cat("\n## Combined alignment width:", unique(width(combined_aln)), "(may exceed 148 due to added gaps)\n")

# Ambiguity + NUMT check on gapless sequences
ungapped <- DECIPHER::RemoveGaps(combined_aln)

amb_letters <- c("N","R","Y","S","W","K","M","B","D","H","V")
amb_tbl <- as.data.frame(letterFrequency(ungapped, amb_letters))
amb_tbl$seqname <- names(ungapped)
amb_tbl$ambigs <- rowSums(amb_tbl[amb_letters])

# Vertebrate mito translation (robust to versions)
gc <- tryCatch(
  Biostrings::getGeneticCode("2"),  # NCBI table 2
  error = function(e) Biostrings::getGeneticCode("Vertebrate Mitochondrial", full.search = TRUE)
)
aa <- Biostrings::translate(ungapped, genetic.code = gc, if.fuzzy.codon = "solve")
has_stop <- vcountPattern("*", aa) > 0
names(has_stop) <- names(aa)
cat("Sequences with internal stops (should be none):", length(which(has_stop)), "\n")

# ───────────── 3) Pairwise % identity (gap/ambig tolerant) ─────────────
seq_names <- names(combined_aln)
seq_str   <- as.character(combined_aln)
L <- unique(nchar(seq_str)); stopifnot(length(L) == 1)

mat <- do.call(rbind, strsplit(seq_str, ""))
rownames(mat) <- seq_names
to_upper <- function(x) { x2 <- x; x2[!is.na(x2)] <- toupper(x2[!is.na(x2)]); x2 }

valid_mask <- function(i, j) {
  xi <- mat[i,]; xj <- mat[j,]
  # treat any non-gap character as valid (handles IUPAC codes)
  !(xi %in% c("-", ".", "?", " ")) &
    !(xj %in% c("-", ".", "?", " ")) &
    !is.na(xi) & !is.na(xj)
}

n <- nrow(mat)
pid <- matrix(NA_real_, n, n, dimnames = list(seq_names, seq_names))
valid_counts <- matrix(0L, n, n, dimnames = list(seq_names, seq_names))

for (i in seq_len(n)) {
  for (j in i:n) {
    v <- valid_mask(i, j)
    denom <- sum(v)
    valid_counts[i, j] <- valid_counts[j, i] <- denom
    if (denom > 0) {
      xi <- to_upper(mat[i, v]); xj <- to_upper(mat[j, v])
      matches <- sum(xi == xj)
      p <- matches / denom
      pid[i, j] <- pid[j, i] <- p * 100
    } else {
      pid[i, j] <- pid[j, i] <- NA_real_
    }
  }
}

dist_mat <- 1 - pid/100
diag(dist_mat) <- 0

cat("Min/median valid sites per pair:",
    suppressWarnings(min(valid_counts[upper.tri(valid_counts)], na.rm = TRUE)), "/",
    suppressWarnings(median(valid_counts[upper.tri(valid_counts)], na.rm = TRUE)), "\n")

# ─────────────── 4) Nearest known match per unknown ───────────────
unknown_ix <- which(is_unknown)
known_ix   <- which(is_known)
if (length(known_ix) == 0) stop("No sequences labeled as 'known'. Check cytb_name_labels.csv and your known FASTA.")

nearest_list <- lapply(unknown_ix, function(i) {
  d <- dist_mat[i, known_ix]
  if (all(is.na(d))) {
    return(data.frame(
      unknown = names(combined_aln)[i],
      best_known = NA, dist = NA, pct_identity = NA,
      second_best_dist = NA, second_best_pct_identity = NA,
      gap_to_second_pct = NA,
      ambigs = amb_tbl$ambigs[match(names(combined_aln)[i], amb_tbl$seqname)],
      has_stop = has_stop[match(names(combined_aln)[i], names(has_stop))]
    ))
  }
  ord <- order(d)
  j1 <- ord[1]
  best_d <- d[j1]
  second_d <- if (length(ord) >= 2) d[ord[2]] else NA
  
  best_pct   <- 100 * (1 - best_d)
  second_pct <- if (!is.na(second_d)) 100 * (1 - second_d) else NA
  gap_pct    <- if (!is.na(second_pct)) (best_pct - second_pct) else NA
  
  data.frame(
    unknown = names(combined_aln)[i],
    best_known = names(combined_aln)[known_ix[j1]],
    dist = as.numeric(best_d),
    pct_identity = round(best_pct, 2),
    second_best_dist = as.numeric(second_d),
    second_best_pct_identity = if (!is.na(second_pct)) round(second_pct, 2) else NA_real_,
    gap_to_second_pct = if (!is.na(gap_pct)) round(gap_pct, 2) else NA_real_,
    ambigs = amb_tbl$ambigs[match(names(combined_aln)[i], amb_tbl$seqname)],
    has_stop = has_stop[match(names(combined_aln)[i], names(has_stop))]
  )
})

nearest_tbl <- do.call(rbind, nearest_list)
nearest_tbl <- nearest_tbl[order(-nearest_tbl$pct_identity, nearest_tbl$dist), ]

nearest_tbl$call_strength <- with(nearest_tbl, ifelse(
  !is.na(pct_identity) & pct_identity >= 98 & (is.na(gap_to_second_pct) | gap_to_second_pct >= 1),
  "STRONG",
  ifelse(!is.na(pct_identity) & pct_identity >= 95, "SUGGESTIVE", "WEAK/AMBIGUOUS")
))

cat("\n## Top matches (head):\n")
print(head(nearest_tbl, 10))

write.csv(nearest_tbl, "cytb_unknown_best_matches.csv", row.names = FALSE)
cat("\nWrote: cytb_unknown_best_matches.csv\n")

# ───────────────────────── 5) Quick NJ tree ─────────────────────────
pdf("cytb_nj_tree.pdf", width = 8.5, height = 11)
tree <- njs(as.dist(dist_mat))
plot(tree, cex = 0.6)
unk_tips <- match(names(combined_aln)[unknown_ix], tree$tip.label)
tiplabels(pch = 19, tip = unk_tips, col = "red", cex = 0.6)
legend("topleft", legend = c("Unknown"), pch = 19, col = "red", bty = "n", cex = 0.8)
dev.off()
cat("Wrote: cytb_nj_tree.pdf\n")

# ───────────────────── 6) Visual alignment check ─────────────────────
try({
  DECIPHER::BrowseSeqs(combined_aln, highlight = which(is_unknown))
}, silent = TRUE)

cat("\nOpen 'cytb_name_labels.csv' to confirm labels, ",
    "'cytb_unknown_best_matches.csv' for IDs, and 'cytb_nj_tree.pdf' to visualize clustering.\n", sep = "")
#-------------------------
  
  
  
  library(ape)
library(Biostrings)

# Read nearest-match results
nearest <- read.csv("cytb_unknown_best_matches.csv", stringsAsFactors = FALSE)

# Extract "Genus species" from best_known headers like "MF135770.1 Tadarida brasiliensis voucher ..."
extract_species <- function(x) {
  m <- regexec("^[^ ]+\\s+([A-Z][a-z]+\\s+[a-z]+)", x)
  sapply(regmatches(x, m), function(v) if (length(v) >= 2) v[2] else NA_character_)
}
nearest$pred_species <- extract_species(nearest$best_known)

# Map unknown IDs to their indices in the combined alignment / tree
id_in_tree <- match(nearest$unknown, names(combined_aln))

# 1) Relabel tips with predicted species and % identity
new_labels <- tree$tip.label
ok <- !is.na(id_in_tree)
new_labels[id_in_tree[ok]] <- paste0(
  new_labels[id_in_tree[ok]], " [",
  nearest$pred_species[ok], " ",
  nearest$pct_identity[ok], "%]"
)
tree_annot <- tree
tree_annot$tip.label <- new_labels

pdf("cytb_nj_tree_annotated.pdf", width = 8.5, height = 11)
plot(tree_annot, cex = 0.55)
title("NJ tree with unknowns relabeled by predicted species (%)")
dev.off()

# 2) Color unknowns by predicted species
sp_levels <- sort(unique(na.omit(nearest$pred_species)))
cols <- setNames(rainbow(length(sp_levels)), sp_levels)   # pick any palette you like
tip_cols <- rep("black", length(tree$tip.label))
tip_cols[id_in_tree[ok]] <- cols[nearest$pred_species[ok]]

pdf("cytb_nj_tree_colored.pdf", width = 8.5, height = 11)
plot(tree, tip.color = tip_cols, cex = 0.55)
legend("topleft", legend = names(cols), col = cols, pch = 19, cex = 0.6, bty = "n")
title("Unknowns colored by predicted species")
dev.off()

# 3) Representatives-only tree (keep all references + best unknown per species)
rep_ix <- tapply(seq_len(nrow(nearest))[ok], nearest$pred_species[ok],
                 function(ii) ii[ which.max(nearest$pct_identity[ii]) ])
rep_ix <- unlist(rep_ix, use.names = FALSE)

keep_names <- c(
  names(combined_aln)[is_known],
  nearest$unknown[rep_ix]
)

tree_small <- keep.tip(tree, keep_names)

pdf("cytb_nj_tree_representatives.pdf", width = 8.5, height = 11)
plot(tree_small, cex = 0.7)
title("References + one best unknown per predicted species")
dev.off()

