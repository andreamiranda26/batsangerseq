# ==================== BAT CYTB — 148-bp & ~300-bp (NJ + bootstrap) ====================
# Robust: align unknowns + refs together with AlignSeqs(); strip all '-' gaps first
# Outputs:
#   cytb_unknown_best_matches148.csv
#   cytb_unknown_best_matches300.csv
#   cytb_tree_full_148.{pdf,png}
#   cytb_tree_full_300.{pdf,png}
#   plots_margin_hist.{pdf,png}
#   plots_identity_vs_margin.{pdf,png}
# ================================================================================

suppressPackageStartupMessages({
  library(sangerseqR)
  library(DECIPHER)
  library(Biostrings)
  library(ape)
  library(phytools)
  library(ggplot2)
  library(dplyr)
})

# ---- 1) Parameters ------------------------------------------------------------
folder           <- "data/cytb"                 # .ab1 files
trim_len_raw     <- 170                         # for 148-bp
align_trim_start <- 23                          # drop first 22 bases
len_mini         <- 148                         # 148-bp window

trim_len_long    <- 550                         # for ~300-bp
len_frame        <- 300                         # ~300-bp window

known_seq_file   <- "cytb_known_batseqs.fasta"  # references (multi-FASTA; unaligned OK)
outgroup_fa      <- "cytb_outgroup.fasta"       # single sequence (optional)
outgroup_hint    <- "^Glossophaga soricina$"    # for rooting if present

set.seed(1)

# ---- 2) Helpers ---------------------------------------------------------------
remove_gaps <- function(dna) DNAStringSet(gsub("-", "", as.character(dna), fixed = TRUE))

extract_species <- function(x) {
  sapply(x, function(xx) {
    m1 <- regexec("^[^ ]+\\s+([A-Z][a-z]+\\s+[a-z]+)", xx); g1 <- regmatches(xx, m1)[[1]]
    if (length(g1) >= 2) return(g1[2])
    m2 <- regexec("^([A-Z][a-z]+\\s+[a-z]+)", xx); g2 <- regmatches(xx, m2)[[1]]
    if (length(g2) >= 2) return(g2[2])
    NA_character_
  })
}

SCI_TO_COMMON <- c(
  "Eptesicus fuscus"          = "Big brown bat",
  "Lasionycteris noctivagans" = "Silver-haired bat",
  "Lasiurus borealis"         = "Eastern red bat",
  "Lasiurus cinereus"         = "Hoary bat",
  "Lasiurus intermedius"      = "Northern yellow bat",
  "Lasiurus seminolus"        = "Seminole bat",
  "Myotis austroriparius"     = "Southeastern myotis",
  "Myotis grisescens"         = "Gray bat",
  "Myotis lucifugus"          = "Little brown bat",
  "Myotis septentrionalis"    = "Northern long-eared bat",
  "Nycticeius humeralis"      = "Evening bat",
  # Historical vs current genus — both map to the same common name:
  "Pipistrellus subflavus"    = "Tri-colored bat",
  "Perimyotis subflavus"      = "Tri-colored bat",
  "Tadarida brasiliensis"     = "Mexican free-tailed bat",
  "Glossophaga soricina"      = "Pallas’s long-tongued bat"
)
to_common <- function(species_vec) {
  out <- unname(SCI_TO_COMMON[species_vec]); out[is.na(out)] <- species_vec; out
}

pid_matrix_from_alignment <- function(aln) {
  bases <- c("A","C","G","T")
  nm <- names(aln)
  M <- do.call(rbind, strsplit(as.character(aln), ""))
  M <- toupper(M); rownames(M) <- nm
  n <- nrow(M)
  pid <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  valid_counts <- matrix(0L, n, n, dimnames = list(nm, nm))
  if (n >= 2) for (i in 1:(n-1)) for (j in (i+1):n) {
    mask <- (M[i,] %in% bases) & (M[j,] %in% bases)
    comp <- sum(mask); valid_counts[i,j] <- valid_counts[j,i] <- comp
    if (comp > 0) pid[i,j] <- pid[j,i] <- 100 * sum(M[i,mask]==M[j,mask]) / comp
  }
  diag(pid) <- 100; diag(valid_counts) <- ncol(M)
  list(pid=pid, valid=valid_counts)
}

call_from_pid_vec <- function(x){
  out <- rep("WEAK/AMBIGUOUS", length(x))
  out[is.na(x)] <- "UNRESOLVED"
  out[!is.na(x) & x >= 80] <- "SUGGESTIVE"
  out[!is.na(x) & x >= 90] <- "STRONG"
  out
}

nearest_table <- function(pid, valid, unknown_ids, known_ids) {
  seq_names <- rownames(pid)
  dist_mat <- 1 - pid/100; dist_mat[is.na(dist_mat)] <- 1; diag(dist_mat) <- 0
  unknown_ids <- intersect(unknown_ids, seq_names)
  known_ids   <- intersect(known_ids,   seq_names)
  do.call(rbind, lapply(unknown_ids, function(u){
    drow <- dist_mat[u, known_ids]; vrow <- valid[u, known_ids]
    ok <- !is.na(drow) & vrow > 0
    if (!any(ok)) return(data.frame(unknown=u,best_known=NA,dist=NA,pct_identity=NA,
                                    second_best_dist=NA,second_best_pct_identity=NA,
                                    gap_to_second_pct=NA, stringsAsFactors=FALSE))
    ord <- order(drow[ok]); k_ok <- names(drow[ok])
    best <- k_ok[ord[1]]; d1 <- drow[ok][ord[1]]; pi1 <- (1-d1)*100
    if (length(ord) >= 2) { d2 <- drow[ok][ord[2]]; pi2 <- (1-d2)*100; gap <- pi1 - pi2
    } else { d2 <- NA; pi2 <- NA; gap <- NA }
    data.frame(unknown=u, best_known=best,
               dist=round(d1,5), pct_identity=round(pi1,2),
               second_best_dist=ifelse(is.na(d2),NA,round(d2,5)),
               second_best_pct_identity=ifelse(is.na(pi2),NA,round(pi2,2)),
               gap_to_second_pct=ifelse(is.na(gap),NA,round(gap,2)),
               stringsAsFactors=FALSE)
  }))
}

bootstrap_nj_on_alignment <- function(aln, tree, B = 1000, seed = 1L) {
  set.seed(seed)
  M <- toupper(do.call(rbind, strsplit(as.character(aln), "")))
  rownames(M) <- names(aln)
  bases <- c("A","C","G","T")
  build_tree_boot <- function(cols) {
    X <- M[, cols, drop = FALSE]
    n  <- nrow(X); rn <- rownames(X)
    pid <- matrix(NA_real_, n, n, dimnames = list(rn, rn))
    for (i in 1:(n-1)) for (j in (i+1):n) {
      mask <- (X[i,] %in% bases) & (X[j,] %in% bases)
      comp <- sum(mask)
      if (comp > 0) pid[i,j] <- pid[j,i] <- 100 * sum(X[i,mask]==X[j,mask]) / comp
    }
    diag(pid) <- 100
    D <- 1 - pid/100; D[is.na(D)] <- 1; diag(D) <- 0
    ape::njs(as.dist(D))
  }
  boot_fun <- function(dat) {
    cols <- sample(seq_len(ncol(M)), replace = TRUE)
    build_tree_boot(cols)
  }
  ape::boot.phylo(tree, M, FUN = boot_fun, B = B)
}

plot_tree_common <- function(tree, is_unknown, out_stub, refs_common = NULL) {
  tr <- tree
  labs <- tr$tip.label
  if (!is.null(refs_common)) {
    labs_new <- labs
    labs_new[!is_unknown] <- refs_common[match(labs[!is_unknown], names(refs_common))]
    labs_new[is.na(labs_new)] <- labs[is.na(labs_new)]
    tr$tip.label <- labs_new
  }
  H <- max(ape::node.depth.edgelength(tr)); extend <- 1.18
  for (dev in c("pdf","png")) {
    if (dev=="pdf") pdf(paste0(out_stub,".pdf"), width=12, height=16, family="Helvetica")
    else             png(paste0(out_stub,".png"),  width=2200, height=3000, res=300)
    par(mar=c(5,4,2,1))
    plot(tr, type="phylogram", direction="rightwards",
         tip.color=ifelse(is_unknown, "black", "grey35"),
         edge.color="grey60", cex=max(0.45, min(0.9, 95/length(tr$tip.label))),
         label.offset=0.006, no.margin=FALSE,
         x.lim=c(0, H*extend))
    bs <- attr(tr, "boot")
    if (!is.null(bs)) {
      bs_show <- ifelse(!is.na(bs) & bs >= 70, as.character(round(bs)), "")
      try(nodelabels(text = bs_show, frame = "none", cex = 0.6, col = "grey30"), silent = TRUE)
    }
    legend("topleft", bty="n", cex=0.9,
           legend=c("Unknown sample","Reference"),
           pch=c(19,1), col=c("black","grey40"))
    dev.off()
  }
}

# ---- 3) References ------------------------------------------------------------
stopifnot("Known FASTA not found" = file.exists(known_seq_file))
knowns_raw <- readDNAStringSet(known_seq_file)
knowns_sel <- remove_gaps(knowns_raw)  # ensure no gaps up front

wanted_acc <- c(
  "AF376835.1","MF038479.1","EU350026.1","EU350025.1","KC747682.1",
  "KP341708.1","KP341709.1","KP341731.1","KP341713.1","KC747687.1",
  "KP341748.1","KP341753.1","KP341751.1","AM261885.1","AM261892.1",
  "OM160895.1","OM160889.1","DQ503551.1","AM262335.1",
  "OP157144.1","KC747697.1","AJ504449.1","MF135779.1","MF135770.1"
)
all_acc <- sub("\\s.*$", "", names(knowns_sel))
keep_idx <- match(wanted_acc[wanted_acc %in% all_acc], all_acc)
knowns_sel <- knowns_sel[keep_idx]
drop_ix <- grep("(^KC747677\\.1\\b)|Choeronycteris", names(knowns_sel), ignore.case = TRUE)
if (length(drop_ix)) knowns_sel <- knowns_sel[-drop_ix]

# Map reference label -> species (for common names)
ref_species_all <- extract_species(names(knowns_sel))
refs_common_map <- setNames(to_common(ref_species_all), names(knowns_sel))

# ---- 4) Unknowns: make 148-bp and ~300-bp windows --------------------------------
ab1_files <- list.files(folder, pattern="\\.ab1$", full.names=TRUE)
stopifnot("No .ab1 files found" = length(ab1_files) > 0)

make_window <- function(files, hard_trim, start, width) {
  lst <- vector("list", length(files))
  nm  <- tools::file_path_sans_ext(basename(files))
  for (i in seq_along(files)) {
    abif <- read.abif(files[i]); s <- sangerseq(abif)
    raw  <- toupper(as.character(primarySeq(s)))
    lst[[i]] <- DNAString(substr(raw, 1, min(hard_trim, nchar(raw))))
  }
  unk_raw <- DNAStringSet(lst); names(unk_raw) <- nm
  # Align unknowns to stabilize insertions/deletions, then crop window;
  # this step introduces gaps, which we'll remove later before the joint alignment.
  unk_aln <- DECIPHER::AlignSeqs(unk_raw, verbose = FALSE)
  Biostrings::subseq(unk_aln, start = start, width = width)
}

unknowns_148 <- make_window(ab1_files, trim_len_raw,  align_trim_start, len_mini)
unknowns_300 <- make_window(ab1_files, trim_len_long, align_trim_start, len_frame)

# ---- 5) STRIP GAPS before joint alignment (this fixes your error) --------------
unknowns_148_ng <- remove_gaps(unknowns_148)
unknowns_300_ng <- remove_gaps(unknowns_300)
knowns_sel_ng   <- remove_gaps(knowns_sel)   # redundant but safe

# ---- 6) Build joint alignments (unknowns + refs) --------------------------------
combined_148 <- DECIPHER::AlignSeqs(c(unknowns_148_ng, knowns_sel_ng), verbose = FALSE)
combined_300 <- DECIPHER::AlignSeqs(c(unknowns_300_ng, knowns_sel_ng), verbose = FALSE)

unk_names_148 <- names(unknowns_148_ng)
ref_names_148 <- names(knowns_sel_ng)
unk_names_300 <- names(unknowns_300_ng)
ref_names_300 <- names(knowns_sel_ng)

# ---- 7) Nearest matches (CSV) ---------------------------------------------------
pm_148 <- pid_matrix_from_alignment(combined_148)
nearest_148 <- nearest_table(pm_148$pid, pm_148$valid,
                             unknown_ids = unk_names_148,
                             known_ids   = ref_names_148)
nearest_148$pred_species  <- extract_species(nearest_148$best_known)
nearest_148$call_strength <- call_from_pid_vec(nearest_148$pct_identity)
write.csv(nearest_148, "cytb_unknown_best_matches148.csv", row.names = FALSE)

pm_300 <- pid_matrix_from_alignment(combined_300)
nearest_300 <- nearest_table(pm_300$pid, pm_300$valid,
                             unknown_ids = unk_names_300,
                             known_ids   = ref_names_300)
nearest_300$pred_species  <- extract_species(nearest_300$best_known)
nearest_300$call_strength <- call_from_pid_vec(nearest_300$pct_identity)
write.csv(nearest_300, "cytb_unknown_best_matches300.csv", row.names = FALSE)

# ---- 8) NJ trees + bootstrap (optional outgroup) --------------------------------
build_tree_with_optional_outgroup <- function(combined_aln, out_stub){
  aln_for_tree <- combined_aln
  if (file.exists(outgroup_fa)) {
    og <- readDNAStringSet(outgroup_fa)
    if (length(og) >= 1) {
      names(og)[1] <- {
        nm <- names(og)[1]
        nm2 <- sub("^.*([A-Z][a-z]+)[ _]([a-z]+).*$", "\\1 \\2", nm)
        if (!grepl("^[A-Z][a-z]+\\s+[a-z]+$", nm2)) "Glossophaga soricina" else nm2
      }
      og <- remove_gaps(og)
      # Realign everything together for plotting
      aln_for_tree <- DECIPHER::AlignSeqs(c(remove_gaps(aln_for_tree), og), verbose = FALSE)
    }
  }
  pm <- pid_matrix_from_alignment(aln_for_tree)
  D  <- 1 - pm$pid/100; D[is.na(D)] <- 1; diag(D) <- 0
  tr <- ape::njs(as.dist(D))
  if (any(grepl(outgroup_hint, tr$tip.label))) {
    hit <- tr$tip.label[grepl(outgroup_hint, tr$tip.label)][1]
    tr <- tryCatch(ape::root(tr, outgroup=hit, resolve.root=TRUE), error=function(e) tr)
  }
  tr <- ape::ladderize(tr, right = TRUE)
  bs <- bootstrap_nj_on_alignment(aln_for_tree, tr, B = 1000)
  attr(tr, "boot") <- bs
  list(tree = tr)
}

# 148 tree
tree148 <- build_tree_with_optional_outgroup(combined_148, out_stub = "cytb_tree_full_148")$tree
is_unknown_plot_148 <- tree148$tip.label %in% unk_names_148
plot_tree_common(tree148, is_unknown_plot_148, out_stub = "cytb_tree_full_148",
                 refs_common = setNames(to_common(extract_species(names(knowns_sel_ng))), names(knowns_sel_ng)))

# 300 tree
tree300 <- build_tree_with_optional_outgroup(combined_300, out_stub = "cytb_tree_full_300")$tree
is_unknown_plot_300 <- tree300$tip.label %in% unk_names_300
plot_tree_common(tree300, is_unknown_plot_300, out_stub = "cytb_tree_full_300",
                 refs_common = setNames(to_common(extract_species(names(knowns_sel_ng))), names(knowns_sel_ng)))

# ---- 9) Identity vs margin plots ----------------------------------------------
d148 <- nearest_148 %>%
  transmute(sample = unknown,
            frame  = "148-bp mini-barcode",
            pct_identity,
            margin = gap_to_second_pct,
            tier   = call_from_pid_vec(pct_identity))

d300 <- nearest_300 %>%
  transmute(sample = unknown,
            frame  = "~300-bp window",
            pct_identity,
            margin = gap_to_second_pct,
            tier   = call_from_pid_vec(pct_identity))

dd <- bind_rows(d148, d300)

p_hist <- ggplot(dd %>% filter(!is.na(margin)), aes(x = margin)) +
  geom_histogram(bins = 25, color = "white") +
  facet_wrap(~ frame, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 1.0, linetype = "dotted") +
  labs(title = "Margin (best − second-best) distributions",
       x = "Margin to second-best (percentage points)",
       y = "Count") +
  theme_classic(base_size = 12)
ggsave("plots_margin_hist.pdf", p_hist, width = 8.5, height = 4.4)
ggsave("plots_margin_hist.png", p_hist, width = 8.5, height = 4.4, dpi = 300)

p_scatter <- ggplot(dd %>% filter(!is.na(margin)),
                    aes(x = pct_identity, y = margin, color = tier)) +
  geom_point(size = 2, alpha = 0.9) +
  facet_wrap(~ frame, ncol = 2) +
  scale_color_manual(values = c(
    "STRONG" = "#1b9e77",
    "SUGGESTIVE" = "#7570b3",
    "WEAK/AMBIGUOUS" = "#d95f02",
    "UNRESOLVED" = "grey50"
  ), drop = FALSE) +
  geom_vline(xintercept = 80, linetype = "dashed") +
  geom_vline(xintercept = 90, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 1.0, linetype = "dotted") +
  labs(title = "Percent identity vs margin",
       x = "Percent identity (%)",
       y = "Margin to second-best (percentage points)",
       color = "Confidence") +
  theme_classic(base_size = 12)
ggsave("plots_identity_vs_margin.pdf", p_scatter, width = 9.5, height = 4.8)
ggsave("plots_identity_vs_margin.png", p_scatter, width = 9.5, height = 4.8, dpi = 300)

# ==================== Paired per-sample comparison: 148 vs ~300 ====================

# Helper to get an object if present, else load from CSV
get_df <- function(obj_name, csv_path) {
  if (exists(obj_name, inherits = TRUE)) {
    get(obj_name, inherits = TRUE)
  } else {
    stopifnot("Missing file: {csv_path}" = file.exists(csv_path))
    read.csv(csv_path, stringsAsFactors = FALSE)
  }
}

# Ensure rubric function is available
if (!exists("call_from_pid_vec")) {
  call_from_pid_vec <- function(x){
    out <- rep("WEAK/AMBIGUOUS", length(x))
    out[is.na(x)] <- "UNRESOLVED"
    out[!is.na(x) & x >= 80] <- "SUGGESTIVE"
    out[!is.na(x) & x >= 90] <- "STRONG"
    out
  }
}

# Load results (from memory if available, else from disk)
df148 <- get_df("nearest_148", "cytb_unknown_best_matches148.csv")
df300 <- get_df("nearest_300", "cytb_unknown_best_matches300.csv")

# Be robust to presence/absence of pred_species
if (!"pred_species" %in% names(df148)) {
  df148$pred_species <- extract_species(df148$best_known)
}
if (!"pred_species" %in% names(df300)) {
  df300$pred_species <- extract_species(df300$best_known)
}

# Recompute tiers (just to be sure)
df148$tier <- call_from_pid_vec(df148$pct_identity)
df300$tier <- call_from_pid_vec(df300$pct_identity)

# Pair by sample id (the 'unknown' column is the sample ID)
paired <- dplyr::inner_join(
  df148 %>% dplyr::transmute(sample = unknown,
                             species_148 = pred_species,
                             pct_148 = pct_identity,
                             margin_148 = gap_to_second_pct,
                             tier_148 = tier),
  df300 %>% dplyr::transmute(sample = unknown,
                             species_300 = pred_species,
                             pct_300 = pct_identity,
                             margin_300 = gap_to_second_pct,
                             tier_300 = tier),
  by = "sample"
) %>%
  dplyr::mutate(
    delta_pct    = pct_300 - pct_148,
    delta_margin = margin_300 - margin_148,
    tier_rank_148 = match(tier_148, c("WEAK/AMBIGUOUS","SUGGESTIVE","STRONG")),
    tier_rank_300 = match(tier_300, c("WEAK/AMBIGUOUS","SUGGESTIVE","STRONG")),
    delta_tier    = tier_rank_300 - tier_rank_148,
    tier_change   = dplyr::case_when(
      is.na(delta_tier)      ~ "Unknown",
      delta_tier > 0         ~ "Improved",
      delta_tier < 0         ~ "Declined",
      TRUE                   ~ "Unchanged"
    ),
    species_switched = dplyr::case_when(
      is.na(species_148) | is.na(species_300) ~ NA,
      TRUE                                    ~ species_148 != species_300
    ),
    near_tie_148 = ifelse(is.na(margin_148), NA, margin_148 < 1)  # <1 pp margin in 148-bp run
  )

# Write the table
write.csv(paired, "paired_comparison_148_vs_300.csv", row.names = FALSE)

# Brief console summary
n <- nrow(paired)
more_margin <- sum(paired$delta_margin > 0, na.rm = TRUE)
less_margin <- sum(paired$delta_margin < 0, na.rm = TRUE)
same_margin <- sum(paired$delta_margin == 0, na.rm = TRUE)

tab_tier <- table(paired$tier_change, useNA = "no")
improved <- unname(tab_tier["Improved"]); if (is.na(improved)) improved <- 0
declined <- unname(tab_tier["Declined"]); if (is.na(declined)) declined <- 0
unchanged <- unname(tab_tier["Unchanged"]); if (is.na(unchanged)) unchanged <- 0

switched <- sum(paired$species_switched %in% TRUE, na.rm = TRUE)
switched_with_gain <- sum(paired$species_switched %in% TRUE & paired$delta_pct > 0, na.rm = TRUE)
switched_near_tie_148 <- sum(paired$species_switched %in% TRUE & (paired$near_tie_148 %in% TRUE), na.rm = TRUE)

cat(sprintf(
  "\nPAIRED SUMMARY (n = %d):\n", n))
cat(sprintf(" • Margin larger in ~300 for %d/%d samples (%.1f%%); smaller for %d; unchanged for %d.\n",
            more_margin, n, ifelse(n>0, 100*more_margin/n, NA), less_margin, same_margin))
cat(sprintf(" • Confidence tier changes — Improved: %d, Declined: %d, Unchanged: %d.\n",
            improved, declined, unchanged))
cat(sprintf(" • Species switches: %d/%d samples; of those, %d had higher %% identity in ~300; %d were near-ties (<1 pp) in the 148-bp run.\n\n",
            switched, n, switched_with_gain, switched_near_tie_148))

cat(" - paired_comparison_148_vs_300.csv\n")


cat("\nDone. Outputs:\n",
    " - cytb_unknown_best_matches148.csv\n",
    " - cytb_unknown_best_matches300.csv\n",
    " - cytb_tree_full_148.{pdf,png}\n",
    " - cytb_tree_full_300.{pdf,png}\n",
    " - plots_margin_hist.{pdf,png}\n",
    " - plots_identity_vs_margin.{pdf,png}\n", sep="")
