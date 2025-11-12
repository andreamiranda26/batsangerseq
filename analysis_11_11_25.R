# ==================== BAT CYTB — 148-bp & ~300-bp (NJ + bootstrap) ====================
# Clean, journal-ready outputs:
#  - cytb_unknown_best_matches148.csv
#  - cytb_unknown_best_matches300.csv
#  - paired_comparison_148_vs_300.csv
#  - cytb_tree_full_148.{pdf,png}  (unknowns colored by predicted species, refs gray, common names)
#  - cytb_tree_full_300.{pdf,png}  (same styling)
#  - plots_margin_hist.{pdf,png}
#  - plots_identity_vs_margin.{pdf,png}  (only vline at 90%)
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

# Scientific -> Common (includes both Perimyotis/Pipistrellus)
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
  "Perimyotis subflavus"      = "Tri-colored bat",
  "Pipistrellus subflavus"    = "Tri-colored bat",
  "Tadarida brasiliensis"     = "Mexican free-tailed bat",
  "Glossophaga soricina"      = "Pallas’s long-tongued bat"
)
to_common <- function(species_vec) {
  out <- unname(SCI_TO_COMMON[species_vec])
  out[is.na(out)] <- species_vec
  out
}

# Distinct color per *common-name* species for unknowns (12+ distinct hues)
COMMON_PALETTE <- c(
  "Big brown bat"            = "#1b9e77",
  "Silver-haired bat"        = "#d95f02",
  "Eastern red bat"          = "#7570b3",
  "Hoary bat"                = "#e7298a",
  "Northern yellow bat"      = "#66a61e",
  "Seminole bat"             = "#e6ab02",
  "Southeastern myotis"      = "#a6761d",
  "Gray bat"                 = "#666666",
  "Little brown bat"         = "#1f78b4",
  "Northern long-eared bat"  = "#b2df8a",
  "Evening bat"              = "#33a02c",
  "Tri-colored bat"          = "#fb9a99",
  "Mexican free-tailed bat"  = "#cab2d6"
  # (Outgroup not colored here—references stay gray)
)

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

# Plot: unknown tips colored by predicted species (common names); refs gray with common names
plot_tree_colored <- function(tree, unknown_names, ref_common_map, unknown_to_common, out_stub) {
  tr <- ape::ladderize(tree, right = TRUE)

  # Build tip labels: refs -> common names; unknowns -> keep sample IDs
  labs <- tr$tip.label
  is_unknown <- labs %in% unknown_names
  labs_new <- labs
  labs_new[!is_unknown] <- ref_common_map[match(labs[!is_unknown], names(ref_common_map))]
  labs_new[is.na(labs_new)] <- labs[is.na(labs_new)]
  tr$tip.label <- labs_new

  # Colors: unknowns colored by predicted species (common); refs gray
  tip_cols <- rep("grey35", length(labs))
  names(tip_cols) <- tr$tip.label
  # Map unknown sample ID (original label before relabel) to species-common
  # We need a vector aligned to current order: create lookup by original name
  orig_lbl <- labs
  for (i in seq_along(orig_lbl)) {
    if (orig_lbl[i] %in% unknown_names) {
      sp_common <- unknown_to_common[orig_lbl[i]]
      if (!is.na(sp_common) && sp_common %in% names(COMMON_PALETTE)) {
        tip_cols[i] <- COMMON_PALETTE[[sp_common]]
      }
    }
  }

  # Legend: unique species present among unknowns
  unk_sp_present <- unique(unknown_to_common[intersect(names(unknown_to_common), unknown_names)])
  unk_sp_present <- unk_sp_present[!is.na(unk_sp_present)]
  leg_cols <- COMMON_PALETTE[unk_sp_present]

  # Plot with big right margin and legend outside the panel
  H <- max(ape::node.depth.edgelength(tr)); extend <- 1.18
  for (dev in c("pdf","png")) {
    if (dev=="pdf") pdf(paste0(out_stub,".pdf"), width=12, height=16, family="Helvetica")
    else             png(paste0(out_stub,".png"),  width=2200, height=3000, res=300)

    par(mar=c(5,4,2,14), xpd=NA)  # extra right margin for legend
    plot(tr, type="phylogram", direction="rightwards",
         tip.color=tip_cols, edge.color="grey60",
         cex=max(0.45, min(0.9, 95/length(tr$tip.label))),
         label.offset=0.006, no.margin=FALSE,
         x.lim=c(0, H*extend))

    # Draw node bootstraps if present (>=70)
    bs <- attr(tr, "boot")
    if (!is.null(bs)) {
      bs_show <- ifelse(!is.na(bs) & bs >= 70, as.character(round(bs)), "")
      try(nodelabels(text = bs_show, frame = "none", cex = 0.6, col = "grey30"), silent = TRUE)
    }

    # Legends (species and tip types) placed in right outer margin
    usr <- par("usr")
    x_leg <- usr[2] + (usr[2]-usr[1])*0.06
    y_top <- usr[4]
    if (length(unk_sp_present)) {
      legend(x_leg, y_top, title="Unknowns: predicted species",
             legend=unk_sp_present, col=leg_cols, pch=19,
             bty="n", cex=0.9, xjust=0, yjust=1)
    }
    legend(x_leg, y_top - (usr[4]-usr[3])*0.25, title="Tip type",
           legend=c("Unknown sample","Reference"),
           pch=c(19,1), col=c("black","grey35"),
           bty="n", cex=0.9, xjust=0, yjust=1)

    dev.off()
  }
}

# ---- 3) References (Eumops glaucinus EXCLUDED) ----------------------------------
stopifnot("Known FASTA not found" = file.exists(known_seq_file))
knowns_raw <- readDNAStringSet(known_seq_file)
knowns_sel <- remove_gaps(knowns_raw)  # ensure no gaps

# Alabama set (no Eumops glaucinus; also drop any Choeronycteris just in case)
wanted_acc <- c(
  # Eptesicus fuscus
  "AF376835.1","MF038479.1",
  # Lasionycteris noctivagans
  "KC747682.1",
  # Lasiurus borealis
  "KP341708.1","KP341709.1",
  # Lasiurus cinereus
  "KP341731.1","KP341713.1",
  # Lasiurus intermedius
  "KC747687.1","KP341748.1",
  # Lasiurus seminolus
  "KP341753.1","KP341751.1",
  # Myotis austroriparius
  "AM261885.1",
  # Myotis grisescens
  "AM261892.1",
  # Myotis lucifugus
  "OM160895.1","OM160889.1",
  # Myotis septentrionalis
  "DQ503551.1","AM262335.1",
  # Nycticeius humeralis
  "OP157144.1","KC747697.1",
  # Perimyotis/Pipistrellus subflavus
  "AJ504449.1",
  # Tadarida brasiliensis
  "MF135779.1","MF135770.1"
)

all_acc <- sub("\\s.*$", "", names(knowns_sel))
keep_idx <- match(wanted_acc[wanted_acc %in% all_acc], all_acc)
knowns_sel <- knowns_sel[keep_idx]

# Safety drops
drop_ix <- grep("(Choeronycteris|^KC747677\\.1\\b|Eumops\\s+glaucinus|EU350026\\.1|EU350025\\.1)", 
                names(knowns_sel), ignore.case = TRUE)
if (length(drop_ix)) knowns_sel <- knowns_sel[-drop_ix]

# Map reference label -> species (common names)
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
  unk_aln <- DECIPHER::AlignSeqs(unk_raw, verbose = FALSE)  # stabilize indels
  Biostrings::subseq(unk_aln, start = start, width = width)
}

unknowns_148 <- make_window(ab1_files, trim_len_raw,  align_trim_start, len_mini)
unknowns_300 <- make_window(ab1_files, trim_len_long, align_trim_start, len_frame)

# Strip gaps before the joint alignment
unknowns_148_ng <- remove_gaps(unknowns_148)
unknowns_300_ng <- remove_gaps(unknowns_300)
knowns_sel_ng   <- remove_gaps(knowns_sel)

# ---- 5) Joint alignments (unknowns + refs) --------------------------------------
combined_148 <- DECIPHER::AlignSeqs(c(unknowns_148_ng, knowns_sel_ng), verbose = FALSE)
combined_300 <- DECIPHER::AlignSeqs(c(unknowns_300_ng, knowns_sel_ng), verbose = FALSE)

unk_names_148 <- names(unknowns_148_ng)
ref_names_148 <- names(knowns_sel_ng)
unk_names_300 <- names(unknowns_300_ng)
ref_names_300 <- names(knowns_sel_ng)

# ---- 6) Nearest matches (CSV) ---------------------------------------------------
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

# For coloring: map unknown sample -> predicted species *common name*
unknown_to_common_148 <- setNames(to_common(nearest_148$pred_species), nearest_148$unknown)
unknown_to_common_300 <- setNames(to_common(nearest_300$pred_species), nearest_300$unknown)

# ---- 7) NJ trees + bootstrap (optional outgroup) --------------------------------
build_tree_with_optional_outgroup <- function(combined_aln){
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
  tr
}

# Build trees
tree148 <- build_tree_with_optional_outgroup(combined_148)
tree300 <- build_tree_with_optional_outgroup(combined_300)

# Plot trees (unknowns colored by predicted species; refs gray, labeled by common names)
plot_tree_colored(tree148, unk_names_148, refs_common_map, unknown_to_common_148,
                  out_stub = "cytb_tree_full_148")
plot_tree_colored(tree300, unk_names_300, refs_common_map, unknown_to_common_300,
                  out_stub = "cytb_tree_full_300")

# ---- 8) Identity vs margin plots ----------------------------------------------
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
  labs(title = "Margin (best - second-best) distributions",
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
  geom_vline(xintercept = 90, linetype = "dashed") +  # only 90%
  labs(title = "Percent identity vs margin",
       x = "Percent identity (%)",
       y = "Margin to second-best (percentage points)",
       color = "Confidence") +
  theme_classic(base_size = 12)
ggsave("plots_identity_vs_margin.pdf", p_scatter, width = 9.5, height = 4.8)
ggsave("plots_identity_vs_margin.png", p_scatter, width = 9.5, height = 4.8, dpi = 300)

# ---- 9) Paired per-sample comparison: 148 vs ~300 --------------------------------
get_df <- function(obj_name, csv_path) {
  if (exists(obj_name, inherits = TRUE)) get(obj_name, inherits = TRUE)
  else { stopifnot("Missing file: {csv_path}" = file.exists(csv_path))
         read.csv(csv_path, stringsAsFactors = FALSE) }
}

df148 <- get_df("nearest_148", "cytb_unknown_best_matches148.csv")
df300 <- get_df("nearest_300", "cytb_unknown_best_matches300.csv")

if (!"pred_species" %in% names(df148)) df148$pred_species <- extract_species(df148$best_known)
if (!"pred_species" %in% names(df300)) df300$pred_species <- extract_species(df300$best_known)
df148$tier <- call_from_pid_vec(df148$pct_identity)
df300$tier <- call_from_pid_vec(df300$pct_identity)

paired <- dplyr::inner_join(
  df148 %>% transmute(sample = unknown,
                      species_148 = to_common(pred_species),
                      pct_148 = pct_identity,
                      margin_148 = gap_to_second_pct,
                      tier_148 = tier),
  df300 %>% transmute(sample = unknown,
                      species_300 = to_common(pred_species),
                      pct_300 = pct_identity,
                      margin_300 = gap_to_second_pct,
                      tier_300 = tier),
  by = "sample"
) %>%
  mutate(
    delta_pct    = pct_300 - pct_148,
    delta_margin = margin_300 - margin_148,
    tier_rank_148 = match(tier_148, c("WEAK/AMBIGUOUS","SUGGESTIVE","STRONG")),
    tier_rank_300 = match(tier_300, c("WEAK/AMBIGUOUS","SUGGESTIVE","STRONG")),
    delta_tier    = tier_rank_300 - tier_rank_148,
    tier_change   = case_when(
      is.na(delta_tier) ~ "Unknown",
      delta_tier > 0    ~ "Improved",
      delta_tier < 0    ~ "Declined",
      TRUE              ~ "Unchanged"
    ),
    species_switched = ifelse(is.na(species_148) | is.na(species_300), NA, species_148 != species_300),
    near_tie_148 = ifelse(is.na(margin_148), NA, margin_148 < 1)
  )

write.csv(paired, "paired_comparison_148_vs_300.csv", row.names = FALSE)


# ---- 10) Text-aligned summaries & reproduction checks -------------------------
num_fmt <- function(x, d=2) ifelse(is.na(x), NA, format(round(x, d), nsmall=d))
range_str <- function(v) paste0(num_fmt(min(v, na.rm=TRUE)), "–", num_fmt(max(v, na.rm=TRUE)))
med_str   <- function(v) num_fmt(stats::median(v, na.rm=TRUE))

# Ensure tiers exist
if (!"tier" %in% names(df148)) df148$tier <- call_from_pid_vec(df148$pct_identity)
if (!"tier" %in% names(df300)) df300$tier <- call_from_pid_vec(df300$pct_identity)

# Helper: species composition (common names)
species_comp <- function(df) {
  sp <- to_common(df$pred_species)
  sort(table(sp), decreasing = TRUE)
}

sink("results_reproduction_checks.txt")

cat("\n=== 148-bp mini-barcode ===\n")
cat("N unknowns with best match:", nrow(df148), "\n")
cat("Percent identity range:", range_str(df148$pct_identity), "\n")
cat("Percent identity median:", med_str(df148$pct_identity), "\n")
cat("Margin (best-second) range:", range_str(df148$gap_to_second_pct), "\n")
cat("Margin median:", med_str(df148$gap_to_second_pct), "\n")
tab148 <- table(df148$tier, useNA = "no")
cat("Confidence counts:", paste(paste(names(tab148), as.integer(tab148), sep="="), collapse=", "), "\n")
cat("Species composition (common names):\n"); print(species_comp(df148))

cat("\n=== ~300-bp reference-anchored frame ===\n")
cat("N unknowns with best match:", nrow(df300), "\n")
cat("Percent identity range:", range_str(df300$pct_identity), "\n")
cat("Percent identity median:", med_str(df300$pct_identity), "\n")
cat("Margin (best-second) range:", range_str(df300$gap_to_second_pct), "\n")
cat("Margin median:", med_str(df300$gap_to_second_pct), "\n")
tab300 <- table(df300$tier, useNA = "no")
cat("Confidence counts:", paste(paste(names(tab300), as.integer(tab300), sep="="), collapse=", "), "\n")
cat("Species composition (common names):\n"); print(species_comp(df300))

cat("\n=== Paired per-sample comparisons (148 vs ~300) ===\n")
n <- nrow(paired)
cat("Total paired samples:", n, "\n")
cat("Larger margin in ~300:", sum(paired$delta_margin > 0, na.rm=TRUE), "of", n, "\n")
cat("Median per-sample Δmargin (pp):", med_str(paired$delta_margin), "\n")

tab_change <- table(paired$tier_change, useNA="no")
improved <- unname(tab_change["Improved"]); if (is.na(improved)) improved <- 0
declined <- unname(tab_change["Declined"]); if (is.na(declined)) declined <- 0
unchanged <- unname(tab_change["Unchanged"]); if (is.na(unchanged)) unchanged <- 0
cat("Confidence tier changes — Improved:", improved, 
    "Declined:", declined, "Unchanged:", unchanged, "\n")

switched <- paired %>% filter(species_switched %in% TRUE)
cat("Samples with top-species switch:", nrow(switched), "\n")
if (nrow(switched)) {
  cat("  with higher % identity in ~300:", sum(switched$delta_pct > 0, na.rm=TRUE), "\n")
  cat("  that were near-ties (<1 pp) in 148:", sum(switched$near_tie_148 %in% TRUE, na.rm=TRUE), "\n")
}

cat("\n=== Near-ties in 148 (<1 pp) that gained separation in ~300 ===\n")
near_ties <- paired %>% 
  filter(!is.na(margin_148), margin_148 < 1, !is.na(margin_300), margin_300 > margin_148) %>%
  arrange(sample)

cat("Count:", nrow(near_ties), "\n")
if (nrow(near_ties)) {
  pretty_rows <- near_ties %>%
    transmute(
      sample,
      species_148,
      species_300,
      margin_148 = num_fmt(margin_148, 2),
      margin_300 = num_fmt(margin_300, 2),
      pct_148 = num_fmt(pct_148, 2),
      pct_300 = num_fmt(pct_300, 2)
    )
  print(pretty_rows, row.names = FALSE)
}

sink()

cat("\nWrote a human-readable check file: results_reproduction_checks.txt\n")


# Console summary
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

cat("\nDone. Outputs:\n",
    " - cytb_unknown_best_matches148.csv\n",
    " - cytb_unknown_best_matches300.csv\n",
    " - paired_comparison_148_vs_300.csv\n",
    " - cytb_tree_full_148.{pdf,png}\n",
    " - cytb_tree_full_300.{pdf,png}\n",
    " - plots_margin_hist.{pdf,png}\n",
    " - plots_identity_vs_margin.{pdf,png}\n", sep="")

