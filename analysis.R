# ==================== BAT CYTB (148-bp) — AL accessions + external outgroup ====================

# ---- 0) Setup -------------------------------------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Bats")  # <- adjust if needed

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("sangerseqR","DECIPHER","Biostrings","ape")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("sangerseqR","DECIPHER","Biostrings")) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}
library(sangerseqR); library(DECIPHER); library(Biostrings); library(ape)

# ---- 1) Parameters --------------------------------------------------------------------------
folder           <- "data/cytb"                 # folder with .ab1 chromatograms
trim_len_raw     <- 170                         # hard trim of raw reads
align_trim_start <- 23                          # drop first 22 bases (primer/poly-A)
len_mini         <- 148                         # 148-bp mini-barcode
known_seq_file   <- "cytb_known_batseqs.fasta"  # your reference panel (multi-FASTA)
outgroup_fa  <- "cytb_outgroup.fasta"  # your A. jamaicensis cytb FASTA (single sequence)
# Force rooting by name pattern:
outgroup_hint <- "^Glossophaga soricina$" 


# ---- 2) Utilities ----------------------------------------------------------------------------
remove_gaps <- function(dna) DNAStringSet(gsub("-", "", as.character(dna), fixed = TRUE))

# robust species extractor: handles headers that start with accession or directly with "Genus species"
extract_species <- function(x) {
  sapply(x, function(xx) {
    m1 <- regexec("^[^ ]+\\s+([A-Z][a-z]+\\s+[a-z]+)", xx); g1 <- regmatches(xx, m1)[[1]]
    if (length(g1) >= 2) return(g1[2])
    m2 <- regexec("^([A-Z][a-z]+\\s+[a-z]+)", xx); g2 <- regmatches(xx, m2)[[1]]
    if (length(g2) >= 2) return(g2[2])
    NA_character_
  })
}
abbr_species <- function(s) {
  ifelse(grepl("^[A-Z][a-z]+\\s+[a-z]+$", s),
         sub("^([A-Z])[a-z]+\\s+([a-z]+)$", "\\1. \\2", s), s)
}

# Pairwise-deletion % identity from an alignment
pid_matrix_from_alignment <- function(aln) {
  bases <- c("A","C","G","T")
  M <- do.call(rbind, strsplit(as.character(aln), ""))
  M <- toupper(M)
  n <- nrow(M)
  pid <- matrix(NA_real_, n, n, dimnames = list(rownames(M), rownames(M)))
  valid_counts <- matrix(0L, n, n, dimnames = dimnames(pid))
  for (i in 1:(n-1)) for (j in (i+1):n) {
    mask <- (M[i,] %in% bases) & (M[j,] %in% bases)
    comp <- sum(mask)
    valid_counts[i,j] <- valid_counts[j,i] <- comp
    if (comp > 0) {
      matches <- sum(M[i, mask] == M[j, mask])
      pid[i,j] <- pid[j,i] <- (matches/comp)*100
    }
  }
  diag(pid) <- 100
  diag(valid_counts) <- ncol(M)
  list(pid = pid, valid = valid_counts)
}

# Confidence rubric (advisor’s request)
call_from_pid <- function(pid) {
  if (is.na(pid)) return("UNRESOLVED")
  if (pid >= 90) "STRONG" else if (pid >= 80) "SUGGESTIVE" else "WEAK/AMBIGUOUS"
}

nearest_table <- function(pid, valid, is_unknown) {
  dist_mat <- 1 - pid/100
  dist_mat[is.na(dist_mat)] <- 1
  diag(dist_mat) <- 0
  seq_names <- rownames(dist_mat)
  unknown_ids <- if (is.logical(is_unknown)) seq_names[is_unknown] else seq_names[seq_names %in% is_unknown]
  known_ids   <- setdiff(seq_names, unknown_ids)
  res <- lapply(unknown_ids, function(u) {
    drow <- dist_mat[u, known_ids]; vrow <- valid[u, known_ids]
    ok   <- !is.na(drow) & vrow > 0
    if (!any(ok)) return(data.frame(unknown=u, best_known=NA, dist=NA, pct_identity=NA,
                                    second_best_dist=NA, second_best_pct_identity=NA,
                                    gap_to_second_pct=NA, stringsAsFactors=FALSE))
    ord <- order(drow[ok], decreasing = FALSE); k_ok <- names(drow[ok])
    best   <- k_ok[ord[1]]; best_d <- drow[ok][ord[1]]; best_pid <- (1 - best_d)*100
    if (length(ord) >= 2) { second_d <- drow[ok][ord[2]]; second_pid <- (1 - second_d)*100; gap <- best_pid - second_pid
    } else { second_d <- NA; second_pid <- NA; gap <- NA }
    data.frame(
      unknown = u, best_known = best,
      dist = round(best_d, 5), pct_identity = round(best_pid, 2),
      second_best_dist = ifelse(is.na(second_d), NA, round(second_d,5)),
      second_best_pct_identity = ifelse(is.na(second_pid), NA, round(second_pid,2)),
      gap_to_second_pct = ifelse(is.na(gap), NA, round(gap,2)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

# ---- 3) Read unknown chromatograms; align & window to 148 bp --------------------------------
ab1_files <- list.files(folder, pattern = "\\.ab1$", full.names = TRUE)
stopifnot("No .ab1 files found" = length(ab1_files) > 0)

seqs_out <- setNames(vector("list", length(ab1_files)),
                     tools::file_path_sans_ext(basename(ab1_files)))
for (i in seq_along(ab1_files)) {
  abif <- read.abif(ab1_files[i])
  s <- sangerseq(abif)
  raw_seq <- toupper(as.character(primarySeq(s)))
  seqs_out[[i]] <- DNAString(substr(raw_seq, 1, trim_len_raw))
}
unknowns_raw <- DNAStringSet(seqs_out)
writeXStringSet(unknowns_raw, "trimmed_170.fasta")

unknowns_aln <- DECIPHER::AlignSeqs(unknowns_raw, verbose = FALSE)
unknowns_148 <- Biostrings::subseq(unknowns_aln, start = align_trim_start, width = len_mini)

# ---- 4) Load references; keep EXACT Alabama accession list (C. mexicana REMOVED) ----
stopifnot("Known FASTA not found" = file.exists(known_seq_file))
knowns_raw <- readDNAStringSet(known_seq_file)
knowns_raw <- DNAStringSet(gsub("-", "", as.character(knowns_raw), fixed = TRUE))  # remove gaps

# Accessions to include (EXCLUDES KC747677.1 Choeronycteris mexicana)
wanted_acc <- c(
  # Eptesicus fuscus
  "AF376835.1","MF038479.1",
  # Eumops glaucinus
  "EU350026.1","EU350025.1",
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
  # Pipistrellus (Perimyotis) subflavus
  "AJ504449.1",
  # Tadarida brasiliensis
  "MF135779.1","MF135770.1"
)

# Subset in your specified order
all_acc <- sub("\\s.*$", "", names(knowns_raw))
present <- wanted_acc %in% all_acc
if (!all(present)) {
  warning("These requested accessions were NOT found in ", known_seq_file, ":\n  ",
          paste(wanted_acc[!present], collapse = ", "))
}
keep_idx   <- match(wanted_acc[present], all_acc)
knowns_sel <- knowns_raw[keep_idx]

# Safety: if any Choeronycteris sneaks in by mistake, drop it
drop_ix <- grep("(^KC747677\\.1\\b)|Choeronycteris", names(knowns_sel), ignore.case = TRUE)
if (length(drop_ix)) knowns_sel <- knowns_sel[-drop_ix]

# Align the selected references
knowns_AL_aln <- DECIPHER::AlignSeqs(knowns_sel, verbose = FALSE)


# ---- 5) Combine unknowns (148) + references --------------------------------------------------
if ("AlignProfiles" %in% getNamespaceExports("DECIPHER")) {
  combined_148_AL <- tryCatch(
    DECIPHER::AlignProfiles(unknowns_148, knowns_AL_aln),
    error = function(e) DECIPHER::AlignSeqs(c(unknowns_148, knowns_sel), verbose = FALSE)
  )
} else {
  combined_148_AL <- DECIPHER::AlignSeqs(c(unknowns_148, knowns_sel), verbose = FALSE)
}
is_unknown_AL <- names(combined_148_AL) %in% names(unknowns_148)

# ---- 6) Nearest matches (148 bp) -------------------------------------------------------------
pm_148_AL <- pid_matrix_from_alignment(combined_148_AL)
nearest_148_AL <- nearest_table(pm_148_AL$pid, pm_148_AL$valid, is_unknown_AL)
nearest_148_AL$pred_species  <- extract_species(nearest_148_AL$best_known)
nearest_148_AL$call_strength <- vapply(nearest_148_AL$pct_identity, call_from_pid, character(1))
write.csv(nearest_148_AL, "cytb_unknown_best_matches148.csv", row.names = FALSE)
cat("Wrote: cytb_unknown_best_matches148.csv\n")

# ==============================================================================================
# Rectangular plotting helpers (FULL + REDUCED) with legends bottom-right & outgroup
# ==============================================================================================
plot_rect_full_148 <- function(tree, D, nearest_df, is_unknown_vec,
                               tip_labels_original, out_stub,
                               outgroup_hint = NULL,
                               direction = "rightwards",
                               label_refs = TRUE) {
  extract_species <- function(x){
    sapply(x, function(xx) {
      m1 <- regexec("^[^ ]+\\s+([A-Z][a-z]+\\s+[a-z]+)", xx); g1 <- regmatches(xx, m1)[[1]]
      if (length(g1) >= 2) return(g1[2])
      m2 <- regexec("^([A-Z][a-z]+\\s+[a-z]+)", xx); g2 <- regmatches(xx, m2)[[1]]
      if (length(g2) >= 2) return(g2[2])
      NA_character_
    })
  }
  abbr_species <- function(s){
    ifelse(grepl("^[A-Z][a-z]+\\s+[a-z]+$", s),
           sub("^([A-Z])[a-z]+\\s+([a-z]+)$","\\1. \\2", s), s)
  }
  root_with_outgroup <- function(tree, D, tip_labels_original, is_unknown_vec, outgroup_hint=NULL){
    known_labels <- tip_labels_original[!is_unknown_vec]
    chosen_label <- NULL
    if (!is.null(outgroup_hint)) {
      hits <- known_labels[grepl(outgroup_hint, known_labels, ignore.case=TRUE)]
      if (length(hits) >= 1) chosen_label <- hits[1]
    }
    if (is.null(chosen_label)) {
      Dm <- as.matrix(D)
      mu <- rowMeans(Dm[known_labels, known_labels, drop=FALSE], na.rm=TRUE)
      chosen_label <- names(which.max(mu))
    }
    rooted <- tryCatch(ape::root(tree, outgroup=chosen_label, resolve.root=TRUE),
                       error=function(e) tree)
    list(tree=rooted, outgroup_label=chosen_label)
  }
  
  pred_map <- setNames(extract_species(nearest_df$best_known), nearest_df$unknown)
  ref_species_map <- setNames(extract_species(tip_labels_original), tip_labels_original)
  species_by_name <- ref_species_map; species_by_name[names(pred_map)] <- pred_map[names(pred_map)]
  species_by_name[is.na(species_by_name)] <- "Unassigned"
  
  rooted <- root_with_outgroup(tree, D, tip_labels_original, is_unknown_vec, outgroup_hint)
  tr <- rooted$tree; out_label <- rooted$outgroup_label
  
  is_unknown_by_name <- setNames(is_unknown_vec, tip_labels_original)
  unk_flags <- unname(ifelse(is.na(is_unknown_by_name[tr$tip.label]),
                             FALSE, is_unknown_by_name[tr$tip.label]))
  species_tips <- unname(ifelse(is.na(species_by_name[tr$tip.label]),
                                "Unassigned", species_by_name[tr$tip.label]))
  
  sp_levels <- sort(unique(na.omit(species_tips[unk_flags])))
  okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73",
                 "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
  pal <- rep(okabe_ito, length.out=length(sp_levels))
  cols_unknown <- setNames(pal, sp_levels)
  tip_cols <- ifelse(unk_flags, cols_unknown[species_tips], "grey55")
  
  lab_refs <- if (label_refs) abbr_species(species_tips) else rep("", length(species_tips))
  lab_new  <- ifelse(unk_flags, tr$tip.label, lab_refs)
  
  # outgroup index BEFORE relabeling
  out_tip_idx <- match(out_label, tr$tip.label)
  
  # spacing & limits
  H <- max(ape::node.depth.edgelength(tr))
  extend <- 1.18
  xlim <- if (direction %in% c("rightwards","leftwards")) c(0, H*extend) else NULL
  ylim <- if (direction %in% c("upwards","downwards")) c(0, H*extend) else NULL
  cex_auto    <- max(0.38, min(0.85, 90/length(tr$tip.label)))
  offset_auto <- 0.006 + 0.0012*log1p(length(tr$tip.label))
  
  for (dev in c("pdf","png")) {
    if (dev=="pdf") pdf(paste0(out_stub,".pdf"), width=18, height=26, family="Helvetica")
    else             png(paste0(out_stub,".png"),  width=3200, height=4600, res=300)
    
    par(mar=c(8,8,6,26), xpd=NA)  # big right margin
    
    tr_plot <- tr; tr_plot$tip.label <- lab_new
    plot(tr_plot, type="phylogram", direction=direction,
         tip.color=tip_cols, edge.color="grey80",
         cex=cex_auto, label.offset=offset_auto,
         no.margin=FALSE, use.edge.length=TRUE, adj=0,
         x.lim=xlim, y.lim=ylim)
    
    unk_tips <- which(unk_flags); ref_tips <- which(!unk_flags)
    if (length(unk_tips)) tiplabels(pch=19, tip=unk_tips, col=tip_cols[unk_tips], cex=0.7)
    if (length(ref_tips)) tiplabels(pch=1,  tip=ref_tips, col="grey40", cex=0.6)
    if (!is.na(out_tip_idx)) tiplabels(pch=17, tip=out_tip_idx, col="black", cex=1.0)
    
    # legends in lower-right margin (explicit coords)
    usr <- par("usr")
    x_leg <- usr[2] + (usr[2]-usr[1])*0.12
    y_leg_species <- usr[3] + (usr[4]-usr[3])*0.26
    y_leg_tip     <- usr[3] + (usr[4]-usr[3])*0.08
    
    if (length(sp_levels)) {
      legend(x_leg, y_leg_species, title="Unknowns: predicted species",
             legend=names(cols_unknown), col=cols_unknown, pch=19,
             bty="n", cex=0.95, xjust=0, yjust=0)
    }
    legend(x_leg, y_leg_tip, title="Tip type",
           legend=c("Unknown sample","Reference","Outgroup"),
           pch=c(19,1,17), col=c("black","grey40","black"),
           bty="n", cex=1.0, xjust=0, yjust=0)
    
    dev.off()
  }
}

plot_rect_reduced_148 <- function(tree, D, nearest_df, is_unknown_vec,
                                  tip_labels_original, out_stub,
                                  outgroup_hint = NULL,
                                  direction = "rightwards") {
  extract_species <- function(x){
    sapply(x, function(xx) {
      m1 <- regexec("^[^ ]+\\s+([A-Z][a-z]+\\s+[a-z]+)", xx); g1 <- regmatches(xx, m1)[[1]]
      if (length(g1) >= 2) return(g1[2])
      m2 <- regexec("^([A-Z][a-z]+\\s+[a-z]+)", xx); g2 <- regmatches(xx, m2)[[1]]
      if (length(g2) >= 2) return(g2[2])
      NA_character_
    })
  }
  abbr_species <- function(s){
    ifelse(grepl("^[A-Z][a-z]+\\s+[a-z]+$", s),
           sub("^([A-Z])[a-z]+\\s+([a-z]+)$","\\1. \\2", s), s)
  }
  root_with_outgroup <- function(tree, D, tip_labels_original, is_unknown_vec, outgroup_hint=NULL){
    known_labels <- tip_labels_original[!is_unknown_vec]
    chosen_label <- NULL
    if (!is.null(outgroup_hint)) {
      hits <- known_labels[grepl(outgroup_hint, known_labels, ignore.case=TRUE)]
      if (length(hits) >= 1) chosen_label <- hits[1]
    }
    if (is.null(chosen_label)) {
      Dm <- as.matrix(D)
      mu <- rowMeans(Dm[known_labels, known_labels, drop=FALSE], na.rm=TRUE)
      chosen_label <- names(which.max(mu))
    }
    rooted <- tryCatch(ape::root(tree, outgroup=chosen_label, resolve.root=TRUE),
                       error=function(e) tree)
    list(tree=rooted, outgroup_label=chosen_label)
  }
  
  # one reference per species (pick first by name)
  ref_species_map <- setNames(extract_species(tip_labels_original), tip_labels_original)
  ref_names <- tip_labels_original[!is_unknown_vec]
  by_sp <- split(ref_names, ref_species_map[ref_names])
  rep_refs <- unique(na.omit(vapply(by_sp, function(v) v[1], character(1))))
  
  keep_names <- c(tip_labels_original[is_unknown_vec], rep_refs)
  # ensure outgroup kept if present in labels
  rooted_tmp <- root_with_outgroup(tree, D, tip_labels_original, is_unknown_vec, outgroup_hint)
  out_label  <- rooted_tmp$outgroup_label
  if (!is.null(out_label) && length(out_label) && !(out_label %in% keep_names)) {
    keep_names <- c(keep_names, out_label)
  }
  
  small <- keep.tip(tree, keep_names)
  Dm <- as.matrix(D); D_small <- as.dist(Dm[keep_names, keep_names])
  
  rooted <- root_with_outgroup(njs(D_small), Dm[keep_names, keep_names], keep_names,
                               is_unknown_vec = small$tip.label %in% tip_labels_original[is_unknown_vec],
                               outgroup_hint = outgroup_hint)
  tr <- rooted$tree; out_label <- rooted$outgroup_label
  
  # flags/colors/labels
  is_unknown_small <- tr$tip.label %in% tip_labels_original[is_unknown_vec]
  pred_map <- setNames(extract_species(nearest_df$best_known), nearest_df$unknown)
  species_all <- setNames(extract_species(keep_names), keep_names)
  species_all[names(pred_map)] <- pred_map[names(pred_map)]
  species_tips <- species_all[tr$tip.label]; species_tips[is.na(species_tips)] <- "Unassigned"
  
  sp_levels <- sort(unique(na.omit(species_tips[is_unknown_small])))
  okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73",
                 "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
  pal <- rep(okabe_ito, length.out=length(sp_levels))
  cols_unknown <- setNames(pal, sp_levels)
  tip_cols <- ifelse(is_unknown_small, cols_unknown[ species_tips ], "grey55")
  
  lab_new <- ifelse(is_unknown_small, tr$tip.label, abbr_species(species_tips))
  out_tip_idx <- match(out_label, tr$tip.label)
  
  H <- max(ape::node.depth.edgelength(tr))
  extend <- 1.18
  xlim <- if (direction %in% c("rightwards","leftwards")) c(0, H*extend) else NULL
  ylim <- if (direction %in% c("upwards","downwards")) c(0, H*extend) else NULL
  cex_auto    <- max(0.52, min(1.0, 95/length(tr$tip.label)))
  offset_auto <- 0.007 + 0.0012*log1p(length(tr$tip.label))
  
  for (dev in c("pdf","png")) {
    if (dev=="pdf") pdf(paste0(out_stub,".pdf"), width=16, height=22, family="Helvetica")
    else             png(paste0(out_stub,".png"),  width=2800, height=4000, res=300)
    par(mar=c(8,8,6,22), xpd=NA)
    
    tr_plot <- tr; tr_plot$tip.label <- lab_new
    plot(tr_plot, type="phylogram", direction=direction,
         tip.color=tip_cols, edge.color="grey80",
         cex=cex_auto, label.offset=offset_auto,
         no.margin=FALSE, use.edge.length=TRUE, adj=0,
         x.lim=xlim, y.lim=ylim)
    
    unk_tips <- which(is_unknown_small); ref_tips <- which(!is_unknown_small)
    if (length(unk_tips)) tiplabels(pch=19, tip=unk_tips, col=tip_cols[unk_tips], cex=0.8)
    if (length(ref_tips)) tiplabels(pch=1,  tip=ref_tips, col="grey40", cex=0.75)
    if (!is.na(out_tip_idx)) tiplabels(pch=17, tip=out_tip_idx, col="black", cex=1.0)
    
    usr <- par("usr")
    x_leg <- usr[2] + (usr[2]-usr[1])*0.10
    y_leg_species <- usr[3] + (usr[4]-usr[3])*0.26
    y_leg_tip     <- usr[3] + (usr[4]-usr[3])*0.08
    if (length(sp_levels)) {
      legend(x_leg, y_leg_species, title="Unknowns: predicted species",
             legend=names(cols_unknown), col=cols_unknown, pch=19,
             bty="n", cex=0.95, xjust=0, yjust=0)
    }
    legend(x_leg, y_leg_tip, title="Tip type",
           legend=c("Unknown sample","Reference","Outgroup"),
           pch=c(19,1,17), col=c("black","grey40","black"),
           bty="n", cex=1.0, xjust=0, yjust=0)
    
    dev.off()
  }
}

# ---- 7) Build trees (add external outgroup for plotting, if present) -------------------------
D148_AL <- 1 - pm_148_AL$pid/100; D148_AL[is.na(D148_AL)] <- 1; diag(D148_AL) <- 0
tree_148_AL <- ape::njs(as.dist(D148_AL))

if (file.exists(outgroup_fa)) {
  og_raw <- readDNAStringSet(outgroup_fa) |> remove_gaps()
  
  if (length(og_raw) == 0) {
    warning("Outgroup FASTA exists but contains 0 sequences; proceeding without outgroup.")
    # Fall back to trees rooted within the Alabama references
    plot_rect_full_148(
      tree_148_AL, D148_AL, nearest_148_AL, is_unknown_AL,
      tip_labels_original = names(combined_148_AL),
      out_stub = "cytb_rect_full_148_AL",
      outgroup_hint = NULL, direction = "rightwards", label_refs = TRUE
    )
    plot_rect_reduced_148(
      tree_148_AL, D148_AL, nearest_148_AL, is_unknown_AL,
      tip_labels_original = names(combined_148_AL),
      out_stub = "cytb_rect_reduced_148_AL",
      outgroup_hint = NULL, direction = "rightwards"
    )
    
  } else {
    # DO NOT call AlignSeqs() on the outgroup alone (needs >=2 seqs).
    # Instead, align the outgroup together with your existing alignment.
    if ("AlignProfiles" %in% getNamespaceExports("DECIPHER")) {
      combined_for_tree <- tryCatch(
        DECIPHER::AlignProfiles(combined_148_AL, og_raw),
        error = function(e) DECIPHER::AlignSeqs(c(combined_148_AL, og_raw), verbose = FALSE)
      )
    } else {
      combined_for_tree <- DECIPHER::AlignSeqs(c(combined_148_AL, og_raw), verbose = FALSE)
    }
    
    # Flags and names for plotting
    is_unknown_tree <- c(is_unknown_AL, rep(FALSE, length(og_raw)))
    tip_names_tree  <- c(names(combined_148_AL), names(og_raw))
    # ensure names are carried through (Align* usually preserves them)
    if (!identical(names(combined_for_tree), tip_names_tree)) {
      names(combined_for_tree) <- tip_names_tree
    }
    
    # Distances & NJ
    pm_tree <- pid_matrix_from_alignment(combined_for_tree)
    D_tree  <- 1 - pm_tree$pid/100; D_tree[is.na(D_tree)] <- 1; diag(D_tree) <- 0
    tree_for_plot <- ape::njs(as.dist(D_tree))
    
    # FULL tree (all refs)
    plot_rect_full_148(
      tree_for_plot, D_tree, nearest_148_AL, is_unknown_tree,
      tip_labels_original = tip_names_tree,
      out_stub = "cytb_rect_full_148_AL_with_external_outgroup",
      outgroup_hint = NULL, direction = "rightwards", label_refs = TRUE
    )
    
    # REDUCED tree (all unknowns + one ref per species)
    plot_rect_reduced_148(
      tree_for_plot, D_tree, nearest_148_AL, is_unknown_tree,
      tip_labels_original = tip_names_tree,
      out_stub = "cytb_rect_reduced_148_AL_with_external_outgroup",
      outgroup_hint = NULL, direction = "rightwards"
    )
  }
  
} else {
  message("NOTE: '", outgroup_fa, "' not found. Trees will be rooted within the AL references.")
  plot_rect_full_148(
    tree_148_AL, D148_AL, nearest_148_AL, is_unknown_AL,
    tip_labels_original = names(combined_148_AL),
    out_stub = "cytb_rect_full_148_AL",
    outgroup_hint = NULL, direction = "rightwards", label_refs = TRUE
  )
  plot_rect_reduced_148(
    tree_148_AL, D148_AL, nearest_148_AL, is_unknown_AL,
    tip_labels_original = names(combined_148_AL),
    out_stub = "cytb_rect_reduced_148_AL",
    outgroup_hint = NULL, direction = "rightwards"
  )
}


cat("\nDone. Outputs:\n",
    " - trimmed_170.fasta\n",
    " - cytb_unknown_best_matches148.csv\n",
    " - cytb_rect_full_148_AL_with_external_outgroup.{pdf,png} (if outgroup present)\n",
    " - cytb_rect_reduced_148_AL_with_external_outgroup.{pdf,png} (if outgroup present)\n", sep = "")
# =================================================================================================

# ---- Plots: margins histogram + identity vs margin scatter ----------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(ggplot2); library(dplyr)

# Build plotting data (148 always present)
d148 <- nearest_148_AL %>%
  transmute(sample = unknown,
            frame  = "148-bp mini-barcode",
            pct_identity,
            margin = gap_to_second_pct,
            tier   = call_strength)

df_list <- list(d148)

# If you also generated a 300-bp CSV earlier, we’ll include it for side-by-side facets
if (file.exists("cytb_unknown_best_matches300.csv")) {
  n300 <- read.csv("cytb_unknown_best_matches300.csv", stringsAsFactors = FALSE)
  if (all(c("unknown","pct_identity","gap_to_second_pct","call_strength") %in% names(n300))) {
    d300 <- n300 %>%
      transmute(sample = unknown,
                frame  = "~300-bp reference-anchored frame",
                pct_identity,
                margin = gap_to_second_pct,
                tier   = call_strength)
    df_list[[length(df_list) + 1]] <- d300
  }
}

dd <- bind_rows(df_list)

# (i) Histogram of margins
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

# (ii) Percent identity vs margin scatter
p_scatter <- ggplot(dd %>% filter(!is.na(margin)),
                    aes(x = pct_identity, y = margin, color = tier)) +
  geom_point(size = 2, alpha = 0.85) +
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
       color = "Confidence (≥90% = STRONG)") +
  theme_classic(base_size = 12)

ggsave("plots_identity_vs_margin.pdf", p_scatter, width = 9.5, height = 4.8)
ggsave("plots_identity_vs_margin.png", p_scatter, width = 9.5, height = 4.8, dpi = 300)

# Quick console summary (optional)
print(dd %>% group_by(frame, tier) %>% summarise(n = n(), .groups = "drop"))

# A tibble: 3 × 3
# frame               tier               n
# <chr>               <chr>          <int>
#   1 148-bp mini-barcode STRONG          39
# 2 148-bp mini-barcode SUGGESTIVE         1
# 3 148-bp mini-barcode WEAK/AMBIGUOUS    11
