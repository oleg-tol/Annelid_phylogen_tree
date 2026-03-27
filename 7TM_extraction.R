#!/usr/bin/env Rscript
library(Biostrings)
library(seqinr)

work_dir   <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/"
input_file <- file.path(work_dir, "opsin_newA3.fa")
output_file <- file.path(work_dir, "listA_filtered_add_7tm_extracted4.fa")

aln     <- readAAStringSet(input_file)
aln_mat <- as.matrix(aln)
cat("Sequences:", length(aln), "| Alignment length:", width(aln)[1], "\n\n")

# ---- CONSERVED OPSIN MOTIFS ----
opsin_motifs <- list(
  N_term_conserved = list(pattern = "G[GAVLI].{0,3}[FYWH]", flexibility = "medium"),
  TM1_Asn          = list(pattern = "N[LIVMF].{2,4}[LIVMF]", flexibility = "medium"),
  TM2_Asp          = list(pattern = "[DE][LIVMF].{2,4}[LIVMF]", flexibility = "high"),
  DRY_motif        = list(pattern = "[DE].{0,2}[RK].{0,2}[YFW]", flexibility = "low"),
  TM4_Trp          = list(pattern = "W.{2,4}[GAVLI]", flexibility = "medium"),
  CWxP_motif       = list(pattern = "C.{0,1}W.{1,2}P", flexibility = "low"),
  FxxCWxP_motif    = list(pattern = "F.{2}C.W.{1,2}P", flexibility = "medium"),
  NPxxY_motif      = list(pattern = "NP.{2}[YF]", flexibility = "low"),
  NPxxY_variant    = list(pattern = "[ND]P.{2}[YF]", flexibility = "medium")
)

# ---- MOTIF SEARCH ----
find_motif_positions <- function(aln_mat, pattern, min_conservation = 0.3) {
  consensus <- apply(aln_mat, 2, function(col) {
    col <- col[col != "-"]
    if (!length(col)) return("-")
    names(which.max(table(col)))
  })
  conservation <- apply(aln_mat, 2, function(col) {
    col <- col[col != "-"]
    if (!length(col)) return(0)
    max(table(col)) / length(col)
  })
  consensus_seq <- paste(consensus, collapse = "")
  m <- gregexpr(pattern, consensus_seq, perl = TRUE)[[1]]
  if (m[1] == -1) return(NULL)
  
  results <- Filter(Negate(is.null), lapply(seq_along(m), function(i) {
    pos <- m[i]; end <- pos + attr(m, "match.length")[i] - 1
    cons <- mean(conservation[pos:end])
    if (cons < min_conservation) return(NULL)
    list(start = pos, end = end, conservation = cons,
         sequence = substr(consensus_seq, pos, end))
  }))
  if (!length(results)) NULL else results
}

min_cons_map <- c(low = 0.4, medium = 0.3, high = 0.2)

motif_results <- lapply(names(opsin_motifs), function(nm) {
  m   <- opsin_motifs[[nm]]
  res <- find_motif_positions(aln_mat, m$pattern, min_cons_map[m$flexibility])
  cat(nm, "->", if (is.null(res)) "NOT FOUND" else
    paste(sapply(res, function(x) sprintf("%d-%d (%.2f)", x$start, x$end, x$conservation)),
          collapse = ", "), "\n")
  res
})
names(motif_results) <- names(opsin_motifs)

# ---- DEFINE 7TM BOUNDARIES ----
dry_matches <- motif_results[["DRY_motif"]]
npy_matches <- motif_results[["NPxxY_motif"]] %||% motif_results[["NPxxY_variant"]]

`%||%` <- function(a, b) if (!is.null(a) && length(a)) a else b

if (!is.null(dry_matches) && !is.null(npy_matches)) {
  pairs <- do.call(rbind, lapply(dry_matches, function(d)
    lapply(npy_matches, function(n) {
      sp <- n$start - d$start
      if (sp < 150 || sp > 300) return(NULL)
      list(dry_pos = d$start, npy_pos = n$start, spacing = sp,
           score = (d$conservation + n$conservation) / 2 - abs(sp - 200) / 100)
    })))
  pairs <- Filter(Negate(is.null), pairs)
  
  if (length(pairs)) {
    best     <- pairs[[which.max(sapply(pairs, `[[`, "score"))]]
    dry_pos  <- best$dry_pos
    npy_pos  <- best$npy_pos
    cat(sprintf("\nDRY@%d + NPxxY@%d (spacing %d)\n", dry_pos, npy_pos, best$spacing))
  } else {
    dry_pos <- dry_matches[[which.max(sapply(dry_matches, `[[`, "conservation"))]]$start
    npy_pos <- npy_matches[[which.max(sapply(npy_matches, `[[`, "conservation"))]]$start
    cat("WARNING: No valid DRY-NPxxY pair; using most conserved matches\n")
  }
} else {
  stop("ERROR: Could not find both DRY and NPxxY motifs.")
}

start_pos <- max(1, dry_pos - 120)
end_pos   <- min(ncol(aln_mat), npy_pos + 30)
cat(sprintf("7TM region: %d-%d (%d positions)\n\n", start_pos, end_pos, end_pos - start_pos + 1))

# ---- EXTRACT AND CLEAN ----
aln_7tm   <- aln_mat[, start_pos:end_pos]
gap_prop  <- colMeans(aln_7tm == "-")
aln_7tm   <- aln_7tm[, gap_prop < 0.8]
cat("After gap-column removal:", ncol(aln_7tm), "positions\n")

# ---- QUALITY CHECK ----
set.seed(42)
ids <- replicate(min(100, nrow(aln_7tm) * (nrow(aln_7tm) - 1) / 2),
                 sample(nrow(aln_7tm), 2), simplify = FALSE)
ident <- sapply(ids, function(i) {
  v <- aln_7tm[i[1], ] != "-" & aln_7tm[i[2], ] != "-"
  if (!any(v)) return(0)
  mean(aln_7tm[i[1], v] == aln_7tm[i[2], v])
})
cat(sprintf("Pairwise identity: %.1f%% (%.1f-%.1f%%)\n\n",
            mean(ident)*100, min(ident)*100, max(ident)*100))

# ---- SAVE ----
out <- AAStringSet(apply(aln_7tm, 1, paste, collapse = ""))
names(out) <- rownames(aln_7tm)
writeXStringSet(out, output_file)
cat("Saved:", length(out), "sequences,", width(out)[1], "positions ->", output_file, "\n")