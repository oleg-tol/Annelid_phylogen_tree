library(ape)
library(Biostrings)

tree_dir       <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/"
alignment_file <- file.path(tree_dir, "opsin_newA4.fa")
ml_tree_file   <- file.path(tree_dir, "opsin_newA4.fa.treefile")
meta_file      <- file.path(tree_dir, "opsins_listA_Tree_Labels_updated.csv")

alignment <- readAAStringSet(alignment_file)
meta      <- read.csv(meta_file, stringsAsFactors = FALSE)
phy_ml    <- read.tree(ml_tree_file)

outgroup <- grep("^Anthozoan", phy_ml$tip.label, value = TRUE)
if (length(outgroup)) phy_ml <- root(phy_ml, outgroup = outgroup, resolve.root = TRUE)

labeled_alignment <- function(alignment, tree_tips, metadata, output_file) {
  seq_ids    <- sub("/.*$", "", tree_tips)
  idx        <- match(seq_ids, metadata$Seq_ID)
  labels     <- gsub(" ", "_", ifelse(is.na(metadata$Tree_labels[idx]), seq_ids, metadata$Tree_labels[idx]))
  
  align_idx  <- match(seq_ids, sub("/.*$", "", names(alignment)))
  valid      <- !is.na(align_idx)
  if (any(!valid)) cat("WARNING:", sum(!valid), "tree tips missing from alignment\n")
  
  out        <- alignment[align_idx[valid]]
  names(out) <- labels[valid]
  writeXStringSet(out, output_file, format = "fasta")
  cat(length(out), "sequences ->", output_file, "\n")
  invisible(labels[valid])
}

labeled_alignment(alignment, phy_ml$tip.label,
                  meta, file.path(tree_dir, "opsin_alignment_labeled_full.fa"))

node72_tips <- phy_ml$tip.label[Descendants(phy_ml, 72, type = "tips")[[1]]]
labeled_alignment(alignment, node72_tips,
                  meta, file.path(tree_dir, "opsin_alignment_labeled_node72.fa"))