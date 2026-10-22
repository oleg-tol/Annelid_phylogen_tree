library(ape)
library(phangorn)
library(Biostrings)

tree_dir       <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/revision_capitella/"
alignment_file <- file.path(tree_dir, "opsin_newA4_rev2.fa")
ml_tree_file   <- file.path(tree_dir, "opsin_newA4_rev2.treefile")
meta_file      <- file.path(tree_dir, "opsins_listA_Tree_Labels_rev.csv")

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

# r-opsin subtree, found by tip label (MRCA of Schistosoma r-opsins and Mfu r-opsin1)
tip_lab     <- meta$Tree_labels[match(sub("/.*$", "", phy_ml$tip.label), meta$Seq_ID)]
ropsin_node <- getMRCA(phy_ml, phy_ml$tip.label[grepl("Schistosoma|Malacoceros fuliginosus r-opsin1", tip_lab)])
ropsin_tips <- phy_ml$tip.label[Descendants(phy_ml, ropsin_node, type = "tips")[[1]]]
labeled_alignment(alignment, ropsin_tips,
                  meta, file.path(tree_dir, "opsin_alignment_labeled_ropsins.fa"))
