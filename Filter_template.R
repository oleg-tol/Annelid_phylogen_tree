#!/usr/bin/env Rscript
library(Biostrings)

work_dir      <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/"
template_file <- file.path(work_dir, "opsin_newA3.fa")
source_file   <- file.path(work_dir, "opsins_all.fasta")
output_file   <- file.path(work_dir, "listA_filtered4.fa")

template_ids <- names(readAAStringSet(template_file))
source_seqs  <- readAAStringSet(source_file)
source_ids   <- names(source_seqs)

cat("Template:", length(template_ids), "IDs | Source:", length(source_seqs), "sequences\n")

matches <- source_ids %in% template_ids
if (!any(matches)) stop("No matching IDs found — check ID formatting in both files.")

filtered <- source_seqs[matches]
cat("Matched:", length(filtered), "sequences\n")

missing <- template_ids[!template_ids %in% source_ids]
if (length(missing)) cat("WARNING:", length(missing), "template IDs not in source:\n",
                         paste(" ", missing, collapse = "\n"), "\n")

writeXStringSet(filtered, output_file)
cat("Saved ->", output_file, "\n")
cat(paste(sprintf("  %s: %d aa", names(filtered), width(filtered)), collapse = "\n"), "\n")