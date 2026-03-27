library(Biostrings)

work_dir        <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/"
input_file      <- file.path(work_dir, "opsin_newA4.fa")
output_filtered <- file.path(work_dir, "opsin_newA4_80pct.fa")

aln     <- readAAStringSet(input_file)
aln_mat <- as.matrix(aln)
cat("Sequences:", length(aln), "| Length:", width(aln)[1], "\n")

occupancy <- colMeans(!(aln_mat %in% c("-", "X", "?")))
keep_cols <- occupancy >= 0.80
cat("Columns kept:", sum(keep_cols), "| Removed:", sum(!keep_cols), "\n")

out <- AAStringSet(apply(aln_mat[, keep_cols], 1, paste, collapse = ""))
names(out) <- rownames(aln_mat)
writeXStringSet(out, output_filtered)
cat("Saved:", width(out)[1], "positions ->", output_filtered, "\n")