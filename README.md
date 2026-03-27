# Annelid_phylogen_tree

Phylogenetic analysis of annelid opsins cover alignment pre-processing, 7TM domain extraction, and visualisation of maximum likelihood / Bayesian consensus trees with dual support values.

## Modules

| File | Description |
|------|-------------|
| `Filter_template.R` | Filter a large FASTA by retaining only sequences whose IDs match a template FASTA; reports missing IDs |
| `Cleaning_sparse_columns.R` | Remove poorly occupied alignment columns (default ≥80% occupancy); output ready for re-alignment with MAFFT |
| `7TM_extraction.R` | Extract the 7-transmembrane domain from a full-length opsin alignment using conserved GPCR motifs (DRY/ERY and NPxxY as primary anchors, CWxP and TM1–TM4 motifs as secondary); removes gap-rich columns and reports pairwise identity and conservation statistics |
| `Seq_alignment_renaming.R` | Replace sequence IDs with tree labels from metadata CSV; reorder sequences to match ML tree structure; write Jalview-ready FASTA for the full tree and a focal subtree |
| `Phylogeny_visualisation.R` | Load IQ-TREE ML and MrBayes consensus trees; parse SH-aLRT, UFBoot, and posterior probability support; match internal nodes between trees by descendant tip sets; reconstruct ancestral taxon assignments; produce figures |

## Pipeline

```
opsins_all.fasta
  → Filter_template.R          →  listA_filtered.fa   (sequences matching template IDs)

listA_filtered.fa
  → Cleaning_sparse_columns.R  →  opsin_*_filtered.fa  (columns ≥80% occupied)
                                   [re-align externally with MAFFT E-INS-i]

opsin_alignment.fa
  → 7TM_extraction.R           →  opsin_7tm_extracted.fa  (7TM domain only, gappy cols removed)
                                   [tree inference with IQ-TREE and MrBayes run externally]

opsin_alignment.fa  +  ML tree  +  metadata CSV
  → Seq_alignment_renaming.R   →  opsin_alignment_labeled_full.fa
                                   opsin_alignment_labeled_node72.fa  (Jalview-ready)

ML treefile  +  MrBayes .con.tre  +  metadata CSV
  → Phylogeny_visualisation.R  →  Fig_Opsin_phylogeny_ML_MB.pdf/.png/.svg
                                       
```

## Key parameters

| Parameter | Value | Location |
|-----------|-------|----------|
| Column occupancy threshold | ≥80% non-gap | `Cleaning_sparse_columns.R` |
| 7TM motif conservation minimum | 0.4 (low-flexibility), 0.3 (medium), 0.2 (high) | `7TM_extraction.R` |
| 7TM upstream buffer from DRY | 120 positions | `7TM_extraction.R` |
| 7TM downstream buffer from NPxxY | 30 positions | `7TM_extraction.R` |
| Gap column removal threshold | >80% gaps | `7TM_extraction.R` |
| ML support threshold (SH-aLRT) | ≥75 | `Phylogeny_visualisation.R` |
| ML support threshold (UFBoot) | ≥85 | `Phylogeny_visualisation.R` |
| Bayesian support threshold (PP) | ≥95% | `Phylogeny_visualisation.R` |
| Node matching strategy | Exact descendant tip set identity | `Phylogeny_visualisation.R` |
| Outgroup for rooting | `Anthozoan*` tip labels | `Phylogeny_visualisation.R` |
| Focal clade (annelid-focused view) | Node 72 subtree | `Phylogeny_visualisation.R` |

## Requirements

```r
# Bioconductor
BiocManager::install(c("Biostrings"))

# CRAN
install.packages(c("ape", "treeio", "ggtree", "ggplot2",
                   "dplyr", "tibble", "tidytree", "phangorn", "seqinr"))
```

External tools: [IQ-TREE2](http://www.iqtree.org/), [MrBayes](https://nbisweden.github.io/MrBayes/), [MAFFT](https://mafft.cbrc.jp/alignment/software/)

Edit the `work_dir`, `tree_dir`, and file path variables at the top of each script before running.
