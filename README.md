# Annelid_phylogen_tree

Phylogenetic analysis of annelid opsins: alignment pre-processing, 7TM domain extraction, maximum-likelihood / Bayesian trees with dual support values, a constrained-topology (AU) test, and analysis of residues distinguishing r-opsin1 from r-opsin3.

## Modules

| File | Description |
|------|-------------|
| `Filter_template.R` | Filter a large FASTA by retaining only sequences whose IDs match a template FASTA; reports missing IDs |
| `Cleaning_sparse_columns.R` | Remove poorly occupied alignment columns (default ≥80% occupancy) |
| `7TM_extraction.R` | Extract the 7-transmembrane domain from a full-length opsin alignment using conserved GPCR motifs (DRY/ERY and NPxxY as primary anchors); removes columns with >80% gaps; reports pairwise identity |
| `Seq_alignment_renaming.R` | Replace sequence IDs with tree labels from the metadata CSV; order sequences as in the ML tree; write Jalview-ready FASTA for the full tree and the r-opsin subtree |
| `Phylogeny_visualisation.R` | Load IQ-TREE ML and MrBayes consensus trees; parse SH-aLRT, UFBoot and posterior probabilities; transfer PP to ML nodes by matching bipartitions; print support of key nodes; produce the supplementary full tree and the main r-opsin figure |
| `Diagnostic_residues.py` | Identify alignment positions distinguishing r-opsin1 from r-opsin3 orthologs, classify states as ancestral (shared with mollusc rhodopsins) or lineage-specific, and plot them with *M. fuliginosus* and *O. fusiformis* highlighted (Sup. Fig. 3) |
| `Revision_steps.sh` | Commands for the revision: addition of two *C. teleta* r-opsins, removal of redundant *H. robusta* entries, ML and Bayesian re-inference, AU test, residue figure |

## Pipeline

```
opsins_all.fasta
  → Filter_template.R          →  listA_filtered.fa
                                   [align with MAFFT E-INS-i]
full-length alignment
  → 7TM_extraction.R           →  7TM block (columns with >80% gaps removed)
                                   [re-align with MAFFT E-INS-i → opsin_newA4.fa, 52 seqs × 319 positions]

Revision (Revision_steps.sh)
opsin_newA4.fa
  → MAFFT --add --keeplength   →  + C. teleta r-opsin2/3 (AXB38864.1, AXB38865.1)
  → remove redundant H. robusta entries (XP_009022296.1, AID66633.1)
                               →  opsin_newA4_rev2.fa  (52 seqs × 319 positions)
  → IQ-TREE 2.4.0  LG+F+R8, UFBoot 2000, SH-aLRT 1000
  → MrBayes 3.2.7a  LG+Γ6, 2 runs × 4 chains, 1.5×10⁷ generations, 25% burn-in
  → IQ-TREE -g + AU test       →  O. fusiformis CAH1779377 constrained as sister to r-opsin1 + r-opsin3

ML treefile + MrBayes .con.tre + metadata CSV
  → Phylogeny_visualisation.R  →  Fig_Opsin_phylogeny_supplement.pdf/.png/.svg  (full tree)
                                   Fig_Main_ropsins.pdf/.png/.svg               (r-opsin part)
opsin_newA4_rev2.fa + metadata CSV
  → Diagnostic_residues.py     →  SupFig3_diagnostic.png/.svg + _table.tsv
```

## Key parameters

| Parameter | Value | Location |
|-----------|-------|----------|
| 7TM upstream / downstream buffer | 120 from DRY / 30 from NPxxY | `7TM_extraction.R` |
| Gap column removal threshold | >80% gaps | `7TM_extraction.R` |
| ML support threshold | SH-aLRT ≥75 and UFBoot ≥85 | `Phylogeny_visualisation.R` |
| Bayesian support threshold | PP ≥0.95 | `Phylogeny_visualisation.R` |
| Node matching ML ↔ MB | Identical bipartition (clade or its complement) | `Phylogeny_visualisation.R` |
| Outgroup for rooting | `Anthozoan*` tips (*Nematostella vectensis*) | `Phylogeny_visualisation.R` |
| Main-figure clade | MRCA of all r-opsins, found by tip label (not node number) | `Phylogeny_visualisation.R` |
| Diagnostic position | r-opsin1 and r-opsin3 residue sets do not overlap (one exception allowed) | `Diagnostic_residues.py` |
| Ancestral state | residue shared with mollusc rhodopsins | `Diagnostic_residues.py` |

## Requirements

```r
BiocManager::install("Biostrings")
install.packages(c("ape", "treeio", "ggtree", "ggplot2",
                   "dplyr", "tibble", "tidytree", "phangorn", "seqinr"))
```

Python ≥3.8 with `matplotlib`.

External tools: [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.490, [IQ-TREE 2](http://www.iqtree.org/) v2.4.0, [MrBayes](https://nbisweden.github.io/MrBayes/) v3.2.7a

Edit the `work_dir`, `tree_dir`, and file path variables at the top of each script before running.Annelid_phylogen_tree

Phylogenetic analysis of annelid opsins cover alignment pre-processing, 7TM domain extraction, and visualisation of maximum likelihood / Bayesian consensus trees with dual support values.
Modules
File 	Description
Filter_template.R 	Filter a large FASTA by retaining only sequences whose IDs match a template FASTA; reports missing IDs
Cleaning_sparse_columns.R 	Remove poorly occupied alignment columns (default ≥80% occupancy); output ready for re-alignment with MAFFT
7TM_extraction.R 	Extract the 7-transmembrane domain from a full-length opsin alignment using conserved GPCR motifs (DRY/ERY and NPxxY as primary anchors, CWxP and TM1–TM4 motifs as secondary); removes gap-rich columns and reports pairwise identity and conservation statistics
Seq_alignment_renaming.R 	Replace sequence IDs with tree labels from metadata CSV; reorder sequences to match ML tree structure; write Jalview-ready FASTA for the full tree and a focal subtree
Phylogeny_visualisation.R 	Load IQ-TREE ML and MrBayes consensus trees; parse SH-aLRT, UFBoot, and posterior probability support; match internal nodes between trees by descendant tip sets; reconstruct ancestral taxon assignments; produce figures
Pipeline

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
                                       

Key parameters
Parameter 	Value 	Location
Column occupancy threshold 	≥80% non-gap 	Cleaning_sparse_columns.R
7TM motif conservation minimum 	0.4 (low-flexibility), 0.3 (medium), 0.2 (high) 	7TM_extraction.R
7TM upstream buffer from DRY 	120 positions 	7TM_extraction.R
7TM downstream buffer from NPxxY 	30 positions 	7TM_extraction.R
Gap column removal threshold 	>80% gaps 	7TM_extraction.R
ML support threshold (SH-aLRT) 	≥75 	Phylogeny_visualisation.R
ML support threshold (UFBoot) 	≥85 	Phylogeny_visualisation.R
Bayesian support threshold (PP) 	≥95% 	Phylogeny_visualisation.R
Node matching strategy 	Exact descendant tip set identity 	Phylogeny_visualisation.R
Outgroup for rooting 	Anthozoan* tip labels 	Phylogeny_visualisation.R
Focal clade (annelid-focused view) 	Node 72 subtree 	Phylogeny_visualisation.R
Requirements

# Bioconductor
BiocManager::install(c("Biostrings"))

# CRAN
install.packages(c("ape", "treeio", "ggtree", "ggplot2",
                   "dplyr", "tibble", "tidytree", "phangorn", "seqinr"))

External tools: IQ-TREE2, MrBayes, MAFFT

Edit the work_dir, tree_dir, and file path variables at the top of each script before running.

External tools: [IQ-TREE2](http://www.iqtree.org/), [MrBayes](https://nbisweden.github.io/MrBayes/), [MAFFT](https://mafft.cbrc.jp/alignment/software/)

Edit the `work_dir`, `tree_dir`, and file path variables at the top of each script before running.
