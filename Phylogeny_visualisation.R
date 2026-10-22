library(ape); library(treeio); library(ggtree); library(ggplot2)
library(dplyr); library(tibble); library(tidytree); library(phangorn)

tree_dir      <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/revision_capitella/"
ml_tree_file  <- file.path(tree_dir, "opsin_newA4_rev2.treefile")
mb_tree_file  <- file.path(tree_dir, "opsin_newA4_rev2_mb.con.tre")
meta_file     <- file.path(tree_dir, "opsins_listA_Tree_Labels_rev.csv")

meta      <- read.csv(meta_file, stringsAsFactors = FALSE)
phy_ml    <- read.tree(ml_tree_file)
mb_treeio <- read.mrbayes(mb_tree_file)
phy_mb    <- mb_treeio@phylo

# Reroot on Anthozoa outgroup
outgroup <- grep("^Anthozoan", phy_ml$tip.label, value = TRUE)
# edgelabel = TRUE keeps SH-aLRT/UFBoot attached to the correct bipartitions after rerooting
if (length(outgroup)) phy_ml <- root(phy_ml, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)

n_tip  <- length(phy_ml$tip.label)
n_node <- phy_ml$Nnode

# ---- PARSE SUPPORT VALUES ----
support_df <- data.frame(node = 1:(n_tip + n_node),
                         SH = NA_real_, UF = NA_real_, PP = NA_real_)

if (length(phy_ml$node.label)) {
  for (k in seq_along(phy_ml$node.label)) {
    lbl <- phy_ml$node.label[k]
    if (!is.na(lbl) && lbl != "") {
      parts <- strsplit(lbl, "/")[[1]]
      if (length(parts) >= 2) {
        support_df$SH[n_tip + k] <- suppressWarnings(as.numeric(parts[1]))
        support_df$UF[n_tip + k] <- suppressWarnings(as.numeric(parts[2]))
      }
    }
  }
}

# Extract MB posterior probabilities, indexed by node number
n_tips_mb    <- length(phy_mb$tip.label)
mb_d         <- mb_treeio@data
mb_pp_values <- rep(NA_real_, n_tips_mb + phy_mb$Nnode)
int_rows     <- as.integer(mb_d$node) > n_tips_mb & !is.na(mb_d$prob)
mb_pp_values[as.integer(mb_d$node[int_rows])] <- as.numeric(mb_d$prob[int_rows]) * 100

# Transfer PP to ML nodes by matching bipartitions (clade or its complement,
# because the MrBayes consensus is rooted differently)
ml_tips_clean <- sub("/.*$", "", phy_ml$tip.label)
mb_tips_all   <- phy_mb$tip.label

for (ml_node in (n_tip + 1):(n_tip + n_node)) {
  ml_desc <- sort(ml_tips_clean[Descendants(phy_ml, ml_node, type = "tips")[[1]]])
  for (mb_idx in seq_len(phy_mb$Nnode)) {
    mb_node <- n_tips_mb + mb_idx
    mb_desc <- sort(phy_mb$tip.label[Descendants(phy_mb, mb_node, type = "tips")[[1]]])
    if ((setequal(ml_desc, mb_desc) || setequal(ml_desc, setdiff(mb_tips_all, mb_desc))) &&
        !is.na(mb_pp_values[mb_node])) {
      support_df$PP[ml_node] <- mb_pp_values[mb_node]; break
    }
  }
}
cat("PP matched:", sum(!is.na(support_df$PP)), "nodes\n")

# ---- HELPER FUNCTIONS ----
assign_taxa_by_species <- function(sp, taxa) {
  if (is.na(sp) || sp == "") return(taxa)
  s <- tolower(sp)
  if (grepl("owenia|malacoceros|capitella|platynereis|perionyx|helobdella", s)) return("Annelida")
  if (grepl("sepia|leptochiton|terebratalia|lingula|magallana", s))            return("Mollusca")
  if (grepl("schmidtea|prostheceraeus|schistosoma|dugesia", s))                return("Platyhelminthes")
  if (grepl("tricellaria", s))                                                  return("Bryozoa")
  if (grepl("tribolium|apis", s))                                               return("Arthropoda")
  if (grepl("nematostella|hydra", s))                                           return("Cnidaria")
  if (grepl("danio|xenopus|gallus|mus|homo", s))                               return("Chordata")
  taxa
}

consolidate_taxa <- function(taxa_vec, species_vec) {
  mapply(function(t, s) {
    if (t == "Invertebrate" || is.na(t) || t == "") t <- assign_taxa_by_species(s, t)
    switch(t,
           "Vertebrate"          = , "Chordata"              = "Chordata",
           "Annelid"             = , "Mfu"                   = "Annelida",
           "Platyhelmintes"      =                             "Platyhelminthes",
           "Artropode"           = , "Arthropod"             = "Arthropoda",
           "Mollusc"             = , "BrachiopodaMolluscs"   = "Mollusca",
           "Anthozoan"           =                             "Cnidaria",
           t)
  }, taxa_vec, species_vec, USE.NAMES = FALSE)
}

build_tip_metadata <- function(tree_obj, metadata) {
  tip_df         <- data.frame(label = tree_obj@phylo$tip.label,
                               node  = seq_len(length(tree_obj@phylo$tip.label)),
                               stringsAsFactors = FALSE)
  tip_df$Seq_ID  <- sub("/.*$", "", tip_df$label)
  idx            <- match(tip_df$Seq_ID, metadata$Seq_ID)
  tip_df$Species <- metadata$Species[idx]
  tip_df$Taxa    <- consolidate_taxa(
    replace(metadata$Taxa[idx], is.na(metadata$Taxa[idx]), "Unknown"),
    tip_df$Species)
  tip_df$Opsin_Type <- {
    ot <- replace(metadata$Opsin_Type[idx], is.na(metadata$Opsin_Type[idx]), "unknown")
    ot[ot == "otheropsins"] <- "opsins"
    ot[ot == "ropsins"]     <- "r-opsin"
    ot[ot == "copsins"]     <- "c-opsin"
    ot[ot == "xenopsins"]   <- "xenopsin"
    ot
  }
  tip_df$is_Mfu       <- grepl("^Mfu_", tip_df$Seq_ID)
  tip_df$display_label <- {
    dl <- metadata$Tree_labels[idx]
    dl[is.na(dl)] <- tip_df$Seq_ID[is.na(dl)]
    dl
  }
  tip_df
}

reconstruct_ancestral_taxa <- function(tree_obj, tip_df) {
  phy     <- tree_obj@phylo
  n_tip   <- length(phy$tip.label)
  all_taxa <- rep(NA_character_, n_tip + phy$Nnode)
  all_taxa[1:n_tip] <- setNames(tip_df$Taxa, tip_df$label)[phy$tip.label]
  post_nodes <- rev(unique(rev(phy$edge[postorder(phy), 1])))   # children before parents
  for (node in post_nodes) {
    ct <- all_taxa[phy$edge[phy$edge[, 1] == node, 2]]
    ct <- ct[!is.na(ct)]
    if (!length(ct)) next
    u <- unique(ct)
    if (length(u) == 1) { all_taxa[node] <- u; next }
    tc <- table(ct)
    all_taxa[node] <- if (max(tc) / length(ct) >= 0.75) names(which.max(tc)) else "Mixed"
  }
  all_taxa
}

# ---- BUILD TREE DATA ----
tr       <- as.treedata(phy_ml)
tr@data  <- as_tibble(support_df)
tip_df   <- build_tip_metadata(tr, meta)

complete_data <- data.frame(node = 1:(n_tip + n_node),
                            Taxa = reconstruct_ancestral_taxa(tr, tip_df),
                            stringsAsFactors = FALSE) %>%
  left_join(tip_df[, c("node", "Opsin_Type", "is_Mfu", "display_label")], by = "node")

tr@data <- tr@data %>%
  dplyr::select(node, SH, UF, PP) %>%
  left_join(complete_data, by = "node") %>%
  mutate(
    isTip = node <= n_tip,
    support_category = case_when(
      node == n_tip + 1                                                                 ~ "none",  # root
      !isTip & !is.na(SH) & !is.na(UF) & !is.na(PP) & SH >= 75 & UF >= 85 & PP >= 95 ~ "both_high",
      !isTip & !is.na(SH) & !is.na(UF) & !is.na(PP) & SH >= 75 & UF >= 85             ~ "ml_only",
      !isTip & !is.na(PP) & PP >= 95                                                    ~ "mb_only",
      TRUE ~ "none"
    )
  )

cat("Support summary:\n"); print(table(tr@data$support_category[!tr@data$isTip]))

# ---- NODES BY NAME (robust to renumbering) ----
lab_of   <- tip_df$display_label[match(tr@phylo$tip.label, tip_df$label)]
mrca_lab <- function(pattern) getMRCA(tr@phylo, tr@phylo$tip.label[grepl(pattern, lab_of)])
node_view <- mrca_lab("Schistosoma|Malacoceros fuliginosus r-opsin1")   # all r-opsins (main figure)
node_moll <- mrca_lab("rhodops")                                         # Mollusca r-opsins
node_plat <- mrca_lab("Schistosoma")                                     # Platyhelminthes r-opsins

# Support of key nodes (reported in the text); Seq_IDs as in the metadata CSV
r1 <- c("Annelid_ropsins_1", "Annelid_ropsins_8", "Mfu_opsins_4")                  # Pdu, Pex, Mfu r-opsin1
r3 <- c("Annelid_ropsins_2", "Annelid_ropsins_3", "Mfu_opsins_3", "Annelid_ropsins_6",
        "Annelid_ropsins_7", "Annelid_ropsins_9", "Annelid_ropsins_11",
        "Annelid_ropsins_13", "Annelid_ropsins_14", "Annelid_ropsins_15")           # r-opsin3 subclade
ow <- "Annelid_ropsins_16"                                                          # Owenia CAH1779377
cb_core <- c("Annelid_ropsins_4", "Annelid_ropsins_5", "Mfu_opsins_2", "Annelid_ropsins_23", "Annelid_ropsins_24")
key_nodes <- list(
  "Clade A"           = c(r1, ow, r3),
  "Owenia + r-opsin1" = c(r1, ow),
  "Clade B"           = c(cb_core, "Annelid_ropsins_17", "Annelid_ropsins_18", "Annelid_ropsins_19"),
  "Clade B core"      = cb_core,
  "Mfu r4 + Ct r2/r3" = c("Mfu_opsins_2", "Annelid_ropsins_23", "Annelid_ropsins_24"),
  "Ct r2 + r3"        = c("Annelid_ropsins_23", "Annelid_ropsins_24"))
for (k in names(key_nodes)) {
  hit <- ml_tips_clean %in% key_nodes[[k]]
  nd  <- getMRCA(tr@phylo, tr@phylo$tip.label[hit])
  r   <- tr@data[tr@data$node == nd, ]
  cat(sprintf("  %-18s tips %2d/%2d  SH-aLRT %5.1f  UFBoot %3.0f  PP %5.1f\n", k,
              length(Descendants(tr@phylo, nd, type = "tips")[[1]]), sum(hit), r$SH, r$UF, r$PP))
}

# ---- COLOR / SHAPE SCHEMES ----
taxa_colors <- c(Chordata = "#c2a5cf", Annelida = "#762a83", Arthropoda = "#35978f",
                 Mollusca = "#2166ac", Platyhelminthes = "#8c510a", Cnidaria = "#dfc27d",
                 Bryozoa = "#a6611a", Mixed = "gray60")
opsin_shapes <- c("c-opsin" = 16, "r-opsin" = 17, "xenopsin" = 15, "opsins" = 18)

# ---- SHARED LAYER BUILDER ----
support_layers <- function() {
  list(
    geom_point2(aes(subset = (!isTip & support_category == "both_high")),
                size = 2.0, shape = 21, fill = "black", color = "black", stroke = 1.5),
    geom_point2(aes(subset = (!isTip & support_category == "ml_only")),
                size = 2.5, shape = 21, fill = "#1b9e77", color = "black", stroke = 0.8),
    geom_point2(aes(subset = (!isTip & support_category == "mb_only")),
                size = 2.5, shape = 21, fill = "#d95f02", color = "black", stroke = 0.8)
  )
}

# ---- SUPPLEMENT FIGURE: FULL TREE ----
tree_data <- fortify(tr)
tree_data$show_symbol <- tree_data$isTip &
  !is.na(tree_data$Opsin_Type) &
  tree_data$Opsin_Type %in% c("c-opsin", "r-opsin", "xenopsin", "opsins")

mx <- max(tree_data$x, na.rm = TRUE)

p_supp <- ggtree(tr, layout = "rectangular", ladderize = TRUE) %<+% tree_data +
  aes(color = Taxa) +
  geom_tree(linewidth = 1.0) +
  scale_color_manual(values = taxa_colors, name = "Taxon", na.value = "gray70") +
  support_layers() +
  geom_tippoint(aes(subset = show_symbol, shape = Opsin_Type, fill = Taxa),
                size = 2.6, color = "black", stroke = 0.6) +
  scale_shape_manual(values = opsin_shapes, name = "Opsin type",
                     breaks = c("c-opsin", "r-opsin", "xenopsin", "opsins")) +
  scale_fill_manual(values = taxa_colors, guide = "none") +
  geom_tiplab(aes(subset = (is.na(is_Mfu) | !is_Mfu), label = display_label),
              size = 3, color = "gray20", fontface = "italic", offset = 0.1) +
  geom_tiplab(aes(subset = (!is.na(is_Mfu) & is_Mfu), label = display_label),
              size = 3.2, color = "#b2182b", fontface = "bold.italic", offset = 0.1) +
  geom_treescale(x = 0, y = -3, width = 0.5, fontsize = 3.5, linesize = 0.75, offset = 1) +
  annotate("text",  x = mx * 0.05, y = max(tree_data$y, na.rm=TRUE) * c(0.98, 0.95, 0.92, 0.89),
           label = c("Node support:", "ML & MB (UFBoot ≥85%, SH-aLRT ≥75%, PP ≥95%)",
                     "ML only (UFBoot ≥85%, SH-aLRT ≥75%)", "MB only (PP ≥95%)"),
           hjust = 0, size = c(3.5, 3, 3, 3), color = "gray20",
           fontface = c("bold", "plain", "plain", "plain")) +
  annotate("point", x = mx * 0.05,
           y = max(tree_data$y, na.rm=TRUE) * c(0.95, 0.92, 0.89),
           size = c(2.0, 2.5, 2.5), shape = 21,
           fill = c("black", "#1b9e77", "#d95f02"), color = "black",
           stroke = c(1.5, 0.8, 0.8)) +
  theme_tree2(bgcolor = "white") +
  theme(legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        panel.grid = element_blank(),
        plot.margin = margin(8, 8, 8, 8, "mm")) +
  xlim(NA, mx * 1.5)

for (fmt in c("pdf", "png", "svg"))
  ggsave(paste0("Fig_Opsin_phylogeny_supplement.", fmt), plot = p_supp, path = tree_dir,
         width = 13, height = 16, units = "in",
         device = if (fmt == "pdf") cairo_pdf else fmt,
         dpi = 300)

# ---- MAIN FIGURE: r-opsin part of the same tree, outgroups collapsed ----
p_main <- ggtree(tr, layout = "rectangular", ladderize = TRUE) %<+% tree_data + aes(color = Taxa)
p_main <- viewClade(p_main, node = node_view, xmax_adjust = 1.7)
p_main <- collapse(p_main, node = node_moll, mode = "max", fill = taxa_colors["Mollusca"],
                   color = taxa_colors["Mollusca"], alpha = 0.5, linewidth = 1.2)
p_main <- collapse(p_main, node = node_plat, mode = "max", fill = taxa_colors["Platyhelminthes"],
                   color = taxa_colors["Platyhelminthes"], alpha = 0.5, linewidth = 1.2)

d_main <- p_main$data
in_view <- d_main$node %in% c(node_view, Descendants(tr@phylo, node_view, type = "all"))
x0 <- d_main$x[d_main$node == node_view]
y_top <- max(d_main$y[in_view]); y_bot <- min(d_main$y[in_view])

p_main <- p_main +
  geom_tree(linewidth = 1.0) +
  scale_color_manual(values = taxa_colors, name = "Taxon", na.value = "gray70",
                     breaks = c("Annelida", "Mollusca", "Platyhelminthes")) +
  geom_cladelab(node = node_moll, label = "Mollusca r-opsins", color = taxa_colors["Mollusca"],
                fontsize = 3.5, offset = 0.2, barsize = 0, fontface = "bold") +
  geom_cladelab(node = node_plat, label = "Platyhelminthes r-opsins", color = taxa_colors["Platyhelminthes"],
                fontsize = 3.5, offset = 0.6, barsize = 0, fontface = "bold") +
  support_layers() +
  geom_tippoint(aes(subset = show_symbol, shape = Opsin_Type, fill = Taxa),
                size = 2.6, color = "black", stroke = 0.6) +
  scale_shape_manual(values = opsin_shapes, name = "Opsin type", breaks = "r-opsin") +
  scale_fill_manual(values = taxa_colors, guide = "none") +
  geom_tiplab(aes(subset = (is.na(is_Mfu) | !is_Mfu), label = display_label),
              size = 3, color = "gray20", fontface = "italic", offset = 0.05) +
  geom_tiplab(aes(subset = (!is.na(is_Mfu) & is_Mfu), label = display_label),
              size = 3.2, color = "#b2182b", fontface = "bold.italic", offset = 0.05) +
  geom_treescale(x = x0, y = y_bot - 0.8, width = 0.5, fontsize = 3.5, linesize = 0.75, offset = 0.3) +
  annotate("text", x = x0 + 0.05, y = y_top - 0:3,
           label = c("Node support:", "ML & MB (SH-aLRT \u226575%, UFBoot \u226585%, PP \u22650.95)",
                     "ML only", "MB only"),
           hjust = 0, size = c(3.5, 3, 3, 3), color = "gray20",
           fontface = c("bold", "plain", "plain", "plain")) +
  annotate("point", x = x0 + 0.02, y = y_top - 1:3, size = c(2.0, 2.5, 2.5), shape = 21,
           fill = c("black", "#1b9e77", "#d95f02"), color = "black", stroke = c(1.5, 0.8, 0.8)) +
  theme_tree(bgcolor = "white") +
  theme(legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = margin(8, 50, 8, 8, "mm"),
        plot.clip = "off")

for (fmt in c("pdf", "png", "svg"))
  ggsave(paste0("Fig_Main_ropsins.", fmt), plot = p_main, path = tree_dir,
         width = 18, height = 11, units = "in",
         device = if (fmt == "pdf") cairo_pdf else fmt,
         dpi = 300)library(ape); library(treeio); library(ggtree); library(ggplot2)
library(dplyr); library(tibble); library(tidytree); library(phangorn)

tree_dir      <- "/work/oleg/Anna_data/Tree_01.02.25/Last_update2/"
ml_tree_file  <- file.path(tree_dir, "opsin_newA4.fa.treefile")
mb_tree_file  <- file.path(tree_dir, "opsin_7tm_mb.con.tre")
meta_file     <- file.path(tree_dir, "opsins_listA_Tree_Labels_updated.csv")

meta      <- read.csv(meta_file, stringsAsFactors = FALSE)
phy_ml    <- read.tree(ml_tree_file)
mb_treeio <- read.mrbayes(mb_tree_file)
phy_mb    <- mb_treeio@phylo

# Reroot on Anthozoa outgroup
outgroup <- grep("^Anthozoan", phy_ml$tip.label, value = TRUE)
if (length(outgroup)) phy_ml <- root(phy_ml, outgroup = outgroup, resolve.root = TRUE)

n_tip  <- length(phy_ml$tip.label)
n_node <- phy_ml$Nnode

# ---- PARSE SUPPORT VALUES ----
support_df <- data.frame(node = 1:(n_tip + n_node),
                         SH = NA_real_, UF = NA_real_, PP = NA_real_)

if (length(phy_ml$node.label)) {
  for (k in seq_along(phy_ml$node.label)) {
    lbl <- phy_ml$node.label[k]
    if (!is.na(lbl) && lbl != "") {
      parts <- strsplit(lbl, "/")[[1]]
      if (length(parts) >= 2) {
        support_df$SH[n_tip + k] <- suppressWarnings(as.numeric(parts[1]))
        support_df$UF[n_tip + k] <- suppressWarnings(as.numeric(parts[2]))
      }
    }
  }
}

# Extract and transfer MB posterior probabilities by matching tip sets
n_tips_mb    <- length(phy_mb$tip.label)
mb_pp_values <- ifelse(seq_along(mb_treeio@data$prob) > n_tips_mb,
                       as.numeric(mb_treeio@data$prob) * 100, NA_real_)
ml_tips_clean <- sub("/.*$", "", phy_ml$tip.label)

for (ml_node in (n_tip + 1):(n_tip + n_node)) {
  ml_desc <- sort(ml_tips_clean[Descendants(phy_ml, ml_node, type = "tips")[[1]]])
  for (mb_idx in seq_len(phy_mb$Nnode)) {
    mb_node <- n_tips_mb + mb_idx
    mb_desc <- sort(phy_mb$tip.label[Descendants(phy_mb, mb_node, type = "tips")[[1]]])
    if (setequal(ml_desc, mb_desc) && !is.na(mb_pp_values[mb_node])) {
      support_df$PP[ml_node] <- mb_pp_values[mb_node]; break
    }
  }
}
cat("PP matched:", sum(!is.na(support_df$PP)), "nodes\n")

# ---- HELPER FUNCTIONS ----
assign_taxa_by_species <- function(sp, taxa) {
  if (is.na(sp) || sp == "") return(taxa)
  s <- tolower(sp)
  if (grepl("owenia|malacoceros|capitella|platynereis|perionyx|helobdella", s)) return("Annelida")
  if (grepl("sepia|leptochiton|terebratalia|lingula|magallana", s))            return("Mollusca")
  if (grepl("schmidtea|prostheceraeus|schistosoma|dugesia", s))                return("Platyhelminthes")
  if (grepl("tricellaria", s))                                                  return("Bryozoa")
  if (grepl("tribolium|apis", s))                                               return("Arthropoda")
  if (grepl("nematostella|hydra", s))                                           return("Cnidaria")
  if (grepl("danio|xenopus|gallus|mus|homo", s))                               return("Chordata")
  taxa
}

consolidate_taxa <- function(taxa_vec, species_vec) {
  mapply(function(t, s) {
    if (t == "Invertebrate" || is.na(t) || t == "") t <- assign_taxa_by_species(s, t)
    switch(t,
           "Vertebrate"          = , "Chordata"              = "Chordata",
           "Annelid"             = , "Mfu"                   = "Annelida",
           "Platyhelmintes"      =                             "Platyhelminthes",
           "Artropode"           = , "Arthropod"             = "Arthropoda",
           "Mollusc"             = , "BrachiopodaMolluscs"   = "Mollusca",
           "Anthozoan"           =                             "Cnidaria",
           t)
  }, taxa_vec, species_vec, USE.NAMES = FALSE)
}

build_tip_metadata <- function(tree_obj, metadata) {
  tip_df         <- data.frame(label = tree_obj@phylo$tip.label,
                               node  = seq_len(length(tree_obj@phylo$tip.label)),
                               stringsAsFactors = FALSE)
  tip_df$Seq_ID  <- sub("/.*$", "", tip_df$label)
  idx            <- match(tip_df$Seq_ID, metadata$Seq_ID)
  tip_df$Species <- metadata$Species[idx]
  tip_df$Taxa    <- consolidate_taxa(
    replace(metadata$Taxa[idx], is.na(metadata$Taxa[idx]), "Unknown"),
    tip_df$Species)
  tip_df$Opsin_Type <- {
    ot <- replace(metadata$Opsin_Type[idx], is.na(metadata$Opsin_Type[idx]), "unknown")
    ot[ot == "otheropsins"] <- "opsins"
    ot[ot == "ropsins"]     <- "r-opsin"
    ot[ot == "copsins"]     <- "c-opsin"
    ot[ot == "xenopsins"]   <- "xenopsin"
    ot
  }
  tip_df$is_Mfu       <- grepl("^Mfu_", tip_df$Seq_ID)
  tip_df$display_label <- {
    dl <- metadata$Tree_labels[idx]
    dl[is.na(dl)] <- tip_df$Seq_ID[is.na(dl)]
    dl
  }
  tip_df
}

reconstruct_ancestral_taxa <- function(tree_obj, tip_df) {
  phy     <- tree_obj@phylo
  n_tip   <- length(phy$tip.label)
  all_taxa <- rep(NA_character_, n_tip + phy$Nnode)
  all_taxa[1:n_tip] <- setNames(tip_df$Taxa, tip_df$label)[phy$tip.label]
  for (node in (n_tip + 1):(n_tip + phy$Nnode)) {
    ct <- all_taxa[phy$edge[phy$edge[, 1] == node, 2]]
    ct <- ct[!is.na(ct)]
    if (!length(ct)) next
    u <- unique(ct)
    if (length(u) == 1) { all_taxa[node] <- u; next }
    tc <- table(ct)
    all_taxa[node] <- if (max(tc) / length(ct) >= 0.75) names(which.max(tc)) else "Mixed"
  }
  all_taxa
}

# ---- BUILD TREE DATA ----
tr       <- as.treedata(phy_ml)
tr@data  <- as_tibble(support_df)
tip_df   <- build_tip_metadata(tr, meta)

complete_data <- data.frame(node = 1:(n_tip + n_node),
                            Taxa = reconstruct_ancestral_taxa(tr, tip_df),
                            stringsAsFactors = FALSE) %>%
  left_join(tip_df[, c("node", "Opsin_Type", "is_Mfu", "display_label")], by = "node")

tr@data <- tr@data %>%
  dplyr::select(node, SH, UF, PP) %>%
  left_join(complete_data, by = "node") %>%
  mutate(
    isTip = node <= n_tip,
    support_category = case_when(
      !isTip & !is.na(SH) & !is.na(UF) & !is.na(PP) & SH >= 75 & UF >= 85 & PP >= 95 ~ "both_high",
      !isTip & !is.na(SH) & !is.na(UF) & SH >= 75 & UF >= 85                          ~ "ml_only",
      !isTip & !is.na(PP) & PP >= 95                                                    ~ "mb_only",
      TRUE ~ "none"
    )
  )

cat("Support summary:\n"); print(table(tr@data$support_category[!tr@data$isTip]))

# ---- COLOR / SHAPE SCHEMES ----
taxa_colors <- c(Chordata = "#c2a5cf", Annelida = "#762a83", Arthropoda = "#35978f",
                 Mollusca = "#2166ac", Platyhelminthes = "#8c510a", Cnidaria = "#dfc27d",
                 Bryozoa = "#a6611a", Mixed = "gray60")
opsin_shapes <- c("c-opsin" = 16, "r-opsin" = 17, "xenopsin" = 15, "opsins" = 18)

# ---- SHARED LAYER BUILDER ----
support_layers <- function() {
  list(
    geom_point2(aes(subset = (!isTip & support_category == "both_high")),
                size = 2.0, shape = 21, fill = "black", color = "black", stroke = 1.5),
    geom_point2(aes(subset = (!isTip & support_category == "ml_only")),
                size = 2.5, shape = 21, fill = "#1b9e77", color = "black", stroke = 0.8),
    geom_point2(aes(subset = (!isTip & support_category == "mb_only")),
                size = 2.5, shape = 21, fill = "#d95f02", color = "black", stroke = 0.8)
  )
}

# ---- SUPPLEMENT FIGURE: FULL TREE ----
tree_data <- fortify(tr)
tree_data$show_symbol <- tree_data$isTip &
  !is.na(tree_data$Opsin_Type) &
  tree_data$Opsin_Type %in% c("c-opsin", "r-opsin", "xenopsin", "opsins")

mx <- max(tree_data$x, na.rm = TRUE)

p_supp <- ggtree(tr, layout = "rectangular", ladderize = TRUE) %<+% tree_data +
  aes(color = Taxa) +
  geom_tree(linewidth = 1.0) +
  scale_color_manual(values = taxa_colors, name = "Taxon", na.value = "gray70") +
  support_layers() +
  geom_tippoint(aes(subset = show_symbol, shape = Opsin_Type, fill = Taxa),
                size = 2.6, color = "black", stroke = 0.6) +
  scale_shape_manual(values = opsin_shapes, name = "Opsin type",
                     breaks = c("c-opsin", "r-opsin", "xenopsin", "opsins")) +
  scale_fill_manual(values = taxa_colors, guide = "none") +
  geom_tiplab(aes(subset = (is.na(is_Mfu) | !is_Mfu), label = display_label),
              size = 3, color = "gray20", fontface = "italic", offset = 0.1) +
  geom_tiplab(aes(subset = (!is.na(is_Mfu) & is_Mfu), label = display_label),
              size = 3.2, color = "#b2182b", fontface = "bold.italic", offset = 0.1) +
  geom_treescale(x = 0, y = -3, width = 0.5, fontsize = 3.5, linesize = 0.75, offset = 1) +
  annotate("text",  x = mx * 0.05, y = max(tree_data$y, na.rm=TRUE) * c(0.98, 0.95, 0.92, 0.89),
           label = c("Node support:", "ML & MB (UFBoot ≥85%, SH-aLRT ≥75%, PP ≥95%)",
                     "ML only (UFBoot ≥85%, SH-aLRT ≥75%)", "MB only (PP ≥95%)"),
           hjust = 0, size = c(3.5, 3, 3, 3), color = "gray20",
           fontface = c("bold", "plain", "plain", "plain")) +
  annotate("point", x = mx * 0.05,
           y = max(tree_data$y, na.rm=TRUE) * c(0.95, 0.92, 0.89),
           size = c(2.0, 2.5, 2.5), shape = 21,
           fill = c("black", "#1b9e77", "#d95f02"), color = "black",
           stroke = c(1.5, 0.8, 0.8)) +
  theme_tree2(bgcolor = "white") +
  theme(legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        panel.grid = element_blank(),
        plot.margin = margin(8, 8, 8, 8, "mm")) +
  xlim(NA, mx * 1.5)

for (fmt in c("pdf", "png", "svg"))
  ggsave(paste0("Fig_Opsin_phylogeny_supplement.", fmt), plot = p_supp, path = tree_dir,
         width = 13, height = 16, units = "in",
         device = if (fmt == "pdf") cairo_pdf else fmt,
         dpi = if (fmt == "png") 600 else NULL)

# ---- MAIN FIGURE: ANNELID-FOCUSED (node 72, collapsed outgroups) ----
collapse_info <- data.frame(
  node  = c(90, 97),
  label = c("Mollusca r-opsins", "Platyhelminthes r-opsins"),
  color = c(taxa_colors["Mollusca"], taxa_colors["Platyhelminthes"]),
  stringsAsFactors = FALSE
)

p_main <- ggtree(tr, layout = "rectangular", ladderize = TRUE) %<+% tr@data +
  aes(color = Taxa)
p_main <- viewClade(p_main, node = 72)

for (i in seq_len(nrow(collapse_info)))
  p_main <- collapse(p_main, node = collapse_info$node[i], mode = "max",
                     fill = collapse_info$color[i], color = collapse_info$color[i],
                     alpha = 0.5, linewidth = 1.2)

p_main <- p_main +
  geom_tree(linewidth = 1.0) +
  scale_color_manual(values = taxa_colors, name = "Taxon", na.value = "gray70") +
  geom_cladelab(node = 90, label = "Mollusca r-opsins",       color = taxa_colors["Mollusca"],
                fontsize = 3.5, offset = 0.2, barsize = 0, fontface = "bold") +
  geom_cladelab(node = 97, label = "Platyhelminthes r-opsins", color = taxa_colors["Platyhelminthes"],
                fontsize = 3.5, offset = 0.2, barsize = 0, fontface = "bold") +
  support_layers() +
  geom_tippoint(aes(subset = !is.na(Opsin_Type) & Opsin_Type %in% names(opsin_shapes),
                    shape = Opsin_Type, fill = Taxa),
                size = 2.6, color = "black", stroke = 0.6) +
  scale_shape_manual(values = opsin_shapes, name = "Opsin type",
                     breaks = c("c-opsin", "r-opsin", "xenopsin", "opsins")) +
  scale_fill_manual(values = taxa_colors, guide = "none") +
  geom_tiplab(aes(subset = (is.na(is_Mfu) | !is_Mfu), label = display_label),
              size = 3, color = "gray20", fontface = "italic", offset = 0.2) +
  geom_tiplab(aes(subset = (!is.na(is_Mfu) & is_Mfu), label = display_label),
              size = 3.2, color = "#b2182b", fontface = "bold.italic", offset = 0.2) +
  geom_treescale(x = 0, y = -3, width = 0.5, fontsize = 3.5, linesize = 0.75, offset = 1) +
  xlim(NA, max(p_main$data$x, na.rm = TRUE) * 1.8) +
  theme_tree2(bgcolor = "white") +
  theme(legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        panel.grid = element_blank(),
        plot.margin = margin(8, 50, 8, 8, "mm"),
        plot.clip = "off")

for (fmt in c("pdf", "png", "svg"))
  ggsave(paste0("Fig_A_Main.", fmt), plot = p_main, path = tree_dir,
         width = 18, height = 11, units = "in",
         device = if (fmt == "pdf") cairo_pdf else fmt,
         dpi = if (fmt == "png") 600 else NULL)
