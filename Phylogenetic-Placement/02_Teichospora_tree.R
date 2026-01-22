#Read in IQTree, improve tip labels, root, color and save as PDF
#Figure 6

#### Load Packages ####
library(ape)
library(ggtree)

##### Make plots ####
# 1–4. Read, map, update labels, identify type strains
tree <- read.tree("Teichospora_LSU_trimmed.fasta.treefile")
label_map <- c(
  "Teichospora_austroafricana_CBS_119330_T"      = "Teichospora austroafricana CBS 119330T",
  "Teichospora_austroafricana_CBS_122674"        = "Teichospora austroafricana CBS 122674",
  "Teichospora_proteae_CBS_122675_T"             = "Teichospora proteae CBS 122675T",
  "Teichospora_claviformis_GKM_1210_T"           = "Teichospora claviformis GKM 1210T",
  "Horn_10"                                      = "Horn 10",
  "Teichospora_kingiae_CPC_29104_T"              = "Teichospora kingiae CPC 29104T",
  "Paulkirkia_arundinis_MFLU_130315_T"           = "Paulkirkia arundinis MFLU 13-0315T",
  "Floricola_striata_JK_5678I"                   = "Floricola striata JK 5678I",
  "Floricola_striata_JK_5603K"                   = "Floricola striata JK 5603K",
  "Pseudomisturatosphaeria_cruciformis_SMH5151T" = "Pseudomisturatosphaeria cruciformis SMH 5151T",
  "Asymmetrispora_mariae_CBS140732T"             = "Asymmetrispora mariae CBS 140732T",
  "Ramusculicola_thailandica_MFLUCC130284T"      = "Ramusculicola thailandica MFLUCC 13-0284T",
  "Misturatosphaeria_radicans_ATCC42522T"        = "Misturatosphaeria radicans ATCC 42522T",
  "Aurantiascoma_minimum_GKM_169N_T"             = "Aurantiascoma minimum GKM 169NT",
  "Aurantiascoma_nephelii_CPC_27539_T"           = "Aurantiascoma nephelii CPC 27539T",
  "Aurantiascoma_quercus_CBS_143396_T"           = "Aurantiascoma quercus CBS 143396T",
  "Pseudoaurantiascoma_kenyensis_GKM1195T"       = "Pseudoaurantiascoma kenyensis GKM 1195T",
  "Decaisnella_formosa_BCC25617"                 = "Decaisnella formosa BCC 25617",
  "Westerdykella_cylindrica_CBS454.72"           = "Westerdykella cylindrica CBS 454.72T",
  "Magnibotryascoma_kunmingense_KUMCC_200254_T"  = "Magnibotryascoma kunmingense KUMCC 20-0254T",
  "Magnibotryascoma_acaciae_CPC_24801_T"         = "Magnibotryascoma acaciae CPC 24801T",
  "Teichospora_grandicipis_CPC_1852"             = "Teichospora grandicipis CPC 1852",
  "Teichospora_grandicipis_CBS_111702_T"         = "Teichospora grandicipis CBS 111702T",
  "Teichospora_trabicola_C_141"                  = "Teichospora trabicola C 141",
  "Teichospora_trabicola_CBS_140730_T"           = "Teichospora trabicola CBS 140730T",
  "Teichospora_trabicola_C_157"                  = "Teichospora trabicola C 157"
)
tree$tip.label   <- unname(label_map[tree$tip.label])
bold_idx         <- grepl("T$", tree$tip.label)
labels_no_T      <- gsub("T$", "", tree$tip.label)
fontface_vec     <- ifelse(bold_idx, "bold", "plain")
tree$tip.label   <- trimws(labels_no_T)

# --- Adjust Ramusculicola branch length ---
target_tip <- "Ramusculicola thailandica MFLUCC 13-0284"
tip_num <- which(tree$tip.label == target_tip)
target_edge_row <- which(tree$edge[,2] == tip_num)
cat("Old Ramusculicola branch length:", tree$edge.length[target_edge_row], "\n")
tree$edge.length[target_edge_row] <- 0.02   # Super short for plotting!
cat("New Ramusculicola branch length:", tree$edge.length[target_edge_row], "\n")

# 5. Root and clean the "Root" label
outgroups   <- c("Westerdykella cylindrica CBS 454.72", "Decaisnella formosa BCC 25617")
tree_rooted <- root(tree, outgroup = outgroups, resolve.root = TRUE)
if (!is.null(tree_rooted$node.label)) {
  tree_rooted$node.label[tree_rooted$node.label == "Root"] <- ""
}

# 6. Add fontface to tip annotations
tipdf <- data.frame(
  label = tree_rooted$tip.label,
  fontface = fontface_vec,
  stringsAsFactors = FALSE
)

# 7. Create base plot (with blue highlight for node 32 if desired)
p <- ggtree(tree_rooted) %<+% tipdf +
  geom_hilight(node = 32, fill = "blue", alpha = 0.18) +
  geom_tiplab(aes(fontface = fontface))

tree_data <- p$data

# 8. Find root coordinates (for the root segment)
root_node <- Ntip(tree_rooted) + 1
root_x    <- tree_data$x[tree_data$node == root_node]
root_y    <- tree_data$y[tree_data$node == root_node]

# 9. Compose final plot
final_plot <- p +
  geom_text(
    data = subset(tree_data, !isTip),
    aes(x = x, y = y, label = label),
    hjust = -0.2, vjust = -0.7, size = 3, na.rm = TRUE
  ) +
  annotate(
    "segment",
    x = root_x, y = root_y,
    xend = root_x - 0.01 * max(tree_data$x),
    yend = root_y,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_treescale(
    x = min(tree_data$x) + 0.0003,
    y = min(tree_data$y) - 1,
    width = 0.01,
    fontsize = 4
  )
# Add extra x padding
final_plot <- final_plot +
  expand_limits(x = max(tree_data$x) + 0.2 * diff(range(tree_data$x)))

# Or try: xlim(min(tree_data$x), max(tree_data$x) + 0.2 * diff(range(tree_data$x)))
#         Width can be tuned as needed.

# Save with ggsave, allowing big plot:
ggsave("teichospora_tree.pdf",
       plot = final_plot,
       width = 20, height = 15, units = "in",
       limitsize = FALSE)
