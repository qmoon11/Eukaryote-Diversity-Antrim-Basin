#This scripts completes the Community analysis of the 18s and ITS Amplicon data

#### Load Packages ####
library(Polychrome)
library(ggtree)
library(ape)
library(VennDiagram)
library(scales)
library(metacoder)
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(ape)
library(picante)
library(tidyverse)
library(ALDEx2)
library(ggsci)
library(RColorBrewer)
library(Polychrome)
library(patchwork)
library(geosphere)
library(purrr)
library(ComplexUpset)
library(ggVennDiagram)
library(pals)
library(metacoder)

setwd("/Users/quinnmoon/Downloads/Antrim_Microbiome")
#### Read in OTU and Taxonomy Tables ####

# PCoA and NMDS grouping, and combining of technical replicates completed in Isolates_taxonomy.R

# ---------- ITS2 DATA ----------

# Read ITS2 OTU table (post-contaminant removal and replicate pooling)
# Contaminated OTUs were identified and removed based on cultured isolates (see Isolates_taxonomy.R)
OTU_table_ITS <- read.csv("ITS2_OTU_table_trimmed.csv")

# Remove the Surface Control column (sequenced poorly and excluded from analysis)
OTU_table_ITS <- OTU_table_ITS %>% dplyr::select(-Surface.Control)

# Remove OTUs with no reads across all samples
OTU_table_ITS <- OTU_table_ITS %>%
  dplyr::filter(rowSums(across(where(is.numeric))) > 0)

# Read ITS2 taxonomy file
taxonomy_ITS <- read.csv("ITS2/cleaned_taxonomy_blastn_ITS_final.csv")

# Replace missing or empty taxonomy fields with "Unclassified"
taxonomy_ITS[is.na(taxonomy_ITS)] <- "Unclassified"
taxonomy_ITS[taxonomy_ITS == ""] <- "Unclassified"

# Filter taxonomy to include only OTUs present in the OTU table
taxonomy_ITS <- taxonomy_ITS %>%
  semi_join(OTU_table_ITS, by = "OTU_ID")

# Result: 749 OTUs detected using ITS2 primers after filtering

# ---------- ISOLATE DATA ----------

# Read all cultured isolate OTUs (some match ITS2, some are additional for phylogenetic comparison)
isolates <- read.csv("figures/all_isolates_trimmed.csv")

# 67 fungal OTUs were cultured and sequenced
# An additional 33 OTUs from isolates were included even if not matched in ITS2 amplicons
# Total ITS OTUs considered = 749 (from table) + 33 (extra) = 782

# ---------- 18S DATA ----------

# Read 18S OTU table (animals and land plant OTUs already removed)
OTU_table_18s <- read.csv("18s_OTU_table_trimmed.csv")

# Read 18S taxonomy table
taxonomy_18s <- read.csv("18s/combined_taxonomy_18s.csv") 

# Filter taxonomy to keep only OTUs present in the OTU table
taxonomy_18s <- taxonomy_18s %>%
  dplyr::filter(OTU_ID %in% OTU_table_18s$OTU_ID)

# Result: 342 OTUs detected using 18S primers after filtering


#### OTU Accumulation Curves ####
#### Rarefaction and Species Accumulation Curves ####

# ---------- ITS2 Rarefaction by Sample ----------

# Transpose OTU table so samples are rows and OTUs are columns (as expected by vegan)
otu_mat_ITS <- t(OTU_table_ITS[, -1])  # remove OTU_ID column before transposing

# Generate rarefaction curves for each sample
# Using a fixed step size and sampling depth equal to the minimum library size
rare_list_ITS <- rarecurve(otu_mat_ITS, step = 50, sample = min(rowSums(otu_mat_ITS)), plot = FALSE)
names(rare_list_ITS) <- rownames(otu_mat_ITS)

# Convert rarefaction results to a tidy format for ggplot
rare_df_ITS <- do.call(rbind, lapply(names(rare_list_ITS), function(samp) {
  obs_otus <- rare_list_ITS[[samp]]                   # observed OTUs
  reads <- as.numeric(sub("N", "", names(obs_otus)))  # extract read counts from names
  data.frame(Sample = samp, Reads = reads, Observed_OTUs = obs_otus)
}))

# Define color palette for samples
accumulation_colors <- c(
  Bull = "#332288",
  Conant = "#88CCEE",
  Greg = "#44AA99",
  Horn = "#117733",
  Mortensen = "#DDCC77",
  West = "#CC6677"
)

# Plot rarefaction curves for ITS
acc_ITS_plot <- ggplot(rare_df_ITS, aes(x = Reads, y = Observed_OTUs, color = Sample)) +
  geom_line(size = 1) +
  labs(
    x = "Read Count",
    y = "Observed OTUs (ITS)",
    color = "BNG Well"
  ) +
  scale_color_manual(values = accumulation_colors) +
  theme_minimal() +
  theme(legend.position = "right")

# Results: All six ITS samples clearly show rarefaction saturation

# ---------- 18S Rarefaction by Sample ----------

# Transpose OTU table for 18S
otu_mat_18s <- t(OTU_table_18s[, -1])  # remove OTU_ID

# Remove Surface Control sample (now a row after transposition)
otu_mat_18s <- otu_mat_18s[rownames(otu_mat_18s) != "Surface.Control", ]

# Generate rarefaction curves for 18S
rare_list_18s <- rarecurve(otu_mat_18s, step = 10, sample = min(rowSums(otu_mat_18s)), plot = FALSE)
names(rare_list_18s) <- rownames(otu_mat_18s)

# Convert 18S rarefaction output to tidy format
rare_df_18s <- do.call(rbind, lapply(names(rare_list_18s), function(samp) {
  obs_otus <- rare_list_18s[[samp]]
  reads <- as.numeric(sub("N", "", names(obs_otus)))
  data.frame(Sample = samp, Reads = reads, Observed_OTUs = obs_otus)
}))

# Plot rarefaction curves for 18S
acc_18s_plot <- ggplot(rare_df_18s, aes(x = Reads, y = Observed_OTUs, color = Sample)) +
  geom_line(size = 1) +
  labs(
    x = "Read Count",
    y = "Observed OTUs (18s)",
    color = "BNG Well"
  ) +
  scale_color_manual(values = accumulation_colors) +
  theme_minimal()

# Results: All three 18S samples reach clear rarefaction plateaus

# ---------- ITS2 Species Accumulation (By Sample Count) ----------

# Estimate species richness accumulation curves using 100 random permutations
accum_ITS <- specaccum(otu_mat_ITS, method = "random", permutations = 100)

# Convert to dataframe for plotting
accum_df_ITS <- data.frame(
  Samples = accum_ITS$sites,
  Richness = accum_ITS$richness,
  SD = accum_ITS$sd
)

# Plot species accumulation with confidence shading
acc_overall_ITS <- ggplot(accum_df_ITS, aes(x = Samples, y = Richness)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, fill = "blue") +
  labs(
    x = "Number of BNG Wells",
    y = "Observed OTUs (ITS)"
  ) +
  theme_minimal()

# ---------- 18S Species Accumulation (By Sample Count) ----------

# Run species accumulation analysis for 18S
accum_18s <- specaccum(otu_mat_18s, method = "random", permutations = 100)

# Format output as dataframe
accum_df_18s <- data.frame(
  Samples = accum_18s$sites,
  Richness = accum_18s$richness,
  SD = accum_18s$sd
)

# Plot with shaded error
acc_overall_18s <- ggplot(accum_df_18s, aes(x = Samples, y = Richness)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, fill = "blue") +
  labs(
    x = "Number of BNG Wells",
    y = "Observed OTUs (18s)"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(accum_df_18s$Samples)))) +
  theme_minimal()

# ---------- Combine All Plots into One Figure ----------

# Stack the four plots (18S and ITS, both rarefaction and accumulation)
accplot <- acc_18s_plot / acc_overall_18s / acc_ITS_plot / acc_overall_ITS +
  plot_layout(ncol = 1, heights = c(2, 2, 2, 2)) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )

# Save final combined figure to SVG file
ggsave("accumulation_plots.svg", plot = accplot, width = 8, height = 12, units = "in")



#### Venn Diagram Fungi ITS####
# ----------------- Filter ITS OTUs to Only Fungal Taxa and Generate Venn Diagram ####

# Extract only fungal entries from the ITS taxonomy
Taxonomy_ITS_fungi <- taxonomy_ITS %>%
  dplyr::filter(kingdom == "Fungi")  # Keep only OTUs classified as fungi

# Convert taxonomy to matrix format required by phyloseq::tax_table
# Set OTU_IDs as row names
taxonomy_ITS_mat <- as.matrix(column_to_rownames(Taxonomy_ITS_fungi, var = "OTU_ID"))
tax_table_ITS_phy <- tax_table(taxonomy_ITS_mat)  # Create phyloseq tax_table object

# Filter the ITS OTU table to keep only fungal OTUs based on taxonomy
OTU_table_ITS_fungi <- OTU_table_ITS %>%
  semi_join(Taxonomy_ITS_fungi, by = "OTU_ID")

# Format the OTU table: set OTU IDs as row names and remove the OTU_ID column
rownames(OTU_table_ITS_fungi) <- OTU_table_ITS_fungi[[1]]
OTU_table_ITS_fungi <- OTU_table_ITS_fungi[, -1]

# Convert OTU table to binary presence/absence format (TRUE if OTU > 0 reads)
otu_table_binary <- OTU_table_ITS_fungi > 0

# Generate a list of OTUs present in each site
# The result is a named list where each site has a vector of OTU IDs
otu_list_per_site <- apply(otu_table_binary, 2, function(col) {
  rownames(otu_table_binary)[col]
})

# Define custom colors for each of the 6 sites for plotting
venn_site_colors <- c("red", "blue", "green", "orange", "purple", "brown")
names(venn_site_colors) <- names(otu_list_per_site)[1:6]  # Assign site names to colors

# Plot Venn diagram showing OTU overlaps across sites
venn_plot <- ggVennDiagram(
  otu_list_per_site[1:6],       # use first 6 sites
  label = "count",              # show OTU counts in each intersection
  edge_size = 1.2,
  edge_lty = "solid",
  set_color = venn_site_colors, # custom site colors
  set_name_size = 7             # increase font size of site names
) +
  theme_void() +                # minimal theme for cleaner look
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  labs(fill = "OTU Count")      # customize legend title

# Display plot
venn_plot

# Save Venn diagram to SVG file
ggsave("venn_plot.svg",
       plot = venn_plot,
       width = 10, height = 9, units = "in", device = "svg", limitsize = FALSE)

# ----------------- Identify Shared OTUs -----------------

# Count the number of sites each OTU appears in (row sums of binary table)
venn_otu_site_counts <- rowSums(otu_table_binary)

# Extract OTUs present in all 6 sites
otus_in_6_sites <- names(otu_site_counts[venn_otu_site_counts == 6])

# Extract OTUs present in exactly 5 sites
otus_in_5_sites <- names(otu_site_counts[venn_otu_site_counts == 5])

# OTUs present in all sites:
# "0fac242e0a12de35d46a2b2ac8dde906" — matches Irpex lacteus (100%), also isolated

# OTUs present in 5 of 6 sites:
# "827c3c02a67a4c50e5e10ec4100429f4" — 100% match to Saccharomyces sp. (uncultured)
# "081720e5f9f84e1181553a050cbae7f2" — 92% match to Helotiales sp. (isolated)




#### Stacked Bar Plots for ITS Fungal Community Composition ####

# ---- Load Sample Metadata for ITS Samples ----
metadata_ITS <- read.csv(file = "metadata_ITS.csv")  # Import sample metadata
metadata_phy <- sample_data(column_to_rownames(metadata_ITS, var = "Sample"))  # Set sample IDs as rownames

# ---- Filter ITS Taxonomy to Only Include Fungi ----
# Taxonomy file is assumed already loaded as `taxonomy_ITS`
Taxonomy_ITS_fungi <- taxonomy_ITS %>%
  dplyr::filter(kingdom == "Fungi")  # Keep only fungal OTUs

# Prepare tax_table object from fungal taxonomy
taxonomy_ITS_mat <- as.matrix(column_to_rownames(Taxonomy_ITS_fungi, var = "OTU_ID"))
tax_table_ITS_phy <- tax_table(taxonomy_ITS_mat)

# Filter OTU table to match only fungal OTUs
OTU_table_ITS_fungi <- OTU_table_ITS %>%
  semi_join(Taxonomy_ITS_fungi, by = "OTU_ID")

# Convert OTU table to matrix and make OTU IDs rownames
otu_ITS_mat <- as.matrix(column_to_rownames(OTU_table_ITS_fungi, var = "OTU_ID"))
otu_ITS_phy <- otu_table(otu_ITS_mat, taxa_are_rows = TRUE)  # Indicate OTUs are rows

# ---- Construct Phyloseq Object for ITS Fungi ----
physeq_ITS <- phyloseq(otu_ITS_phy, tax_table_ITS_phy, metadata_phy)  # Combine OTU, taxonomy, and metadata
# Result: 750 OTUs in total, 689 are classified as fungi

# ---- Transform Counts to Relative Abundances ----
physeq_rel_ITS <- transform_sample_counts(physeq_ITS, function(x) x / sum(x))  # Normalize within each sample

# ---- Sample Order for Plotting ----
desired_order <- c("Conant", "Greg", "Mortensen", "Bull", "Horn", "West")  # Desired sample order

# Update sample names and factor levels to enforce order
sample_data(physeq_rel_ITS)$Sample <- rownames(sample_data(physeq_rel_ITS))
sample_data(physeq_rel_ITS)$Sample <- factor(
  sample_data(physeq_rel_ITS)$Sample,
  levels = desired_order
)

# ---- Plot: Relative Abundance at the Phylum Level ----
ITS_phylum <- plot_bar(physeq_rel_ITS, x = "sample_Sample", fill = "phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_brewer(palette = "Set2") +  # Colorblind-friendly palette
  theme_minimal() +
  ylab("Relative Abundance") +
  labs(fill = "Phylum") +
  theme(
    axis.text.x = element_blank(),         # Remove x-axis tick labels
    axis.ticks.x = element_blank(),        # Remove x-axis ticks
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank(),    # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

# Save if needed:
# ggsave("ITS_phylum_stacked_barplot.jpg", plot = ITS_phylum, width = 16, height = 12, dpi = 1000)

# ---- Plot: Relative Abundance at the Class Level ----

# Agglomerate OTUs to the Class level
physeq_ITs_class <- tax_glom(physeq_rel_ITS, taxrank = "class")

# Extract class names and build a custom color palette
class_levels <- unique(tax_table(physeq_ITs_class)[, "class"])
class_levels <- as.character(class_levels[!is.na(class_levels)])  # Remove NA entries

# Generate a named color palette with enough distinct colors for all classes
palette_24_named <- setNames(createPalette(length(class_levels), c("#ffffff", "#000000")),
                             class_levels)

# Build the plot
ITS_class <- plot_bar(physeq_ITs_class, x = "sample_Sample", fill = "class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_24_named) +  # Use custom class color palette
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x labels
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

# Save if needed:
# ggsave("ITS_class_stacked_barplot.jpg", plot = ITS_class, width = 16, height = 12, dpi = 1000)

#### Plot Isolates Stacked Bar ####
# Load isolate taxonomy data and subset relevant taxonomic columns
isolates_data <- read.csv("all_isolates_trimmed.csv", row.names = NULL)
isolates_taxonomy <- isolates_data[, c(1,10:15)]  # OTU_ID and taxonomic ranks

# Set OTU_ID as row names and convert to matrix
taxonomy_mat <- isolates_taxonomy %>%
  column_to_rownames(var = "OTU_ID") %>%
  as.matrix()

# Convert matrix to tax_table and create phyloseq object
tax_table_obj <- tax_table(taxonomy_mat)
physeq_from_taxonomy <- phyloseq(tax_table_obj)

# Extract taxonomy table as data frame
tax_df_isolates <- as.data.frame(tax_table(physeq_from_taxonomy)) %>%
  rownames_to_column("OTU_ID")


# ----- Phylum-Level Stacked Bar Plot -----
# Count number of OTUs per phylum
rank_summaries_phylum <- 
  tax_df_isolates %>%
  group_by(Phylum) %>%
  summarise(OTU_count = n()) %>%
  ungroup() %>%
  arrange(desc(OTU_count))

# Plot phylum distribution as stacked bar
Isolates_phylum <- ggplot(rank_summaries_phylum, aes(x = "", y = OTU_count, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_phylum) +
  theme_minimal() +
  ylab("") +
  xlab("") +
  labs(fill = "Phylum") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
  )


# ----- Class-Level Stacked Bar Plot -----
# Count number of OTUs per class
rank_summaries_class <- 
  tax_df_isolates %>%
  group_by(Class) %>%
  summarise(OTU_count = n()) %>%
  ungroup() %>%
  arrange(desc(OTU_count))

# Plot class distribution
Isolates_class <- ggplot(rank_summaries_class, aes(x = "", y = OTU_count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_class) +
  theme_minimal() +
  ylab("") +
  xlab("") +
  labs(fill = "Class") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

Isolates_class


# ----- Order-Level Stacked Bar Plot -----
# Count OTUs by order and generate custom color palette
rank_summaries_order <- 
  tax_df_isolates %>%
  group_by(Order) %>%
  summarise(OTU_count = n()) %>%
  ungroup() %>%
  arrange(desc(OTU_count))

rank_summaries_order$Order <- as.character(rank_summaries_order$Order)
unique_orders <- setdiff(unique(rank_summaries_order$Order), "Unclassified")
n_orders <- length(unique_orders)

# Generate color palette for orders, with "Unclassified" as black
full_palette_alphabet <- alphabet()
clean_palette_order <- full_palette_alphabet[full_palette_alphabet != "#191919"]
palette_colors_order <- clean_palette_order[1:n_orders]
palette_orders <- setNames(palette_colors_order, unique_orders)
palette_orders["Unclassified"] <- "black"

rank_summaries_order$Order <- factor(rank_summaries_order$Order, levels = c(unique_orders, "Unclassified"))

# Plot order distribution
Isolates_order <- ggplot(rank_summaries_order, aes(x = "", y = OTU_count, fill = Order)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_orders, guide = guide_legend(ncol = 2)) +
  theme_minimal() +
  ylab("Isolate OTU Count") +
  xlab("") +
  labs(fill = "Order") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.y = element_text(size = 25),
    axis.text.y = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

Isolates_order


# ----- Family-Level Stacked Bar Plot -----
# Count OTUs by family and assign custom colors
rank_summaries_family <- 
  tax_df_isolates %>%
  group_by(Family) %>%
  summarise(OTU_count = n()) %>%
  ungroup() %>%
  arrange(desc(OTU_count))

unique_families <- setdiff(unique(rank_summaries_family$Family), "Unclassified")
n_families <- length(unique_families)

# Use Polychrome palette, ensuring enough unique colors
full_palette_polychrome <- polychrome()
clean_palette_family <- full_palette_polychrome[full_palette_polychrome != "#191919"]
if (n_families > length(clean_palette_family)) {
  stop("Not enough distinct colors after removing #191919 from Polychrome palette.")
}

palette_colors_families <- clean_palette_family[1:n_families]
palette_family <- setNames(palette_colors_families, unique_families)
palette_family["Unclassified"] <- "black"

rank_summaries_family$Family <- factor(rank_summaries_family$Family, levels = c(unique_families, "Unclassified"))

# Plot family distribution
Isolates_family <- ggplot(rank_summaries_family, aes(x = "", y = OTU_count, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_family, guide = guide_legend(ncol = 3)) +
  theme_minimal() +
  ylab("") +
  xlab("") +
  labs(fill = "Family") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )


# ----- Genus-Level Stacked Bar Plot -----
# Count OTUs per genus and assign colors
rank_summaries_genus <- 
  tax_df_isolates %>%
  group_by(Genus) %>%
  summarise(OTU_count = n()) %>%
  ungroup() %>%
  arrange(desc(OTU_count))

unique_genera <- setdiff(unique(rank_summaries_genus$Genus), "Unclassified")
n_genera <- length(unique_genera)

full_palette_polychrome <- polychrome()
clean_palette_genus <- full_palette_polychrome[full_palette_polychrome != "#191919"]

palette_colors_genus <- clean_palette_genus[1:n_genera]
palette_genus <- setNames(palette_colors_genus, unique_genera)
palette_genus["Unclassified"] <- "black"

rank_summaries_genus$Genus <- factor(rank_summaries_genus$Genus, levels = c(unique_genera, "Unclassified"))

# Plot genus distribution
Isolates_genus <- ggplot(rank_summaries_genus, aes(x = "", y = OTU_count, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_genus, guide = guide_legend(ncol = 3)) +
  theme_minimal() +
  ylab("") +
  xlab("") +
  labs(fill = "Genus") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

Isolates_genus


# ----- Combine All Taxonomic Plots -----
# Assemble all stacked bar plots vertically using patchwork
combined_stacked_bar_plot_isolates <- Isolates_phylum / Isolates_class / Isolates_order / Isolates_family / Isolates_genus +
  plot_layout(ncol = 1, heights = rep(1, 5)) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold", family = "Source Sans Pro")  # Optional plot-level styling
    )
  )

combined_stacked_bar_plot_isolates

# Save full composite figure to PDF
ggsave("combined_stacked_bar_plot_isolates.pdf",
       plot = combined_stacked_bar_plot_isolates,
       width = 15, height = 40, units = "in", device = "pdf", limitsize = FALSE)




#### Stacked Bar Plots: 18S rRNA Data ####

# ----- Prepare Metadata for 18S Analysis ####
# Use ITS metadata for 18S samples (assumed compatible structure)
metadata_18s <- metadata_ITS

# Convert metadata to phyloseq sample_data format
metadata_phy_18s <- sample_data(column_to_rownames(metadata_18s, var = "Sample"))

# Remove surface control column from OTU table
OTU_table_18s <- OTU_table_18s %>% dplyr::select(-Surface.Control)


# ----- Filter Fungal Taxa for Stacked Bar Plot ####
# Subset taxonomy to only fungal OTUs (based on Subdivision.Kingdom == "Fungi")
taxonomy_18s_fungi <- taxonomy_18s %>%
  dplyr::filter(Subdivision.Kingdom == "Fungi")
# Result: 95 fungal OTUs detected with 18S primers

# Convert taxonomy to matrix and set OTU IDs as rownames
taxonomy_18s_fungi_mat <- as.matrix(column_to_rownames(taxonomy_18s_fungi , var = "OTU_ID"))
tax_table_18s_fungi_phy <- tax_table(taxonomy_18s_fungi_mat)

# Filter OTU table to only fungal OTUs
OTU_table_18s_fungi <- OTU_table_18s %>%
  semi_join(taxonomy_18s_fungi, by = "OTU_ID")

# Convert filtered OTU table to matrix and create phyloseq OTU object
otu_18s_fungi_mat <- as.matrix(column_to_rownames(OTU_table_18s_fungi, var = "OTU_ID"))
otu_18s_fungi_phy <- otu_table(otu_18s_fungi_mat, taxa_are_rows = TRUE)

# Create phyloseq object for fungal data
physeq_18s_fungi <- phyloseq(otu_18s_fungi_phy, tax_table_18s_fungi_phy, metadata_phy_18s)

# Transform to relative abundance
physeq_rel_18s_fungi <- transform_sample_counts(physeq_18s_fungi, function(x) x / sum(x))

# Ensure sample factor levels follow desired order
sample_data(physeq_rel_18s_fungi)$Sample <- rownames(sample_data(physeq_rel_18s_fungi))
sample_data(physeq_rel_18s_fungi)$Sample <- factor(
  sample_data(physeq_rel_18s_fungi)$Sample,
  levels = desired_order
)


# ----- Plot Fungal Phylum-Level Bar Plot ####
phylum_18s <- plot_bar(physeq_rel_18s_fungi, x = "sample_Sample", fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Phylum") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

#ggsave("18s_fungi_phylum_stacked_barplot.jpg", plot = phylum_18s, width = 16, height = 12, dpi = 1000)


# ----- Plot Fungal Class-Level Bar Plot ####
# Agglomerate OTUs to Class level
physeq_18s_class <- tax_glom(physeq_rel_18s_fungi, taxrank = "Class")

# Generate unique class names and define palette
class_levels_18s <- unique(tax_table(physeq_18s_class)[, "Class"])
class_levels_18s <- as.character(class_levels_18s[!is.na(class_levels_18s)])
palette_18s_fungi <- setNames(createPalette(length(class_levels_18s), c("#ffffff", "#000000")),
                              class_levels_18s)

# Plot fungal class-level composition
class_18s <- plot_bar(physeq_18s_class, x = "sample_Sample", fill = "Class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_18s_fungi) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

#ggsave("18s_fungi_class_stacked_barplot.jpg", plot = class_18s, width = 16, height = 12, dpi = 1000)


# ----- Filter Microeukaryotes (Non-Fungal Eukaryotes) ####
taxonomy_18s_microeuk <- taxonomy_18s %>%
  dplyr::filter(Subdivision.Kingdom != "Fungi")

# Convert to matrix and tax_table
taxonomy_18s_euk_mat <- as.matrix(column_to_rownames(taxonomy_18s_microeuk , var = "OTU_ID"))
tax_table_18s_euk_phy <- tax_table(taxonomy_18s_euk_mat)

# Filter OTU table to microeukaryotes only
OTU_table_18s_microeuk <- OTU_table_18s %>%
  semi_join(taxonomy_18s_microeuk, by = "OTU_ID")

# Convert OTU table to matrix and build OTU object
otu_18s_euk_mat <- as.matrix(column_to_rownames(OTU_table_18s_microeuk, var = "OTU_ID"))
otu_18s_euk_phy <- otu_table(otu_18s_euk_mat, taxa_are_rows = TRUE)

# Create phyloseq object for microeukaryotes
physeq_18s_euk <- phyloseq(otu_18s_euk_phy, tax_table_18s_euk_phy, metadata_phy_18s)

# Transform to relative abundance
physeq_rel_18s_euk <- transform_sample_counts(physeq_18s_euk , function(x) x / sum(x))


# ----- Plot Microeukaryotic Division-Level Bar Plot ####
# Ensure sample order
sample_data(physeq_rel_18s_euk)$Sample <- rownames(sample_data(physeq_rel_18s_euk))
sample_data(physeq_rel_18s_euk)$Sample <- factor(
  sample_data(physeq_rel_18s_euk)$Sample,
  levels = desired_order
)

# Plot Division-level bar plot
euk_Division <- plot_bar(physeq_rel_18s_euk, x = "sample_Sample", fill = "Division") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Division") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

#ggsave("18s_microeuk_division_stacked_barplot.jpg", plot = euk_Division, width = 16, height = 12, dpi = 1000)


# ----- Plot Microeukaryotic Class-Level Bar Plot ####
# Agglomerate to Class level
physeq_18s_euk_class <- tax_glom(physeq_rel_18s_euk, taxrank = "Class")

# Get unique classes and assign colors
class_levels_18s_euk <- unique(tax_table(physeq_18s_euk_class)[, "Class"])
class_levels_18s_euk <- as.character(class_levels_18s_euk[!is.na(class_levels_18s_euk)])

palette_18s_euk <- setNames(createPalette(length(class_levels_18s_euk), c("#ffffff", "#000000")),
                            class_levels_18s_euk)

# Plot microeukaryotic class composition
euk_class <- plot_bar(physeq_18s_euk_class, x = "sample_Sample", fill = "Class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_18s_euk) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

#ggsave("18s_microeuk_class_stacked_barplot.jpg", plot = euk_class, width = 16, height = 12, dpi = 1000)



#### Merge All Stacked Bar Plots ####

# ----- Define consistent color palettes for taxa -----
# Extract unique phyla from ITS, 18s eukaryotes, and 18s fungi; remove NA
all_phyla <- unique(c(
  as.character(tax_table(physeq_rel_ITS)[, "phylum"]),
  as.character(tax_table(physeq_rel_18s_euk)[, "Phylum"]),
  as.character(tax_table(physeq_rel_18s_fungi)[, "Phylum"])
))
all_phyla <- all_phyla[!is.na(all_phyla)]
# Assign colors from Set2 palette for phyla
palette_phylum <- setNames(
  brewer.pal(n = 6, name = "Set2"),
  all_phyla
)

# Extract unique class names from ITS, 18s eukaryotes, and fungi; remove NA
all_class <- unique(c(
  as.character(tax_table(physeq_rel_ITS)[, "class"]),
  as.character(tax_table(physeq_rel_18s_euk)[, "Class"]),
  as.character(tax_table(physeq_rel_18s_fungi)[, "Class"])
))
all_class <- all_class[!is.na(all_class)]
# Create class palette seeded with black
palette_class <- setNames(
  createPalette(length(all_class), seedcolors = c("#000000")),
  all_class
)
# Set "Unclassified" to black if present
if ("Unclassified" %in% all_class) {
  palette_class["Unclassified"] <- "#000000"
}

# Extract unique divisions from 18s eukaryotes; remove NA
all_division <- unique(c(
  as.character(tax_table(physeq_rel_18s_euk)[, "Division"])
))
all_division <- all_division[!is.na(all_division)]
# Create division palette seeded with black
palette_division <- setNames(
  createPalette(length(all_division), seedcolors = c("#000000")),
  all_division
)
# Set "Unclassified" to black if present
if ("Unclassified" %in% all_division) {
  palette_division["Unclassified"] <- "#000000"
}


# ----- Generate stacked bar plots for ITS datasets -----
# ITS phylum plot with no x-axis labels for compactness
plot_ITS_phyla <- plot_bar(physeq_rel_ITS, x = "sample_Sample", fill = "phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  scale_fill_manual(values = palette_phylum) +
  theme_minimal() +
  ylab("Relative Abundance") +
  labs(fill = "Phylum") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ITS class plot with rotated x-axis labels and full axis details
plot_ITS_class <- plot_bar(physeq_ITs_class, x = "sample_Sample", fill = "class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_class) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )


# ----- Generate stacked bar plots for 18s fungi datasets -----
# Fungi phylum plot without x-axis labels
plot_fungi_18s_phyla <- plot_bar(physeq_rel_18s_fungi, x = "sample_Sample", fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_phylum) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Phylum") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Fungi class plot without x-axis labels
plot_fungi_18s_class <- plot_bar(physeq_18s_class, x = "sample_Sample", fill = "Class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_class) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )


# ----- Generate stacked bar plots for 18s eukaryote datasets -----
# Division plot without x-axis labels
plot_euk_18s_division <- plot_bar(physeq_rel_18s_euk, x = "sample_Sample", fill = "Division") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_division) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Division") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Eukaryote class plot with rotated x-axis labels
plot_euk_18s_class <- plot_bar(physeq_18s_euk_class, x = "sample_Sample", fill = "Class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = palette_class) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  xlab("BNG Well") +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )


# ----- Combine and save stacked bar plots -----
combined_stacked_bar_plot_ITS <- plot_ITS_phyla / plot_ITS_class +
  plot_layout(ncol = 1, heights = c(1, 1)) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(", 
    tag_suffix = ")",
    theme = theme(
      plot.tag.text = element_text(size = 16, face = "bold")
    )
  )

ggsave("combined_stacked_bar_plot_ITS.pdf",
       plot = combined_stacked_bar_plot_ITS,
       width = 15, height = 15, units = "in", device = "pdf", limitsize = FALSE)

combined_stacked_bar_plot_18s <- plot_fungi_18s_phyla / plot_fungi_18s_class / plot_euk_18s_division / plot_euk_18s_class  +
  plot_layout(ncol = 1, heights = c(1,1,1,1))

ggsave("combined_stacked_bar_plot_18s.pdf",
       plot = combined_stacked_bar_plot_18s,
       width = 15, height = 26, units = "in", device = "pdf", limitsize = FALSE)






#### Calculate Chao Richness ####

# Filter taxonomy to include only fungi from ITS taxonomy data
Taxonomy_ITS_fungi <- taxonomy_ITS %>%
  dplyr::filter(kingdom == "Fungi")

# Subset ITS OTU table to include only fungal OTUs using semi_join
OTU_table_ITS_fungi <- OTU_table_ITS %>%
  semi_join(Taxonomy_ITS_fungi, by = "OTU_ID")

# Remove OTU_ID column to keep only abundance data for diversity calculation
OTU_table_ITS_fungi_diversity <- OTU_table_ITS_fungi %>%
  dplyr::select(-OTU_ID)

# Transpose the OTU table so rows are samples and columns are OTUs (required format)
OTU_table_ITS_fungi_diversity <- t(OTU_table_ITS_fungi_diversity)

# Calculate Chao1 richness estimates for each sample
chao1_values_ITS <- estimateR(OTU_table_ITS_fungi_diversity)["S.chao1", ]

# Convert Chao1 vector to data frame with sample names
chao1_df <- data.frame(Sample = rownames(OTU_table_ITS_fungi_diversity), Chao1 = chao1_values_ITS)

# Merge the Chao1 richness estimates with sample metadata for downstream analysis
metadata_ITS <- merge(chao1_df, metadata_ITS, by = "Sample")

#### Calculate Shannon Diversity ####

# Display the OTU table used for diversity calculations (fungal ITS only)
OTU_table_ITS_fungi_diversity

# Determine minimum sequencing depth across all samples for rarefaction
min_depth <- min(rowSums(OTU_table_ITS_fungi_diversity))

# Rarefy OTU counts to the minimum sequencing depth to standardize sample sizes
shannon_otu_rarefied <- rrarefy(OTU_table_ITS_fungi_diversity, sample = min_depth)

# Calculate Shannon diversity index on the rarefied OTU table
shannon <- diversity(shannon_otu_rarefied, index = "shannon")

# Append Shannon diversity values to the metadata by matching sample names
metadata_ITS$Shannon <- shannon[match(metadata_ITS$Sample, names(shannon))]

#### Calculate Evenness ####

# Calculate observed richness (number of OTUs per sample)
richness_ITS <- specnumber(OTU_table_ITS_fungi_diversity)

# Compute Pielou's Evenness as Shannon diversity divided by log of richness
evenness_ITS <- shannon / log(richness_ITS)

# Add Evenness values to the metadata by matching sample names
metadata_ITS$Evenness <- evenness_ITS[match(metadata_ITS$Sample, rownames(OTU_table_ITS_fungi_diversity))]

#### Phylogenetic Dispersion ####

# Load OTU table and taxonomy data
pd_otu <- read.csv("ITS2_OTU_table_trimmed.csv", row.names = 1)
pd_tax <- read.csv("ITS2/cleaned_taxonomy_blastn_ITS_final.csv", row.names = 1)

# Filter taxonomy to keep only fungi OTUs
pd_tax <- pd_tax %>% dplyr::filter(kingdom == "Fungi")

# Subset OTU table to include only OTUs present in filtered taxonomy
pd_otu <- pd_otu[rownames(pd_otu) %in% rownames(pd_tax), ]
pd_tax <- pd_tax[rownames(pd_tax) %in% rownames(pd_otu), ]

# Replace NA or empty taxonomy cells with "Unclassified"
pd_tax[is.na(pd_tax)] <- "Unclassified"
pd_tax[pd_tax == ""] <- "Unclassified"

# Export taxonomy for tree construction (Perl script input)
write.table(pd_tax,
            file = "taxonomy_for_tree.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = FALSE)

# Read phylogenetic tree generated externally (newick format)
ITS_tree <- read.tree("tree.nwk")

# Transpose OTU table to samples as rows for phylogenetic analysis
pd_otu_t <- t(pd_otu)

# Prune tree and OTU table to ensure they contain the same taxa
pruned_result <- prune.sample(pd_otu_t, ITS_tree)
pruned_tree <- pruned_result
taxa_in_tree <- pruned_tree$tip.label

# Subset OTU table to match pruned tree taxa
pruned_otu <- pd_otu_t[, taxa_in_tree, drop = FALSE]

# Fix branch lengths to uniform value (if necessary)
pruned_tree$edge.length <- rep(60, length(pruned_tree$edge.length))

# Calculate phylogenetic diversity (PD) per sample
pd_results <- pd(pruned_otu, pruned_tree) %>%
  dplyr::mutate(Sample = rownames(.)) %>%
  dplyr::select(Sample, PD)

# Merge PD results with metadata
metadata_ITS <- merge(pd_results, metadata_ITS, by = "Sample")

# Calculate standardized effect size of PD (ses.pd) using taxa.labels null model
sespd_results <- ses.pd(pd_otu_t, ITS_tree, null.model = "taxa.labels", runs = 999)

# Format ses.pd results for downstream use
sespd_df <- data.frame(
  Sample = rownames(pd_otu_t),
  SES_PD_Z = sespd_results$pd.obs.z,
  SES_PD_P = sespd_results$pd.obs.p,
  fungal_OTUs = sespd_results$ntaxa
)

# Calculate phylogenetic distance matrix from tree
phydist <- cophenetic(ITS_tree)

# Calculate standardized effect size of mean pairwise distance (ses.mpd)
ses.mpd.result <- ses.mpd(
  pd_otu_t,
  phydist,
  null.model = "taxa.labels",
  abundance.weighted = FALSE,
  runs = 999
)
ses.mpd_df <- data.frame(
  Sample = rownames(pd_otu_t),
  SES_MPD_Z = ses.mpd.result$mpd.obs.z,
  SES_MPD_P = ses.mpd.result$mpd.obs.p
)

# Calculate standardized effect size of mean nearest taxon distance (ses.mntd)
ses.mntd.result <- ses.mntd(
  pd_otu_t, 
  phydist, 
  null.model = "taxa.labels",
  abundance.weighted = FALSE, 
  runs = 999
)
ses.mntd_df <- data.frame(
  Sample = rownames(pd_otu_t),
  SES_MNTD_Z = ses.mntd.result$mntd.obs.z,
  SES_MNTD_P = ses.mntd.result$mntd.obs.p
)

# Combine all phylogenetic dispersion metrics by sample
combined_pd_results <- pd_results %>%
  left_join(sespd_df, by = "Sample") %>%
  left_join(ses.mpd_df, by = "Sample") %>%
  left_join(ses.mntd_df, by = "Sample")

# Merge combined phylogenetic metrics with metadata
metadata_ITS <- merge(combined_pd_results, metadata_ITS, by = "Sample")

#### Genetic Dissimilarity ####

# Load percent identity data generated from BLASTn results (final_blastn_taxonmy_ITS.R)
ITS_per_ident <- read.csv("ITS2/highest_%identity_per_otu.csv")

# Filter percent identity data to keep only fungal OTUs present in taxonomy
ITS_per_ident <- ITS_per_ident %>%
  semi_join(Taxonomy_ITS_fungi, by = "OTU_ID")

# Merge filtered OTU abundance table with percent identity data by OTU_ID
ITS_per_ident <- merge(OTU_table_ITS_fungi, ITS_per_ident, by = "OTU_ID")

# Identify sample columns in merged table (exclude OTU_ID and percent_identity)
site_columns <- setdiff(names(ITS_per_ident), c("OTU_ID", "percent_identity"))

# Calculate unweighted average percent identity per sample (only OTUs present)
avg_pid_unweighted <- sapply(site_columns, function(site) {
  present <- ITS_per_ident[[site]] > 0                     # Identify OTUs present in this sample
  mean(ITS_per_ident$percent_identity[present], na.rm = TRUE)  # Mean percent identity of present OTUs
})

# Convert unweighted averages to data frame with sample names
avg_pid_unweighted_df <- data.frame(
  site = names(avg_pid_unweighted),
  unweighted_avg_percent_identity = as.numeric(avg_pid_unweighted)
)

# Rename column for consistency and merge into metadata
names(avg_pid_unweighted_df)[names(avg_pid_unweighted_df) == "site"] <- "Sample"
metadata_ITS <- merge(metadata_ITS, avg_pid_unweighted_df, by = "Sample", all.x = TRUE)

# Calculate abundance-weighted average percent identity per sample
avg_pid_weighted <- sapply(site_columns, function(site) {
  abundances <- ITS_per_ident[[site]]
  total_abundance <- sum(abundances, na.rm = TRUE)
  
  # Return NA if no abundance; else weighted mean of percent identity
  if (total_abundance == 0) {
    NA
  } else {
    weighted.mean(ITS_per_ident$percent_identity, w = abundances, na.rm = TRUE)
  }
})

# Convert weighted averages to data frame and rename column
avg_pid_weighted_df <- data.frame(
  site = names(avg_pid_weighted),
  weighted_avg_percent_identity = as.numeric(avg_pid_weighted)
)

names(avg_pid_weighted_df)[names(avg_pid_weighted_df) == "site"] <- "Sample"

# Merge weighted averages into metadata
metadata_ITS <- merge(metadata_ITS, avg_pid_weighted_df, by = "Sample", all.x = TRUE)

# Convert similarity metrics to dissimilarity by subtracting from 100
metadata_ITS$unweighted_avg_percent_identity <- 100 - metadata_ITS$unweighted_avg_percent_identity
metadata_ITS$weighted_avg_percent_identity <- 100 - metadata_ITS$weighted_avg_percent_identity



#### Range per OTU ####
# Calculate maximum distance (km) per OTU based on sites where OTU is present

# Make sure metadata is keyed by site/sample
# Move OTU_ID column to rownames if present
if ("OTU_ID" %in% colnames(OTU_table_ITS_fungi)) {
  rownames(OTU_table_ITS_fungi) <- OTU_table_ITS_fungi$OTU_ID
  OTU_table_ITS_fungi$OTU_ID <- NULL
}

# Extract coordinates (latitude, longitude) for each sample
coords <- metadata_ITS[, c("Sample", "latitude", "longitude")]
rownames(coords) <- coords$Sample
coords$Sample <- NULL

# Function to calculate max range for a given set of sites
calc_max_range <- function(sites_present) {
  if(length(sites_present) < 2) {
    return(0)  # range is zero if only one site
  }
  site_coords <- coords[sites_present, ]
  dist_matrix <- distm(site_coords[, c("longitude", "latitude")], fun = distHaversine) / 1000  # meters to km
  return(max(dist_matrix))
}

# For each OTU, identify sites where it is present
otu_sites_list <- apply(OTU_table_ITS_fungi, 1, function(abund) {
  sites <- colnames(OTU_table_ITS_fungi)[abund > 0]
  return(sites)
})

# Calculate maximum geographic range per OTU
otu_ranges <- sapply(otu_sites_list, calc_max_range)

# Initialize vectors to store unweighted and weighted range per site
unweighted_range_per_site <- numeric(ncol(OTU_table_ITS_fungi))
weighted_range_per_site <- numeric(ncol(OTU_table_ITS_fungi))

sites <- colnames(OTU_table_ITS_fungi)

# Calculate unweighted and abundance-weighted average range per site
for (i in seq_along(sites)) {
  site <- sites[i]
  
  # Abundance of all OTUs at this site
  abundances <- OTU_table_ITS_fungi[, site]
  
  # Identify OTUs present at this site
  present_otus <- abundances > 0
  
  if (sum(present_otus) == 0) {
    unweighted_range_per_site[i] <- NA
    weighted_range_per_site[i] <- NA
    next
  }
  
  # Unweighted mean of OTU ranges present at the site
  unweighted_range_per_site[i] <- mean(otu_ranges[present_otus], na.rm = TRUE)
  
  # Weighted mean of OTU ranges weighted by their abundances at the site
  weighted_range_per_site[i] <- weighted.mean(otu_ranges[present_otus], w = abundances[present_otus], na.rm = TRUE)
}

# Create data frame with unweighted and weighted range per site
range_df <- data.frame(
  site = sites,
  unweighted_range_km = unweighted_range_per_site,
  weighted_range_km = weighted_range_per_site
)

# Rename for consistency before merging with metadata
names(range_df)[names(range_df) == "site"] <- "Sample"

# Merge range data into metadata
metadata_ITS <- merge(metadata_ITS, range_df, by = "Sample", all.x = TRUE)

#### Plot Richness, Evenness, Shannon, Genetic Similarity, Range, PD, MPD vs Depth ####

# Panel 1: Richness vs Depth
richness_plot <- ggplot(metadata_ITS, aes(x = depth, y = Chao1)) +
  geom_point(color = "steelblue", size = 3) +                      # Scatter plot points
  geom_smooth(method = "lm", se = FALSE, color = "black",          # Linear regression trend line without confidence interval
              linetype = "dashed") +
  labs(y = "Richness") +                                            # Y-axis label
  theme_minimal() +                                                 # Minimal theme for clean look
  theme(
    panel.grid.major = element_blank(),                            # Remove major grid lines
    panel.grid.minor = element_blank(),                            # Remove minor grid lines
    axis.line = element_line(color = "black", size = 0.5),         # Add axis lines in black
    axis.title.x = element_blank()                                 # Remove x-axis title
  )

# Panel 1: Evenness vs Depth
evenness_plot <- ggplot(metadata_ITS, aes(x = depth, y = Evenness)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "Evenness") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Panel 1: Shannon Diversity vs Depth
beta_plot <- ggplot(metadata_ITS, aes(x = depth, y = shannon)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(x = "Sample Depth (m)", y = "Shannon Diversity") +          # X and Y axis labels
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Panel 2: Phylogenetic Dispersion (SES_MNTD_Z) vs Depth
pd_plot <- ggplot(metadata_ITS, aes(x = depth, y = SES_MNTD_Z)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "Phylogenetic Dispersion") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Panel 2: Genetic Dissimilarity (weighted average percent identity) vs Depth
gd_plot <- ggplot(metadata_ITS, aes(x = depth, y = weighted_avg_percent_identity)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "Dissimilarity to Surface (%)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Panel 2: Range (weighted range in km) vs Depth
range_plot <- ggplot(metadata_ITS, aes(x = depth, y = weighted_range_km)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(x = "Sample Depth (m)", y = "Range (km)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Combine all six plots into a 2-column by 3-row layout with tags (A)-(F)
combined_ITS_plot <- (richness_plot | pd_plot) /
  (evenness_plot | gd_plot) /
  (beta_plot | range_plot) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold")  # Customize tag font size and style
    )
  )

# Print the combined plot to the output device
combined_ITS_plot

# Save the combined plot as an SVG file with specified dimensions
ggsave("combined_ITS_plot.svg",
       plot = combined_ITS_plot,
       width = 15, height = 15, units = "in", device = "svg", limitsize = FALSE)

#### Metacoder ITS ####
# Remove OTUs from the taxonomy table that belong to the phylum "Anthophyta"
taxonomy_ITS_meta <- taxonomy_ITS %>%
  filter(phylum != "Anthophyta")

# Keep only OTUs in the OTU table that are still present in the filtered taxonomy table
OTU_table_ITS_meta <- OTU_table_ITS %>%
  semi_join(taxonomy_ITS_meta, by = "OTU_ID")

# Add a "Cultured" column to the OTU table:
# OTUs found in the isolate dataset are marked as 1 (cultured), others as 0
OTU_table_ITS_meta$Cultured <- ifelse(OTU_table_ITS_meta$OTU_ID %in% isolates_data$OTU_ID, 1, 0)

# Identify OTUs that are in the isolates but missing from the OTU table
missing_otus_ITS <- setdiff(isolates_data$OTU_ID, OTU_table_ITS_meta$OTU_ID)

# Create a new data frame for those missing OTUs with zero counts and Cultured = 1
new_rows_ITS <- data.frame(
  OTU_ID = missing_otus_ITS,
  matrix(0, nrow = length(missing_otus_ITS), ncol = 6),  # assumes 6 count/sample columns
  Cultured = 1
)

# Assign column names to match the OTU table so the rows can be combined
colnames(new_rows_ITS) <- colnames(OTU_table_ITS_meta)

# Append the new cultured-only rows to the main OTU table
OTU_table_ITS_meta <- rbind(OTU_table_ITS_meta, new_rows_ITS)

# Identify OTUs in isolates that are missing from the taxonomy table
missing_otus_ITS_tax <- setdiff(isolates_data$OTU_ID, taxonomy_ITS_meta$OTU_ID)

# Create empty rows for these OTUs to populate with taxonomy info later
new_rows_ITS_tax <- data.frame(
  OTU_ID = missing_otus_ITS_tax,
  matrix(0, nrow = length(missing_otus_ITS_tax), ncol = 6)  # assuming 7 columns total
)

# Assign the same column names as the taxonomy table
colnames(new_rows_ITS_tax) <- colnames(taxonomy_ITS_meta)

# Make a copy of the isolate metadata and rename relevant columns to match taxonomy format
isolates_data_meta <- isolates_data
isolates_data_meta <- isolates_data_meta %>%
  rename(
    phylum = Phylum,
    class = Class,
    order = Order,
    family = Family,
    genus = Genus,
  )

# Fill in the taxonomy fields for each missing OTU using isolate metadata
for (i in seq_len(nrow(new_rows_ITS_tax))) {
  otu <- new_rows_ITS_tax$OTU_ID[i]
  
  # Find the corresponding row in the isolate metadata
  match_row <- isolates_data_meta[isolates_data_meta$OTU_ID == otu, ]
  
  # Copy over only the matching taxonomy columns
  common_cols <- intersect(colnames(new_rows_ITS_tax), colnames(match_row))
  
  new_rows_ITS_tax[i, common_cols] <- match_row[1, common_cols]
}

# Add the new rows with filled-in taxonomy to the main taxonomy table
taxonomy_ITS_meta <- rbind(taxonomy_ITS_meta, new_rows_ITS_tax)

taxonomy_ITS_meta <- taxonomy_ITS_meta %>%
  dplyr::select(OTU_ID, everything()) %>%      # reorder columns if needed (adjust if your first column is different)
  mutate(domain = "Eukarya") %>%               # add new column
  dplyr::select(OTU_ID, domain, everything()) # place domain as second column


# Compute minimum sample depth (library size), excluding OTU_ID and Cultured columns
min_depth_ITS_meta <- min(colSums(OTU_table_ITS_meta %>% dplyr::select(-OTU_ID, -Cultured)))


# Normalize sample counts to the minimum library size (proportional scaling)
OTU_table_ITS_meta_rarified <- OTU_table_ITS_meta %>%
  mutate(across(-c(OTU_ID, Cultured), ~ (.x / sum(.x)) * min_depth_ITS_meta))


merged_table_ITS <- left_join(OTU_table_ITS_meta_rarified, taxonomy_ITS_meta, by = "OTU_ID")

# Add semicolon-separated taxonomy string by uniting taxonomy columns 9 through 13
merged_table_ITS <- merged_table_ITS %>%
  tidyr::unite("taxonomy", 9:13, sep = ";", remove = FALSE)


# Convert to taxmap format for metacoder analysis
obj_ITS <- parse_tax_data(merged_table_ITS,
                          class_cols = "taxonomy", # The column in the input table containing taxonomy string
                          class_sep = ";")         # Separator between taxonomic ranks

# Calculate taxon abundance per sample
obj_ITS$data$tax_abund <- calc_taxon_abund(obj_ITS, "tax_data",
                                           cols = c("Bull", "Mortensen", "West", "Horn", "Conant", "Greg"))

# Calculate total abundance across all samples for each taxon
obj_ITS$data$tax_abund$total_abund <- rowSums(obj_ITS$data$tax_abund[, c("Bull", "Mortensen", "West", "Horn", "Conant", "Greg")])

# 1. Calculate culture abundance for cultured taxa only
obj_ITS$data$tax_abund$culture_abund <- calc_taxon_abund(obj_ITS, "tax_data", cols = "Cultured")$Cultured


# Generate heat tree plot showing abundance and cultured status
heat_tree(obj_ITS,
          node_label = taxon_names,
          node_size  = obj_ITS$data$tax_abund$total_abund,
          node_color = obj_ITS$data$tax_abund$culture_abund,
          node_color_range = c("gray", "yellow", "red"),
          node_color_axis_label = "Culture Abundance",
          node_size_axis_label  = "Amplicon Abundance",
          edge_size_axis_label  = "Number of OTUs",
          edge_size = n_obs,
          node_size_range = c(0.009, 0.03),
          layout = "davidson-harel",              # primary layout algorithm
          initial_layout = "reingold-tilford", 
          output_file = "meta_coder_ITS_colored.pdf")




#### Metacoder 18s ####

# Read taxonomy data for 18S OTUs; large animals and land plants already removed
meta_18s_tax <- read.csv("taxonomy_18s.csv")

# Join OTU counts table with taxonomy data by OTU_ID
merged_table_18s <- left_join(OTU_table_18s, meta_18s_tax, by = "OTU_ID")

# Add semicolon-separated taxonomy string by uniting taxonomy columns 6 through 11
merged_table_18s <- merged_table_18s %>%
  tidyr::unite("taxonomy", 6:11, sep = ";", remove = FALSE)

# Drop the 'Surface.Control' column, presumably a control sample
merged_table_18s <- merged_table_18s %>% dplyr::select(-Surface.Control)

# Normalize read counts across three samples:
# Step 1: Calculate total reads per sample (columns 2 to 4)
sample_totals_18s <- colSums(merged_table_18s[, 2:4])

# Step 2: Find the smallest library size among samples
min_depth_18s <- min(sample_totals_18s)

# Step 3: Normalize each sample to the smallest depth by proportional scaling
normalized_counts18s <- sweep(merged_table_18s[, 2:4], 2, sample_totals_18s, FUN = "/") * min_depth_18s

# Step 4: Replace original counts with normalized counts in the data frame
merged_table_18s[, 2:4] <- normalized_counts18s

# Convert the merged data frame to metacoder taxmap format, specifying taxonomy column and separator
obj_18s <- parse_tax_data(merged_table_18s,
                          class_cols = "taxonomy",  # Column containing taxonomy string
                          class_sep = ";")          # Taxonomic rank separator

# Calculate taxon abundance across samples Bull, Mortensen, and West
obj_18s$data$tax_abund <- calc_taxon_abund(obj_18s, "tax_data",
                                           cols = c("Bull", "Mortensen", "West"))

# Calculate total abundance for each taxon across the three samples
obj_18s$data$tax_abund$total_abund <- rowSums(obj_18s$data$tax_abund[, c("Bull", "Mortensen", "West")])

# Generate and save a heat tree plot visualizing taxon abundance
heat_tree(obj_18s,
          node_label = taxon_names,                         # Labels for nodes are taxon names
          node_size = obj_18s$data$tax_abund$total_abund,  # Node size scaled by total abundance
          node_color = obj_18s$data$tax_abund$total_abund, # Node color scaled by total abundance
          node_color_axis_label = "Total abundance",
          edge_size_axis_label = "Number of OTUs",
          edge_size = n_obs,                                # Edge size proportional to number of OTUs
          node_size_range = c(0.009, 0.03),                 # Range of node sizes
          layout = "davidson-harel",                        # Primary layout algorithm
          initial_layout = "reingold-tilford",              # Initial layout for optimization
          output_file = "meta_coder_18s.pdf")               # Save plot directly to file

#### Rank Abundance Plot ####

# Set row names of rarified OTU table to OTU IDs (first column), then remove the OTU ID column
rownames(OTU_table_ITS_meta_rarified) <- OTU_table_ITS_meta_rarified[[1]]  
OTU_table_ITS_meta_rarified <- OTU_table_ITS_meta_rarified[, -1]  

# Sum abundances across all samples for each OTU
rank_abun_sums <- rowSums(OTU_table_ITS_meta_rarified)

# Sort OTUs by abundance in descending order
otu_ranked <- sort(rank_abun_sums, decreasing = TRUE)

# Get cultured OTU IDs from isolate data
cultured_otus <- isolates_data$OTU_ID

# Names of ranked OTUs
otu_names_rank_abund <- names(otu_ranked)

# Logical vector indicating if each OTU is cultured
is_cultured_rank_abund <- otu_names_rank_abund %in% cultured_otus

# Create a data frame with ranks and abundances
rank_abundance_df <- data.frame(
  Rank = seq_along(otu_ranked),
  Abundance = otu_ranked
)

# Plot total rank abundance (log scale)
rank_abundance_total <- ggplot(rank_abundance_df, aes(x = Rank, y = Abundance)) +
  geom_line(color = "red") +
  geom_point(shape = 1, size = 1, color = "black") +  # open circles
  scale_y_log10() +
  labs(x = "OTU Rank (ITS2)",
       y = "Abundance (log scale)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),                    # remove grid
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # black border
    axis.line = element_line(color = "black"),       # black axis lines
    axis.ticks = element_line(color = "black")       # black ticks
  )

# Save plot to PDF
ggsave("rank_abundance_total.pdf", plot = rank_abundance_total, width = 14, height = 10)

# Plot only top 50 OTUs
otu_top50 <- otu_ranked[1:51]
otu_names_top50 <- names(otu_top50)

# Cultured status for top 50 OTUs
is_cultured_top50 <- otu_names_top50 %in% cultured_otus

# Data frame with cultured status as a factor for fill color
rank_abundance_top50_df <- data.frame(
  Rank = seq_along(otu_top50),
  Abundance = otu_top50,
  Cultured = factor(ifelse(is_cultured_top50, "Cultured", "Uncultured"))
)

# Plot top 50 with filled circles for cultured, open circles for uncultured
rank_abundance_top50 <- ggplot(rank_abundance_top50_df, aes(x = Rank, y = Abundance)) +
  geom_line(color = "red") +
  geom_point(aes(fill = Cultured), shape = 21, color = "black", size = 2) +  # filled circles with black border
  scale_fill_manual(values = c("Cultured" = "black", "Uncultured" = "white")) +
  scale_y_log10() +
  labs(x = "OTU Rank (ITS2)",
       y = "Abundance (log scale)",
       fill = "") +  # no legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.95, 0.95),   # position legend top right inside plot
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha('white', 0.6), color = "black")
  )

# Display the plot
rank_abundance_top50

# Save the plot
ggsave("rank_abundance_top50_plot.pdf", plot = rank_abundance_top50, width = 14, height = 10)

# Add taxonomy info (phylum) by joining with fungal taxonomy data
rank_abundance_top50_df <- rownames_to_column(rank_abundance_top50_df, var = "OTU_ID")

rank_abundance_top50_df <- rank_abundance_top50_df %>%
  left_join(Taxonomy_ITS_fungi %>% dplyr::select(OTU_ID, phylum), by = "OTU_ID")

# Calculate dissimilarity from percent identity (convert percent to dissimilarity)
ITS_per_ident_diss <- ITS_per_ident %>%
  mutate(Dissimilarity = 100 - percent_identity)

# Join dissimilarity data
rank_abundance_top50_df <- rank_abundance_top50_df %>%
  left_join(ITS_per_ident_diss %>% dplyr::select(OTU_ID, Dissimilarity), by = "OTU_ID")

# Manual labels for top 10 OTUs with formatted text (italic, plain)
manual_labels <- data.frame(
  Rank = 1:10,
  Label = c(
    "italic('Trametes') * ' sp.'",
    "plain('Herpotrichiellaceae sp.')",
    "italic('Naganishia') * ' sp.'",
    "plain('Helotiales sp.')",
    "plain('Didymellaceae sp.')",
    "italic('Alternaria') * ' sp.'",
    "italic('Alternaria') * ' sp.'",
    "italic('Aureobasidium') * ' sp.'",
    "italic('Sporobolomyces') * ' sp.'",
    "italic('Alternaria') * ' sp.'"
  )
)

# Final plot with taxonomy, cultured status, and dissimilarity mapped
rank_abund_manuscript <- ggplot(rank_abundance_top50_df, aes(x = Rank, y = Abundance)) +
  geom_line(color = "black") +
  geom_point(
    aes(color = phylum, shape = Cultured, size = Dissimilarity),
    stroke = 0.3
  ) +
  scale_shape_manual(values = c("Uncultured" = 16, "Cultured" = 17)) +
  scale_size_continuous(range = c(2, 6)) +
  scale_y_log10() +
  labs(
    x = "OTU Rank (ITS2)",
    y = "Abundance (log scale)",
    color = "Phylum",
    shape = "Culture Status",
    size = "Dissimilarity (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha('white', 0.6), color = "black"),
    legend.key.size = unit(0.8, "lines")
  ) +
  geom_text_repel(
    data = manual_labels %>%
      mutate(Abundance = rank_abundance_top50_df$Abundance[Rank]),
    aes(x = Rank, y = Abundance, label = Label),
    parse = TRUE,
    size = 3,
    segment.color = "black",   # line color for label connectors
    nudge_y = 0.1,             # nudge labels upward
    ylim = c(NA, Inf),         # restrict label movement upward only
    box.padding = 0.7,
    force = 1.2,
    max.overlaps = Inf
  )

# Save manuscript-ready plot
ggsave("rank_abund_manuscript.pdf", plot = rank_abund_manuscript, width = 12, height = 8)




