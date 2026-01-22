#This scripts completes the Community analysis of the 16s, 18s and ITS Amplicon data, as well as the Isolates data. And, is used
# to generate the staked bar plots, metacoder trees, isolates, abundance, and diversity figures. (Figs. 3,4,5,7 and Supp Figs 1,2,4,5,6)

#### Load Packages ####
library(Polychrome)
library(scales)
library(broom)
library(ggrastr)
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
library(eulerr)
library(jpeg)         
library(grid) 
library(ggpattern)

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



# ---------- 16s DATA ----------
#OTU table
OTU_table_16s <- read.csv("16s/16s_OTU_table_trimmed.csv") # 16S OTU table

#taxonomy table (based on 0.7 default confidence)
taxonomy_16s <- read_tsv("16s/outputs/taxonomy_16s.tsv") %>%
  dplyr::rename(OTU_ID = `Feature ID`)

# Filter taxonomy to keep only OTUs present in the OTU table
taxonomy_16s <- taxonomy_16s %>%
  dplyr::filter(OTU_ID %in% OTU_table_16s$OTU_ID)

#parse taxon string
taxonomy_16_clean <- taxonomy_16s %>%
  # Split taxon into individual taxonomic ranks
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ", fill = "right") %>%
  
  # Remove prefixes like "d__", "p__", etc.
  mutate(across(everything(), ~ sub("^[a-z]__*", "", .))) %>%
  
  # Replace empty strings or NAs with "unclassified"
  mutate(across(everything(), ~ ifelse(is.na(.) | . == "", "Unclassified", .)))


# ---------- Metadata ----------
metadata_AS <- read.csv(file = "metadata_AS.csv") 

#### Rarefaction and Species Accumulation Curves (Supp Fig 1)####

# ---------- ITS2 Rarefaction by Sample ----------

# Transpose OTU table so samples are rows and OTUs are columns
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
    x = NULL,
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
    x = NULL,
    y = "Observed OTUs (18s)",
    color = "BNG Well"
  ) +
  scale_color_manual(values = accumulation_colors) +
  theme_minimal()

# Results: All three 18S samples reach clear rarefaction plateaus

# 
# ---------- 16S Rarefaction by Sample ----------
# Transpose OTU table so samples are rows and OTUs are columns (as expected by vegan)
otu_mat_16s <- t(OTU_table_16s[, -1])  # remove OTU_ID column before transposing

# Generate rarefaction curves for each sample
# Using a fixed step size and sampling depth equal to the minimum library size
rare_list_16s <- rarecurve(otu_mat_16s, step = 50, sample = min(rowSums(otu_mat_16s)), plot = FALSE)
names(rare_list_16s) <- rownames(otu_mat_16s)

# Convert rarefaction results to a tidy format for ggplot
rare_df_16s <- do.call(rbind, lapply(names(rare_list_16s), function(samp) {
  obs_otus <- rare_list_16s[[samp]]                   # observed OTUs
  reads <- as.numeric(sub("N", "", names(obs_otus)))  # extract read counts from names
  data.frame(Sample = samp, Reads = reads, Observed_OTUs = obs_otus)
}))

# Define color palette for samples
accumulation_colors_16s <- c(
  Bull = "#332288",
  Conant = "#88CCEE",
  Greg = "#44AA99",
  Horn = "#117733",
  Mortensen = "#DDCC77",
  West = "#CC6677"
)

# Plot rarefaction curves for 16S
acc_16s_plot <- ggplot(rare_df_16s, aes(x = Reads, y = Observed_OTUs, color = Sample)) +
  geom_line(size = 1) +
  labs(
    x = "Read Count",
    y = "Observed ASVs (16S)",
    color = "BNG Well"
  ) +
  scale_color_manual(values = accumulation_colors_16s) +
  theme_minimal() +
  theme(legend.position = "right")


# ---------- ITS2 Species Accumulation (By Sample Count) ----------
# Run species accumulation analysis for ITS2
accum_ITS <- specaccum(otu_mat_ITS, method = "random", permutations = 100)
# Format output as dataframe
accum_df_ITS <- data.frame(
  Samples = accum_ITS$sites,
  Richness = accum_ITS$richness,
  SD = accum_ITS$sd
)
# Plot with shaded error
acc_overall_ITS <- ggplot(accum_df_ITS, aes(x = Samples, y = Richness)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, fill = "blue") +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(accum_df_ITS$Samples)))) +
  theme_minimal()

#acc_overall_ITS
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
    x = NULL,
    y = NULL
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(accum_df_18s$Samples)))) +
  theme_minimal()

# ---------- 16s Species Accumulation (By Sample Count) ----------
# Run species accumulation analysis for 16S
accum_16s <- specaccum(otu_mat_16s, method = "random", permutations = 100)
# Format output as dataframe
accum_df_16s <- data.frame(
  Samples = accum_16s$sites,
  Richness = accum_16s$richness,
  SD = accum_16s$sd
)

# Plot with shaded error
acc_overall_16s <- ggplot(accum_df_16s, aes(x = Samples, y = Richness)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, fill = "blue") +
  labs(
    x = "Number of BNG Wells",
    y = NULL
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(accum_df_16s$Samples)))) +
  theme_minimal()





# ---------- Combine All Plots into One Figure ----------

# Stack the four plots (18S and ITS, both rarefaction and accumulation)

accplot <- (
  acc_18s_plot + acc_overall_18s +
    acc_ITS_plot + acc_overall_ITS +
    acc_16s_plot + acc_overall_16s
) +
  plot_layout(ncol = 2, nrow = 3, byrow = TRUE) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 12, face = "bold")  # bigger size & bold here, applied to all plots
  )

accplot
# Save final combined figure to SVG file
#ggsave("accumulation_plots.pdf", plot = accplot, width = 14, height = 10)





#### Stacked Bar Plots (Supp Fig 5)####
#### ITS Fungal Community Composition ####

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

# ---- Plot: Relative Abundance at the Class Level ----

# Agglomerate OTUs to the Class level
physeq_ITs_class <- tax_glom(physeq_rel_ITS, taxrank = "class")

# Extract class names and build a custom color palette
class_levels <- unique(tax_table(physeq_ITs_class)[, "class"])
class_levels <- as.character(class_levels[!is.na(class_levels)])  # Remove NA entries

# Generate a named color palette with enough distinct colors for all classes
palette_24_named <- setNames(createPalette(length(class_levels), c("#ffffff", "#000000")),
                             class_levels)


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






#### Stacked Bar Plots: 16s ####
# ----- Prepare Metadata for 16s Analysis ####
metadata_16s <-read.csv(file = "metadata_AS.csv") 

# Convert metadata to phyloseq sample_data format
metadata_phy_16s <- sample_data(column_to_rownames(metadata_16s, var = "Sample"))

# ----- Filter 16s ASVs ####
#Drop the 1 ASV that was not classified at Domain level
taxonomy_16s_phy <- taxonomy_16_clean %>%
  dplyr::filter(Domain != "Unclassified")

# ----- Create 16s phy seq object ####
# Prepare tax_table object from fungal taxonomy
taxonomy_16s_mat <- as.matrix(column_to_rownames(taxonomy_16s_phy, var = "OTU_ID"))
tax_table_16s_phy <- tax_table(taxonomy_16s_mat)

# Filter OTU table to match only ASVs in taxonomy file
OTU_table_16s_phy <- OTU_table_16s %>%
  semi_join(taxonomy_16s_phy, by = "OTU_ID")

# Convert OTU table to matrix and make OTU IDs rownames
otu_16s_mat <- as.matrix(column_to_rownames(OTU_table_16s_phy, var = "OTU_ID"))
otu_16s_phy <- otu_table(otu_16s_mat, taxa_are_rows = TRUE)  # Indicate OTUs are rows


# ---- Construct Phyloseq Object for 16s ----
physeq_16s <- phyloseq(otu_16s_phy, tax_table_16s_phy, metadata_phy_16s)  # Combine OTU, taxonomy, and metadata

# ---- Transform Counts to Relative Abundances ----
physeq_rel_16s <- transform_sample_counts(physeq_16s, function(x) x / sum(x))  # Normalize within each sample

# ---- Sample Order for Plotting ----
desired_order <- c("Conant", "Greg", "Mortensen", "Bull", "Horn", "West")  # Desired sample order

# Update sample names and factor levels to enforce order
sample_data(physeq_rel_16s)$Sample <- rownames(sample_data(physeq_rel_16s))
sample_data(physeq_rel_16s)$Sample <- factor(
  sample_data(physeq_rel_16s)$Sample,
  levels = desired_order
)
# ---- Create "Other" for less than 1% phyla ----
# Transform physeq object to data frame
df_16s <- psmelt(physeq_rel_16s)  # psmelt() flattens phyloseq into long-form

# Sum relative abundances across all samples
low_phyla_16s <- df_16s %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  mutate(rel_abund = total_abundance / sum(total_abundance)) %>%
  filter(rel_abund < 0.01) %>%
  pull(Phylum)

df_16s <- df_16s %>%
  dplyr::mutate(Phylum = ifelse(Phylum %in% low_phyla_16s, "Other", Phylum))

#Assign colors to other and unclassified
# 1. Define archaeal phyla
archaea_phyla <- c("Halobacterota", "Euryarchaeota")

# 2. Get phylum levels, put "Other" and "Unclassified" last
phyla_levels <- unique(df_16s$Phylum)
phyla_levels <- setdiff(phyla_levels, c("Other", "Unclassified"))
phyla_levels <- sort(phyla_levels)
phyla_levels <- c(phyla_levels, "Other", "Unclassified")

# 3. Create display labels with * for Archaea
phyla_display <- phyla_levels
phyla_display[phyla_display %in% archaea_phyla] <- paste0(phyla_display[phyla_display %in% archaea_phyla], "*")

# 4. Generate a distinct palette for non-gray/black
palette_main <- createPalette(length(phyla_levels) - 2, 
                              seedcolors = c("#FF0000", "#00FF00", "#0000FF"))

# 5. Assign colors to display names
phyla_colors <- setNames(palette_main, phyla_display[1:(length(phyla_levels) - 2)])
phyla_colors["Other"] <- "gray70"
phyla_colors["Unclassified"] <- "black"

# 6. Set factor with original levels, but new labels
df_16s$Phylum <- factor(df_16s$Phylum, levels = phyla_levels, labels = phyla_display)


# ---- Plot: Relative Abundance at the Phylum Level ----
phylum_16s_plot <- ggplot(df_16s, aes(x = sample_Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  scale_fill_manual(values = phyla_colors, guide = guide_legend(ncol = 2)) +
  theme_minimal() +
  ylab("Relative Abundance") +
  labs(
    fill = "Bacterial and Archaeal Phylum (16s)",
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

phylum_16s_plot
# ---- Plot: Relative Abundance at the Class Level ----
# ---- Collapse low abundance classes ----
low_classes_16s <- df_16s %>%
  group_by(Class) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  mutate(rel_abund = total_abundance / sum(total_abundance)) %>%
  filter(rel_abund < 0.01) %>%
  pull(Class)

df_16s <- df_16s %>%
  dplyr::mutate(Class = ifelse(Class %in% low_classes_16s, "Other", Class))

# ---- Assign colors and annotate archaeal classes ----
# 2. Get class levels, put "Other" and "Unclassified" last
class_levels <- unique(df_16s$Class)
class_levels <- setdiff(class_levels, c("Other", "Unclassified"))
class_levels <- sort(class_levels)
class_levels <- c(class_levels, "Other", "Unclassified")


# 1. Define archaeal classes (edit as needed based on your taxonomy)
archaea_classes <- c("Methanosarcinia", "Methanomicrobia", "Methanobacteria")

# 3. Create display labels with * for Archaea
class_display <- class_levels
class_display[class_display %in% archaea_classes] <- paste0(class_display[class_display %in% archaea_classes], "*")

# 4. Generate distinct color palette
palette_main <- createPalette(length(class_levels) - 2, 
                              seedcolors = c("#FF0000", "#00FF00", "#0000FF"))

# 5. Assign colors to display names
class_colors <- setNames(palette_main, class_display[1:(length(class_levels) - 2)])
class_colors["Other"] <- "gray70"
class_colors["Unclassified"] <- "black"

# 6. Set factor in df_16s with display labels
df_16s$Class <- factor(df_16s$Class, levels = class_levels, labels = class_display)


# ---- Draw Plot ----
class_16s_plot <- ggplot(df_16s, aes(x = sample_Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  scale_fill_manual(values = class_colors, guide = guide_legend(ncol = 2)) +
  theme_minimal() +
  ylab("") +
  labs(fill = "Bacterial and Archaeal Class (16s)") +
  theme(
    axis.text.x = element_blank(),                     # Remove x-axis labels
    axis.ticks.x = element_blank(),                    # Remove x-axis ticks
    axis.title.x = element_blank(),                    # Remove x-axis title
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )
class_16s_plot 
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
  labs(fill = "Fungal Phylum (ITS)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ITS class plot with rotated x-axis labels and full axis details
plot_ITS_class <- plot_bar(physeq_ITs_class, x = "sample_Sample", fill = "class") +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +  # Set linewidth to 0 to match
  scale_fill_manual(values = palette_class) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("") +
  labs(fill = "Fungal Class (ITS)") +
  theme(
    axis.text.x = element_blank(),                       # No x-axis text
    axis.ticks.x = element_blank(),                      # No x-axis ticks
    axis.title.x = element_blank(),                      # No x-axis title
    axis.text.y = element_text(size = 16, color = "black"),  # Big y-axis numbers
    axis.ticks.y = element_line(),                       # Y-axis ticks
    axis.title.y = element_text(size = 18),              # Big y-axis label
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),              # Legends match
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )



# ----- Generate stacked bar plots for 18s fungi datasets -----

# Create a data frame for the missing samples. Add * to empty sites
asterisk_df <- data.frame(
  sample_Sample = c("Conant", "Greg", "Horn"),
  y = 0.45,  # Adjust this for positioning
  label = "*"
)


# Fungi phylum plot without x-axis labels
plot_fungi_18s_phyla <- plot_bar(physeq_rel_18s_fungi, x = "sample_Sample", fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  geom_text(
    data = asterisk_df,
    aes(x = sample_Sample, y = y, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 12,  # Increase for better visibility
    vjust = 0
  ) +
  scale_fill_manual(values = palette_phylum) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  labs(fill = "Fungal Phylum (18s)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Fungi class plot without x-axis labels
plot_fungi_18s_class <- plot_bar(physeq_18s_class, x = "sample_Sample", fill = "Class") +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  geom_text(
    data = asterisk_df,
    aes(x = sample_Sample, y = y, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 12,  # Increase for better visibility
    vjust = 0
  ) +
  scale_fill_manual(values = palette_class) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("") +
  labs(fill = "Fungal Class (18s)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ----- Generate stacked bar plots for 18s eukaryote datasets -----
# Division plot without x-axis labels
plot_euk_18s_division <- plot_bar(physeq_rel_18s_euk, x = "sample_Sample", fill = "Division") +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  geom_text(
    data = asterisk_df,
    aes(x = sample_Sample, y = y, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 12,
    vjust = 0
  ) +
  scale_fill_manual(values = palette_division) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("Relative Abundance") +
  labs(
    x = "BNG Well",
    fill = "Microeukaryotic Division (18s)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"), # bigger x-axis labels
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 22),                                         # bigger x-axis title
    axis.text.y = element_text(size = 16, color = "black"),                         # bigger y-axis numbers
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),                                         # bigger y-axis title
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )


# plot euk class
plot_euk_18s_class <- plot_bar(physeq_18s_euk_class, x = "sample_Sample", fill = "Class") +
  geom_bar(stat = "identity", position = "stack", color = NA, linewidth = 0) +
  geom_text(
    data = asterisk_df,
    aes(x = sample_Sample, y = y, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 12,
    vjust = 0
  ) +
  scale_fill_manual(values = palette_class) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  ylab("") +
  labs(
    x = "BNG Well",
    fill = "Microeukaryotic Class (18s)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"), # bigger x-axis labels
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 22),                                         # bigger x-axis title
    axis.text.y = element_text(size = 16, color = "black"),                         # bigger y-axis numbers
    axis.ticks.y = element_line(),
    axis.title.y = element_text(size = 18),                                         # bigger y-axis title
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )
plot_euk_18s_class

# ----- Combine and save stacked bar plots -----
  

#Combine all into single plot
total_combined_bars <- plot_ITS_phyla + plot_ITS_class +
  phylum_16s_plot + class_16s_plot +
  plot_fungi_18s_phyla + plot_fungi_18s_class +
  plot_euk_18s_division + plot_euk_18s_class +
  plot_layout(ncol = 2, nrow = 4) +
  plot_annotation(
    tag_levels = list(c("(A) 1.", "2.",
                        "(B) 1.", "2.",
                        "(C) 1.", "2.",
                        "    3.", "4.")),  # <-- Closed properly
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(plot.tag = element_text(size = 35, face = "bold"))  # <-- Custom tag style
  ) &
  theme(
    plot.tag = element_text(size = 30, face = "bold"),  # Applies to all plots
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("total_combined_bars.png",
       plot = total_combined_bars,
       width = 30, height = 26, units = "in",
       dpi = 600,
       device = "png",
       limitsize = FALSE)















#### Community Analysis ####
#----Calculate Chao Richness -----
# ----- 16s ####
### --- TOTAL 16S RICHNESS (ALL TAXA) --- ###

# Transpose OTU table so samples are rows, OTUs are columns
# Then remove first row (often taxonomy row or header)
OTU_table_16s_diversity <- t(OTU_table_16s_phy)[-1, ]

# Convert all values to numeric (needed for estimateR)
OTU_table_16s_diversity_numeric <- apply(OTU_table_16s_diversity, 2, as.numeric)

# Calculate Chao1 richness using vegan::estimateR()
chao1_total <- estimateR(OTU_table_16s_diversity_numeric)["S.chao1", ]

# Create dataframe with Chao1 values and sample names
chao1_df_total <- data.frame(
  Sample = rownames(OTU_table_16s_diversity),
  Chao1_16S_total = chao1_total
)


### --- BACTERIAL RICHNESS --- ###

# Filter taxonomy table to include only OTUs classified as Bacteria
taxonomy_bacteria <- taxonomy_16s_phy %>%
  filter(Domain == "Bacteria")

# Keep only bacterial OTUs from OTU table
OTU_table_16s_bacteria <- OTU_table_16s_phy %>%
  semi_join(taxonomy_bacteria, by = "OTU_ID")

# Transpose and clean bacterial OTU table
OTU_table_16s_bact_div <- t(OTU_table_16s_bacteria)[-1, ]
OTU_table_16s_bact_numeric <- apply(OTU_table_16s_bact_div, 2, as.numeric)

# Calculate Chao1 richness for bacteria
chao1_bact <- estimateR(OTU_table_16s_bact_numeric)["S.chao1", ]

# Create dataframe for merging
chao1_df_bact <- data.frame(
  Sample = rownames(OTU_table_16s_bact_div),
  Chao1_16S_bacteria = chao1_bact
)


### --- ARCHAEAL RICHNESS --- ###

# Filter taxonomy to only include Archaea
taxonomy_archaea <- taxonomy_16s_phy %>%
  filter(Domain == "Archaea")

# Keep only archaeal OTUs from OTU table
OTU_table_16s_archaea <- OTU_table_16s_phy %>%
  semi_join(taxonomy_archaea, by = "OTU_ID")

# Transpose and clean archaeal OTU table
OTU_table_16s_arch_div <- t(OTU_table_16s_archaea)[-1, ]
OTU_table_16s_arch_numeric <- apply(OTU_table_16s_arch_div, 2, as.numeric)

# Calculate Chao1 richness for archaea
chao1_arch <- estimateR(OTU_table_16s_arch_numeric)["S.chao1", ]

# Create dataframe for merging
chao1_df_arch <- data.frame(
  Sample = rownames(OTU_table_16s_arch_div),
  Chao1_16S_archaea = chao1_arch
)


### --- MERGE ALL CHAO1 RICHNESS VALUES INTO METADATA --- ###

# Merge total, bacterial, and archaeal richness into metadata
metadata_AS <- metadata_AS %>%
  left_join(chao1_df_total, by = "Sample") %>%
  left_join(chao1_df_bact, by = "Sample") %>%
  left_join(chao1_df_arch, by = "Sample")

#Results: 5297 total bacteria ASVs. 801 ASV for archea

# ----- ITS ####
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

# Rename Chao1 column in chao1_df to Chao1_ITS
names(chao1_df)[names(chao1_df) == "Chao1"] <- "Chao1_fungi"

#Add to metadata
metadata_AS <- merge(metadata_AS, chao1_df, by = "Sample", all.x = TRUE)

# Merge the Chao1 richness estimates with sample metadata for downstream analysis
#metadata_ITS <- merge(chao1_df, metadata_ITS, by = "Sample")


#comined plot
# Calculate scale factor for rescaling bacteria and total 16S richness
scale_factor <- max(metadata_AS$Chao1_16S_archaea, na.rm = TRUE) / max(metadata_AS$Chao1_16S_bacteria, na.rm = TRUE)

dual_richness_plot <- ggplot(metadata_AS, aes(x = Depth..m.)) + 
  
  # Archaea (primary y-axis)
  geom_line(aes(y = Chao1_16S_archaea, color = factor("Archaea", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 1.2) +
  geom_point(aes(y = Chao1_16S_archaea, color = factor("Archaea", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 2) +
  
  # Fungi (primary y-axis)
  geom_line(aes(y = Chao1_fungi, color = factor("Fungi", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 1.2) +
  geom_point(aes(y = Chao1_fungi, color = factor("Fungi", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 2) +
  
  # Bacteria (secondary y-axis, rescaled)
  geom_line(aes(y = Chao1_16S_bacteria * scale_factor, color = factor("Bacteria", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 1.2) +
  geom_point(aes(y = Chao1_16S_bacteria * scale_factor, color = factor("Bacteria", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 2) +
  
  # Total 16S (Bacteria + Archaea)
  geom_line(aes(y = Chao1_16S_total * scale_factor, color = factor("Bacteria and Archaea", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 1.2) +
  geom_point(aes(y = Chao1_16S_total * scale_factor, color = factor("Bacteria and Archaea", levels = c("Fungi", "Archaea", "Bacteria", "Bacteria and Archaea"))), size = 2) +
  
  # Y axis and secondary y-axis
  scale_y_continuous(
    name = "Fungal OTU Richness\nArchaeal ASV Richness",
    sec.axis = sec_axis(~ . / scale_factor, name = "Bacterial ASV Richness\nBacterial and Archaeal ASV Richness")
  ) +
  
  # Legend colors and title
  scale_color_manual(
    values = c(
      "Fungi" = "forestgreen",
      "Archaea" = "firebrick",
      "Bacteria" = "steelblue",
      "Bacteria and Archaea" = "purple"
    )
  ) +
  
  labs(
    x = "Depth (m)",
    color = "Taxonomic Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.title.y.right = element_text(size = 14)
  )


# Print plot
print(dual_richness_plot)





#test for significance between richness trends
cor.test(metadata_AS$Chao1_16S_total, metadata_AS$Chao1_fungi, method = "spearman")
cor.test(metadata_AS$Chao1_16S_bacteria, metadata_AS$Chao1_fungi, method = "spearman")
cor.test(metadata_AS$Chao1_16S_archaea, metadata_AS$Chao1_fungi, method = "spearman")
#result: fungal richness is significantly positively correlated with total 16S abundance, but not with bacterial or archaeal richness alone.

# cor_values and p_values contain the correlations and p-values respectively





#---- Shannon -----
# --- TOTAL 16S SHANNON DIVERSITY ---
OTU_table_16s_diversity_numeric <- apply(OTU_table_16s_diversity, 2, as.numeric)

shannon_total <- diversity(OTU_table_16s_diversity_numeric, index = "shannon")

shannon_df_total <- data.frame(
  Sample = rownames(OTU_table_16s_diversity),
  Shannon_16S_total = shannon_total
)


# --- BACTERIAL SHANNON DIVERSITY ---
OTU_table_16s_bact_numeric <- apply(OTU_table_16s_bact_div, 2, as.numeric)

shannon_bact <- diversity(OTU_table_16s_bact_numeric, index = "shannon")

shannon_df_bact <- data.frame(
  Sample = rownames(OTU_table_16s_bact_div),
  Shannon_16S_bacteria = shannon_bact
)


# --- ARCHAEAL SHANNON DIVERSITY ---
OTU_table_16s_arch_numeric <- apply(OTU_table_16s_arch_div, 2, as.numeric)

shannon_arch <- diversity(OTU_table_16s_arch_numeric, index = "shannon")

shannon_df_arch <- data.frame(
  Sample = rownames(OTU_table_16s_arch_div),
  Shannon_16S_archaea = shannon_arch
)


# --- FUNGAL SHANNON DIVERSITY ---
OTU_table_ITS_fungi_numeric <- apply(OTU_table_ITS_fungi_diversity, 2, as.numeric)

shannon_fungi <- diversity(OTU_table_ITS_fungi_numeric, index = "shannon")

shannon_df_fungi <- data.frame(
  Sample = rownames(OTU_table_ITS_fungi_diversity),
  Shannon_fungi = shannon_fungi
)

# Merge all Shannon diversity data frames into metadata_AS
metadata_AS <- metadata_AS %>%
  left_join(shannon_df_total, by = "Sample") %>%
  left_join(shannon_df_bact, by = "Sample") %>%
  left_join(shannon_df_arch, by = "Sample") %>%
  left_join(shannon_df_fungi, by = "Sample")



shannon_plot <- ggplot(metadata_AS, aes(x = Depth..m.)) + 
  
  # Archaea
  geom_line(aes(y = Shannon_16S_archaea, color = "Archaea"), size = 1.2) +
  geom_point(aes(y = Shannon_16S_archaea, color = "Archaea"), size = 2) +
  
  # Fungi
  geom_line(aes(y = Shannon_fungi, color = "Fungi"), size = 1.2) +
  geom_point(aes(y = Shannon_fungi, color = "Fungi"), size = 2) +
  
  # Bacteria
  geom_line(aes(y = Shannon_16S_bacteria, color = "Bacteria"), size = 1.2) +
  geom_point(aes(y = Shannon_16S_bacteria, color = "Bacteria"), size = 2) +
  
  # Bacteria and Archaea combined line
  geom_line(aes(y = Shannon_16S_total, color = "Bacteria and Archaea"), size = 1.2) +
  geom_point(aes(y = Shannon_16S_total, color = "Bacteria and Archaea"), size = 2) +
  
  # Y axis label
  scale_y_continuous(name = "Shannon") +
  
  # Colors for groups
  scale_color_manual(values = c(
    "Bacteria" = "steelblue",
    "Archaea" = "firebrick",
    "Fungi" = "forestgreen",
    "Bacteria and Archaea" = "purple"
  )) +
  
  # Labels and theme
  labs(x = "Depth (m)", color = "Taxonomic Group") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 14)
  )

# Print plot
print(shannon_plot)








##### PIE #####
# Function to calculate Hurlbert's PIE for a single sample vector
calc_hurlbert_pie <- function(abundances) {
  N <- sum(abundances)
  if (N <= 1) {
    return(NA)  # Undefined if N <= 1
  }
  p <- abundances / N
  pie <- (N / (N - 1)) * (1 - sum(p^2))
  return(pie)
}


# --- TOTAL 16S Hurlbert's PIE ---
OTU_table_16s_diversity_numeric <- apply(OTU_table_16s_diversity, 2, as.numeric)

pie_total <- apply(OTU_table_16s_diversity_numeric, 1, calc_hurlbert_pie)

pie_df_total <- data.frame(
  Sample = rownames(OTU_table_16s_diversity),
  Hurlbert_PIE_16S_total = pie_total
)


# --- BACTERIAL Hurlbert's PIE ---
OTU_table_16s_bact_numeric <- apply(OTU_table_16s_bact_div, 2, as.numeric)

pie_bact <- apply(OTU_table_16s_bact_numeric, 1, calc_hurlbert_pie)

pie_df_bact <- data.frame(
  Sample = rownames(OTU_table_16s_bact_div),
  Hurlbert_PIE_16S_bacteria = pie_bact
)


# --- ARCHAEAL Hurlbert's PIE ---
OTU_table_16s_arch_numeric <- apply(OTU_table_16s_arch_div, 2, as.numeric)

pie_arch <- apply(OTU_table_16s_arch_numeric, 1, calc_hurlbert_pie)

pie_df_arch <- data.frame(
  Sample = rownames(OTU_table_16s_arch_div),
  Hurlbert_PIE_16S_archaea = pie_arch
)


# --- FUNGAL Hurlbert's PIE ---
OTU_table_ITS_fungi_numeric <- apply(OTU_table_ITS_fungi_diversity, 2, as.numeric)

pie_fungi <- apply(OTU_table_ITS_fungi_numeric, 1, calc_hurlbert_pie)

pie_df_fungi <- data.frame(
  Sample = rownames(OTU_table_ITS_fungi_diversity),
  Hurlbert_PIE_fungi = pie_fungi
)


# Merge all into metadata_AS
metadata_AS <- metadata_AS %>%
  left_join(pie_df_total, by = "Sample") %>%
  left_join(pie_df_bact, by = "Sample") %>%
  left_join(pie_df_arch, by = "Sample") %>%
  left_join(pie_df_fungi, by = "Sample")





evenness_plot <- ggplot(metadata_AS, aes(x = Depth..m.)) + 
  
  # Archaea
  geom_line(aes(y = Hurlbert_PIE_16S_archaea, color = "Archaea"), size = 1.2) +
  geom_point(aes(y = Hurlbert_PIE_16S_archaea, color = "Archaea"), size = 2) +
  
  # Fungi
  geom_line(aes(y = Hurlbert_PIE_fungi, color = "Fungi"), size = 1.2) +
  geom_point(aes(y = Hurlbert_PIE_fungi, color = "Fungi"), size = 2) +
  
  # Bacteria
  geom_line(aes(y = Hurlbert_PIE_16S_bacteria, color = "Bacteria"), size = 1.2) +
  geom_point(aes(y = Hurlbert_PIE_16S_bacteria, color = "Bacteria"), size = 2) +
  
  # Combined Bacteria and Archaea
  geom_line(aes(y = Hurlbert_PIE_16S_total, color = "Bacteria and Archaea"), size = 1.2) +
  geom_point(aes(y = Hurlbert_PIE_16S_total, color = "Bacteria and Archaea"), size = 2) +
  
  # Y axis label
  scale_y_continuous(name = "Evenness (PIE)") +
  
  # Colors for groups
  scale_color_manual(values = c(
    "Bacteria" = "steelblue",
    "Archaea" = "firebrick",
    "Fungi" = "forestgreen",
    "Bacteria and Archaea" = "purple"
  )) +
  
  # Labels and theme
  labs(x = "Depth (m)", color = "Taxonomic Group") +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 14)
  )

print(evenness_plot)



























#----Distance Decay -----

# ------------------------------------
# PREP: Metadata & Coordinates
# ------------------------------------

rownames(metadata_AS) <- metadata_AS$Sample
coords <- metadata_AS[, c("longitude", "latitude")]

# ------------------------------------
# 1. FUNGI (ITS)
# ------------------------------------

# Prepare & rarefy ITS
otu_counts_ITS <- OTU_table_ITS_fungi[, -1]
rownames(otu_counts_ITS) <- OTU_table_ITS_fungi[[1]]
otu_counts_mat_ITS <- as.matrix(otu_counts_ITS)
mode(otu_counts_mat_ITS) <- "numeric"
otu_counts_mat_ITS_t <- t(otu_counts_mat_ITS)
otu_rarefied_ITS_t <- rrarefy(otu_counts_mat_ITS_t, sample = min(rowSums(otu_counts_mat_ITS_t)))
otu_rarefied_ITS <- t(otu_rarefied_ITS_t)

# Bray-Curtis similarity
bc_sim_ITS <- 1 - as.matrix(vegdist(t(otu_rarefied_ITS), method = "bray"))
geo_dist_ITS <- distm(coords[colnames(otu_rarefied_ITS), ], fun = distHaversine) / 1000
pairs_ITS <- combn(colnames(otu_rarefied_ITS), 2, simplify = FALSE)

# Ensure dimnames are set
samples_ITS <- colnames(otu_rarefied_ITS)
rownames(bc_sim_ITS) <- samples_ITS
colnames(bc_sim_ITS) <- samples_ITS

rownames(geo_dist_ITS) <- samples_ITS
colnames(geo_dist_ITS) <- samples_ITS


decay_df_ITS <- do.call(rbind, lapply(pairs_ITS, function(pair) {
  i <- pair[1]; j <- pair[2]
  data.frame(Sample1 = i, Sample2 = j,
             Geo_Distance_km = geo_dist_ITS[i, j],
             Bray_Curtis_Similarity = bc_sim_ITS[i, j],
             Group = "Fungi")
}))


# ------------------------------------
# 2. TOTAL 16S (Rarefy full table first)
# ------------------------------------
otu_counts_16S <- OTU_table_16s_phy[, -1]
rownames(otu_counts_16S) <- OTU_table_16s_phy[[1]]
otu_mat_16S <- as.matrix(otu_counts_16S)
mode(otu_mat_16S) <- "numeric"
otu_mat_16S_t <- t(otu_mat_16S)
otu_rarefied_16S_t <- rrarefy(otu_mat_16S_t, sample = min(rowSums(otu_mat_16S_t)))
otu_rarefied_16S <- t(otu_rarefied_16S_t)

samples_total <- colnames(otu_rarefied_16S)

bc_sim_total <- 1 - as.matrix(vegdist(t(otu_rarefied_16S), method = "bray"))
rownames(bc_sim_total) <- samples_total
colnames(bc_sim_total) <- samples_total

geo_dist_total <- distm(coords[samples_total, ], fun = distHaversine) / 1000
rownames(geo_dist_total) <- samples_total
colnames(geo_dist_total) <- samples_total

pairs_total <- combn(samples_total, 2, simplify = FALSE)

decay_df_total <- do.call(rbind, lapply(pairs_total, function(pair) {
  i <- pair[1]; j <- pair[2]
  data.frame(Sample1 = i, Sample2 = j,
             Geo_Distance_km = geo_dist_total[i, j],
             Bray_Curtis_Similarity = bc_sim_total[i, j],
             Group = "Total_16S")
}))

# ------------------------------------
# 3. BACTERIA (subset from rarefied total 16S)
# ------------------------------------
bact_otus <- taxonomy_bacteria %>%
  pull(OTU_ID) %>%
  intersect(rownames(otu_rarefied_16S))

otu_bact <- otu_rarefied_16S[bact_otus, ]
samples_bact <- colnames(otu_bact)

bc_sim_bact <- 1 - as.matrix(vegdist(t(otu_bact), method = "bray"))
rownames(bc_sim_bact) <- samples_bact
colnames(bc_sim_bact) <- samples_bact

geo_dist_bact <- distm(coords[samples_bact, ], fun = distHaversine) / 1000
rownames(geo_dist_bact) <- samples_bact
colnames(geo_dist_bact) <- samples_bact

pairs_bact <- combn(samples_bact, 2, simplify = FALSE)

decay_df_bact <- do.call(rbind, lapply(pairs_bact, function(pair) {
  i <- pair[1]; j <- pair[2]
  data.frame(Sample1 = i, Sample2 = j,
             Geo_Distance_km = geo_dist_bact[i, j],
             Bray_Curtis_Similarity = bc_sim_bact[i, j],
             Group = "Bacteria")
}))

# ------------------------------------
# 4. ARCHAEA (subset from rarefied total 16S)
# ------------------------------------
arch_otus <- taxonomy_archaea %>%
  pull(OTU_ID) %>%
  intersect(rownames(otu_rarefied_16S))

otu_arch <- otu_rarefied_16S[arch_otus, ]
samples_arch <- colnames(otu_arch)

bc_sim_arch <- 1 - as.matrix(vegdist(t(otu_arch), method = "bray"))
rownames(bc_sim_arch) <- samples_arch
colnames(bc_sim_arch) <- samples_arch

geo_dist_arch <- distm(coords[samples_arch, ], fun = distHaversine) / 1000
rownames(geo_dist_arch) <- samples_arch
colnames(geo_dist_arch) <- samples_arch

pairs_arch <- combn(samples_arch, 2, simplify = FALSE)

decay_df_arch <- do.call(rbind, lapply(pairs_arch, function(pair) {
  i <- pair[1]; j <- pair[2]
  data.frame(Sample1 = i, Sample2 = j,
             Geo_Distance_km = geo_dist_arch[i, j],
             Bray_Curtis_Similarity = bc_sim_arch[i, j],
             Group = "Archaea")
}))

# ------------------------------------
# COMBINE ALL DECAY DATAFRAMES
# ------------------------------------
decay_all <- bind_rows(decay_df_total, decay_df_bact, decay_df_arch, decay_df_ITS)

decay_all$Group <- recode(decay_all$Group, "Total_16S" = "Bacteria and Archaea")



combined_decay <- ggplot(decay_all, aes(x = Geo_Distance_km, y = Bray_Curtis_Similarity, color = Group, fill = Group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linetype = "solid",
    size = 1,
    alpha = 0.15
  ) +
  labs(
    x = "Geographic Distance (km)",
    y = "Community Similarity",
    color = "Taxonomic Group",
    fill = "Taxonomic Group"
  ) +
  scale_color_manual(values = c(
    "Bacteria and Archaea" = "purple",
    "Bacteria" = "steelblue",
    "Archaea" = "firebrick",
    "Fungi" = "forestgreen"
  )) +
  scale_fill_manual(values = c(
    "Bacteria and Archaea" = alpha("purple", 0.15),
    "Bacteria" = alpha("steelblue", 0.15),
    "Archaea" = alpha("firebrick", 0.15),
    "Fungi" = alpha("forestgreen", 0.15)
  )) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank()
  )





#linear models of decay

# Group-wise linear models
slopes_by_group <- decay_all %>%
  group_by(Group) %>%
  do(tidy(lm(Bray_Curtis_Similarity ~ Geo_Distance_km, data = .))) %>%
  filter(term == "Geo_Distance_km")  # Only keep the slope term

print(slopes_by_group)

# Group-wise linear models
slopes_by_group <- decay_all %>%
  group_by(Group) %>%
  do(tidy(lm(Bray_Curtis_Similarity ~ Geo_Distance_km, data = .))) %>%
  filter(term == "Geo_Distance_km")  # Only keep the slope term

print(slopes_by_group)








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



#### 16s vs ITS richness ####
ITS_vs_16s_rich <- ggplot(metadata_AS, aes(x = Chao1_16S_total, y = Chao1_fungi)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "solid", size = 1.2) +
  labs(
    x = "Bacterial + Archaeal ASV Richness",
    y = "Fungal OTU Richness",
    title = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.y = element_text(color = "forestgreen", size = 15),
    axis.title.x = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    legend.position = "none"
  )

ITS_vs_16s_rich
#----Combined plot -----
# Get ranges
decay_df <- merge(decay_df_ITS, decay_df_16s, 
                  by = c("Sample1", "Sample2", "Geo_Distance_km"),
                  suffixes = c("_ITS", "_16S"))

range_ITS <- range(decay_df$Bray_Curtis_Similarity_ITS, na.rm = TRUE)
range_16S <- range(decay_df$Bray_Curtis_Similarity_16S, na.rm = TRUE)

decay_df$Bray_Curtis_16S_rescaled <- scales::rescale(decay_df$Bray_Curtis_Similarity_16S,
                                                     to = range_ITS,
                                                     from = range_16S)


ggplot(decay_df, aes(x = Geo_Distance_km)) +
  # ITS similarity points and line with blue shadow
  geom_point(aes(y = Bray_Curtis_Similarity_ITS), color = "blue", alpha = 0.6, size = 2) +
  geom_smooth(aes(y = Bray_Curtis_Similarity_ITS), method = "lm", color = "blue", linetype = "dashed", fill = "lightblue", alpha = 0.3) +
  
  # 16S similarity points and line with red shadow (rescaled)
  geom_point(aes(y = Bray_Curtis_16S_rescaled), color = "red", alpha = 0.6, size = 2) +
  geom_smooth(aes(y = Bray_Curtis_16S_rescaled), method = "lm", color = "red", linetype = "dotted", fill = "salmon", alpha = 0.3) +
  
  # Primary and secondary y-axis
  scale_y_continuous(
    name = "Bray-Curtis Similarity (ITS)",
    sec.axis = sec_axis(
      trans = ~ scales::rescale(., to = range_16S, from = range_ITS),
      name = "Bray-Curtis Similarity (16S)"
    )
  ) +
  xlab("Geographic Distance (km)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red")
  )





# ITS linear model
lm_ITS <- lm(Bray_Curtis_Similarity_ITS ~ Geo_Distance_km, data = decay_df)
summary(lm_ITS)

# 16S linear model
lm_16S <- lm(Bray_Curtis_Similarity_16S ~ Geo_Distance_km, data = decay_df)
summary(lm_16S)


# Spearman correlations
spearman_ITS <- cor.test(decay_df$Geo_Distance_km, decay_df$Bray_Curtis_Similarity_ITS, method = "spearman")
spearman_16S <- cor.test(decay_df$Geo_Distance_km, decay_df$Bray_Curtis_Similarity_16S, method = "spearman")



##### Combined Community Plot #####
#make plots without legends
dual_richness_plot_nolegend <- dual_richness_plot + theme(legend.position = "none")
shannon_plot_nolegend      <- shannon_plot + theme(legend.position = "none")
evenness_plot_nolegend     <- evenness_plot + theme(legend.position = "none")
combined_decay_nolegend    <- combined_decay + theme(legend.position = "none")

# Extract legend (from one plot is enough)
legend_grob <- cowplot::get_legend(
  dual_richness_plot + theme(legend.position = "bottom")
)

# Build the layout
combined_community_plots <-
  (
    (dual_richness_plot_nolegend + ITS_vs_16s_rich) /
      (evenness_plot_nolegend + shannon_plot_nolegend) /
      combined_decay_nolegend
  ) /
  wrap_elements(legend_grob) +
  plot_annotation(
    tag_levels = list(c("(A)", "(B)", "(C)", "(D)", "(E)")),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.tag = element_text(size = 30, face = "bold")
    )
  ) &
  theme(plot.title.position = "plot")

# Specify patchwork layout heights
combined_community_plots <- combined_community_plots +
  plot_layout(heights = c(1, 1, 1, 0.1))

# Print or save
print(combined_community_plots)



#ggsave("combined_community_plots.pdf", plot = combined_community_plots, width = 17, height = 12)


  
#### Test Relationships between variables ####


cor.test(metadata_AS$Total.Microbial.Abundance..cells.mL., metadata_AS$NPDOC, method = "spearman", conf.level = 0.95)


#### Metacoder/ Core Microbiome ####
#### Metacoder ITS ####
# only keep fungi
taxonomy_ITS_meta <- taxonomy_ITS %>%
  filter(kingdom == "Fungi")

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
heat_tree_ITS <- heat_tree(obj_ITS,
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
          initial_layout = "reingold-tilford")


#output_file = "meta_coder_ITS_fungi.pdf")




#Let us remake the plot but only showing fungi (excludes the algae)

#### Metacoder 18s ####

# Read taxonomy data for 18S OTUs; large animals and land plants already removed
meta_18s_tax <- read.csv("taxonomy_18s.csv")

# Join OTU counts table with taxonomy data by OTU_ID
merged_table_18s <- left_join(OTU_table_18s, meta_18s_tax, by = "OTU_ID")

# Add semicolon-separated taxonomy string by uniting taxonomy columns 6 through 11
merged_table_18s <- merged_table_18s %>%
  tidyr::unite("taxonomy", 6:11, sep = ";", remove = FALSE)

# Drop the 'Surface.Control' column, a control sample
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
heat_tree_18s <- heat_tree(obj_18s,
          node_label = taxon_names,                         # Labels for nodes are taxon names
          node_size = obj_18s$data$tax_abund$total_abund,  # Node size scaled by total abundance
          node_color = obj_18s$data$tax_abund$total_abund, # Node color scaled by total abundance
          node_color_axis_label = "Total abundance",
          edge_size_axis_label = "Number of OTUs",
          edge_size = n_obs,                                # Edge size proportional to number of OTUs
          node_size_range = c(0.009, 0.03),                 # Range of node sizes
          layout = "davidson-harel",                        # Primary layout algorithm
          initial_layout = "reingold-tilford",  
          edge_size_range = c(0.0005, 0.013))    # Initial layout for optimization
   
# output_file = "meta_coder_18s.pdf")               # Save plot directly to file

#### Metacoder 16s ####
#read in taxonomy and OTU data for 16S OTUs. Merge files
merged_table_16s <- left_join(OTU_table_16s_phy, taxonomy_16s_phy, by = "OTU_ID")

#Add in LUCA to connect bacteria and archea
merged_table_16s <- cbind(
  merged_table_16s[, 1:7],
  LUCA = "LUCA",
  merged_table_16s[, 8:ncol(merged_table_16s)]
)


#drop any rows with merged_table_16s$Domain == Unassigned
merged_table_16s <- merged_table_16s %>%
  filter(Domain != "Unassigned") %>%  # Remove rows where Domain is "Unassigned"
  filter(!is.na(Domain))               # Also remove rows with NA in Domain


# Add semicolon-separated taxonomy string by uniting taxonomy columns 6 through 11
merged_table_16s <- merged_table_16s %>%
  tidyr::unite("taxonomy", 8:11, sep = ";", remove = FALSE)

# Normalize read counts across three samples:
# Step 1: Calculate total reads per sample (columns 2 to 4)
sample_totals_16s <- colSums(merged_table_16s[, 2:7])

# Step 2: Find the smallest library size among samples
min_depth_16s <- min(sample_totals_16s)

# Step 3: Normalize each sample to the smallest depth by proportional scaling
normalized_counts16s <- sweep(merged_table_16s[, 2:7], 2, sample_totals_16s, FUN = "/") * min_depth_16s

# Step 4: Replace original counts with normalized counts in the data frame
merged_table_16s[, 2:7] <- normalized_counts16s

# Convert the merged data frame to metacoder taxmap format, specifying taxonomy column and separator
obj_16s <- parse_tax_data(merged_table_16s,
                          class_cols = "taxonomy",  # Column containing taxonomy string
                          class_sep = ";")          # Taxonomic rank separator

# Calculate taxon abundance across samples Bull, Mortensen, and West
obj_16s$data$tax_abund <- calc_taxon_abund(obj_16s, "tax_data",
                                           cols = c("Bull", "Mortensen", "West", "Horn", "Conant", "Greg"))

# Calculate total abundance for each taxon across the three samples
obj_16s$data$tax_abund$total_abund <- rowSums(obj_16s$data$tax_abund[, c("Bull", "Mortensen", "West", "Horn", "Conant", "Greg")])

#remove LUCA label
custom_labels <- taxon_names(obj_16s)
custom_labels[custom_labels == "LUCA"] <- NA  # or "" to blank it
heat_tree(obj_16s, node_label = custom_labels)


# Generate and save a heat tree plot visualizing taxon abundance
heat_tree_16s <- heat_tree(obj_16s,
          node_label = custom_labels,                         # Labels for nodes are taxon names
          node_size = obj_16s$data$tax_abund$total_abund,  # Node size scaled by total abundance
          node_color = obj_16s$data$tax_abund$total_abund, # Node color scaled by total abundance
          node_color_axis_label = "Total abundance",
          edge_size_axis_label = "Number of ASVs",
          edge_size = n_obs,                                # Edge size proportional to number of OTUs
          node_size_range = c(0.009, 0.03),                 # Range of node sizes
          initial_layout = "fr", layout = "da",
          edge_size_range = c(0.0005, 0.013))

#output_file = "meta_coder_16s.pdf")


# ----------------- Merge Meta Coder Plots -----------------

# Combine the three heat tree plots into a single layout
combined_meta_plot <- (
  heat_tree_ITS /
    heat_tree_16s /
    heat_tree_18s
) +
  plot_layout(ncol = 1, heights = c(5, 5, 5)) +  # Equal tallness
  plot_annotation(
    tag_levels = 'A',
    tag_suffix = ")",
    tag_prefix = "("
  ) &
  theme(
    plot.tag = element_text(size = 18, face = "bold")
  )

combined_meta_plot


#save as pdf
ggsave("combined_meta_plot.pdf",
       plot = combined_meta_plot,
       width = 10, height = 18, units = "in", device = "pdf", limitsize = FALSE)











#### Isolates ####
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
#ggsave("venn_plot.svg",
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




#### Venn Isolate vs Total OTUs (need to do) ####
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
  ylab("Isolate OTU Count") +
  xlab("") +
  labs(fill = "Phylum") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
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
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
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
  ylab("") +
  xlab("") +
  labs(fill = "Order") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
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
  scale_fill_manual(values = palette_family, guide = guide_legend(ncol = 2)) +
  theme_minimal() +
  ylab("Isolate OTU Count") +
  xlab("") +
  labs(fill = "Family") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
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
  scale_fill_manual(values = palette_genus, guide = guide_legend(ncol = 2)) +
  theme_minimal() +
  ylab("") +
  xlab("") +
  labs(fill = "Genus") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, color = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )


# ----- Combine All Taxonomic Plots -----
# Assemble all stacked bar plots vertically using patchwork
combined_stacked_bar_plot_isolates <-
  (Isolates_phylum | Isolates_class | Isolates_order) /
  (Isolates_family | Isolates_genus) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 20, face = "bold")
  )

combined_stacked_bar_plot_isolates

# Save full composite figure to PDF
ggsave("combined_stacked_bar_plot_isolates.pdf",
       plot = combined_stacked_bar_plot_isolates,
       width = 30, height = 20, units = "in", device = "pdf", limitsize = FALSE)

# ----------------- Percent Similarity -----------------

# 1. Create breaks and labels for similarity bins (74% to 100% by 1%)
breaks <- seq(74, 100, by = 1)
labels <- paste(head(breaks, -1), tail(breaks, -1), sep = "–")
descending_labels <- rev(labels)  # descending order for x-axis

# 2. Convert percent_identity_Unite to numeric (if needed)
isolates$percent_identity_Unite <- as.numeric(as.character(isolates$percent_identity_Unite))

# 3. Bin similarity values using cut()
isolates$similarity_bin <- cut(
  isolates$percent_identity_Unite,
  breaks = breaks,
  labels = labels,
  include.lowest = TRUE,
  right = FALSE
)

# 4. Set similarity_bin as a factor with descending levels for plotting order
isolates$similarity_bin <- factor(isolates$similarity_bin, levels = descending_labels)

# 5. Summarize counts per bin, ensure all bins present
bin_counts <- isolates %>%
  count(similarity_bin, name = "count") %>%
  complete(similarity_bin = descending_labels, fill = list(count = 0)) %>%
  mutate(
    similarity_bin = factor(similarity_bin, levels = descending_labels),
    bin_start = as.numeric(sub("–.*", "", as.character(similarity_bin))),
    group = ifelse(bin_start < 98, "Novel Candidate Taxa", "Known Taxa"),
    asterisk = ifelse(bin_start < 95 & count > 0, "*", NA)  # Flag bins for asterisk
  ) %>%
  arrange(similarity_bin)

# 6. Select every 3rd label for x-axis breaks to reduce clutter
x_breaks <- descending_labels[seq(1, length(descending_labels), by = 3)]

# 7. Create dummy data frame for asterisk legend entry
asterisk_legend <- data.frame(
  similarity_bin = factor(NA, levels = descending_labels),
  count = 0,
  group = "* Phylogenetic Inference"
)

#plot
isolates_similarity <- ggplot(bin_counts, aes(x = similarity_bin, y = count, fill = group)) +
  geom_col() +
  geom_text(aes(label = asterisk), vjust = -0.5, size = 5, fontface = "bold", na.rm = TRUE) +
  geom_point(data = asterisk_legend, aes(x = similarity_bin, y = count, fill = group),
             shape = 22, size = 5, color = "black", alpha = 0, show.legend = TRUE, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "Novel Candidate Taxa" = "red",
      "Known Taxa" = "steelblue",
      "* Phylogenetic Inference" = "white"
    ),
    breaks = c("Novel Candidate Taxa", "* Phylogenetic Inference"),
    name = NULL
  ) +
  scale_x_discrete(drop = FALSE, breaks = x_breaks) +
  labs(
    x = "Similarity to UNITE Reference (%)",
    y = "Isolate OTU Count"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # full box
    axis.line = element_blank(),  # Remove axis lines to avoid doubling
    axis.ticks.length.x = unit(0, "pt"),  # Remove x axis ticks to avoid confusion with border
    legend.position = c(0.99, 0.99),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha('white', 0.6), color = "NA", size = 0),
    legend.key.size = unit(1.2, "lines"),
    legend.text = element_text(face = "bold")
  ) +
  guides(
    fill = guide_legend(order = 1)
  )

print(isolates_similarity)



# ----- Combine All Isolate Plots -----
#read in jpeg version of photos of isolate images
isolate_photos <- readJPEG("manuscript/isolate_photo.jpg")     # Read JPEG image
isolate_photos <- rasterGrob(isolate_photos, interpolate = TRUE)  # Convert to grob

#read in jpeg version of isolates venn (minor edits in adobe illustrator)
isolate_amplicon_venn <- readJPEG("manuscript/updates_isoaltes_amplicon_venn.jpg")     # Read JPEG image
isolate_amplicon_venn <- rasterGrob(isolate_amplicon_venn, interpolate = TRUE)  # Convert to grob


# isolates in amplicon venn (These values were made in isolates_taxonomy_finalOTU_tables.R)
# Define sets with meaningful names
venn_data <- euler(c(
  "ITS Amplicon" = 742,        # only in ITS
  "Isolated" = 33,            # only in Isolated
  "ITS Amplicon&Isolated" = 34  # in both
))

# Plot with custom labels
isolates_venn <- plot(venn_data,
     fills = list(fill = c("skyblue", "salmon"), alpha = 0.6),
     labels = list(font = 2),
     edges = TRUE,
     quantities = TRUE)

#plot all together
# Wrap all grobs
isolate_photos <- wrap_elements(full = isolate_photos)
isolates_venn <- wrap_elements(full = isolates_venn)
isolate_amplicon_venn <- wrap_elements(full = isolate_amplicon_venn)


total_isolates_plot <- 
  (Isolates_phylum | Isolates_class | Isolates_order) /         # Row 1
  (Isolates_family | Isolates_genus) /                         # Row 2
  (isolates_venn | isolate_photos) /                           # Row 3: B and C side-by-side
  (isolates_similarity | isolate_amplicon_venn) +              # Row 4: D and E side-by-side
  plot_layout(
    nrow = 4,
    heights = c(1, 1, 0.9, 1.2)   # Last row 1.5 times taller
  ) +
  plot_annotation(
    tag_levels = list(c(
      "(A) 1.", "2.", "3.",        # Row 1 tags
      "4.", "5.",                  # Row 2 tags
      "(B)", "(C)",                # Row 3 tags
      "(D)", "(E)"                 # Row 4 tags
    )),
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(plot.tag = element_text(size = 30, face = "bold"))
  ) &
  theme(
    plot.tag = element_text(size = 30, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )


#save as pdf
ggsave("total_isolates_plot.pdf",
       plot = total_isolates_plot,
       width = 23, height = 36, units = "in", device = "pdf", limitsize = FALSE)





#### Abundance ####

# ----------------- Calculate ratios and their SEs -----------------
metadata_AS <- metadata_AS %>%
  rowwise() %>%
  mutate(
    # Vectors for triplicates as list-columns
    fungi_qpcr_vec    = list(c(qPCR_fungi_1, qPCR_fungi_2, qPCR_fungi_3)),
    bacteria_qpcr_vec = list(c(qPCR_bacteria_1, qPCR_bacteria_2, qPCR_bacteria_3)),
    
    fungi_mean   = mean(unlist(fungi_qpcr_vec), na.rm = TRUE),
    bacteria_mean= mean(unlist(bacteria_qpcr_vec), na.rm = TRUE),
    fungi_se     = sd(unlist(fungi_qpcr_vec), na.rm = TRUE) / sqrt(sum(!is.na(unlist(fungi_qpcr_vec)))),
    bacteria_se  = sd(unlist(bacteria_qpcr_vec), na.rm = TRUE) / sqrt(sum(!is.na(unlist(bacteria_qpcr_vec)))),
    
    # 1. Cell ratio by qPCR (mean of triplicates)
    Cell_ratio_qpcr = fungi_mean / bacteria_mean,
    # 2. Cell ratio by enumeration
    Cell_ratio_enum = Fungal_Abundance.cells.mL. / Total.Microbial.Abundance..cells.mL.,
    # 3. Standard error for the qPCR ratio (propagate error)
    ster_ratio_qpcr = abs(Cell_ratio_qpcr) * sqrt((fungi_se / fungi_mean)^2 + (bacteria_se / bacteria_mean)^2),
    
    # 4. Carbon ratio by qPCR, scaled by cell size
    carbon_ratio_qpcr = (fungi_mean * 6469) / (bacteria_mean * 10),
    # 5. SE for the carbon ratio (propagate error)
    ster_carbon_qpcr = abs(carbon_ratio_qpcr) * sqrt((fungi_se / fungi_mean)^2 + (bacteria_se / bacteria_mean)^2)
  ) %>%
  ungroup() %>%
  select(
    -fungi_qpcr_vec, -bacteria_qpcr_vec, -fungi_mean, -bacteria_mean, -fungi_se, -bacteria_se
  )

# ----------------- Prepare long data for plotting -----------------
plotdata_both <- metadata_AS %>%
  mutate(ster_ratio_enum = NA_real_) %>%
  select(Sample, Cell_ratio_enum, ster_ratio_enum, Cell_ratio_qpcr, ster_ratio_qpcr) %>%
  pivot_longer(
    cols = c(Cell_ratio_enum, Cell_ratio_qpcr, ster_ratio_enum, ster_ratio_qpcr),
    names_to = c("Metric", "Method"),
    names_pattern = "(Cell_ratio|ster_ratio)_(enum|qpcr)",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value,
    values_fn = mean          # ensure no duplicate rows per (Sample, Method)
  ) %>%
  mutate(
    Method = recode(Method, qpcr = "qPCR cell ratio", enum = "Enumeration cell ratio"),
    Cell_ratio = as.numeric(Cell_ratio),
    ster_ratio = as.numeric(ster_ratio)
  ) %>%
  filter(!is.na(Cell_ratio))

plotdata_both$Method <- factor(
  plotdata_both$Method,
  levels = c("qPCR cell ratio", "Enumeration cell ratio")
)
plotdata_both$Sample <- factor(plotdata_both$Sample, levels = site_order)

## Shared axis limits for both plots, so they're visually comparable


# ----------------- Enumeration plot (no SE bar, just dot) -----------------

enum_df   <- filter(plotdata_both, Method == "Enumeration cell ratio")
mean_enum <- mean(enum_df$Cell_ratio, na.rm = TRUE)
med_enum  <- median(enum_df$Cell_ratio, na.rm = TRUE)

enum_plot <- ggplot(enum_df, aes(x = Cell_ratio, y = Sample)) +
  geom_point(size = 5, color = "black") +
  geom_vline(xintercept = mean_enum, color = "red", linetype = "dashed", size = 1.3) +
  geom_vline(xintercept = med_enum, color = "red", linetype = "solid", size = 1.3) +
  labs(
    title = NULL,
    x = "Fungal: Total Microbial Cell (Enumeration)",
    y = "BNG Well"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.12, 0.05))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x     = element_text(size = 18, color = "black"),
    axis.title.y     = element_text(size = 18, color = "black"),
    axis.text.x      = element_text(size = 16, color = "black"),
    axis.text.y      = element_text(size = 16, color = "black"),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    legend.key.size  = grid::unit(1.2, "cm"),
    plot.margin      = margin(5, 5, 5, 5)
  )
#enum_plot

# > med_enum [1] 0.06986329 > mean_enum [1] 0.07040717
# ----------------- qPCR plot (with SE error bars) -----------------
# Filter only rows with a valid cell ratio to show dots for every site with data
qpcr_df <- plotdata_both %>%
  filter(Method == "qPCR cell ratio") %>%
  filter(!is.na(Cell_ratio))

# Calculate mean and median across all available ratios for the line overlays
mean_qpcr <- mean(qpcr_df$Cell_ratio, na.rm = TRUE)
med_qpcr  <- median(qpcr_df$Cell_ratio, na.rm = TRUE)

# Determine a suitable axis limit to cover all ratios and their error bars
max_x <- max(qpcr_df$Cell_ratio + ifelse(is.na(qpcr_df$ster_ratio), 0, qpcr_df$ster_ratio), na.rm = TRUE) * 1.1

# Plot: point for each ratio, error bar only if SE is present, mean & median lines
qpcr_plot <- ggplot(qpcr_df, aes(x = Cell_ratio, y = Sample)) +
  geom_point(size = 5, color = "black") +
  geom_errorbarh(
    data = subset(qpcr_df, !is.na(ster_ratio)),
    aes(xmin = Cell_ratio - ster_ratio, xmax = Cell_ratio + ster_ratio),
    height = 0.2, color = "black", linewidth = 1.1
  ) +
  geom_vline(xintercept = mean_qpcr, color = "red", linetype = "dashed", size = 1.3) +
  geom_vline(xintercept = med_qpcr, color = "red", linetype = "solid", size = 1.3) +
  labs(
    title = NULL,
    x = "Fungal: Bacterial Cell (qPCR)",
    y = "BNG Well"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.12, 0.05))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x     = element_text(size = 18, color = "black"),
    axis.title.y     = element_text(size = 18, color = "black"),
    axis.text.x      = element_text(size = 16, color = "black"),
    axis.text.y      = element_text(size = 16, color = "black"),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    legend.key.size  = grid::unit(1.2, "cm"),
    plot.margin      = margin(5, 5, 5, 5)
  )

qpcr_plot

# > med_qpcr [1] 0.0003274277 > mean_qpcr [1] 0.0004708453
# ----------------- carbon ratio -----------------
plotdata_C_qpcr <- metadata_AS %>%
  select(Sample, carbon_ratio_qpcr, ster_carbon_qpcr) %>%
  filter(!is.na(carbon_ratio_qpcr)) %>%
  mutate(Sample = factor(Sample, levels = site_order))

mean_c_qpcr   <- mean(plotdata_C_qpcr$carbon_ratio_qpcr, na.rm = TRUE)
median_c_qpcr <- median(plotdata_C_qpcr$carbon_ratio_qpcr, na.rm = TRUE)

max_x <- max(plotdata_C_qpcr$carbon_ratio_qpcr + ifelse(is.na(plotdata_C_qpcr$ster_carbon_qpcr), 0, plotdata_C_qpcr$ster_carbon_qpcr), na.rm = TRUE) * 1.1

carbon_qpcr_plot <- ggplot(plotdata_C_qpcr, aes(x = carbon_ratio_qpcr, y = Sample)) +
  geom_point(size = 5, color = "black") +
  geom_errorbarh(
    aes(xmin = carbon_ratio_qpcr - ster_carbon_qpcr,
        xmax = carbon_ratio_qpcr + ster_carbon_qpcr),
    height = 0.2, color = "black", linewidth = 1.1,
    na.rm = TRUE       # so only plots error bars if SE is present
  ) +
  geom_vline(xintercept = mean_c_qpcr, color = "red",   linetype = "dashed", size = 1.3) +
  geom_vline(xintercept = median_c_qpcr, color = "red", linetype = "solid", size = 1.3) +
  labs(
    x = "Fungal: Bacterial Carbon",
    y = "BNG Well"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.12, 0.05)),
    limits = c(0, max_x)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x     = element_text(size = 18, color = "black"),
    axis.title.y     = element_text(size = 18, color = "black"),
    axis.text.x      = element_text(size = 16, color = "black"),
    axis.text.y      = element_text(size = 16, color = "black"),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    legend.key.size  = grid::unit(1.2, "cm"),
    plot.margin      = margin(5, 5, 5, 5)
  )

carbon_qpcr_plot


# > median_c_qpcr [1] 0.211813 > mean_c_qpcr [1] 0.3045898

# ----------------- Depth vs fungal vs total abundance -----------------

# Prepare dataframe with clear names
enum_abundane_df <- metadata_AS %>%
  select(Depth = Depth..m., 
         Fungal = Fungal_Abundance.cells.mL., 
         Total = Total.Microbial.Abundance..cells.mL.)

# Scaling factor
scale_factor <- max(enum_abundane_df$Fungal, na.rm = TRUE) / max(enum_abundane_df$Total, na.rm = TRUE)

depth_fungi_total_dual <- ggplot(enum_abundane_df, aes(x = Depth)) +
  # Fungal abundance (primary y-axis, green)
  geom_line(aes(y = Fungal), color = "forestgreen", size = 1.2, linetype = "dashed") +
  geom_point(aes(y = Fungal), color = "forestgreen", size = 3) +
  
  # Total abundance (secondary y-axis, black, scaled)
  geom_line(aes(y = Total * scale_factor), color = "black", size = 1.2, linetype = "dashed") +
  geom_point(aes(y = Total * scale_factor), color = "black", size = 3) +
  
  scale_y_continuous(
    name = "Fungal Abundance (cells/mL)",
    labels = scales::scientific,
    sec.axis = sec_axis(~./scale_factor, name = "Total Abundance (cells/mL)", labels = scales::scientific)
  ) +
  scale_x_continuous(name = "Depth Below Surface (m)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x     = element_text(size = 18, color = "black"),
    axis.title.y.left  = element_text(color = "forestgreen", size = 18),
    axis.title.y.right = element_text(color = "black", size = 18),
    axis.text.x      = element_text(size = 16, color = "black"),
    axis.text.y.left  = element_text(size = 16, color = "forestgreen"),
    axis.text.y.right = element_text(size = 16, color = "black"),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(5, 5, 5, 5)
  )

depth_fungi_total_dual



# ----------------- NPDOC vs fungal vs total abundance -----------------
npdoc_plot <- ggplot() +
  # Fungal (primary y-axis, green, dashed)
  geom_point(data = fungi_enum, aes(x = NPDOC, y = Abundance),
             color = "forestgreen", size = 3) +
  geom_smooth(
    data = fungi_enum,
    aes(x = NPDOC, y = Abundance),
    method = "lm",
    color = "forestgreen",
    linetype = "dashed",
    size = 1.2,
    se = FALSE
  ) +
  # Total (secondary y-axis, black, solid)
  geom_point(data = total_enum, aes(x = NPDOC, y = Abundance * scale_factor),
             color = "black", size = 3) +
  geom_smooth(
    data = total_enum,
    aes(x = NPDOC, y = Abundance * scale_factor),
    method = "lm",
    color = "black",
    linetype = "solid",
    size = 1.2,
    se = FALSE
  ) +
  scale_y_continuous(
    name = "Fungal Enumeration (cells/mL)",
    labels = scales::scientific,
    sec.axis = sec_axis(~./scale_factor, name = "Total Enumeration (cells/mL)", labels = scales::scientific)
  ) +
  scale_x_continuous(name = "NPDOC (mg/L)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line        = element_blank(),
    legend.position  = "none",
    axis.title.y.left  = element_text(color = "forestgreen", size = 18),
    axis.title.y.right = element_text(color = "black", size = 18),
    axis.title.x      = element_text(size = 18, color = "black"),
    axis.text.x       = element_text(size = 16, color = "black"),
    axis.text.y.left  = element_text(size = 16, color = "forestgreen"),
    axis.text.y.right = element_text(size = 16, color = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_blank(),
    plot.margin       = margin(5, 5, 5, 5)
  )

npdoc_plot

ggsave("npdoc_plot_lm.pdf",
       plot = npdoc_plot,
       width = 16, height = 8, units = "in", device = "pdf", limitsize = FALSE)

# Fungal linear model
fungal_lm <- lm(Abundance ~ NPDOC, data = fungi_enum)
# Get R²
fungal_r2 <- summary(fungal_lm)$r.squared
cat("Fungal Enumeration R²:", fungal_r2, "\n")


# Total linear model
total_lm <- lm(Abundance ~ NPDOC, data = total_enum)
# Get R²
total_r2 <- summary(total_lm)$r.squared
cat("Total Enumeration R²:", total_r2, "\n")


# Fungal abundance vs NPDOC
spearman_fungal <- cor.test(metadata_AS$Fungal_Abundance.cells.mL., metadata_AS$NPDOC, method = "spearman", use = "complete.obs")

# Total abundance vs NPDOC
spearman_total <- cor.test(metadata_AS$Total.Microbial.Abundance..cells.mL., metadata_AS$NPDOC, method = "spearman", use = "complete.obs")

# View results
spearman_fungal #p-value = 0.6583,rho 0.2571429 
spearman_total #p-value = 0.05833,rho 0.8285714 

spearman_fungi_total <- cor.test(
  metadata_AS$Fungal_Abundance.cells.mL.,
  metadata_AS$Total.Microbial.Abundance..cells.mL.,
  method = "spearman",
  use = "complete.obs"
)

# View results
spearman_fungi_total


# ----------------- combine abundance plots -----------------
#Combine plots
combined_plot_abundance <- (depth_fungi_total_dual | npdoc_plot) / 
  (enum_plot | qpcr_plot | plot_spacer()) /   # Big empty right space on C row!
  (carbon_qpcr_plot | plot_spacer()) +        # D only on bottom left
  plot_layout(
    nrow = 4,
    heights = c(1, 1, 1, 1)   # Adjust as you like
  ) +
  plot_annotation(
    tag_levels = list(c("(A)", "(B)", "(C) 1.", "(C) 2.", "(D)")),
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(
      plot.tag = element_text(size = 30, face = "bold")
    )
  ) &
  theme(
    plot.tag = element_text(size = 30, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("combined_plot_abundance.pdf",
       plot = combined_plot_abundance,
       width = 20, height = 16, units = "in", device = "pdf", limitsize = FALSE)




#### Supplemental ####
#### Supplemental ITS Richness####
# --- Prepare and merge data ---

# Get read counts per sample from ITS OTU table
ITS_read_counts <- data.frame(
  Sample = colnames(OTU_table_ITS[, -1]),          # assumes first column is OTU_ID
  Read_Count = colSums(OTU_table_ITS[, -1])
)

# Merge read counts with fungal richness metadata
ITS_data <- merge(
  ITS_read_counts,
  metadata_AS[, c("Sample", "Chao1_fungi")],
  by = "Sample"
)

# --- Calculate Spearman's rho and p-value ---
spearman_1 <- cor.test(metadata_AS$Chao1_16S_total, metadata_AS$Chao1_fungi, method = "spearman")
rho_1 <- spearman_1$estimate
pval_1 <- spearman_1$p.value

spearman_2 <- cor.test(ITS_data$Read_Count, ITS_data$Chao1_fungi, method = "spearman")
rho_2 <- spearman_2$estimate
pval_2 <- spearman_2$p.value

# --- Annotate plots ---
ITS_vs_16s_rich <- ggplot(metadata_AS, aes(x = Chao1_16S_total, y = Chao1_fungi)) +
  geom_point(size = 5, color = "purple") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "Total 16S (Bacterial + Archaeal) ASV Richness",
    y = "Fungal OTU (ITS) Richness",
    subtitle = paste0("Spearman's \u03C1 = ", round(rho_1, 3), ", p = ", signif(pval_1, 3))
  ) +
  theme_minimal(base_size = 17)

ITSrich_reads <- ggplot(ITS_data, aes(x = Read_Count, y = Chao1_fungi)) +
  geom_point(size = 5, color = "darkorchid") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(
    x = "ITS Read Count per Sample",
    y = "Fungal OTU (ITS) Richness",
    subtitle = paste0("Spearman's \u03C1 = ", round(rho_2, 3), ", p = ", signif(pval_2, 3))
  ) +
  theme_minimal(base_size = 17)

# --- Stack plots and add tags ---
combined_plots_supp <- ITS_vs_16s_rich / ITSrich_reads +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(plot.tag = element_text(size = 18, face = "bold")) # optional formatting for subplot tags

# Display
combined_plots_supp

# --- Save combined figure ---
ggsave("combined_plots_supp_ITS_richness.pdf",
       plot = combined_plots_supp,
       width = 10, height = 12, units = "in", device = "pdf", limitsize = FALSE)

# --- Print stats to console ---
cat("\n(A) ITS vs 16S richness: Spearman's rho =", round(rho_1, 3), ", p =", signif(pval_1, 3), "\n")
cat("(B) ITS read count vs fungal richness: Spearman's rho =", round(rho_2, 3), ", p =", signif(pval_2, 3), "\n")







#### Supplemental Total vs Fungal Enumeration ####

total_vs_fun_enum <- ggplot(metadata_AS, aes(x = Total.Microbial.Abundance..cells.mL., y = Fungal_Abundance.cells.mL.)) +
  geom_point(size = 3, color = "darkred") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  labs(
    x = "Total Enumeration (cells/ml or log10 scale)",
    y = "Fungal Enumeration (cells/ml or log10 scale)",
    title = "Total vs. Fungal Enumeration"
  ) +
  theme_minimal(base_size = 15)

total_vs_fun_enum
ggsave2("total_vs_fun_enum.pdf",
        plot = total_vs_fun_enum,
        width = 8, height = 6, units = "in", device = "pdf", limitsize = FALSE)



ggplot(metadata_AS, aes(x = Cell_ratio_qpcr, y = Depth..m.)) +
  geom_point(size = 3, color = "navy") +
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
  scale_y_reverse(name = "Depth Below Surface (m)") +  # Shallow at top, deep at bottom
  labs(
    x = "Cell Ratio (qPCR; e.g. Fungi/Bacteria)",
    title = "Cell Ratio (qPCR) vs Depth"
  ) +
  theme_minimal(base_size = 16)

ggsave2("cell_ratio_vs_depth.pdf",
        plot = last_plot(),
        width = 8, height = 6, units = "in", device = "pdf", limitsize = FALSE)




res1 <- cor.test(metadata_AS$Chao1_16S_total, metadata_AS$Chao1_fungi, method = "spearman")
cat("Total 16S richness vs Fungal richness:\n")
cat(sprintf("  Spearman rho = %.3f, p-value = %.3g\n", res1$estimate, res1$p.value))



res2 <- cor.test(metadata_AS$Cell_ratio_qpcr, metadata_AS$Depth..m., method = "spearman")
cat("ratio vs depth:\n")
cat(sprintf("  Spearman rho = %.3f, p-value = %.3g\n", res3$estimate, res3$p.value))



res3 <- cor.test(metadata_AS$Total.Microbial.Abundance..cells.mL., metadata_AS$Fungal_Abundance.cells.mL., method = "spearman")
cat("Bacterial enumeration vs Fungal enumeration:\n")
cat(sprintf("  Spearman rho = %.3f, p-value = %.3g\n", res3$estimate, res3$p.value))






# Load your data
metadata <- read.csv("metadata_AS.csv", stringsAsFactors = FALSE)

# Select numeric columns
numdata <- metadata_AS[sapply(metadata_AS, is.numeric)]

# Get all unique variable pairs
pair_mat <- combn(names(numdata), 2, simplify = TRUE)

# Create an empty results data frame
res_list <- list()

for(i in seq_len(ncol(pair_mat))) {
  var1 <- pair_mat[1, i]
  var2 <- pair_mat[2, i]
  x <- numdata[[var1]]
  y <- numdata[[var2]]
  keep <- complete.cases(x, y)
  N <- sum(keep)
  if (N >= 3) {
    ct <- cor.test(x[keep], y[keep], method = "spearman")
    res_list[[i]] <- data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Spearman_rho = unname(ct$estimate),
      p_value = ct$p.value,
      N = N
    )
  } else {
    res_list[[i]] <- data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Spearman_rho = NA,
      p_value = NA,
      N = N
    )
  }
}

# Combine into a single data frame
spearman_table <- do.call(rbind, res_list)

# Print as nice table
View(spearman_table)

# Optional: view as spreadsheet-style popup
# View(spearman_table)

# Optional: save to CSV
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













