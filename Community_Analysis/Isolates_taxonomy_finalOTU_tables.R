#This script reads in the results from blasting the ITS2 amplicons against a database of sanger isolates. For isolates that match an amplicon,
#The amplicon taxonomy is used for the isolate taxonomy as well. For isolate sequences not matched to an amplicon sequence, they are binned at 0.98 vsearch 
# and blasted against the Unite database to assign taxonomy (in seperate slurm script). These files unifies all of the isolate sequences and their taxonomy into a single table.


#### Load Packages ####
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(ggplot2)
library(tibble)
library(compositions)
library(ALDEx2)
library(Biostrings)
library(VennDiagram)
library(grid)
library(metacoder)
library(ggrepel)
library(systemfonts)
library(readr)
library(eulerr)
setwd("/Users/quinnmoon/Downloads/Antrim_Microbiome/")


#### Matched Isolated ####
#Read in the filtered blastn results (amplicon against sanger isolate database)
column_names <- c(
  "OTU_ID", "Isolate_ID", "percent_identity", "alignment_length", "mismatches", 
  "gap_opens",
  "e_value", "bit_score"
)

filtered_blast_results <- read_tsv("filtered_hits_150bp.tsv", col_names = column_names)

#check to see if any isolates match multiple OTUs
duplicated_otus <- unique(filtered_blast_results$Isolate_ID[duplicated(filtered_blast_results$Isolate_ID)])
print(duplicated_otus)

#Results: "AS_Bull_09_ITS1_500bp"      "AS_Bull_04-ITS1_300bp"      "AS_Mortensen_08_ITS1_516bp".
#Allow only 1 best match. Manually adjust .csv file.

#Read back in trimmed file. 
filtered_blast_results_trimmed <- read_csv("filtered_hits_150bp_trimmed.csv", col_names = column_names)

#Add well information as a column
filtered_blast_results_trimmed <- filtered_blast_results_trimmed %>%
  mutate(well_name = str_split(Isolate_ID, "_", simplify = TRUE)[, 2])

#drop columns 3-8
filtered_blast_results_trimmed <- filtered_blast_results_trimmed %>%
  select(OTU_ID, Isolate_ID, well_name)

#Now lets collapse the Isolate_IDs into a single string per OTU_ID and create a presence/absence matrix for wells.
# Step 1: Collapse all Isolate_IDs into a single string per OTU_ID
collapsed <- filtered_blast_results_trimmed %>%
  group_by(OTU_ID) %>%
  summarise(Isolates = paste(Isolate_ID, collapse = ", "), .groups = "drop")

# Step 2: Create 1/0 presence matrix for wells
presence_matrix <- filtered_blast_results_trimmed %>%
  distinct(OTU_ID, well_name) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = well_name,
    values_from = present,
    values_fill = 0
  )

# Step 3: Join collapsed Isolate_IDs and presence matrix
final_blast <- collapsed %>%
  dplyr::left_join(presence_matrix, by = "OTU_ID")

#Results: 42 isolate OTUs at 98%. 8 of those OTUs have a representative that was cultured in the controls. 
#Lets read in the taxonomy and add it to the final_blast table.
taxonomy <- read.csv("ITS2/final_taxonomy_ITS.csv", stringsAsFactors = FALSE) 

taxonomy_subset <- dplyr::select(taxonomy, OTU_ID, percent_identity, kingdom, phylum, class, order, family, genus)

taxonomy_subset <- taxonomy_subset %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | . == "", "Unclassified", .)))

# Join the taxonomy information with the final_blast results
blast_with_taxonomy <- final_blast %>%
  left_join(taxonomy_subset, by = "OTU_ID")

#update column name to be clear where the percent identity comes from
blast_with_taxonomy <- blast_with_taxonomy %>%
  dplyr::rename(percent_identity_Unite = percent_identity_unite)

#update col names to be more clear
blast_with_taxonomy <- blast_with_taxonomy %>%
  dplyr::rename(
    Horn_Isolate = Horn,
    Mortensen_Isolate = Mortensen,
    Control_Isolate = Control,
    West_Isolate = West,
    Shooks_Isolate = Shooks,
    Bull_Isolate = Bull
  )

#Now lets add on the amplicon OTU table to see where the isolates appear in the amplicon data
otu_table <- read_tsv("ITS2/ITS2_OTU_table.tsv", skip = 1)

colnames(otu_table)[1] <- "OTU_ID"

#Based on NMDS, we can pool technical replicates together. 
replicate_groups <- list(
  Bull = c("Bull-A", "Bull-B", "Bull-C"),
  Conant = c("Conant-A", "Conant-B", "Conant-C"),
  Greg = c("Greg-A", "Greg-B", "Greg-C"),
  Horn = c("Horn-Elliott"),
  Mortensen = c("Mortensen-A", "Mortensen-B", "Mortensen-C"),
  West = c("West-A", "West-B", "West-C"),
  Control = c("Sand-Control")
)

# Start with just the OTU_ID column
combined_df <- otu_table %>%
  dplyr::select("OTU_ID")

# Loop over each group in replicate_groups and sum their columns
for (group_name in names(replicate_groups)) {
  replicate_cols <- replicate_groups[[group_name]]
  
  # Add the new combined column by summing across the replicates
  combined_df[[group_name]] <- otu_table %>%
    dplyr::select(all_of(replicate_cols)) %>%
    rowSums(na.rm = TRUE)
}


# Step 1: Calculate the rank abundance for each OTU in each well
rank_abundance_matched <- combined_df %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "Well",
    values_to = "Abundance"
  ) %>%
  group_by(Well) %>%
  mutate(
    total_OTUs = sum(Abundance > 0),
    Rank = if_else(
      Abundance > 0,
      rank(-Abundance, ties.method = "min"),
      NA_real_
    ),
    RankString = if_else(
      !is.na(Rank),
      paste0(Rank, "/", total_OTUs),
      "0"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(OTU_ID, Well, RankString) %>%
  pivot_wider(
    names_from = Well,
    values_from = RankString,
    names_glue = "{Well}_Amplicon",
    values_fill = "0"
  )

#Merge ranks to final_df
final_df_with_ranks <- blast_with_taxonomy %>%
  left_join(rank_abundance_matched, by = "OTU_ID")

write_csv(final_df_with_ranks, "matched_isolates.csv")


#### Unmatched Isolates ####
#Lets bin the remaining isolate sanger sequences into OTUs at 98% identity and assign taxonomy.
#compare original fasta and pull out isolates that matches amplicon data. 
fasta_all <- readDNAStringSet("Antrim6-24_Isolates.fasta")

identified_isolates <- unique(filtered_blast_results$Isolate_ID)

unidentified_fasta <- fasta_all[!names(fasta_all) %in% identified_isolates]

writeXStringSet(unidentified_fasta, filepath = "unmatched_isolates.fasta")
#Results: 56 isolate sequences do not match the amplicon data. Lets bin these at 98% identity and assign taxonomy using custom Unite 10 best hits. 
#Sequences binned using vsearch and taxonomy with blastn of Unite
column_names <- c(
  "OTU_ID", "subject_id", "percent_identity", "alignment_length", "mismatches", 
  "gap_opens", "query_start", "query_end", "subject_start", "subject_end", 
  "e_value", "bit_score", "taxonomic_classification"
)

# Read the BLAST result
blast_results_unmatched <- read_tsv("blastn_unidentified.tsv", col_names = column_names)

#Parse out taxonomic information
#Parse out taxonomic classification into seperate coloumns
blast_results_unmatched <- blast_results_unmatched %>%
  # Extract only the taxonomy portion starting at k__
  mutate(taxonomy_string = str_extract(taxonomic_classification, "k__.*")) %>%
  
  # Separate into columns by taxonomic rank
  separate(taxonomy_string,
           into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right", remove = FALSE) %>%
  
  # Remove k__, p__, etc. prefixes
  mutate(across(kingdom:species, ~str_remove(., "^[a-z]__")))

#Write the file out for manual inspection as used before to assign taxonomy to ITS2 amplicons
write.csv(blast_results_unmatched, "blastn_unidentified_taxonomy.csv", row.names = FALSE)
#Results: 37 OTUs from unidentified isolates.
concensus_unmatched <- read.csv("blastn_unidentified_taxonomy_concensus.csv")

#fill in empty strings with NA and replace with "Unclassified"
concensus_unmatched[concensus_unmatched == ""] <- NA  # Convert empty strings to NA (optional)
concensus_unmatched[is.na(concensus_unmatched)] <- "Unclassified"

#Lets add percent match to Unite database
min_alignment_length <- 175

# Apply logic using the variable
highest_identity_per_otu <- blast_results_unmatched %>%
  group_by(OTU_ID) %>%
  summarise(
    percent_identity = if (any(alignment_length > min_alignment_length)) {
      max(percent_identity[alignment_length > min_alignment_length])
    } else {
      mean(percent_identity)
    },
    .groups = "drop"
  )

taxonomy_unmatched <- concensus_unmatched %>%
  left_join(highest_identity_per_otu, by = "OTU_ID")

#Let us check whoch isolates belong to each isolate OTU
unmatched_clusters <- unmatched_clusters %>%
  mutate(
    all_isolates = if_else(
      is.na(Member_Sequences) | Member_Sequences == "",
      Representative_Sequence,
      paste(Representative_Sequence, Member_Sequences, sep = ", ")
    )
  )

#drop columns 2 and 3
unmatched_clusters <- unmatched_clusters[, -c(2, 3)]

#lets merge the taxonomy with the unmatched clusters
taxonomy_unmatched <- taxonomy_unmatched %>%
  left_join(unmatched_clusters, by = "OTU_ID")

# Step 1: Split the comma-separated sequences into rows
unmatched_long <- unmatched_clusters %>%
  separate_rows(all_isolates, sep = ",\\s*")


# Step 2: Extract the second element after splitting by "_"
unmatched_long <- unmatched_long %>%
  mutate(well_name = str_split(all_isolates, "_", simplify = TRUE)[, 2])

# Step 3: Collapse all isolate IDs per OTU
collapsed_unmatched <- unmatched_long %>%
  group_by(OTU_ID) %>%
  summarise(Isolates = paste(all_isolates, collapse = ", "), .groups = "drop")

# Step 4: Create presence/absence matrix by well
presence_matrix_unmatched <- unmatched_long %>%
  distinct(OTU_ID, well_name) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = well_name,
    values_from = present,
    values_fill = 0
  )

# Step 5: Merge collapsed isolate names and presence matrix
final_unmatched <- collapsed_unmatched %>%
  left_join(presence_matrix_unmatched, by = "OTU_ID")


final_unmatched <- final_unmatched %>%
  dplyr::rename(
    Horn_Isolate = Horn,
    Mortensen_Isolate = Mortensen,
    Control_Isolate = Control,
    West_Isolate = West,
    Shooks_Isolate = Shooks,
    Bull_Isolate = Bull
  )

#Results: Additional 37 isolate OTUs at 0.98 that did not match the amplicon data. 4 of these OTUs appear in controls and can be discarded later. 

#merge taxonomy_unmatched to final_unmatched
unmatched_isolates_pa_taxonomy <- final_unmatched %>%
  left_join(taxonomy_unmatched, by = "OTU_ID")

#clean up final df and write to .csv
unmatched_isolates_pa_taxonomy <- unmatched_isolates_pa_taxonomy %>%
  dplyr::select(-all_isolates)

unmatched_isolates_pa_taxonomy <- unmatched_isolates_pa_taxonomy %>%
  dplyr::rename(
    percent_identity_Unite = percent_identity,
  )

#write to .csv
write_csv(unmatched_isolates_pa_taxonomy, "unmatched_isolates.csv")

#### Venn Diagram Isolates ####
# Filter OTUs: exclude those appearing in control isolates
otus_filtered <- all_isolates %>%
  filter(appears_in_control_isolates != "Yes") %>%
  pull(OTU_ID) %>%
  unique()

# Cultured OTUs (filtered)
cultured_otus <- otus_filtered

# Amplicon OTUs from combined_df
amplicon_otus <- combined_df$OTU_ID %>%
  unique()

# Create counts
only_cultured <- setdiff(cultured_otus, amplicon_otus)
only_amplicon <- setdiff(amplicon_otus, cultured_otus)
both <- intersect(cultured_otus, amplicon_otus)

# Prepare data
venn_data <- euler(c(
  Cultured = length(only_cultured),
  "ITS2 Amplicon" = length(only_amplicon),
  "Cultured&ITS2 Amplicon" = length(both)
))

# Save to JPG
jpeg("venn_cultured_amplicon.jpg", width = 2500, height = 1700, res = 300)

# Create plot
plot(
  venn_data,
  fills = list(fill = c("skyblue", "salmon"), alpha = 0.5),
  edges = TRUE,
  labels = list(font = 2),
  quantities = list(type = "counts", font = 1, col = "black"),
  main = "Proportional Overlap Between Cultured and ITS2 Amplicon OTUs"
)

# Close the file
dev.off()



#### NMDS Grouping ITS####
#Read in OTU table and 
ITS_otu_table <- read_tsv("ITS2/ITS2_OTU_table.tsv", skip = 1) %>%
  dplyr::rename(OTU_ID = `#OTU ID`)

# Get OTUs that appear in controls
control_otus <- all_isolates %>%
  filter(appears_in_control_isolates == "Yes") %>%
  pull(OTU_ID) %>%
  unique()

#filter out control OTUs from isolates
otu_table_ordination_ITS <- ITS_otu_table %>%
  filter(!OTU_ID %in% control_otus)

#check reads per sample
colSums(otu_table_ordination_ITS[, -1])

#filter out samples with less than 230 reads
otu_table_ordination_ITS <- otu_table_ordination_ITS %>%
  {
    sample_cols <- names(.)[-1]
    keep_cols <- sample_cols[colSums(.[, sample_cols], na.rm = TRUE) > 300]
    keep_cols <- union(keep_cols, "Sand-Control")  # Always include Sand-Control
    dplyr::select(., OTU_ID, all_of(keep_cols))
  }

#Samples removed from downstream analysis=Greg-B,Greg-C,Horn-Concentrate,Mortensen-Elliott West-B
#13 samples, across 6 sites included. Also have surface control sample. 35887 total reads across all samples.

#drop any OTUs that no longer occur in table
otu_nmds_data_filtered_ITS <- otu_table_ordination_ITS %>%
  dplyr::filter(rowSums(dplyr::select(., -OTU_ID)) > 0)

# 2. Convert to relative abundance (proportions per sample)
# 1. Extract OTU_IDs
otu_ids_ITS <- otu_nmds_data_filtered_ITS[[1]]

# 2. Convert counts to relative abundance (excluding OTU_ID)
rel_abund_matrix_ITS <- sweep(
  otu_nmds_data_filtered_ITS[ , -1],           # exclude OTU_ID
  2,
  colSums(otu_nmds_data_filtered_ITS[ , -1]),
  FUN = "/"
)

# 3. Combine back with OTU_IDs
ITS_rel_abund <- cbind(OTU_ID = otu_ids_ITS, rel_abund_matrix_ITS)

otu_nmds_matrix_ITS <- ITS_rel_abund %>%
  column_to_rownames(var = "OTU_ID")
#transpose matrix
otu_nmds_matrix_ITS_t <- t(otu_nmds_matrix_ITS)
# Run NMDS on the transposed matrix
nmds_result_ITS <- metaMDS(otu_nmds_matrix_ITS_t, distance = "bray", k = 2, trymax = 100)


# 2. Extract NMDS points and add grouping info
nmds_points_ITS <- as.data.frame(nmds_result_ITS$points)
nmds_points_ITS$Sample <- rownames(nmds_points_ITS)

# Extract group name (e.g., "Bull" from "Bull-A")
nmds_points_ITS$Group <- sub("-[A-Za-z]+$", "", nmds_points_ITS$Sample)

# 3. Function to get convex hull points for each group
find_hull_ITS <- function(df) df[chull(df$MDS1, df$MDS2), ]

hulls <- nmds_points_ITS %>%
  group_by(Group) %>%
  do(find_hull_ITS(.))
# Identify groups with exactly two samples
groups_two_points <- nmds_points_ITS %>%
  group_by(Group) %>%
  filter(n() == 2) %>%
  ungroup()

# Your plot with polygons, points, labels, and lines connecting two-point groups
nmds_points_ITS <- nmds_points_ITS %>%
  mutate(Sample = case_when(
    Sample == "Horn-Elliott" ~ "Horn-E",
    Sample == "Sand-Control" ~ "Surface Control",
    TRUE ~ Sample
  ))



# Rename "Sand" to "Surface Control" in all relevant data frames
nmds_points_ITS <- nmds_points_ITS %>%
  mutate(Group = ifelse(Group == "Sand", "Surface Control", Group))

hulls <- hulls %>%
  mutate(Group = ifelse(Group == "Sand", "Surface Control", Group))

groups_two_points <- groups_two_points %>%
  mutate(Group = ifelse(Group == "Sand", "Surface Control", Group))

# Then create the plot with updated group names
nmds_rel_ITS <- ggplot(nmds_points_ITS, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = Group),
               alpha = 0.2, color = NA) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample),
                  size = 4, max.overlaps = 20, box.padding = 0.4) +
  geom_line(data = groups_two_points,
            aes(x = MDS1, y = MDS2, group = Group, color = Group),
            size = 2, alpha = 0.4, show.legend = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) +
  labs(color = "BNG Well", fill = "BNG Well")


# Assuming your ggplot object is the last plot you made
ggsave("nmds_replicates_rel_ITS.svg", plot = nmds_rel_ITS, width = 12, height = 10, units = "in")
ggsave("nmds_replicates_rel_ITS.jpg", plot = nmds_rel_ITS, width = 12, height = 10, units = "in", dpi = 2000)

#Results: Strong Technical replicate grouping using NMDS for ITS with relative abundances.








#### NMDS Grouping 18s ####
#Read in OTU table and 
otu_table_18s <- read_tsv("18s/18s_OTU_table.tsv") 

#Fix column names
# Step 1: Remove current column names (make them empty or NA)
colnames(otu_table_18s) <- rep(NA, ncol(otu_table_18s))  # or colnames(df) <- NULL

# Step 2: Set first row as column names
colnames(otu_table_18s) <- as.character(otu_table_18s[1, ])

# Step 3: Remove the first row (now used as column names)
otu_table_18s <- otu_table_18s[-1, ]
otu_table_18s[, -1] <- sapply(otu_table_18s[, -1], as.numeric)
# Step 4: Rename the first column to "OTU_ID"
colnames(otu_table_18s)[1] <- "OTU_ID"


#check reads per sample
colSums(otu_table_18s[, -1])

#Lets filter out land plants (division Streptophyta) and Craniata and Arthropoda animal classes from the 18s data.
#load taxonomy table
taxonomy_18s <- read.csv("18s/combined_taxonomy_18s.csv", stringsAsFactors = FALSE)

# Step 1: Filter OTUs from the taxonomy table that are NOT contamination
otus_to_keep_18s <- taxonomy_18s %>%
  filter(Division != "Streptophyta") %>%
  filter(!(Class %in% c("Craniata", "Arthropoda"))) %>%
  pull(OTU_ID)

# Step 2: Subset the OTU table to keep only those OTUs
cleaned_18s <- otu_table_18s_filtered %>%
  filter(OTU_ID %in% otus_to_keep_18s)

#filter out samples with less than 230 reads
cleaned_18s <- cleaned_18s %>%
  {
    keep_cols <- names(.)[-1][colSums(.[, -1], na.rm = TRUE) > 320]
    dplyr:: select(., OTU_ID, all_of(keep_cols))
  }

#drop any OTUs that no longer occur in table
cleaned_18s <- cleaned_18s %>%
  dplyr::filter(rowSums(dplyr::select(., -OTU_ID)) > 0)

colSums(cleaned_18s[, -1]) # 393 OTUs, 12 samples

#Result:This removed about 75% of the reads. A pinus OTU was most common OTU. 
#But it removed very little of the OTUs. 342 OTUs remain of the original 393. 
#12 non control samples remain. They cover Wells: Bull, Conant, Mortensen, West


# 2. Convert to relative abundance (proportions per sample)
# 1. Extract OTU_IDs
otu_ids_18s <- cleaned_18s[[1]]

# 2. Convert counts to relative abundance (excluding OTU_ID)
rel_abund_matrix_18s <- sweep(
  cleaned_18s[ , -1],           # exclude OTU_ID
  2,
  colSums(cleaned_18s[ , -1]),
  FUN = "/"
)


# 3. Combine back with OTU_IDs
rel_abund_18s <- cbind(OTU_ID = otu_ids_18s, rel_abund_matrix_18s)

otu_nmds_matrix_18s <- rel_abund_18s %>%
  column_to_rownames(var = "OTU_ID")
#transpose matrix
otu_nmds_matrix_18s_t <- t(otu_nmds_matrix_18s)
# Run NMDS on the transposed matrix
nmds_result_18s <- metaMDS(otu_nmds_matrix_18s_t, distance = "bray", k = 2, trymax = 100)

# 2. Extract NMDS points and add grouping info
nmds_points_18s <- as.data.frame(nmds_result_18s$points)
nmds_points_18s$Sample <- rownames(nmds_points_18s)

# Extract group name (e.g., "Bull" from "Bull-A")
nmds_points_18s$Group <- sub("-[A-Za-z]+$", "", nmds_points_18s$Sample)

# 3. Function to get convex hull points for each group
find_hull_18s <- function(df) df[chull(df$MDS1, df$MDS2), ]

hulls_18s <- nmds_points_18s %>%
  group_by(Group) %>%
  do(find_hull_18s(.))
# Identify groups with exactly two samples
groups_two_points_18s <- nmds_points_18s %>%
  group_by(Group) %>%
  filter(n() == 2) %>%
  ungroup()

# Your plot with polygons, points, labels, and lines connecting two-point groups
nmds_points_18s <- nmds_points_18s %>%
  mutate(Sample = case_when(
    Sample == "Horn-Elliott" ~ "Horn-E",
    Sample == "Sand-Control" ~ "Surface Control",
    TRUE ~ Sample
  ))

# Rename "Sand" to "Surface Control" in all relevant data frames
nmds_points_18s <- nmds_points_18s %>%
  mutate(Group = ifelse(Group == "Sand", "Surface Control", Group))

hulls_18s <- hulls_18s %>%
  mutate(Group = ifelse(Group == "Sand", "Surface Control", Group))

groups_two_points_18s <- groups_two_points_18s %>%
  mutate(Group = ifelse(Group == "Sand", "Surface Control", Group))

# Then create the plot with updated group names
nmds_rel_18s <- ggplot(nmds_points_18s, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_polygon(data = hulls_18s, aes(x = MDS1, y = MDS2, fill = Group),
               alpha = 0.2, color = NA) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample),
                  size = 4, max.overlaps = 20, box.padding = 0.4) +
  geom_line(data = groups_two_points_18s,
            aes(x = MDS1, y = MDS2, group = Group, color = Group),
            size = 2, alpha = 0.4, show.legend = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) +
  labs(color = "BNG Well", fill = "BNG Well")

# Assuming your ggplot object is the last plot you made
ggsave("nmds_replicates_rel_18a.svg", plot = nmds_rel_18s, width = 12, height = 10, units = "in")
ggsave("nmds_replicates_rel_18s.jpg", plot = nmds_rel_18s, width = 12, height = 10, units = "in", dpi = 2000)

#Results: Strong Technical replicate grouping using NMDS for 18s with relative abundances.







#### NMDS Grouping 16s ####
#Read in OTU table and
otu_table_16s <- read_tsv("16s/outputs/16s_OTU_table.tsv") %>%
  dplyr::rename(OTU_ID = `#OTU ID`)

colSums(otu_table_16s[, -1])

#filter out samples with less than 5,000 reads (Only Horn-Concentrate has less than 5,000 reads)
filtered_16s <- otu_table_16s %>%
  {
    keep_cols <- names(.)[-1][colSums(.[, -1], na.rm = TRUE) > 5000]
    dplyr::select(., OTU_ID, all_of(keep_cols))
  }

# Exclude specific samples from the dataset
trimmed_16s <- filtered_16s %>%
  dplyr::select(-`Sand-Control`, -`Mortensen-Elliott`) 

#convert to relative abundance
# Extract OTU IDs
otu_ids_16s <- trimmed_16s[[1]]

# Convert counts to relative abundances (excluding OTU_ID)
rel_abund_matrix_16s <- sweep(
  trimmed_16s[, -1],
  2,
  colSums(trimmed_16s[, -1]),
  FUN = "/"
)

# Combine back with OTU_ID
rel_abund_16s <- cbind(OTU_ID = otu_ids_16s, rel_abund_matrix_16s)

# Convert OTU_ID to rownames and transpose
otu_nmds_matrix_16s <- rel_abund_16s %>%
  column_to_rownames(var = "OTU_ID")

otu_nmds_matrix_16s_t <- t(otu_nmds_matrix_16s)

# Run NMDS

nmds_result_16s <- metaMDS(otu_nmds_matrix_16s_t, distance = "bray", k = 2, trymax = 100)

nmds_points_16s <- as.data.frame(nmds_result_16s$points)
nmds_points_16s$Sample <- rownames(nmds_points_16s)

# Extract group (e.g. "Bull" from "Bull-A")
nmds_points_16s$Group <- sub("-[A-Za-z]+$", "", nmds_points_16s$Sample)

# Convex hull function
find_hull_16s <- function(df) df[chull(df$MDS1, df$MDS2), ]

# Find hulls
hulls_16s <- nmds_points_16s %>%
  group_by(Group) %>%
  do(find_hull_16s(.))

# Identify 2-point groups
groups_two_points_16s <- nmds_points_16s %>%
  group_by(Group) %>%
  filter(n() == 2) %>%
  ungroup()

# Sample-level adjustments
nmds_points_16s <- nmds_points_16s %>%
  mutate(Sample = case_when(
    Sample == "Horn-Elliott" ~ "Horn-E",
    TRUE ~ Sample
  ))

# Group-level renaming
nmds_rel_16s <- ggplot(nmds_points_16s, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_polygon(data = hulls_16s, aes(x = MDS1, y = MDS2, fill = Group),
               alpha = 0.2, color = NA) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample),
                  size = 4, max.overlaps = 20, box.padding = 0.4) +
  geom_line(data = groups_two_points_16s,
            aes(x = MDS1, y = MDS2, group = Group, color = Group),
            size = 2, alpha = 0.4, show.legend = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  ) +
  labs(color = "BNG Well", fill = "BNG Well")

#nmds_rel_16s

ggsave("nmds_plot_16s.pdf", nmds_rel_16s, width = 10, height = 8)


#### Output final tables ####
#Lets make the final OTU table to drop OTUs isolated as controls and low quality samples
# Filter out OTUs that are present in the control isolates
otu_table_clean <- otu_table %>%
  dplyr::filter(!OTU_ID %in% otus_in_control)



otu_table_drop <- otu_table_clean %>%
  dplyr::select(-`Mortensen-Elliott`, -`Horn-Concentrate`)
#dropping these samples because they have extremely low read counts (below 300)


replicate_groups <- list(
  Bull = c("Bull-A", "Bull-B", "Bull-C"),
  Conant = c("Conant-A", "Conant-B", "Conant-C"),
  Greg = c("Greg-A", "Greg-B", "Greg-C"),
  Horn = c("Horn-Elliott"),
  Mortensen = c("Mortensen-A", "Mortensen-B", "Mortensen-C"),
  West = c("West-A", "West-B", "West-C"),
  Control = c("Sand-Control")
)


# Start with just the OTU_ID column
combined_df <- otu_table %>%
  select("OTU_ID")

# Loop over each group in replicate_groups and sum their columns
for (group_name in names(replicate_groups)) {
  replicate_cols <- replicate_groups[[group_name]]
  
  # Add the new combined column by summing across the replicates
  combined_df[[group_name]] <- otu_table %>%
    select(all_of(replicate_cols)) %>%
    rowSums(na.rm = TRUE)
}

#save as csv
write.csv(combined_df, "ITS2_OTU_table_unrare_control.csv", row.names = FALSE)

#Let us drop the control here because also has very few reads. And empty rows
otu_table_drop_all <- combined_df %>%
  select( -`Control`) %>%
  filter(rowSums(select(., -OTU_ID)) > 0)

write.csv(otu_table_drop_all, "ITS2_OTU_table_unrare.csv", row.names = FALSE)

#Rarefied OTU table
# Step 1: Remove OTU_ID column and save it separately
otu_counts <- otu_table_drop_all %>%
  column_to_rownames("OTU_ID")

# Step 2: Find the minimum sample depth
min_depth <- min(rowSums(t(otu_counts)))

# Step 3: Rarefy
otu_rarefied <- t(rrarefy(t(otu_counts), sample = min_depth))

# Step 4: Turn back into a data frame with OTU_ID
otu_rarefied_df <- as.data.frame(otu_rarefied) %>%
  rownames_to_column("OTU_ID")

write.csv(otu_rarefied_df, "ITS2_OTU_table_rare.csv", row.names = FALSE)


#Lets do CLR using ALDEx2 as recommended in "Microbiome differential abundance methods produce different results across 38 datasets" Nat Comm

# Assume OTU_ID is in the first column
otu_counts <- otu_table_drop_all %>%
  column_to_rownames("OTU_ID")

# Transpose so ALDEx2 has OTUs as rows and samples as columns
otu_counts_t <- as.data.frame(t(otu_counts))

# Perform CLR transformation (no group needed if you're not testing)
otu_clr_aldex <- aldex.clr(otu_counts_t, mc.samples = 128, denom = "all", verbose = FALSE)

# Extract just the CLR matrix
clr_matrix <- otu_clr_aldex@clr

# Transpose back to OTUs as columns, samples as rows
clr_matrix_t <- t(clr_matrix)

# Convert to data frame and add OTU_ID as column
otu_clr_df <- as.data.frame(clr_matrix_t) %>%
  rownames_to_column("Sample")



otu_counts <- otu_table_drop_all %>%
  column_to_rownames("OTU_ID") %>%
  t()

jpeg("rarefaction_curves_base.jpg", width = 8, height = 6, units = "in", res = 300)

rarecurve(otu_counts, 
          step = 100, 
          sample = min(rowSums(otu_counts)), 
          col = rainbow(nrow(otu_counts)), 
          cex = 0.6, 
          label = TRUE,
          xlab = "Sequencing Depth (Reads)",
          ylab = "Observed OTUs",
          main = "Rarefaction Curves for Each Pooled Sample")

dev.off()


#### Write out final OTU Tables and drop control OTUs ####
#The final df for matched isolates is final_df_with_ranks (written out to csv)
#final df for umatched is unmatched_isolates_pa_taxonomy (written to csv)

#Lets now make a complete df of matched and unmatched isolate OTUs
all_isolates <- bind_rows(final_df_with_ranks, unmatched_isolates_pa_taxonomy)

#Add a column that indicates if the isolate appears in the amplicon data at all
all_isolates <- all_isolates %>%
  mutate(
    appears_in_amplicon = if_any(ends_with("_Amplicon"), ~ .x > 0)
  )

all_isolates$appears_in_amplicon[is.na(all_isolates$appears_in_amplicon)] <- FALSE

#Add a column that indicates if the isolate appears in the isolate controls at all
all_isolates <- all_isolates %>%
  mutate(
    appears_in_control_isolates = if_else(
      !is.na(Control_Isolate) & Control_Isolate > 0,
      "Yes",
      "No"
    )
  )

#write to .csv
write_csv(all_isolates, "all_isolates.csv")


#Lets make and write out a final table of the non contaminated isolates. 
otus_filtered_isolates <- all_isolates %>%
  filter(appears_in_control_isolates != "Yes")

#This will serve as final OTU table for isolates
write_csv(otus_filtered_isolates, "all_isolates_trimmed.csv")


#Results: 79 total isolate OTUs at 0.98. 12 of these OTUs appear in the controls and can be dropped from the final OTU table.
# 42 of these OTUs appear in the amplicon data.37 do not match amplicon data at 0.98. 
#12 of the isolate OTUs appear in the controls (8 matched, 4 unmatched). 
#Total= 36 isolate OTUs that are in amplicon and not contamination. 33 isolates that do not match amplicon data and are not contamination. 69 total isolate OTUs that are not contamination.

#Lets check to see which isolate OTUs appear to be undescribed species.
otus_below_98 <- all_isolates %>%
  filter(percent_identity_Unite < 98) %>%
  pull(OTU_ID)

#Results: 14 isolate OTUs have no Unite match over 98% identity.
# [1] "081720e5f9f84e1181553a050cbae7f2" "b58c8c1b809288c1e1404daf2a7a8de3" "d3a69218f40f0b9e4b41ed06a456bb55" "OTU_1"  "OTU_15"                          
#[6] "OTU_17"      "OTU_19"   "OTU_24"    "OTU_28"     "OTU_3"  [11] "OTU_32"   "OTU_33"     "OTU_36"   "OTU_9" 


####ITS2 
#Durig NMDS, we nade (otu_nmds_data_filtered_ITS) that has contaminants removed and samples with less than 230 reads removed.
(otu_nmds_data_filtered_ITS)
colSums(otu_nmds_data_filtered_ITS[, -1]) # 42 OTUs, 13 samples

#Merge tachnical replicates
replicate_groups_ITS <- list(
  Bull = c("Bull-A", "Bull-B", "Bull-C"),
  Conant = c("Conant-A", "Conant-B", "Conant-C"),
  Greg = c("Greg-A"),
  Horn = c("Horn-Elliott"),
  Mortensen = c("Mortensen-A", "Mortensen-B", "Mortensen-C"),
  West = c("West-A","West-C"),
  "Surface Control" = c("Sand-Control")
)
# Start with just the OTU_ID column
pooled_ITS <- otu_nmds_data_filtered_ITS %>%
  dplyr::select("OTU_ID")

# Loop over each group in replicate_groups and sum their columns
for (group_name in names(replicate_groups_ITS)) {
  replicate_cols <- replicate_groups_ITS[[group_name]]
  
  # Add the new combined column by summing across the replicates
  pooled_ITS[[group_name]] <- otu_nmds_data_filtered_ITS %>%
    dplyr::select(all_of(replicate_cols)) %>%
    rowSums(na.rm = TRUE)
}

#Write out OTU table. This will serve as final unrarified OTU table for analysis. 
write_csv(pooled_ITS, "ITS2_OTU_table_trimmed.csv")

sum(rowSums(pooled_ITS[, -1])) 
dim(pooled_ITS) # 42 OTUs, 13 samples
#results: 755 OTUs with ITS2 primers. 35652 total reads included in final OTU ITS table

#### 18s
#Durig NMDS, we nade (cleaned_18s) that has land plants and large animals contaminants removed and samples with less than 400 reads removed.
colSums(cleaned_18s[, -1]) # 342 OTUs, 12 samples
#Merge tachnical replicates
replicate_groups_18s <- list(
  Bull = c("Bull-A"),
  Mortensen = c("Mortensen-A", "Mortensen-B", "Mortensen-C"),
  West = c("West-A"),
  "Surface Control" = c("Sand-Control")
)

# Start with just the OTU_ID column
pooled_18s <- cleaned_18s %>%
  dplyr::select("OTU_ID")
# Loop over each group in replicate_groups and sum their columns
for (group_name in names(replicate_groups_18s)) {
  replicate_cols <- replicate_groups_18s[[group_name]]
  
  # Add the new combined column by summing across the replicates
  pooled_18s[[group_name]] <- cleaned_18s %>%
    dplyr::select(all_of(replicate_cols)) %>%
    rowSums(na.rm = TRUE)
}
View(pooled_18s)
#Write out OTU table. This will serve as final unrarified OTU table for analysis.
write_csv(pooled_18s, "18s_OTU_table_trimmed.csv")
#results: Only 3 of the 6 sites have enough 18s results to proceed with analysis (Bull, Mortensen, West)

####16s
#During NMDS, we nade (trimmed_16s) that has samples with less than 5,000 reads removed, and Horn-Concentrate and Mortensen-Elliott samples removed.
trimmed_16s

#only keep OTUs that have at least 1 read across all samples
trimmed_16s <- trimmed_16s[rowSums(trimmed_16s[ , !colnames(trimmed_16s) %in% "OTU_ID"]) >= 1, ]

dim(trimmed_16s) 
#Results 5043 ASVs, 16 samples


# Separate OTU_IDs from sample data
otu_ids_16s <- trimmed_16s$OTU_ID
otu_data_16s <- trimmed_16s[, !names(trimmed_16s) %in% "OTU_ID"]


#Merge tachnical replicates
replicate_groups_16s <- list(
  Bull = c("Bull-A", "Bull-B", "Bull-C"),
  Conant = c("Conant-A", "Conant-B", "Conant-C"),
  Greg = c("Greg-A", "Greg-B", "Greg-C"),
  Horn = c("Horn-Elliott"),
  Mortensen = c("Mortensen-A", "Mortensen-B", "Mortensen-C"),
  West = c("West-A","West-C", "West-B")
)

# Pool samples by summing replicates
pooled_data <- sapply(replicate_groups_16s, function(samples) {
  rowSums(otu_data_16s[, samples, drop = FALSE])
})

# Combine with OTU_IDs
pooled_trimmed_16s <- cbind(OTU_ID = otu_ids_16s, as.data.frame(pooled_data))

# Write out the final pooled OTU table
write_csv(pooled_trimmed_16s, "16s_OTU_table_trimmed.csv")



