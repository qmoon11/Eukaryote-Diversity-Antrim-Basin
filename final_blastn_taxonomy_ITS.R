#This script loads the blastn results from a custom search against UNITE and parses each taxonomic rank into own columm.
#Next, it calculates the % identity of the best hit for each OTU.

#### Load and Clean blastn taxonomy ####
# Define the column names
column_names <- c(
  "OTU_ID", "subject_id", "percent_identity", "alignment_length", "mismatches", 
  "gap_opens", "query_start", "query_end", "subject_start", "subject_end", 
  "e_value", "bit_score", "taxonomic_classification"
)

# Read the BLAST result
blast_results <- read_tsv("blastn_ITS2.tsv", col_names = column_names)

#Parse out taxonomic classification into seperate coloumns
cleaned_taxonomy <- filtered_blast_results %>%
  # Extract only the taxonomy portion starting at k__
  mutate(taxonomy_string = str_extract(taxonomic_classification, "k__.*")) %>%
  
  # Separate into columns by taxonomic rank
  separate(taxonomy_string,
           into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right", remove = FALSE) %>%
  
  # Remove k__, p__, etc. prefixes
  mutate(across(kingdom:species, ~str_remove(., "^[a-z]__")))

#write out to csv
write.csv(cleaned_taxonomy, "cleaned_taxonomy_blastn_ITS.csv", row.names = FALSE)

#This CV is then inspected according to the blastn parameters provided in 
#"Best practices in metabarcoding of fungi: From experimental design to results", Tedersoo et al., (2022)


#count number of unique OTU_Ids in cleaned_taxonomy
unique_OTU_count <- cleaned_taxonomy %>%
  group_by(OTU_ID) %>%
  summarise(count = n(), .groups = "drop")

unique_OTU_count






#### Calculate highest percent identity #### 

min_alignment_length <- 175

# Apply logic using the variable
highest_identity_per_otu <- blast_results %>%
  group_by(OTU_ID) %>%
  summarise(
    max_pid = if (any(alignment_length > min_alignment_length)) {
      max(percent_identity[alignment_length > min_alignment_length])
    } else {
      mean(percent_identity)
    },
    .groups = "drop"
  ) %>%
  rename(percent_identity = max_pid)

#write out csv
write.csv(highest_identity_per_otu, "highest_%identity_per_otu.csv", row.names = FALSE)

#For OTUs with any results having an alignment over the min_alignment_length parameter, this takes the highest %. 
#If no results over min_alignment_length, then it averages the percent identity of all results for that OTU.

