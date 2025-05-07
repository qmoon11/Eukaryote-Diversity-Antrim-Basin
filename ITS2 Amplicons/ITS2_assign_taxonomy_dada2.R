
library(dada2)
library(Biostrings)
library(R.utils)

fasta_file <- "rep-seqs/dna-sequences.fasta"

unzip_ref_db <- "sh_general_release_dynamic_s_all_19.02.2025.fasta"


# Read in the ASVs from the FASTA file
sequences_ITS <- readDNAStringSet(fasta_file)

# Extract the headers (ASV identifiers) from the FASTA file
headers <- names(sequences_ITS)

taxa <- assignTaxonomy(sequences_ITS,
                       unzip_ref_db,
                       multithread = TRUE,
                       minBoot = 80,           # Set the minimum bootstrap threshold
                       outputBootstraps = TRUE, # Return the bootstrap value  # Use the vector for taxonomic levels
                       verbose = TRUE)

# Convert the taxonomy result to a data frame
taxa_df <- as.data.frame(taxa)

# Add the headers (ASV identifiers) to the taxonomy data frame
taxa_df$ASV <- headers

# Move the row names (sequence identifiers) to a separate column
taxa_df$Sequence_ID <- rownames(taxa_df)


saveRDS(taxa_df, "taxa_df.rds")

# Optionally, save the results to a CSV
write.csv(taxa_df, "ITS_taxonomic_unite_with_bootstrap.csv", row.names = FALSE)
