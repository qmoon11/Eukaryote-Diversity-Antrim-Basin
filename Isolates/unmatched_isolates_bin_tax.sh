#!/bin/bash

#SBATCH -p standard
#SBATCH -J vsearch_blastn_unidentified
#SBATCH -N 1
#SBATCH --account=tyjames0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=qmoon@umich.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#conda install -c bioconda vsearch

# Load paths and modules
export PATH=/scratch/tyjames_root/tyjames0/qmoon/Antrim-Shale-Microbiome/ITS2/ncbi-blast-2.16.0+/bin:$PATH
conda activate vsearch_env

# Create output directories
mkdir -p unmatched_isolates/clustered_otus unmatched_isolates/blast_results

# Step 1: Cluster sequences at 98% and pick longest rep. These are the sequences of isolates that did not match 
# to the ITS2 amplicon data. 
vsearch --cluster_fast unmatched_isolates/unmatched_isolates.fasta \
        --id 0.98 \
        --centroids unmatched_isolates/clustered_otus/rep_seqs_unmatched_98.fasta \
        --uc unmatched_isolates/clustered_otus/clusters.uc \
        --relabel OTU_ \
        --threads 12

#results: 28025 nt in 56 seqs. 37 OTUs at 98%

# Step 2: BLAST representative sequences using Unite database.
blastn \
  -query unmatched_isolates/clustered_otus/rep_seqs_unmatched_98.fasta \
  -db /scratch/tyjames_root/tyjames0/qmoon/Antrim-Shale-Microbiome/ITS2/unite_db \
  -out unmatched_isolates/blast_results/blastn_unidentified.tsv \
  -evalue 1e-20 \
  -word_size 7 \
  -reward 1 \
  -penalty -1 \
  -gapopen 1 \
  -gapextend 2 \
  -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
  -num_threads 12

# Step 3: Convert .uc cluster info to CSV. This will track which isolates belong to each OTU. 
python3 <<EOF
import csv
import os

input_file = "unmatched_isolates/clustered_otus/clusters.uc"
output_file = "unmatched_isolates/clustered_otus/clusters.csv"
clusters = {}

with open(input_file, "r") as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        record_type = parts[0]
        cluster_id = parts[1]
        seq_id = parts[8]
        rep_id = parts[9] if record_type == 'H' else seq_id
        if cluster_id not in clusters:
            clusters[cluster_id] = {"representative": rep_id, "members": []}
        if record_type == "H":
            clusters[cluster_id]["members"].append(seq_id)

with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Cluster_ID", "Representative_Sequence", "Member_Sequences"])
    for cluster_id, cluster_data in clusters.items():
        member_sequences = ", ".join(cluster_data["members"])
        writer.writerow([cluster_id, cluster_data["representative"], member_sequences])

print(f"CSV file has been saved to {output_file}")
EOF