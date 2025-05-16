#!/bin/bash

# The job should run on the standard partition
#SBATCH -p standard

# The name of the job is 
#SBATCH -J uc_to_csv

# The job requires 1 compute node
#SBATCH -N 1

#SBATCH --account=tyjames0

# The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=12

#SBATCH --mem=100gb

# The maximum walltime of the job is 8 days
#SBATCH -t 2-00:00:00

#SBATCH --mail-user=qmoon@umich.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

python <<EOF
import csv

# Define the input and output files
input_file = "clusters.uc"
output_file = "clusters.csv"

# Initialize a dictionary to hold the clusters
clusters = {}

# Read the .uc file
with open(input_file, "r") as infile:
    for line in infile:
        # Skip comment lines starting with '#'
        if line.startswith("#"):
            continue

        # Split the line into components
        parts = line.strip().split("\t")
        record_type = parts[0]  # 'S' for representative, 'H' for member
        cluster_id = parts[1]   # Cluster ID
        seq_id = parts[8]       # Sequence ID

        # If the line is a member (H), we get the representative sequence too
        if record_type == 'H':
            rep_id = parts[9]  # Representative ID
        else:
            rep_id = seq_id    # If it's the representative, use its ID

        # If the cluster ID is not in the dictionary, add it
        if cluster_id not in clusters:
            clusters[cluster_id] = {"representative": rep_id, "members": []}

        # Add the sequence to the list of members for the cluster
        if record_type == "H":
            clusters[cluster_id]["members"].append(seq_id)

# Write the results to a CSV
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Cluster_ID", "Representative_Sequence", "Member_Sequences"])

    # Write each cluster and its sequences
    for cluster_id, cluster_data in clusters.items():
        # Join member sequences with commas for each cluster
        member_sequences = ", ".join(cluster_data["members"])
        writer.writerow([cluster_id, cluster_data["representative"], member_sequences])

print(f"CSV file has been saved to {output_file}")
EOF