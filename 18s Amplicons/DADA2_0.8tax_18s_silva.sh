#!/bin/bash

# The job should run on the standard partition
#SBATCH -p standard

# The name of the job is DADA2 taxonomy 18s silva
#SBATCH -J DADA2_tax_18s_silva

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

# Load R module (adjust if using a specific version/module system)
module load R/4.4.0

# Run the R script
Rscript 18s_taxon_dada2_silva.R
