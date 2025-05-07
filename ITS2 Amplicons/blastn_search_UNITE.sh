#!/bin/bash

# The job should run on the standard partition
#SBATCH -p standard

# The name of the job is 
#SBATCH -J blastn_ITS2

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


# Step 1: # Add BLAST+ to PATH

export PATH=/scratch/tyjames_root/tyjames0/qmoon/Antrim-Shale-Microbiome/ITS2/ncbi-blast-2.16.0+/bin:$PATH



# Step 2: Run BLASTn
blastn \
  -query rep-seqs/dna-sequences.fasta \
  -db unite_db \
  -out blastn_ITS2.tsv \
  -evalue 1e-20 \
  -word_size 7 \
  -reward 1 \
  -penalty -1 \
  -gapopen 1 \
  -gapextend 2 \
  -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
  -num_threads 4

  #Parameters based on Leho best practices barcoding fungi and Tedersoo et al., (2014) Science.