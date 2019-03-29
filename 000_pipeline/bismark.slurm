#!/bin/bash

#SBATCH --time=0-23:59
#SBATCH --mem=20G
#SBATCH --account=def-bruneaua
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

module load r/3.5.0
module load bowtie2
module load samtools

cd $SCRATCH/sophie_breton/Bouvet_Venustaconcha/000_pipeline/

./bismark.R


