#!/bin/bash

#SBATCH --partition=largemem # Job to run on partition 'partname'.(default: 'norm')
#SBATCH --ntasks=1 # Number of task to be run
#SBATCH --cpus-per-task=40 # Number of CPUs required for each task on
#SBATCH --gres=lscratch:400 # Use lscratch for temporary storage
#SBATCH --mem=1000G # Use 1TB
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END,FAIL  # Mail events
#SBATCH --mail-user=huang.lin@nih.gov   # Where to send mail 
#SBATCH --time=1-00:00:00 # Runtime in day-hour:min:sec
#SBATCH --requeue # If an allocated node hangs, the job should be requeued
#SBATCH --output=sim_pair_ancombc2.out # Location of stdout file
#SBATCH --error=sim_pair_ancombc2.err # Location of std class="softBottom"err file

set -e

pwd; hostname; date

echo "start my job"
module load R/4.2.2
echo "finish loading R"

Rscript sim_pair_ancombc2.R
echo "done"
date