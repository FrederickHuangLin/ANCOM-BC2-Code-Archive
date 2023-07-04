#!/bin/bash

#SBATCH --partition=norm # Job to run on partition 'partname'. (default: 'norm')
#SBATCH --ntasks=1 # Number of task to be run
#SBATCH --cpus-per-task=8 # Number of CPUs required for each task
#SBATCH --mem=200g # Memory required for the job. 
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END,FAIL  # Mail events
#SBATCH --mail-user=huang.lin@nih.gov   # Where to send mail 
#SBATCH --time=00-12:00:00 # Runtime in day-hour:min:sec
#SBATCH --requeue # If an allocated node hangs, the job should be requeued
#SBATCH --output=urt_sim_dunn_ancombc.out # Location of stdout file
#SBATCH --error=urt_sim_dunn_ancombc.err # Location of std class="softBottom"err file

pwd; hostname; date

echo "start my job"
module load R/4.2.2
echo "finish loading R"

Rscript urt_sim_dunn_ancombc.R
echo "done"
date