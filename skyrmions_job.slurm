#!/bin/bash
#SBATCH --job-name=running_testcode
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=00:30:00
#SBATCH --mail-user=ac238@rice.edu
#SBATCH --mail-type=ALL
echo "My job ran on:" echo $SLURM_NODELIST
if [[ -d $PROJECTS/mf23/ac238/skyrmions-trial && -w $PROJECTS/mf23/ac238/skyrmions-trial ]]; then
     echo "Conditional satisfied"
     cd $PROJECTS/mf23/ac238/skyrmions-trial
     git checkout dev
     git pull origin dev
     matlab -r "skyrmionstrial; exit"
     echo "Job finished"
     git add -A
     git commit -m "Supercomputer"
     git push origin dev
fi
