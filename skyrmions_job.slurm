#!/bin/bash
#SBATCH --job-name=skyrmion
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
     echo "Trying code: module load"
     module load GCC MATLAB/2023b
     module load GCCcore/12.3.0
     module load Python/3.11.3
     module list
     git checkout nots
     echo "Trying code: git pull origin nots"
     git pull origin nots
     echo "Trying code: matlab -r"
     matlab -r "skyrmionstrial; exit"
     echo "Trying code: git add -A"
     git add -A
     echo "Trying code: git commit"
     git commit -m "57231"
     echo "Trying code: git push origin nots"
     git push origin nots
fi
