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
     module load GCC MATLAB/2023b
     module load GCCcore/12.3.0
     module load Python/3.11.3
     module list
     git checkout nots
     echo "Switched to nots branch"
     git pull origin nots
     echo "Pulled from nots branch"
     matlab -r "skyrmionstrial; exit"
     echo "Ran MATLAB"
     git add -A
     echo "Added all"
     git commit -m "56532"
     echo "Committed"
     git push origin nots
     echo "Pushed to nots"
fi
