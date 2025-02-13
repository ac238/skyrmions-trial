#!/bin/bash
#SBATCH --job-name=10-17_ZSC_scattering_a50_comm
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32g
#SBATCH --time=24:00:00
#SBATCH --mail-user=ac238@rice.edu
#SBATCH --mail-type=ALL
echo "My job ran on:" echo $SLURM_NODELIST
module load GCC
module load GCC MATLAB/2023b
if [[ -d $SHARED_SCRATCH/ac238 && -w $SHARED_SCRATCH/ac238 ]]; then
     echo "Conditional satisfied"
     cd $SHARED_SCRATCH/ac238
     ncpus=1
     framespace=100
     bval=0.13
     stiff=1
     elval=0
     alpha=5.0
     dampval=0.4
     dampfalloff=2
     N=100
     magamp=0.1
     tfinal=500
     matlab -r "skyrmion_passargs('10.17.ZSC.scattering.a5.0.comm', $ncpus, $framespace, $bval, $stiff, $elval, $alpha, $dampval, $dampfalloff, $N, $magamp, $tfinal); exit"
     echo "Job finished"
fi
