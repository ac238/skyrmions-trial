#!/bin/bash
#SBATCH --job-name=25.2.3.ZSCE.bilayer.Vfitscatter.v2
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=32g
#SBATCH --time=24:00:00
#SBATCH --mail-user=ac238@rice.edu
#SBATCH --mail-type=ALL
echo "My job ran on:" echo $SLURM_NODELIST
module load GCC
module load GCC MATLAB/2024a
if [[ -d $SHARED_SCRATCH/ac238 && -w $SHARED_SCRATCH/ac238 ]]; then
     echo "Conditional satisfied"
     cd $SHARED_SCRATCH/ac238
     framespace=100
     bval=0.13
     stiff=1
     elval=1
     alpha1=0
     alpha2=0
     alpha12=25
     time2reach=50
     dampval=0.8
     dampfalloff=2
     N=100
     lambda=5
     magamp=0.25
     magfreq=2.13
     tfinal=150
     Lspace=2
     z1=0
     z2=0
     matlab -r "bilayer_passargs_par('25.2.3.ZSCE.bilayer.Vfitscatter.v2', $framespace, $bval, $stiff, $elval, $alpha1, $alpha2, $alpha12, $time2reach, $dampval, $dampfalloff, $N, $lambda, $magamp, $magfreq, $tfinal, $Lspace, $z1, $z2); exit"
     echo "Job finished"
fi
