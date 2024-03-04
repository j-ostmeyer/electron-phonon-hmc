#!/usr/bin/bash

#SBATCH -p lowpriority,nodes,schaich
#SBATCH -D /users/ostmeyer/volatile2/Organic_Semiconductors/production/NAME/
#SBATCH -J NAME
#SBATCH --output=/users/ostmeyer/volatile2/Organic_Semiconductors/production/NAME/organic.out --error=/users/ostmeyer/volatile2/Organic_Semiconductors/production/NAME/organic.err
#SBATCH -t 1-0:0
#SBATCH -N 1 -n 1 -c 1

module purge
module load libs/intel/2019u5
module load libs/intel-mkl/2019u5/bin
#module load packages/intel-studio-2018
#module load compilers/gcc/9.3.0
#module load libs/openblas/0.3.10/gcc-9.3.0
module load libs/fftw3_double/3.3.9/gcc-5.5.0+openmpi-4.1.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PROC_BIND=close
#export OMP_PLACES=cores

echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
pwd

ulimit -s unlimited
#time ./organic_main input.txt results2.csv corr.csv
#time ./organic_main input.txt results2.csv greens.bin
time ./organic_conv greens.bin smear.txt greens_smeared.csv

exit 0
