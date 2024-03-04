#!/usr/bin/bash

#SBATCH -p lowpriority,nodes,schaich
#SBATCH -D /users/ostmeyer/volatile2/Organic_Semiconductors/production/
#SBATCH -J summarise_zero
#SBATCH --output=/users/ostmeyer/volatile2/Organic_Semiconductors/production/summarise_zero.out --error=/users/ostmeyer/volatile2/Organic_Semiconductors/production/summarise_zero.err
#SBATCH -t 1-0:0
#SBATCH -N 1 -n 1 -c 1

module purge
module load libs/intel/2019u5
module load libs/intel-mkl/2019u5/bin
#module load packages/intel-studio-2018
#module load compilers/gcc/9.3.0
#module load libs/openblas/0.3.10/gcc-9.3.0
#module load libs/fftw3_double/3.3.9/gcc-5.5.0+openmpi-4.1.1

export OMP_NUM_THREADS=1

echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
pwd

#out_name=summary_zero.csv
out_name=summary_zero-zf.csv

rm -f $out_name

nt=64

for dir in L_*_Nt_${nt}_mu0_1_*_hop_*-zf/; do
	res="${dir}/results2.csv"
	if [ ! -e "$res" ]; then continue; fi

	echo $dir

	L=$(echo $dir | cut -d'_' -f2)
	Nt=$(echo $dir | cut -d'_' -f5)
	mu0=$(echo $dir | cut -d'_' -f7)
	beta=$(echo ${dir} | cut -d'_' -f9)
	hop=$(echo ${dir} | cut -d'_' -f11 | awk '{print 0.1*$1}')
	theta=$(echo ${dir%/} | cut -d'_' -f13 | awk '{print 0.3141593*$1}')

	x_avg=$(awk '{print $1}' "$res" | comp-avg -n | cut -f2,4,5)
	x2_avg=$(awk '{print $2}' "$res" | comp-avg -n | cut -f2,4,5)
	n_avg=$(awk '{print $6}' "$res" | comp-avg -n | cut -f2,4,5)
	cond=$(awk '{print $12}' "$res" | comp-avg -n | cut -f2,4,5)
	mob=$(awk '{print $12/$6}' "$res" | comp-avg -n | cut -f2,4,5)
	acc=$(awk '{print $7}' "$res" | comp-avg -n | cut -f2)
	boltz=$(awk '{print $8}' "$res" | comp-avg -n | cut -f2)

	cd $dir
	../../GKSolver/bin/gk_solver -t $((${nt}/2)) -i -c corr.csv -p ../params8.txt > /dev/null
	gk8=$(head -2 rho_final.txt | tail -1 | awk "{print \$2/${nt}, \$3/${nt}}")
	../../GKSolver/bin/gk_solver -t $((${nt}/2)) -i -c corr.csv -p ../params9.txt > /dev/null
	gk9=$(head -2 rho_final.txt | tail -1 | awk "{print \$2/${nt}, \$3/${nt}}")
	cd ..

	echo "$L	$Nt	-0.25	$beta	$hop	$theta	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob	$gk9	$gk8" >> $out_name
	echo "$L	$Nt	0	$beta	$hop	$theta	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob	$gk9	$gk8" >> $out_name
done

exit 0
