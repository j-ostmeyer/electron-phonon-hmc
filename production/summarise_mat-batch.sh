#!/usr/bin/bash

#SBATCH -p schaich
#SBATCH -D /users/ostmeyer/volatile2/Organic_Semiconductors/production/
#SBATCH -J summarise_mat
#SBATCH --output=/users/ostmeyer/volatile2/Organic_Semiconductors/production/summarise_mat.out --error=/users/ostmeyer/volatile2/Organic_Semiconductors/production/summarise_mat.err
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

#out_name=summary_mat.csv
#out_name=summary_mat_z.csv
#out_name=summary_mat_f.csv
#out_name=summary_mat_geom.csv
out_name=summary_mat_abb.csv

rm -f $out_name

nt=48

#for dir in L_*_Nt_*_mat_*/; do
#for dir in L_*_Nt_${nt}_mu0_1_*_mat_*1D/; do
for dir in L_*_Nt_${nt}_mu0_1_*_mat_*abb-zf*/ L_*_Nt_${nt}_mu0_1_*_mat_*abb-sp*/; do
	res="${dir}/results2.csv"
	if [ ! -e "$res" ]; then continue; fi

	echo $dir

	L=$(echo $dir | cut -d'_' -f2)
	Nt=$(echo $dir | cut -d'_' -f5)
	mu0=$(echo $dir | cut -d'_' -f7)
	beta=$(echo ${dir} | cut -d'_' -f9)
	mat=$(echo ${dir%/} | cut -d'_' -f11)

	x_avg=$(awk '{print $1}' "$res" | comp-avg -n | cut -f2,4,5)
	x2_avg=$(awk '{print $2}' "$res" | comp-avg -n | cut -f2,4,5)
	n_avg=$(awk '{print $6}' "$res" | comp-avg -n | cut -f2,4,5)
	cond=$(awk '{print $12}' "$res" | comp-avg -n | cut -f2,4,5)
	mob=$(awk '{print $12/$6}' "$res" | comp-avg -n | cut -f2,4,5)
	acc=$(awk '{print $7}' "$res" | comp-avg -n | cut -f2)
	boltz=$(awk '{print $8}' "$res" | comp-avg -n | cut -f2)

	if [ "$mat" = "rubrene" ] || [ "$mat" = "rubrene-h" ] || [ "$mat" = "rubrene-1D" ]; then
		area=51.84
	elif [ "$mat" = "Rubrene+f" ] || [ "$mat" = "Rubrene+geom" ] || [ "$mat" = "Rubrene+sp" ]; then
		area=51.48
	elif [ "$mat" = "C8-BTBT+f" ] || [ "$mat" = "C8-BTBT+geom" ]; then
		area=22.57
	elif [ "$mat" = "TIPS-PN+f" ] || [ "$mat" = "TIPS-PN+geom" ]; then
		area=28.88
	elif [ "$mat" = "TESADT+f" ] || [ "$mat" = "TESADT+geom" ]; then
		area=49.32
	elif [ "$mat" = "diFTESADT+f" ] || [ "$mat" = "diFTESADT+geom" ]; then
		area=56.93
	elif [ "$mat" = "pentacene" ] || [ "$mat" = "pentacene-h" ]; then
		area=24.28
	elif [ "$mat" = "C10-DNTT" ] || [ "$mat" = "C10-DNTT-h" ]; then
		area=22.83
	elif [ "$mat" = "C10-DNBDT" ] || [ "$mat" = "C10-DNBDT-h" ]; then
		area=23.62
	elif [ "$mat" = "DNTT" ] || [ "$mat" = "DNTT-h" ]; then
		area=47.4
	elif [ "$mat" = "C8-DNTT-C8" ] || [ "$mat" = "C8-DNTT-C8-h" ]; then
		area=47.07
	fi
	area=1

	cd $dir
	../../GKSolver/bin/gk_solver -t $((${nt}/2)) -i -c corr.csv -p ../params8.txt > /dev/null
	gk8=$(head -2 rho_final.txt | tail -1 | awk "{print \$2/${nt}, \$3/${nt}}")
	../../GKSolver/bin/gk_solver -t $((${nt}/2)) -i -c corr.csv -p ../params9.txt > /dev/null
	gk9=$(head -2 rho_final.txt | tail -1 | awk "{print \$2/${nt}, \$3/${nt}}")
	cd ..

	#echo "$L	$Nt	$mu0	$beta	$mat	$area	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob" >> summary_mat.csv
	echo "$L	$Nt	-0.25	$beta	$mat	$area	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob	$gk9	$gk8" >> $out_name
	echo "$L	$Nt	0	$beta	$mat	$area	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob	$gk9	$gk8" >> $out_name
done

exit 0
