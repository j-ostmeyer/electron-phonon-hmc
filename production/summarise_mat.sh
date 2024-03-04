#!/bin/bash

#rm -f summary_mat.csv

#for dir in L_*_Nt_*_mat_*/; do
for dir in L_*_Nt_*_mu0_1_*_mat_*-h/; do
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

	if [ "$mat" = "rubrene-h" ]; then
		area=51.84
	elif [ "$mat" = "pentacene-h" ]; then
		area=24.28
	elif [ "$mat" = "C10-DNTT-h" ]; then
		area=22.83
	elif [ "$mat" = "C10-DNBDT-h" ]; then
		area=23.62
	elif [ "$mat" = "DNTT-h" ]; then
		area=47.4
	elif [ "$mat" = "C8-DNTT-C8-h" ]; then
		area=47.07
	fi

	#echo "$L	$Nt	$mu0	$beta	$mat	$area	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob" >> summary_mat.csv
	echo "$L	$Nt	-0.25	$beta	$mat	$area	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob" >> summary_mat_z.csv
	echo "$L	$Nt	0	$beta	$mat	$area	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob" >> summary_mat_z.csv
done
