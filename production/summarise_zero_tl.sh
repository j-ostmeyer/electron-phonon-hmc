#!/bin/bash

#out_name=summary_mat_abb-tl.csv
out_name=summary_zero_abb-tl.csv

rm -f $out_name

nt=10000

#for dir in L_*_Nt_${nt}_mu0_1_*_mat_*tl*/; do
for dir in L_*_Nt_${nt}_mu0_1_*_hop_*tl*/; do
	res="${dir}/results2.csv"
	corr="${dir}/greens_smeared.csv"
	if [ ! -e "$res" ]; then continue; fi

	echo $dir

	L=$(echo $dir | cut -d'_' -f2)
	Nt=$(echo $dir | cut -d'_' -f5)
	mu0=$(echo $dir | cut -d'_' -f7)
	beta=$(echo ${dir} | cut -d'_' -f9)
	hop=$(echo ${dir} | cut -d'_' -f11 | awk '{print 0.1*$1}')
	theta=$(echo ${dir%/} | cut -d'_' -f13 | awk '{print 0.3141593*$1}')

	n_avg=$(comp-avg -n "$res" | cut -f2)
	cond=$(head -1 "$corr")

	echo "$L	$Nt	-0.25	$beta	$hop	$theta	$n_avg	$cond" >> $out_name
	echo "$L	$Nt	0	$beta	$hop	$theta	$n_avg	$cond" >> $out_name
done
