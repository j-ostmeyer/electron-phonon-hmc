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
	mat=$(echo ${dir%/} | cut -d'_' -f11)

	n_avg=$(comp-avg -n "$res" | cut -f2)
	cond=$(head -1 "$corr")

	echo "$mat	$n_avg	$cond" >> $out_name
done
