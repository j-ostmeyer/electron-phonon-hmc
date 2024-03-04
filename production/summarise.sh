#!/bin/bash

#rm -f summary_angle.csv

for dir in L_*_Nt_*_hop_*/; do
	res="${dir}/results2.csv"
	if [ ! -e "$res" ]; then continue; fi

	L=$(echo $dir | cut -d'_' -f2)
	Nt=$(echo $dir | cut -d'_' -f5)
	mu0=$(echo $dir | cut -d'_' -f7)
	beta=$(echo ${dir} | cut -d'_' -f9)
	hop=$(echo ${dir} | cut -d'_' -f11 | awk '{print 0.1*$1}')
	theta=$(echo ${dir%/} | cut -d'_' -f13 | awk '{print 0.3141593*$1}')

	if [ "$mu0" = "1" ]; then continue; fi

	x_avg=$(awk '{print $1}' "$res" | comp-avg -n | cut -f2,4,5)
	x2_avg=$(awk '{print $2}' "$res" | comp-avg -n | cut -f2,4,5)
	n_avg=$(awk '{print $6}' "$res" | comp-avg -n | cut -f2,4,5)
	cond=$(awk '{print $12}' "$res" | comp-avg -n | cut -f2,4,5)
	mob=$(awk '{print $12/$6}' "$res" | comp-avg -n | cut -f2,4,5)
	acc=$(awk '{print $7}' "$res" | comp-avg -n | cut -f2)
	boltz=$(awk '{print $8}' "$res" | comp-avg -n | cut -f2)

	echo "$L	$Nt	$mu0	$beta	$hop	$theta	$x_avg	$n_avg	$cond	$acc	$boltz	$x2_avg	$mob" >> summary_angle.csv
done
