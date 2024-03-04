#!/bin/bash

fill_input () {
	if [ "$mat" = "rubrene" ]; then
		kappa="0.1407 0.0212 0.0212"
		dKappa="0.282 0.497 0.497"
	elif [ "$mat" = "rubrene-1D" ]; then
		kappa="0.155 0 0"
		dKappa="0.13 0 0"
	elif [ "$mat" = "pentacene" ]; then
		kappa="0.051 -0.0743 0.1306"
		dKappa="0.327 0.330 0.234"
	elif [ "$mat" = "C10-DNTT" ]; then
		kappa="0.1098 -0.0625 -0.0625"
		dKappa="0.21 0.647 0.647"
	elif [ "$mat" = "C10-DNBDT" ]; then
		kappa="0.0768 0.0591 0.0591"
		dKappa="0.38 0.592 0.592"
	elif [ "$mat" = "DNTT" ]; then
		kappa="0.0874 0.0215 -0.1179"
		dKappa="0.186 0.902 0.224"
	elif [ "$mat" = "C8-DNTT-C8" ]; then
		kappa="0.0931 -0.0607 -0.0607"
		dKappa="0.151 0.499 0.499"
	elif [ "$mat" = "Rubrene+f" ]; then
		kappa="0.0961 0.0147 0.0147"
		dKappa="0.2462 0.4213 0.4213"
	elif [ "$mat" = "C8-BTBT+f" ]; then
		kappa="0.0414 -0.0296 -0.0296"
		dKappa="0.2975 0.8662 0.8662"
	elif [ "$mat" = "TIPS-PN+f" ]; then
		kappa="0.0384 0.0011 0.0000"
		dKappa="0.4832 7.7653 0.0000"
	elif [ "$mat" = "TESADT+f" ]; then
		kappa="-0.1049 -0.0360 0.0000"
		dKappa="0.2538 0.4723 0.0000"
	elif [ "$mat" = "diFTESADT+f" ]; then
		kappa="-0.1173 -0.0293 0.0000"
		dKappa="0.2014 0.4599 0.0000"
	elif [ "$mat" = "Rubrene+f-1D" ]; then
		kappa="0.0961 0. 0."
		dKappa="0.2462 0. 0."
	elif [ "$mat" = "diFTESADT+f-1D" ]; then
		kappa="-0.1173 0. 0."
		dKappa="0.2014 0. 0."
	elif [ "$mat" = "Rubrene+geom" ] || [ "$mat" = "Rubrene+sp" ]; then
		kappa="0.0961 0.0147 0.0147"
		dKappa="0.2462 0.4213 0.4213"
		vec_a="7.2 0"
		vec_b="0 14.3"
	elif [ "$mat" = "C8-BTBT+geom" ]; then
		kappa="0.0414 -0.0296 -0.0296"
		dKappa="0.2975 0.8662 0.8662"
		vec_a="6.1 0"
		vec_b="0 7.4"
	elif [ "$mat" = "TIPS-PN+geom" ]; then
		kappa="0.0384 0.0011 0.0000"
		dKappa="0.4832 7.7653 0.0000"
		vec_a="7.5 0"
		vec_b="0.78 7.7"
	elif [ "$mat" = "TESADT+geom" ]; then
		kappa="-0.1049 -0.0360 0.0000"
		dKappa="0.2538 0.4723 0.0000"
		vec_a="13.7 0"
		vec_b="0 7.2"
	elif [ "$mat" = "diFTESADT+geom" ]; then
		kappa="-0.1173 -0.0293 0.0000"
		dKappa="0.2014 0.4599 0.0000"
		vec_a="6.9 0"
		vec_b="1.14 16.5"
	elif [ "$mat" = "Rubrene-abb" ]; then
		kappa="0.0961 0.0147 0.0147"
		dKappa="0.2462 0.4213 0.4213"
		vec_a="7.2 0"
		vec_b="0 14.3"
	elif [ "$mat" = "C8-BTBT-abb" ]; then
		kappa="0.0414 -0.0296 -0.0296"
		dKappa="0.2975 0.8662 0.8662"
		vec_a="6.1 0"
		vec_b="0 7.4"
	elif [ "$mat" = "C10-DNTT-abb" ]; then
		kappa="0.0659 -0.0375 -0.0375"
		dKappa="0.2100 0.6470 0.6470"
		vec_a="6 0"
		vec_b="0 7.6"
	elif [ "$mat" = "C10-DNBDT-abb" ]; then
		kappa="0.0461 0.0355 0.0355"
		dKappa="0.3800 0.5920 0.5920"
		vec_a="6.1 0"
		vec_b="0 7.8"
	elif [ "$mat" = "C8-DNTT-C8-abb" ]; then
		kappa="0.0554 0.0331 0.0331"
		dKappa="0.1528 0.5489 0.5489"
		vec_a="6 0"
		vec_b="0 7.9"
	elif [ "$mat" = "TIPS-PN-abb" ]; then
		kappa="-0.0005 0.0052 0.0062"
		dKappa="0.569 3.25 5.07"
		vec_a="7.7 0"
		vec_b="0 15.6"
	else
		kappa=$(echo "$hop $theta" | awk '{k = 0.1*$1; th = 0.3141593*$2; print k*cos(th), k*sin(th)/sqrt(2), k*sin(th)/sqrt(2);}')
		dKappa=".5 .5 .5"
		vec_a="7.2 0"
		vec_b="0 12.5"
	fi

	#kappa=$(echo "$kappa" | awk '{print -$1, -$2, -$3}')

	#if (( $Nt < 10 )); then
	#	steps=3
	#elif (( $Nt < 20 )); then
	#	steps=2
	#else
	#	steps=1
	#fi
	steps=2

	cat <<EOF > input.txt
$L $L
$Nt
$beta
$kappa
$dKappa
$mu0
0.006
$steps 0
0
10000 1
$vec_a
$vec_b
EOF
}

#for mat in "Rubrene+geom" "C8-BTBT+geom" "TIPS-PN+geom" "TESADT+geom" "diFTESADT+geom"; do
#for mat in "Rubrene-abb" "C8-BTBT-abb" "C10-DNTT-abb" "C10-DNBDT-abb" "C8-DNTT-C8-abb" "TIPS-PN-abb"; do
for L in 15; do
	for Nt in 64; do
		for mu0 in 0; do
		#for mu0 in -0.15 -0.1 -0.05 0; do
			for beta in 40 30 50 20 60; do
			#for beta in 40; do
				#for hop in 1 0.8 0.6 0.4 0.2 0.05; do
				for hop in 1; do
					for theta in {0..10}; do
					#for theta in 0; do
						#name="L_${L}_${L}_Nt_${Nt}_mu0_${mu0}_b_${beta}_hop_${hop}_th_${theta}"
						name="L_${L}_${L}_Nt_${Nt}_mu0_1_b_${beta}_hop_${hop}_th_${theta}-zf"
						#name="L_${L}_${L}_Nt_${Nt}_mu0_${mu0}_b_${beta}_mat_${mat}"
						#name="L_${L}_${L}_Nt_${Nt}_mu0_1_b_${beta}_mat_${mat}-sp"
						rm -rf $name
						mkdir -p $name
						if [ -e "${name}/results2.csv" ]; then continue; fi

						cp ../simulation/organic_main $name
						#cp ../simulation/organic_conv smear.txt $name
						sed "s/NAME/${name}/g" job_script.sh > ${name}/job_script.sh
						cd $name
						fill_input
						sbatch job_script.sh
						cd ..
					done
				done
			done
		done
	done
done
#done
