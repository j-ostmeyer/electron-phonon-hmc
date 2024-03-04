#!/bin/bash

for name in *_hop_*-tl/; do
	cp smear.txt $name
	sed "s/NAME/${name%/}/g" job_script.sh > ${name}/job_script.sh
	cd $name
	rm -f greens_smeared.csv
	sbatch job_script.sh
	cd ..
done
