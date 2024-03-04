#!/bin/bash

mkdir -p all_rho

for dat in */rho_final.txt; do
	dir=$(echo ${dat%/rho_final.txt})

	rsync $dat all_rho/${dir}.txt

done
