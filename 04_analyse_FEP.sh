#!/bin/bash

DATA_DIR=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_data/Data

analyse_FEP () {
/sansom/s137/bioc1535/.local/bin/alchemical_analysis -d $1 -o $1/Analysis_$2 -q "_${2}.xvg" -p "md_" -i 100000 -t 323 -s 2000
}

extract_info () {
:
}

for i in `ls $DATA_DIR`
do
	fep=$DATA_DIR/$i/Data/
	for rep in {1..10}
	do
		mkdir -p $fep/Analysis_$rep/
		while [ `top -n1 | grep alchem | wc -l` -gt 4 ]
		do 
			sleep 5
		done
		analyse_FEP $fep $rep >& $fep/out &
	done
	
done
