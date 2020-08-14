#!/bin/bash

DATA_DIR=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses

analyse_FEP () {
mkdir -p $1/Analysis_$2
if [[ ! -f $1/Analysis_$2/results.txt ]]
then
	/sansom/s137/bioc1535/.local/bin/alchemical_analysis -d $1 -o $1/Analysis_$2 -q "_${2}.xvg" -p "md_" -i 100000 -t 323 -s 2000 
fi
tail -n 1 $1/Analysis_$2/results.txt | awk '{print 93.69-$(NF-2)}' >> $DATA_DIR/FEP_data/energies/$3.txt
}

get_other_values () {
site_dir=$DATA_DIR/PyLipID_poses/$1/lipid_interactions/Interaction_CARD/Binding_Sites_CARD
read -r koff koff_err <<<$(grep -A10 "Binding site $2$" $site_dir/BindingSites_Info_CARD.txt | grep "BS koff Bootstrap" | awk '{print $4" "$6}') 
occ=`grep -A10 "Binding site $2$" $site_dir/BindingSites_Info_CARD.txt | grep "BS Lipid Occupancy" | awk '{print $4}'`
fepval=`tail -n+2 $DATA_DIR/FEP_data/energies/$1_$2_$3.txt | awk '{sum+=$1} END {print (sum/NR)}'`
fep_error=`tail -n+2 $DATA_DIR/FEP_data/energies/$1_$2_$3.txt | awk '{sum+=$1;a[NR]=$1}END{for(i in a)y+=(a[i]-(sum/NR))^2;print sqrt(y/(NR-1))}'`
echo -e "$fepval $fep_error $koff $koff_err $occ 0 " >> $DATA_DIR/FEP_data/energies/comparison.txt 
}

mkdir -p $DATA_DIR/FEP_data/energies
echo "fep fep_error koff koff_error occ" > $DATA_DIR/FEP_data/energies/comparison.txt

for i in `ls $DATA_DIR/FEP_data/Data | grep -v -e old -e longer`
do
	fep=$DATA_DIR/FEP_data/Data/$i/Data
	if [[ -d $fep ]]
	then
		echo $i
		echo $i > $DATA_DIR/FEP_data/energies/$i.txt
		for rep in {1..10}
		do
			while [ `top -n1 | grep alchem | wc -l` -gt 4 ]
			do
				sleep 5
			done	
			analyse_FEP $fep $rep $i >& $fep/out &
		done
		get_other_values `echo $i | cut -f1 -d'_'` `echo $i | cut -f2 -d'_'` `echo $i | cut -f3 -d'_'`
	fi
done
