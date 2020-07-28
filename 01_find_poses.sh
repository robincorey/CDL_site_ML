#!/bin/bash

# data dir
CD=/sansom/s156a/bioc1535/EC_MEMPROT/Data/5us_analyses

mkdir -p $CD/Sites_for_ML

build_system () {
for i in {0..1} #9}
do
	out_dir=$CD/Sites_for_ML/$1/$2/$i
	mkdir -p $out_dir
	rm -f $out_dir/*.*
	mkdir -p $out_dir/out_files
	gmx editconf -f $site_dir/BSid$2_No$i.gro -o $out_dir/pose.gro -d 2 >& $out_dir/out_files/out_edc1
	read -r x y z <<<$(tail -n 1 $out_dir/pose.gro)
	python $CG/insane.py -f $out_dir/pose.gro -o $out_dir/${1}.$2.$i.mem.gro -x $x -y $y -z $z -l POPE:100 -sol W -p $out_dir/temp.top -center >& $out_dir/out_files/out_mem
done
}

# loop through pdbs and extract sites
for pdb in 1FFT
do
	mkdir -p $CD/Sites_for_ML/$pdb
	# this is wrong, of course
	site_dir=$CD/test_PyLipID_sites/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/Binding_Poses
	# gets the range of sites
	total=`(cd $site_dir && ls *gro) | tr -d [[:alpha:]]. | awk -F '_' '{print $1}' | sort -u | wc -l`
	for site in `seq 0 1` #$((total-1))`
	do
		build_system $pdb $site 
	done
done
