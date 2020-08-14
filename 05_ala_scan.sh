#!/bin/bash

DATA_DIR=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
ALA_DIR=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_data/ala_scan
SCRIPTS=/sansom/s156a/bioc1535/CDL_site_ML

combine_data () {
echo -e '"Protein" | "CDHG"' '\n' q | gmx make_ndx -f $1/md_0_1.tpr -o $ALA_DIR/$2/prot_lip.ndx >& $ALA_DIR/$2/out_files/outndx
echo Protein_CDHG | gmx editconf -f  $1/md_0_1.tpr -o  $ALA_DIR/$2/md_0_1.tpr.gro -n $ALA_DIR/$2/prot_lip.ndx >& $ALA_DIR/$2/out_files/outedc
echo Protein_CDHG | gmx trjcat -f $1/Data/md_0_*xtc -b 2000 -n $ALA_DIR/$2/prot_lip.ndx -o $ALA_DIR/$2/all.xtc -cat >& $ALA_DIR/$2/out_files/outtrj
}

get_res_contact () {
python $SCRIPTS/lipid-contact_rc.py $ALA_DIR/$2/md_0_1.tpr.gro $ALA_DIR/$2/all.xtc $ALA_DIR/$2
}

get_binding () {
sed 1d $1/contact.csv | while read line
do 
	read -r res contact <<<$(echo $line | awk -F ',' '{print $1" "$2}')
	if [[ `echo "$contact * 100" | bc | cut -f1 -d .` -gt 19 ]]
	then
		mkdir -p $ALA_DIR/$i/res/res_$res
		echo $res > $ALA_DIR/$i/res/res_$res/$res.txt 
		siteocc=`grep $res[[:alpha:]] $2 | awk '{print $6}'`
		echo $res $contact $siteocc
	fi
done 
}

setup_ala () {
:
}

for i in 1FX8_5_1
do
	pdb=`echo $i | cut -f1 -d '_'`
	site_file=$DATA_DIR/PyLipID_poses/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
	mkdir -p $ALA_DIR/$i/out_files
#	combine_data $DATA_DIR/FEP_data/Data/$i $i
#	get_res_contact $DATA_DIR/FEP_data/Data/$i $i
	get_binding $ALA_DIR/$i
	for resi in `ls $ALA_DIR/$i/res/`
	do
		setup_ala $resi $site_file
	done
done
