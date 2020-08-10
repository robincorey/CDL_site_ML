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

for i in `ls $DATA_DIR/FEP_data/Data | grep -v -e old -e longer`
do
	mkdir -p $ALA_DIR/$i/out_files
	combine_data $DATA_DIR/FEP_data/Data/$i $i
	get_res_contact $DATA_DIR/FEP_data/Data/$i $i
done
