#!/bin/bash
  
# for running locally

SITE=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses

make_load_file () {
echo "
mol new $1/$2/0/frames/1.pdb type pdb
set molID $3
cg_bonds -molid $3 -tpr tprs/$1.bonds.tpr " >> site_info/load_pds.scr
for i in {2..10}
do
        echo "mol addfile $1/$2/0/frames/$i.pdb type pdb" >> site_info/load_pds.scr
done
txt_file=$SITE/Sites_new/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
site_occ=`sed -n "/Binding site $2$/,/^$/p" $txt_file | grep "BS Lipid Occupancy" | awk '{print $4}'`
cutoff=`echo "scale=4; $site_occ / 4" | bc`
resnum=`sed -n "/Binding site $2$/,/^$/p" $txt_file | tail -n+15 | awk -v var=$cutoff '$6>var' | awk '{print $1}' | tr -d [[:alpha:]] | sed ':a;N;$!ba;s/\n/ /g'`
echo -n -e "$1 $2 $site_occ" >> site_info/site_info_$1.txt
set -- "$resnum"
echo -n "mol modselect 1 $3 resid " >> site_info/load_pds.scr
for res in $@
do 
	echo -n "$((res-1)) " >> site_info/load_pds.scr
done
echo "mol modstyle 1 $3 VDW 1.000000 16.000000" >> site_info/load_pds.scr
echo "mol modcolor 1 $3 ResType" >> site_info/load_pds.scr
}

mkdir -p site_info/
mkdir -p tprs/
rm -f site_info/load_pds.scr

count=0
for pdb in 1FFT #1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
        cp /sansom/s156a/bioc1535/EC_MEMPROT/$pdb//md_${pdb}_1.tpr tprs/${pdb}.bonds.tpr
	rm -f site_info/site_info_$pdb.txt
        for site in `ls $SITE/Sites_for_ML/$pdb`
        do
                make_load_file $pdb $site $count
                count=$((count+1))
        done
done
