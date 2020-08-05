#!/bin/bash
  
# makes vmd file for running locally

SITE=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses

make_load_file () {
echo "mol new $1/$2/0/frames/1.pdb type pdb
set molID $3
mol rename $3 {$1 $2}
mol off $3" >> site_info/load_pds.scr
for i in {2..10}
do
        echo "mol addfile $1/$2/0/frames/$i.pdb type pdb" >> site_info/load_pds.scr
done
echo -n "mol modselect 0 $3 name BB
mol modstyle 0 $3 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol modcolor 0 $3 ColorID 6
mol material AOChalky
mol addrep $3
mol modselect 1 $3 resname CARD
mol modstyle 1 $3 VDW 1.500000 12.000000 
mol modcolor 1 $3 Name
mol modmaterial 1 $3 AOEdgy
mol addrep $3
mol modstyle 2 $3 VDW 1.500000 12.000000
mol modcolor 2 $3 ResType
mol modmaterial 2 $3 AOEdgy
mol modselect 2 $3 name BB SC1 SC2 SC3 SC4 SC5 and resid " >> site_info/load_pds.scr
txt_file=$SITE/Sites_new/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
site_occ=`sed -n "/Binding site $2$/,/^$/p" $txt_file | grep "BS Lipid Occupancy" | awk '{print $4}'`
cutoff=`echo "scale=4; $site_occ / 2" | bc`
resnum=`sed -n "/Binding site $2$/,/^$/p" $txt_file | tail -n+15 | awk -v var=$cutoff '$6>var' | awk '{print $1}' | tr -d [[:alpha:]] | sed ':a;N;$!ba;s/\n/ /g'`
echo -n -e "$1 $2 $site_occ" >> site_info/site_info_$1.txt
set -- "$resnum"
for res in $@
do 
	echo -n "$((res-1)) " >> site_info/load_pds.scr
done
echo -e "\n" >> site_info/load_pds.scr
}

mkdir -p site_info/
#mkdir -p tprs/
rm -f site_info/load_pds.scr
echo -e "color Name D green\ncolor Resname GLY red\ncolor Resname ARG blue\ncolor Resname LYS blue\ncolor Resname THR green\ncolor Resname PRO green\ncolor Resname SER green\ncolor Resname PHE orange\ncolor Resname TYR orange\ncolor Resname TRP orange" > site_info/load_pds.scr

count=0
for pdb in 1FFT #1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
#        cp $SITE/Sites_for_ML/$pdb/0/0/em.tpr tprs/${pdb}.bonds.tpr
	rm -f site_info/site_info_$pdb.txt
        for site in `ls $SITE/Sites_for_ML/$pdb`
        do
                make_load_file $pdb $site $count
                count=$((count+1))
        done
done
