#!/bin/bash
  
# makes vmd file for running locally

DATA=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses

make_load_file () {
echo -n "mol new sites/$1.$2.$4.gro type gro
set molID $3
mol rename $3 {$1 $2}
mol off $3
mol modselect 0 $3 name BB SC1 SC2 SC3 SC4 SC5
mol modstyle 0 $3 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol modcolor 0 $3 ColorID 8
mol modmaterial 0 $3 Glass
mol addrep $3
cg_bonds -molid $3 -tpr tprs/bonds_$1.$2.$4.tpr
mol modselect 1 $3 resname CARD
mol modstyle 1 $3 Licorice 0.600000 12.000000 12.000000
mol modcolor 1 $3 Name
mol modmaterial 1 $3 AOEdgy
mol addrep $3
mol modstyle 2 $3 VDW 1.500000 12.000000
mol modcolor 2 $3 ResName
mol modmaterial 2 $3 AOEdgy
mol modselect 2 $3 name BB SC1 SC2 SC3 SC4 SC5 and resid " >> site_info/load_pds.scr
site_occ=`grep site_occ $DATA/FEP/$1/$2/$4/site_specs | tr -d [[:alpha:]]_`
cutoff=`grep cutoff $DATA/FEP/$1/$2/$4/site_specs | tr -d [[:alpha:]]`
resnum=`grep resnum $DATA/FEP/$1/$2/$4/site_specs | tr -d [[:alpha:]]`
echo -n -e "$1 $2 $4 $site_occ" >> site_info/site_info_$1.txt
set -- "$resnum"
for res in $@
do 
	echo -n "$((res-1)) " >> site_info/load_pds.scr
done
echo -e "\n" >> site_info/load_pds.scr
}

mkdir -p site_info/sites/
rm -f site_info/load_pds.scr
echo -e "display rendermode GLSL\ncolor Name D green\ncolor Resname GLY red\ncolor Resname ARG blue\ncolor Resname LYS blue\ncolor Resname THR green\ncolor Resname PRO green\ncolor Resname SER green\ncolor Resname PHE orange\ncolor Resname TYR orange\ncolor Resname TRP orange\ncolor Resname HIS cyan\ncolor Resname LEU yellow\ncolor Resname ILE yellow\ncolor Resname VAL yellow\ncolor Resname ASP pink\ncolor Resname GLU pink\ncolor Resname ASN green\ncolor Resname GLN green " > site_info/load_pds.scr

count=0
module unload gromacs/2019.4-AVX2-GPU
module load gromacs/5.0.2/64
mkdir -p tprs/
mkdir -p sites

while read pdb
do
	rm -f site_info/site_info_$pdb.txt
	for site in `ls $DATA/FEP/$pdb`
        do
		for i in {0..9}
		do
		dir=$DATA/FEP/$pdb/$site/$i
		if [[ -f $dir/eq.gro ]]
		then
			echo $pdb $site $i
			grompp_avx -f /sansom/s137/bioc1535/Desktop/CG_KIT/MDPs/em_memprot.mdp -c $dir/eq.gro -p $dir/topol.top -o tprs/bonds_${pdb}.$site.$i.tpr
			editconf_avx -f $dir/eq.gro -o sites/$pdb.$site.$i.gro 
			make_load_file $pdb $site $count $i
			count=$((count+1))
		fi
		done
        done
rm -f *#* tprs/*#* sites/*#*
done < pdbs.txt
