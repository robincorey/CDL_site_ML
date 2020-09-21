#!/bin/bash
  
# makes vmd file for running locally

DATA=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses

make_load_file () {
echo -n "mol new sites/$1.$2.$4.gro type gro
set molID $3
mol rename $3 {$1 $2}
mol off $3
mol modselect 0 $3 name BB
mol modstyle 0 $3 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol modcolor 0 $3 ColorID 6
mol modmaterial 0 $3 AOChalky
mol addrep $3
cg_bonds -molid $3 -tpr $1.bonds.tpr
mol modselect 1 $3 resname CARD
mol modstyle 1 $3 VDW 1.500000 12.000000 
mol modcolor 1 $3 ColorID 4 
mol modmaterial 1 $3 Transparent
mol addrep $3
mol modstyle 2 $3 VDW 1.500000 12.000000
mol modcolor 2 $3 Resname
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
echo -e "color Name D green\ncolor Resname GLY red\ncolor Resname ARG blue\ncolor Resname LYS blue\ncolor Resname THR green\ncolor Resname PRO green\ncolor Resname SER green\ncolor Resname PHE orange\ncolor Resname TYR orange\ncolor Resname TRP orange\ncolor Resname HIS cyan\ncolor Resname LEU yellow\ncolor Resname ILE yello\ncolor Resname ASP pink\ncolor Resname GLU pink\ncolor Resname ASN green\ncolor Resname GLN green " > site_info/load_pds.scr

count=0
module unload gromacs/2019.4-AVX2-GPU
module load gromacs/5.0.2/64
mkdir -p tprs/
mkdir -p sites
#for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6
for pdb in 5MRW #5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	in=`find /sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP/Ala_scan/$pdb* -name "input.gro" | head -n 1 | awk -F 'input' '{print $1}'`
	out=
	grompp_avx -f /sansom/s137/bioc1535/Desktop/CG_KIT/MDPs/em_memprot.mdp -c $in/input.gro -p $in/topol.top -o tprs/bonds_${pdb}
#	editconf_avx -f $in/input.gro -o sites/sites_${pdb}.pdb
	rm -f site_info/site_info_$pdb.txt
	for site in `ls $DATA/FEP/$pdb`
        do
		for i in {0..9}
		do
		if [[ -f $DATA/FEP/$pdb/$site/$i/eq.gro ]]
		then
			echo $pdb $site $i
			editconf_avx -f $DATA/FEP/$pdb/$site/$i/pose.gro site_info/sites/$pdb.$site.$i.gro 
			#make_load_file $pdb $site $count $i
			count=$((count+1))
		fi
		done
        done
done
