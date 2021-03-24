#!/bin/bash

# script to cycle through sites, extract the correct coords and build new systems for FEP

# define input dirs here
DATA=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
SITES=$DATA/FEP/Sites
SETUP=/sansom/s156a/bioc1535/Ecoli_patch/full_complement/chosen
SCRIPT=/sansom/s156a/bioc1535/CDL_site_ML

module load ubuntu-18/gromacs/2018.6_AVX2_plumed-2.4.4

build_system () {
gmx editconf -f $site_dir/Binding_Poses/BSid$2_No$i.gro -o $4/pose.gro -d 2 >& $4/out_files/out_edc1
read -r x y z <<<$(tail -n 1 $4/pose.gro)
python $CG/insane.py -f $4/pose.gro -o $4/$1.$2.$3.mem.gro -x $x -y $y -z $z -l POPE:100 -sol W -p $4/temp.top -center >& $4/out_files/out_mem
}

plumed_dat () {
# plumed used to keep lipid on site during FEP process
txt_file=$DATA/PyLipID_poses/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
site_occ=`sed -n "/Binding site $2$/,/^$/p" $txt_file | grep "BS Lipid Occupancy" | awk '{print $4}'`
cutoff=`echo "scale=4; $site_occ / 2" | bc`
resnum=`sed -n "/Binding site $2$/,/^$/p" $txt_file | tail -n+15 | awk -v var=$cutoff '$6>var' | awk '{print $1}' | tr -d [[:alpha:]] | sed ':a;N;$!ba;s/\n/ /g'`
file=$1.$2.$3
out=$4
echo -e "site_occ $site_occ\ncutoff $cutoff\nresnum $resnum" > $out/site_specs
rm -f $out/ndxs
cp dist.dat $4/dist.dat
set -- "$resnum"
for res in $@
do
	num=`grep " $((res-1))[[:alpha:]]" $out/$file.mem.gro | grep BB | tr -d '[:alpha:]' | awk '{print $2}'`
	sed "s/#PROT#/#PROT#$num,/g" $out/dist.dat -i
done
# specific for CDL
for bead in GL0 PO1 PO2
do
	num=`grep $bead $out/$file.mem.gro | tr -d '[:alpha:]' | awk '{print $3}'`
	sed "s/#LIP#/#LIP#$num,/g" $out/dist.dat -i
done
sed 's/#LIP#//g' $out/dist.dat -i 
sed 's/#PROT#//g' $out/dist.dat -i
sed 's/,$//g' $out/dist.dat -i 
sed "s-#DIR#-$out-g" $out/dist.dat -i
}

build_topology () {
cat $SETUP/$1/topol.top | grep -v -e POPE -e CARD -e POPG -e W -e TNA -e TCL > $4/topol.top
echo CARD 1 >> $4/topol.top
sed "s-../martini-$SETUP/martini-g" $4/topol.top -i
tail -n 5  $4/temp.top | grep -e POPE -e W | sed 's/^W/WN/g' >> $4/topol.top
cp $SETUP/$1/*itp $4/
gmx grompp -f $CG/MDPs/em_memprot.mdp -c $4/$1.$2.$3.mem.gro -r $4/$1.$2.$3.mem.gro -p $4/topol.top -o $4/preion -maxwarn 2 >& $4/out_files/out_preion
echo WN | gmx genion -s $4/preion.tpr -pname TNA -nname TCL -o $4/ions.gro -p $4/topol.top -neutral -conc 0.15 >& $4/out_files/out_genion
echo -e aBB '\n' q | gmx make_ndx -f $4/ions.gro -o $4/BB.ndx >& $4/out_files/out_ndx1
echo -e '\n#ifdef POSRES\n# include "posres.itp"\n#endif' >> $4/Protein.itp
echo BB | gmx genrestr -f $4/ions.gro -n $4/BB.ndx -o $4/posres.itp >& $4/out_files/out_genr
gmx grompp -f $CG/MDPs/em_memprot.mdp -c $4/ions.gro -r $4/ions.gro -p $4/topol.top -o $4/em -maxwarn 2 >& $4/out_files/out_grompp
gmx mdrun -v -deffnm $4/em -plumed $4/dist.dat >& $4/out_files/out_em
tail -n 1 $4/COLVAR | awk '{print $2}' > $4/dist_from_site
}

equil_system () {
echo -e '"CARD" '"|"' "POPE"' '\n' '"PROTEIN" '"|"' "CARD" '"|"' "POPE"' '\n' '"WN" '"|"' "ION"' '\n' q | gmx make_ndx -f $4/em.gro -o $4/sys.ndx >& $4/out_files/out_ndx2
sed 's/CARD_POPE/LIPID/g' $4/sys.ndx -i
sed 's/WN_ION/SOL_ION/g' $4/sys.ndx -i
sed 's/#UPPER/UPPER/g' $4/dist.dat -i
gmx grompp -f $CG/MDPs/eq_2019.mdp -c $4/em.gro -r $4/em.gro -p $4/topol.top -n $4/sys.ndx -o $4/eq >& $4/out_files/out_grompp2
gmx mdrun -v -deffnm $4/eq -nsteps 50000 -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 -plumed $4/dist.dat >& $4/out_files/out_eq 
}

check_eq () {
out=$4
read -r t0 d0 b0 f0 <<<$(grep " 0.000" $out/COLVAR)
read -r t50 d50 b50 f50 <<<$(grep " 5000" $out/COLVAR)
read -r t100 d100 b100 f100 <<<$(grep " 10000" $out/COLVAR)
}

rm -f site_info/chosen.txt

while read pdb
do
	dir=$DATA/FEP/$pdb
	mkdir -p $dir
	site_dir=$DATA/PyLipID_poses2/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD
	total=`(cd $site_dir/Binding_Poses && ls *gro) | tr -d [[:alpha:]]. | awk -F '_' '{print $1}' | sort -u | wc -l`
	for site in `seq 0 $((total-1))`; do
		occ=`grep -A5 "Binding site ${site}$" $site_dir/BindingSites_Info_CARD.txt | tail -n 1 | awk '{print $4}' | awk -F'.' '{print $1}'`
		if [[ $occ -gt 50 ]] ; then
			for i in {0..9} ; do
				out_dir=$dir/$site/$i
				if ! ls $dir/$site/*/eq.gro 1> /dev/null 2>&1 ; then
					mkdir -p $out_dir/out_files
					if [[ ! -f $out_dir/em.gro ]] ; then
						build_system $pdb $site $i $out_dir
						plumed_dat $pdb $site $i $out_dir
						build_topology $pdb $site $i $out_dir
					fi
					dist_from_site=`awk -F '.' '{print $1}' $out_dir/dist_from_site`
					if [[ $dist_from_site -lt 1 ]] ; then
						echo $pdb $site $i >> site_info/chosen.txt
						equil_system $pdb $site $i $out_dir
					fi
					rm -f $out_dir/*#* $SCRIPT/*step*pdb* *mdp
				fi
			done
		fi
	done
done < pdbs.txt

