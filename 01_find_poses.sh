#!/bin/bash

# data dir
CD=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
SETUP=/sansom/s156a/bioc1535/Ecoli_patch/full_complement/chosen
SCRIPT=/sansom/s156a/bioc1535/CDL_site_ML

# functions
build_system () {
gmx editconf -f $site_dir/Binding_Poses/BSid$2_No$i.gro -o $4/pose.gro -d 2 >& $4/out_files/out_edc1
read -r x y z <<<$(tail -n 1 $4/pose.gro)
python $CG/insane.py -f $4/pose.gro -o $4/$1.$2.$3.mem.gro -x $x -y $y -z $z -l POPE:100 -sol W -p $4/temp.top -center >& $4/out_files/out_mem
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
gmx mdrun -v -deffnm $4/em >& $4/out_files/out_em
}

# this section gets the distance from th card to the site - and filters less than 1 nm
check_site () {
txt_file=$CD/Sites_new/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
site_occ=`sed -n "/Binding site $2$/,/^$/p" $txt_file | grep "BS Lipid Occupancy" | awk '{print $4}'`
cutoff=`echo "scale=4; $site_occ / 2" | bc`
resnum=`sed -n "/Binding site $2$/,/^$/p" $txt_file | tail -n+15 | awk -v var=$cutoff '$6>var' | awk '{print $1}' | tr -d [[:alpha:]] | sed ':a;N;$!ba;s/\n/ /g'`
out=$4
echo -e "site_occ $site_occ\ncutoff $cutoff\nresnum $resnum" > $out/site_specs
rm -f $out/ndxs
set -- "$resnum"
for res in $@
do
        echo -n "ri$((res-1)) " >>  $out/ndxs
done
sed 's/ ri/ | ri/g' $out/ndxs -i
echo "
aGL0 | aPO1 | aPO2
q" >> $out/ndxs
gmx make_ndx -f $out/em.gro -o $out/site.ndx < $out/ndxs >& $out/out_files/out_sitendx
site_ndx=`grep '\[ r_' $out/site.ndx | sed 's/\[//g' | sed 's/\]//g'`
sed "s/$site_ndx/ Site /g" $out/site.ndx -i
card_ndx=`grep '\[ GL0' $out/site.ndx | sed 's/\[//g' | sed 's/\]//g'`
sed "s/$card_ndx/ CARD_HG /g" $out/site.ndx -i
gmx distance -f $out/em.gro -s $out/em.tpr -oxyz $out/site_CL.xvg -select 'cog of group "Site" plus cog of group "CARD_HG"' -rmpbc no -pbc no -n $out/site.ndx >& $out/out_files/out_dist
read -r x y z <<<$(tail -n 1 $out/site_CL.xvg)
echo "scale = 4; sqrt($x^2 + $y^2)" | bc > $out/dist_from_site
}

equil_system () {
echo -e '"CARD" '"|"' "POPE"' '\n' '"PROTEIN" '"|"' "CARD" '"|"' "POPE"' '\n' '"WN" '"|"' "ION"' '\n' q | gmx make_ndx -f $4/em.gro -o $4/sys.ndx >& $4/out_files/out_ndx2
sed 's/CARD_POPE/LIPID/g' $4/sys.ndx -i
sed 's/WN_ION/SOL_ION/g' $4/sys.ndx -i
gmx grompp -f $CG/MDPs/eq_2019.mdp -c $4/em.gro -r $4/em.gro -p $4/topol.top -n $4/sys.ndx -o $4/eq >& $4/out_files/out_grompp2
gmx mdrun -v -deffnm $4/eq -nsteps 500000 -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 >& $4/out_files/out_eq # << needs to be longer
}

get_frames () {
echo -e PROTEIN '\n' CARD | gmx mindist -f $4/eq.xtc -od $4/eq_dist.xvg >& $4/out_files/out_mindist
mkdir -p $4/frames/
for frame in {1..10}
do
	# times need to be tweaked
	echo SYSTEM | gmx trjconv -f $4/eq.xtc -s $4/eq.tpr -dump $((frame*1000)) -o $4/frames/$frame.pdb >& $4/out_files/out_frame_$frame
done
rm -f $4/frames/*#*
}

mkdir -p $CD/Sites_for_ML
rm -f site_info/chosen.txt

# loop through pdbs and extract sites
for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	mkdir -p $CD/Sites_for_ML/$pdb
	# this is wrong, of course
	site_dir=$CD/Sites_new/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD
	# gets the range of sites
	total=`(cd $site_dir/Binding_Poses && ls *gro) | tr -d [[:alpha:]]. | awk -F '_' '{print $1}' | sort -u | wc -l`
	for site in `seq 0 $((total-1))`
	do
		occ=`grep -A5 "Binding site ${site}$" $site_dir/BindingSites_Info_CARD.txt | tail -n 1 | awk '{print $4}' | awk -F'.' '{print $1}'`
		if [[ $occ -gt 50 ]]
		then
			for i in {0..9} 
			do
				#echo starting $pdb $site $i
				out_dir=$CD/Sites_for_ML/$pdb/$site/$i
				mkdir -p $out_dir
				mkdir -p $out_dir/out_files
				#build_system $pdb $site $i $out_dir
				#build_topology $pdb $site $i $out_dir
                                #check_site $pdb $site $i $out_dir
				dist_from_site=`awk -F '.' '{print $1}' $out_dir/dist_from_site`
				if [[ $dist_from_site -lt 1 ]] 
				then
					echo doing $pdb $site $i
					#echo $pdb $site $i >> site_info/chosen.txt
					#equil_system $pdb $site $i $out_dir
					#get_frames $pdb $site $i $out_dir
				fi
				rm -f $out_dir/*#*
				rm -f $SCRIPT/*step*pdb* *mdp
			done
		fi
	done
done

