#!/bin/bash

# data dir
SITE=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
FEP=$SITE/FEP_setup
SCRIPT=/sansom/s156a/bioc1535/CDL_site_ML

# define site extracting protocol
prep_FEP (){
dir=$SITE/FEP_data/$1/$2/$3/EM/
for i in `seq 0 2 30`
do
	build_dir=$dir/EM_$i
	data_dir=$SITE/FEP_data/Data/$1_$2_$3/
	mkdir -p $build_dir/out_files
	mkdir -p $data_dir
	if [ ! -f $build_dir/em_10.gro ]
	then
		sed "s/##INIT##/${i}/g" $FEP/em_FEP.mdp > $build_dir/em_FEP_${i}.mdp
                sed 's/CARD/CDHG/g' $4/topol.top > $4/topol_FEP.top
		sed "/Protein.itp/a #include \"$FEP\/CARD_to_POPC_HG.itp"\" $4/topol_FEP.top -i
		nanumber=`grep '^TNA' $4/topol_FEP.top | awk '{print $2}'`
		sed -e "s/TNA\s\{1,\}$nanumber/TNA $((nanumber-2))/g" $4/topol_FEP.top -i
		sed "/^TNA/a TNAP 2" $4/topol_FEP.top -i
		for rep in {1..10}
		do
			gmx grompp -f $build_dir/em_FEP_${i}.mdp -c $4/frames/eq$rep.pdb -r $4/frames/eq$rep.pdb -p $4/topol_FEP.top -n $4/sys.ndx -o $build_dir/em_$rep.tpr -maxwarn 3 >& $build_dir/out_files/out1
               		gmx mdrun -v -deffnm $build_dir/em_$rep >& $build_dir/out_files/out2
		done
        fi
	if [ ! -f $data_dir/md_${i}_10.tpr ]
        then
		sed "s/##INIT##/${i}/g" $FEP/md_FEP.mdp > $build_dir/md_FEP_${i}.mdp
		for rep in {1..10}
                do
			gmx grompp -f $build_dir/md_FEP_${i}.mdp -c $build_dir/em_$rep.gro -r $build_dir/em_$rep.gro -p $4/topol_FEP.top -n $4/sys.ndx -o $data_dir/md_${i}_$rep.tpr -maxwarn 3 >& $build_dir/out_files/out1
		done
	fi
	rm -f $SCRIPT/*step*pdb*
done
for rep in {1..10}
do
	sed "s/#NAME#/$1_$2_$3/g" $FEP/array.sh | sed "s/#REP#/$rep/g" > $data_dir/array_$rep.sh
done
cp $FEP/run.sh $data_dir/run.sh
}

# loop through PDBs
for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	for site in 0 `ls $SITE/Sites_for_ML/$pdb`
        do
                for i in {0..9}
                do
		build_dir=$SITE/Sites_for_ML/$pdb/$site/$i/
                if [[ -f $build_dir/eq.gro ]]
                then
		#	if [[ ! `ls -d $SITE/FEP_data/Data/${pdb}_${site}_* 2>/dev/null` ]]
                		prep_FEP $pdb $site $i $build_dir
		fi
		done
        done
done
