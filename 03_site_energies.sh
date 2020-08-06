#!/bin/bash

# data dir
SITE=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
FEP=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_setup

# define site extracting protocol
prep_FEP (){
dir=$SITE/FEP_data/$1/$2/$3/EM/
for i in `seq 0 2 30`
do
	em_dir=$dir/EM_$i
	mkdir -p $em_dir/out_files
	if [ ! -f $em_dir/em_10.gro ]
	then
		sed "s/##INIT##/${i}/g" $FEP/em_FEP.mdp > $em_dir/em_FEP_${i}.mdp
                sed 's/CARD/CDHG/g' $4/topol.top > $4/topol_FEP.top
		sed "/Protein.itp/a #include \"$FEP\/CARD_to_POPC_HG.itp"\" $4/topol_FEP.top -i
		nanumber=`grep '^TNA' $4/topol_FEP.top | awk '{print $2}'`
		sed -e "s/TNA\s\{1,\}$nanumber/TNA $((nanumber-2))/g" $4/topol_FEP.top -i
		sed "/^TNA/a TNAP 2" $4/topol_FEP.top -i
		for rep in {1..10}
		do
			gmx grompp -f $em_dir/em_FEP_${i}.mdp -c $4/frames/$rep.pdb -r $4/frames/$rep.pdb -p $4/topol.top -n $4/sys.ndx -o $em_dir/em_$rep.tpr -maxwarn 3 >& $em_dir/out_files/out1
               		gmx mdrun -v -deffnm $em_dir/em_$rep >& $em_dir/out_files/out2
		done
        fi
	md_dir=$dir/EM_$i
	mkdir -p $md_dir/out_files
	data=$SITE/FEP_data/Data/$1_$2_$3/
	mkdir -p $data
	if [ ! -f $data/array_10.sh ]
        then
		sed "s/##INIT##/${i}/g" $FEP/md_FEP.mdp > $em_dir/md_FEP_${i}.mdp
		for rep in {1..10}
                do
			gmx grompp -f $md_dir/md_FEP_${i}.mdp -c $md_dir/em_$rep.gro -r $md_dir/em_$rep.gro -p $4/topol.top -n $4/sys.ndx -o $data/md_${i}_$rep.tpr -maxwarn 3 >& $md_dir/out_files/out1
			sed "s/#NAME#/$1_$2_$3/g" $FEP/array.sh | sed "s/#REP#/$rep/g" > $data/array_$rep.sh
		done
	fi
done
cp $FEP/run.sh $data/run.sh
}

# loop through PDBs
for pdb in 1FFT
do
	for site in 0 # `ls $SITE/Sites_for_ML/$pdb`
        do
                for i in 0 #{0..9}
                do
		build_dir=$SITE/Sites_for_ML/$pdb/$site/$i/
                if [[ -f $build_dir/eq.gro ]]
                then
                	prep_FEP $pdb $site $i $build_dir
		fi
		done
        done
done
