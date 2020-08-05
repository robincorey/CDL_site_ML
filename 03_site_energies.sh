#!/bin/bash

# data dir
SITE=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
FEP=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_setup

# define site extracting protocol
prep_FEP (){
dir=$SITE/FEP_data/$1/$2/$3/EM/
for i in `seq 0 2 0`
do
	em_dir=$dir/EM_$i
	mkdir -p $em_dir/out_files
	if [ ! -f $em_dir/EM_${i}.gro ]
	then
		sed "s/##INIT##/${i}/g" $FEP/em_FEP.mdp > $em_dir/em_FEP_${i}.mdp
                sed 's/CARD/CDHG/g' $4/topol.top > $4/topol_FEP.top
		ls $4/topol_FEP.top
		sed "/Protein.itp/a #include \"$FEP\/CARD_to_POPC_HG.itp"\" $4/topol_FEP.top -i
		nanumber=`grep '^TNA' $4/topol_FEP.top | awk '{print $2}'`
		sed -e "s/TNA\s\{1,\}$nanumber/TNA $((nanumber-2))/g" $4/topol_FEP.top -i
		sed "/^TNA/a TNAP 2" $4/topol_FEP.top -i
		gmx grompp -f em_FEP_${i}.mdp -c $4/eq.gro -r $4/eq.gro -p $4/topol.top -n $4/sys.ndx >& $em_dir/out_files/out1
                #gmx mdrun -v -deffnm EM_${i}
        fi
done
}
	
mkdir -p $SITE/FEP_data

# loop through PDBs
for pdb in 1FFT
do
	for site in `ls $SITE/Sites_for_ML/$pdb`
        do
                for i in {0..9}
                do
		build_dir=$SITE/Sites_for_ML/$pdb/$site/$i/
                if [[ -f $build_dir/eq.gro ]]
                then
                	prep_FEP $pdb $site $i $build_dir
		fi
		exit 0
		done
        done
done
