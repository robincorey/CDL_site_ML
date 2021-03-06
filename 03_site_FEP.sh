#!/bin/bash

# define data dirs here
DATA=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
FEP=$DATA/FEP_setup
SCRIPT=/sansom/s156a/bioc1535/CDL_site_ML

# define site extracting protocol
get_frames () {
echo -e PROTEIN '\n' CARD | gmx mindist -f $4/eq.xtc -od $4/eq_dist.xvg >& $4/out_files/out_mindist
mkdir -p $4/frames/
echo SYSTEM | gmx trjconv -f $4/eq.xtc -s $4/eq -b 50000 -skip 25 -sep -o $4/frames/eq.pdb >& $4/out_files/out_frame_$frame
rm -f $4/frames/*#*
}

prep_FEP (){
dir=$DATA/FEP_data/$1/$2/$3/EM/
for i in `seq 0 2 30`
do
	build_dir=$dir/EM_$i
	data_dir=$DATA/FEP_data/Data/$1_$2_$3/
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
			gmx grompp -f $build_dir/md_FEP_${i}.mdp -c $build_dir/em_$rep.gro -r $build_dir/em_$rep.gro -p $4/topol_FEP.top -n $4/sys.ndx -o $data_dir/md_${i}_$rep.tpr -maxwarn 3 >& $build_dir/out_files/out3
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

while read pdb
do
	for site in 0 `ls $DATA/Sites_for_ML/$pdb`
        do
                for i in {0..9}
                do
		build_dir=$DATA/Sites_for_ML/$pdb/$site/$i/
                if [[ -f $build_dir/eq.gro ]]
                then
        		get_frames $pdb $site $i $build_dir
			prep_FEP $pdb $site $i $build_dir
		fi
		done
        done
done < pdbs.txt
