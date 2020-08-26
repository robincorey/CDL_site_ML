#!/bin/bash
# mutating res instead

DATA=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
FEP=$DATA/FEP/Setup
ALA_DIR=$DATA/FEP/Ala_scan
SCRIPT=/sansom/s156a/bioc1535/CDL_site_ML
CG=/sansom/s137/bioc1535/Desktop/CG_KIT

setup_ala () {
mkdir -p $data
echo $data
cp $coords/topol.top $data/
cp $coords/posres.itp $data/
cp $coords/sys.ndx $data/
cp $coords/em.gro $data/input.gro
if [[ ! -f $data/Protein.itp ]] ; then
	touch $data/Protein.itp
	cat $coords/Protein.itp | while read ;do 
		if [[ "$REPLY" = *$5* ]] && [[ "$REPLY" = *SC* ]] ;then 
			resnum=`echo $REPLY | awk '{print $3}'`
			if  [[ "$resnum" = "$5" ]] ;then
				resname=`echo "$REPLY" | awk '{print $4}'`
				echo "$REPLY" | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "72" Dum "0" "72}' >> $data/Protein.itp 
			else
				echo "$REPLY" >> $data/Protein.itp
			fi
		else
			echo "$REPLY" >> $data/Protein.itp
		fi
	done
fi
}

setup_free () {
mkdir -p $1/out_files
echo $1
sed '/CARD/d' $coords/topol.top > $1/topol.top
cp $data/posres.itp $1/
cp $data/Protein.itp $1/
if [[ ! -f $1/eq.gro ]]; then
gmx editconf -f $coords/em.gro -o $1/input.pdb >& $1/out_files/edc
sed '/CARD/d' $1/input.pdb -i
echo -e '"WN"' '|' '"ION"' '\n' q | gmx make_ndx -f $1/input.pdb -o $1/sys.ndx >& $1/out_files/ndx
sed '0,/POPE/s//LIPID/' $1/sys.ndx -i
sed 's/WN_ION/SOL_ION/g' $1/sys.ndx -i
gmx grompp -f $CG/MDPs/em_memprot.mdp -c $1/input.pdb -r $1/input.pdb -p $1/topol.top -o $1/em -maxwarn 2 >& $1/out_files/g1
gmx mdrun -v -deffnm $1/em >& $1/out_files/m1
gmx grompp -f $CG/MDPs/eq_2019.mdp -c $1/em.gro -r $1/em.gro -p $1/topol.top -o $1/eq -maxwarn 2 -n $1/sys.ndx >& $1/out_files/g2
echo startin eq
gmx mdrun -v -deffnm $1/eq -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 -nsteps 5000000 >& $1/out_files/m2
fi
cp $1/eq.gro $1/input.gro
}

prep_FEP () {
mkdir -p $1/out_files
mkdir -p $1/EM
mkdir -p $1/Data
for i in `seq 0 2 30` ; do
	mdp=$DATA/FEP/Data/$2/$3/$4/EM/EM_$i
	if [[ ! -f $1/EM/em_$i.gro ]]; then
	gmx grompp -f $mdp/em_FEP_${i}.mdp -c $1/input.gro -r $1/input.gro -p $1/topol.top -n $1/sys.ndx -o $1/EM/em_$i.tpr -maxwarn 3 >& $1/out_files/out1
	gmx mdrun -v -deffnm $1/EM/em_$i  >& $1/out_files/emout
	fi
	for run in {1..10} ;do
		echo $run
        	if [ ! -f $1/md_${i}_$run.tpr ] ;then
                        gmx grompp -f $mdp/md_FEP_${i}.mdp -c $1/EM/em_$i.gro -r $1/EM/em_$i.gro -p $1/topol.top -n $1/sys.ndx -o $1/Data/md_${i}_$run.tpr -maxwarn 3 >& $1/out_files/outgrompp
 			#gmx mdrun -v -deffnm $data/md_${i}_$run -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 -nsteps 600000 >& $data/out_files/outmdrun
		fi
	done
done       
for run in {1..10} ; do
        sed "s/#NAME#/$2_$3_$4/g" $FEP/array.sh | sed "s/#REP#/$run/g" > $1/Data/array_$run.sh
done
cp $FEP/run.sh $1/Data/run.sh
}

for site in 5MRW_1_8 ; do
	read -r pdb pose rep <<<$(echo $site | awk -F '_' '{print $1" "$2" "$3}')
	site_file=$DATA/PyLipID_poses/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
	mkdir -p $ALA_DIR/$site/out_files
	mkdir -p $ALA_DIR/$site/res_fep2
	for res in 281 282; do 
		coords=$DATA/FEP/Sites/$pdb/$pose/$rep
		data=$ALA_DIR/$site/res_fep2/$res
		setup_ala $data $pdb $pose $rep $res
		prep_FEP $data $pdb $pose $rep $res
		setup_free ${data}_free $pdb $pose $rep $res
		prep_FEP ${data}_free $pdb $pose $rep $res 

	done
done
