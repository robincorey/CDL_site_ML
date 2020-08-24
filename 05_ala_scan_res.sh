#!/bin/bash
# mutating res instead

DATA=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
FEP=$DATA/FEP/Setup
ALA_DIR=$DATA/FEP/Ala_scan
SCRIPT=/sansom/s156a/bioc1535/CDL_site_ML
CG=/sansom/s137/bioc1535/Desktop/CG_KIT

combine_data () {
echo -e '"Protein" | "CDHG"' '\n' q | gmx make_ndx -f $1/md_0_1.tpr -o $ALA_DIR/$2/prot_lip.ndx >& $ALA_DIR/$2/out_files/outndx
echo Protein_CDHG | gmx editconf -f  $1/md_0_1.tpr -o  $ALA_DIR/$2/md_0_1.tpr.gro -n $ALA_DIR/$2/prot_lip.ndx >& $ALA_DIR/$2/out_files/outedc
echo Protein_CDHG | gmx trjcat -f $1/Data/md_0_*xtc -b 2000 -n $ALA_DIR/$2/prot_lip.ndx -o $ALA_DIR/$2/all.xtc -cat >& $ALA_DIR/$2/out_files/outtrj
}

get_res_contact () {
python $SCRIPT/lipid-contact_rc.py $ALA_DIR/$2/md_0_1.tpr.gro $ALA_DIR/$2/all.xtc $ALA_DIR/$2
}

get_binding () {
echo "res contact occ" > $1/res/res.txt
sed 1d $1/contact.csv | while read line
do
	read -r res contact <<<$(echo $line | awk -F ',' '{print $1" "$2*100}')
	if [[ `echo $contact |cut -f1 -d .` -gt 29 ]]
	then
		siteocc=`grep $res[[:alpha:]] $2 | awk '{print $6}'`
		echo $res $contact $siteocc >> $1/res/res.txt
	fi
done 
}

setup_ala () {
mkdir -p $1/res/$5
echo $1/res/$5
coords=$DATA/FEP/Sites/$2/$3/$4
cp $coords/topol_FEP.top $1/res/$5/
cp $coords/posres.itp $1/res/$5/
if [[ ! -f $1/res/$5/Protein.itp ]] ; then
	touch $1/res/$5/Protein.itp
	cat $coords/Protein.itp | while read ;do 
		if [[ "$REPLY" = *$5* ]] ;then 
			resnum=`echo $REPLY | awk '{print $3}'`
			if  [[ "$resnum" = "$5" ]] && [[ "$REPLY" = *SC* ]] ;then
				t=`echo "$REPLY" | awk '{print $2}'`
				echo "$REPLY" | sed "s/$t/Dum/g" | sed 's/1.0000 ;/0.0000 ;/g' >> $1/res/$5/Protein.itp
			else
				echo "$REPLY" >> $1/res/$5/Protein.itp
			fi
		else
			echo "$REPLY" >> $1/res/$5/Protein.itp
		fi
	done
fi
}

equil_pose () {
mkdir -p $data/out_files
if [[ ! -f $data/eq.gro ]] ; then
	module load ubuntu-18/gromacs/2018.6_AVX2_plumed-2.4.4
	gmx grompp -f $CG/MDPs/eq_2019.mdp -c $coords/em.gro -r $coords/em.gro -p $data/topol_FEP.top -n $coords/sys.ndx -o $data/eq -maxwarn 2 >& $data/out_files/out_grompp2
	gmx mdrun -v -deffnm $data/eq -nsteps 5000000 -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 -plumed $coords/dist.dat >& $data/out_files/out_eq
fi
}

get_frames () {
mkdir -p $data/frames/
echo SYSTEM | gmx trjconv -f $data/eq.xtc -s $data/eq -b 50000 -skip 25 -sep -o $data/frames/eq.pdb >& $data/out_files/out_frame
rm -f $data/frames/*#*
}

prep_FEP () {
mkdir -p $data/out_files
for i in `seq 0 2 30` ; do
	mdp=$DATA/FEP/Data/$2/$3/$4/EM/EM_$i
	for run in {1..10} ;do
		echo $run
        	if [ ! -f $data/md_${i}_$run.tpr ] ;then
                        gmx grompp -f $mdp/md_FEP_${i}.mdp -c $data/frames/eq$run.pdb -r $data/frames/eq$run.pdb -p $data/topol_FEP.top -n $coords/sys.ndx -o $data/md_${i}_$run.tpr -maxwarn 3 >& $data/out_files/outgrompp
 	#		gmx mdrun -v -deffnm $data/md_${i}_$run -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 -nsteps 600000 >& $data/out_files/outmdrun
		fi
	done
done       
for run in {1..10} ; do
        sed "s/#NAME#/$2_$3_$4/g" $FEP/array.sh | sed "s/#REP#/$run/g" > $data/array_$run.sh
done
cp $FEP/run.sh $data/run.sh
}

for site in 5MRW_1_8 ; do
	read -r pdb pose rep <<<$(echo $site | awk -F '_' '{print $1" "$2" "$3}')
	site_file=$DATA/PyLipID_poses/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
	mkdir -p $ALA_DIR/$site/out_files
	mkdir -p $ALA_DIR/$site/res
#	combine_data $DATA/FEP_data/Data/$site $site
#	get_res_contact $DATA/FEP_data/Data/$site $site
#	get_binding $ALA_DIR/$site $site_file
#	sed 1d $ALA_DIR/$site/res/res.txt | while read -r res contact occ
	for res in 268 277 278 283 285 523 ; do # 2
		coords=$DATA/FEP/Sites/$pdb/$pose/$rep
		data=$ALA_DIR/$site/res/$res
		setup_ala $ALA_DIR/$site $pdb $pose $rep $res
		equil_pose $ALA_DIR/$site $pdb $pose $rep $res
		get_frames $ALA_DIR/$site $pdb $pose $rep $res
		prep_FEP $ALA_DIR/$site $pdb $pose $rep $res 
	done
done
