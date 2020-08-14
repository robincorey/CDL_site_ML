#!/bin/bash

DATA_DIR=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
ALA_DIR=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_data/ala_scan
SCRIPTS=/sansom/s156a/bioc1535/CDL_site_ML
CG=/sansom/s137/bioc1535/Desktop/CG_KIT

combine_data () {
echo -e '"Protein" | "CDHG"' '\n' q | gmx make_ndx -f $1/md_0_1.tpr -o $ALA_DIR/$2/prot_lip.ndx >& $ALA_DIR/$2/out_files/outndx
echo Protein_CDHG | gmx editconf -f  $1/md_0_1.tpr -o  $ALA_DIR/$2/md_0_1.tpr.gro -n $ALA_DIR/$2/prot_lip.ndx >& $ALA_DIR/$2/out_files/outedc
echo Protein_CDHG | gmx trjcat -f $1/Data/md_0_*xtc -b 2000 -n $ALA_DIR/$2/prot_lip.ndx -o $ALA_DIR/$2/all.xtc -cat >& $ALA_DIR/$2/out_files/outtrj
}

get_res_contact () {
python $SCRIPTS/lipid-contact_rc.py $ALA_DIR/$2/md_0_1.tpr.gro $ALA_DIR/$2/all.xtc $ALA_DIR/$2
}

get_binding () {
sed 1d $1/contact.csv | while read line
do
	read -r res contact <<<$(echo $line | awk -F ',' '{print $1" "$2*100}')
	if [[ `echo $contact |cut -f1 -d .` -gt 19 ]]
	then
		siteocc=`grep $res[[:alpha:]] $2 | awk '{print $6}'`
		echo $res $contact $siteocc >> $1/res/res.txt
	fi
done 
}

setup_ala () {
mkdir -p $1/res/$5
echo $1/res/$5
#coords=$DATA_DIR/FEP_data/$2/$3/$4/
coords=$DATA_DIR/Sites_for_ML/$2/$3/$4
cp $coords/topol_FEP.top $1/res/$5/
<<'END'
rm -f > $1/res/$5/Protein.itp
touch $1/res/$5/Protein.itp
cat $coords/Protein.itp | while read 
	do 
	if [[ "$REPLY" = *248* ]]
	then 
		resnum=`echo $REPLY | awk '{print $3}'`
		if  [[ "$resnum" = "248" ]] && [[ "$REPLY" = *SC* ]]
			then
			t=`echo "$REPLY" | awk '{print $2}'`
			echo "$REPLY" | sed "s/$t/Dum/g" | sed 's/1.0000 ;/0.0000 ;/g' >> $1/res/$5/Protein.itp
		else
			echo "$REPLY" >> $1/res/$5/Protein.itp
		fi
	else
		echo "$REPLY" >> $1/res/$5/Protein.itp
	fi
done 
END
}

prep_FEP () {
coords=$DATA_DIR/Sites_for_ML/$2/$3/$4
data=$1/res/$5
echo $data
mkdir -p $data/out_files
for i in `seq 0 2 30`
do
	mdp=$DATA_DIR/FEP_data/$2/$3/$4/EM/EM_$i
	for rep in {1..1}
	do
		echo $rep
        	if [ ! -f $data/md_${i}_$rep.gro ]
        	then
                        gmx grompp -f $mdp/md_FEP_${i}.mdp -c $mdp/em_$rep.gro -r $mdp/em_$rep.gro -p $data/topol_FEP.top -n $coords/sys.ndx -o $data/md_${i}_$rep.tpr -maxwarn 3 >& $data/out_files/outgrompp
 			gmx mdrun -v -deffnm $data/md_${i}_$rep -pin on -pinoffset 0 -ntomp 6 -ntmpi 1 -nsteps 600000 >& $data/out_files/outmdrun
		fi
	done
done       
#rm -f $SCRIPT/*step*pdb*
#for rep in {1..10}
#do
  #      sed "s/#NAME#/$1_$2_$3/g" $FEP/array.sh | sed "s/#REP#/$rep/g" > $data_dir/array_$rep.sh
#done
#cp $FEP/run.sh $data_dir/run.sh
}

for i in 1FX8_5_1
do
	read -r pdb pose rep <<<$(echo $i | awk -F '_' '{print $1" "$2" "$3}')
	site_file=$DATA_DIR/PyLipID_poses/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/BindingSites_Info_CARD.txt
	mkdir -p $ALA_DIR/$i/out_files
	mkdir -p $ALA_DIR/$i/res
	echo "res contact occ" > $ALA_DIR/$i/res/res.txt
#	combine_data $DATA_DIR/FEP_data/Data/$i $i
#	get_res_contact $DATA_DIR/FEP_data/Data/$i $i
	get_binding $ALA_DIR/$i $site_file
	sed 1d $ALA_DIR/$i/res/res.txt | while read -r res contact occ
	do
		setup_ala $ALA_DIR/$i $pdb $pose $rep $res
		prep_FEP $ALA_DIR/$i $pdb $pose $rep $res 
		exit 0
	done
done
