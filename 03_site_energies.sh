#!/bin/bash

# data dir
SITE=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses

# define site extracting protocol
prep_FEP (){
echo $1 $2 $3
}
	
# loop through PDBs
for pdb in 1FFT
do
	for site in `ls $SITE/Sites_for_ML/$pdb`
        do
                prep_FEP $pdb $site $count
        done
done
