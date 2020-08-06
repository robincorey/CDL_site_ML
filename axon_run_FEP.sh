#!/bin/bash

for i in `ls /sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_data/Data/`
do
        mkdir -p $i
        cd $i
        if [[ ! -f file.gro ]]
        then
                cp /sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/FEP_data/Data/$i/** .
                sh run.sh
        fi
        cd ../
done
