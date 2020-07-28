Testing gnina installation and setup
Good refs here: 
https://github.com/RMeli/gsoc19
https://github.com/gnina/gnina


To open gnina shell:

export SINGULARITY_BINDPATH="/sansom/s156a/bioc1535"
singularity shell --nv /sansom/s156a/bioc1535/gnina/gnina.simg

Then you're in the shell.

******* 

Approach: Build traning set from subset of data - use actual CG data.

STEP 1: find sites

01_find_poses.sh

STEP 2: site energy

02_site_energies.sh

STEP 3: use sites to train a CNN-based scoring function

03_training_data.sh

