***this is a work in process*** none of these files are intended to be ready for others to use

Intro
====

Building approach to apply ML to CDL-protein binding data, as per CG

Using gnina, as described here:
https://github.com/RMeli/gsoc19
https://github.com/gnina/gnina
https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.6b00740

Note: to open gnina shell need singularity

export SINGULARITY_BINDPATH="/sansom/s156a/bioc1535"
singularity shell --nv /sansom/s156a/bioc1535/gnina/gnina.simg

Approach
====

Build training set from subset of data - use actual CG data.

Scripts
====

STEP 1: find sites

01_find_poses.sh

STEP 2: site energy

02_site_energies.sh

STEP 3: use sites to train a CNN-based scoring function

03_training_data.sh

