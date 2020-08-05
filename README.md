***this is a work in process*** 
none of these files are intended to be ready for others to use

Intro
====

Building approach to apply ML to CDL-protein binding data, as per CG

Using gnina, as described here:
https://github.com/RMeli/gsoc19
https://github.com/gnina/gnina
https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.6b00740

Approach
====

Build training set from subset of data - use actual CG data.

Scripts
====

STEP 1: find sites and build poses. 
Includes a step where the pose is sanity checked.
```
01_find_poses.sh
```

STEP 2:
Visualise the site - essentially builds a VMD vis script
```
02_visualise_site.sh
```

STEP 3: run site energy FEPs
```
03_site_energies.sh
```

STEP 4: prepare data for ML
Gathers the residues in contact with the lipid at lamda 0 and the energy of the site
```
04_prep_ML.sh
```

STEP 5: use sites to train a CNN-based scoring function
```
05_training_data.sh
```
