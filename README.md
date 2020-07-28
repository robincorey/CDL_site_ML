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

STEP 1: find sites and build poses
```
01_find_poses.sh
```

STEP 2: run site energy FEPs

```
02_site_energies.sh
```

STEP 4: prepare data for ML

```
03_prep_ML.sh
```

STEP 5: use sites to train a CNN-based scoring function
```
04_training_data.sh
```
