***this is a work in progress*** 

This is a work in progress, steps 1-5 are now finalised.

Intro
====

Building approach to apply ML to CDL-protein binding data (in coarse grained (CG))

Using gnina, as described here:

https://github.com/RMeli/gsoc19

https://github.com/gnina/gnina

https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.6b00740

Approach
====

Build training set from subset of data - use CG data. Requires Gromacs installation.

Scripts
====

STEP 1: find sites and build poses from existing CG data base. Includes a step where the pose is sanity checked.
```
01_find_poses.sh
```

STEP 2:
Visualise the site - essentially builds a VMD vis script. Requires VMD installation.
```
02_visualise_site.sh
```

STEP 3: run site energy FEPs. Requires Gromacs installation.
```
03_site_FEP.sh
```

Short run script for putting files on local HPC included
```
axon_run_FEP.sh 
```

STEP 4: get FEP data and analyse site FEP values. Requires alchemical analysis installation.
```
04_analyse_FEP.sh
```

STEP 5
Do additional FEP analyses, include Ala scan. Separate version for TRP as virtual sites present. Not essential to workflow.
```
05_ala_scan_res.sh
```

STEP 6: prepare data for ML
Pending trial analysis with CG
```
06_training_data.sh
```

