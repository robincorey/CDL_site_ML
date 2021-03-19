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
-rw-r--r--  1 bioc1535 sansom    0 Jul 29  2020 06_training_data.sh
-rwxr--r--  1 bioc1535 sansom 2.8K Aug 19  2020 03_site_FEP_test_lamda.sh
-rwxr--r--  1 bioc1535 sansom 4.3K Aug 30  2020 05_ala_scan_res.sh
-rwxr--r--  1 bioc1535 sansom 2.6K Sep 22 11:22 02_visualise_site.sh
-rwxr--r--  1 bioc1535 sansom 4.5K Feb  8 09:48 05_ala_scan_res.TRP.sh
-rwxr-xr-x+ 1 bioc1535 sansom 5.0K Mar 19 15:20 01_find_poses.sh
-rwxr--r--+ 1 bioc1535 sansom 2.7K Mar 19 15:20 03_site_FEP.sh
-rwxr--r--+ 1 bioc1535 sansom 1.7K Mar 19 15:20 04_analyse_FEP.sh

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

STEP 3: run site energy FEPs. 
```
03_site_FEP.sh
03_site_FEP_test_lamda.sh
```

Short run script for putting files on local HPC included
```
axon_run_FEP.sh 
```

STEP 4: get FEP data and analyse site FEP values
dependendcy on 'Alchemical analysis'
```
04_analyse_FEP.sh
```

STEP 5
Do additional FEP analyses, include Ala scan. Separate version for TRP as virtual sites present.
```
05_ala_scan_res.sh
05_ala_scan_res.TRP.sh
```

STEP 6: prepare data for ML
Pending trial analysis with CG
```
06_training_data.sh
```

