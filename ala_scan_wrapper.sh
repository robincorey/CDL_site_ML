#!/bin/bash

# wrapper script to control ala-scanning FEPs

# step 1: find poses
01_find_poses.sh

# produce tcl script to open and visualise poses in VMD
02_visualise_site.sh

# run FEP of CDL lipid in each site
03_site_FEP.sh

# analyse site FEP
04_analyse_FEP.sh

# ala scanning FEP of each site
05_ala_scan_res.sh
