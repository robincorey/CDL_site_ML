import numpy as np
import glob

# data dir
data_dir = "/sansom/s156a/bioc1535/EC_MEMPROT/Data/5us_analyses"

# define site extracting protocol
def get_site (pdb):
	
	return site

# loop through PDBs
for pdb in [ '1FFT' ]:
	file_list = glob.glob('%s/Sites_for_ML/%s/*pdb' % (data_dir, pdb)
        for f in file_list:
		site = get_site(pdb, f)

print(site)
