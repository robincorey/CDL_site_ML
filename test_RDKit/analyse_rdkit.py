from rdkit import Chem
import argparse
from rdkit.Chem import AllChem

# script to convert a list of molecules in SMILES format to SDF 
# works, could do with some more arg options

def readSmiles(input):
	suppl = Chem.SmilesMolSupplier(input)
	smiles_list = []
	for i,mol in enumerate(suppl):
		if mol:
			Chem.SanitizeMol(mol)
			smiles_list.append(mol)
	print("Read %s/%s molecules" % (len(smiles_list), len(suppl)))
	return smiles_list

def Smilestomodel(smiles_list):
	all_models = []
	for mol in smiles_list:
		mol_H = Chem.AddHs(mol)
		addbonds = AllChem.EmbedMolecule(mol_H)
		# using MMFF, could also add options for other methods
		AllChem.MMFFOptimizeMolecule(mol_H)
		all_models.append(mol_H)
	return all_models

def write_SDF(all_models):
	outfile = Chem.SDWriter(args.output)
	for mol in all_models:
		outfile.write(mol)
	print("Wrote %s/%s molecules to %s" % (len(all_models), len(smiles_list), args.output))

if __name__ == '__main__':
	# define arguents here. Add more as needed.
	parser = argparse.ArgumentParser()
	group_input = parser.add_argument_group('input')
	group_input.add_argument("-i", "--input", metavar='filename', type=str, required=True, help="name of SMILES file")
	group_output = parser.add_argument_group('outputs')
	group_output.add_argument("-o", "--output", metavar='filename', required=True, type=str, help="name of output SDF file")
	
	# parse arguments
	args = parser.parse_args()

	# read SMILES
	smiles_list = readSmiles(args.input)
	all_models = Smilestomodel(smiles_list)
	write_SDF(all_models)
