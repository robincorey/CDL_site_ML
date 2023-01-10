from rdkit import Chem
import argparse
from rdkit.Chem import AllChem

# script to convert a list of molecules in SMILES format to SDF 

def readSmiles(input):
	suppl = Chem.SmilesMolSupplier(input)
	smiles_list = []
	for i,mol in enumerate(suppl):
		if mol:
			Chem.SanitizeMol(mol)
			smiles_list.append(mol)
	print("Read {}/{} molecules".format(len(smiles_list), len(suppl)))
	return smiles_list

def Smilestomodel(mol):
	mol_H = Chem.AddHs(mol)
	addbonds = AllChem.EmbedMolecule(mol_H)
	# using MMFF, could also add options for other methods
	AllChem.MMFFOptimizeMolecule(mol_H)
	return mol_H

def write_SDF(all_models):
	for mol in all_models:
		

if __name__ == '__main__':
	# define arguents here. Add more as needed.
	parser = argparse.ArgumentParser()
	group_input = parser.add_argument_group('INPUT arguments')
	group_input.add_argument("-i", "--input", metavar='filename', type=str, required=True, help="name of SMILES file")
	
	# parse arguments
	args = parser.parse_args()

	# read SMILES
	smiles_list = readSmiles(args.input)
	for mol in smiles_list:
		all_models = Smilestomodel(mol)
	write_SDF(all_models)
