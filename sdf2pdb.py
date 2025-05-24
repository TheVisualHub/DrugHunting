# sdf2pdb (ver 1.0b) : an Automated, Cutting-Edge Pipeline to convert a multi-state SDFs to separate PDBs
# Turning Big Data into Drug Discoveries
# Required rdkit package installed into your python environment
# Original pipiline of sdf2pdb conversion coding by Valentin Guillaume, 
# Code optimization by Gleb Novikov: split code to functions
#
## How to install rdkit on MAc in the isolate environemnt to prevent any issues:
#python3 -m venv rdkit-env
#source rdkit-env/bin/activate
##pip install --upgrade pip
#pip install rdkit
#
# The Visual Hub - all rights reserved - 2025 (c)

from rdkit import Chem
from rdkit.Chem import AllChem
#from pathlib import Path # used in the first version, changed to os
import os

# basic options: the folders with input and outputs
input_dir = "sdf"
output_dir = "pdb"

def check_output_dir(output_dir):

    # check the output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        print(f"The directory {output_dir} actually exists.")

def converter(input_dir, output_dir):

    # make conversion for each sdf in the input_dir
    for filename in os.listdir(input_dir):
        if filename.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(os.path.join(input_dir, filename))
            for i, mol in enumerate(suppl):
                if mol is None:
                    continue
                
                # initiate magic
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
                AllChem.UFFOptimizeMolecule(mol)
                output_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_{i}.pdb")
                with open(output_path, 'w') as f:
                    f.write(Chem.MolToPDBBlock(mol))

def main():
    check_output_dir(output_dir)
    converter(input_dir, output_dir)

# added in the last revision of the script:
if __name__ == "__main__":
    main()
else:
    print("Ciao bella ciao")