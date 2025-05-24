# sdf2pdb (ver 1.1 delta) : an Automated, Cutting-Edge Pipeline to convert a multi-state SDFs to separate PDBs
# Turning Big Data into Drug Discoveries
# Required rdkit package installed into your python environment
# Original pipiline of sdf2pdb conversion coding by Valentin Guillaume, 
# Code optimization by Gleb Novikov: adapting for any number of input multi-model SDF
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

    # perform conversion:
    # iterate over each sdf in the input_dir
    for filename in os.listdir(input_dir):
        if filename.endswith(".sdf"):
            # i - take the name of the sdf file
            sdf_name = os.path.splitext(filename)[0]
            sdf_path = os.path.join(input_dir, filename)

             # ii- create a subfolder inside the $output_dir for this sdf file
            sdf_output_folder = os.path.join(output_dir, sdf_name)
            os.makedirs(sdf_output_folder, exist_ok=True)

            # iii - perform conversion 
            suppl = Chem.SDMolSupplier(sdf_path)
            # iterate over all structures in the sdf
            for i, mol in enumerate(suppl):
                if mol is None:
                    print(f"⚠️ Skipping molecule {i+1} in {sdf_name} (unreadable)")
                    continue
                
                # Add hydrogens and generate 3D coordinates
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
                AllChem.UFFOptimizeMolecule(mol)

                # Define output pdb filename
                pdb_filename = f"{sdf_name}_{i+1}.pdb"
                pdb_path = os.path.join(sdf_output_folder, pdb_filename)

                # Write pdb file for a given sub-structure
                with open(pdb_path, "w") as f:
                    f.write(Chem.MolToPDBBlock(mol))

            print(f"Processed all structures from the {filename} -> {sdf_output_folder}/")

def superbe():
    check_output_dir(output_dir)
    converter(input_dir, output_dir)

# (!): added in the last revision of the script:
if __name__ == "__main__":
    superbe()
else:
    print("Ciao bella ciao")