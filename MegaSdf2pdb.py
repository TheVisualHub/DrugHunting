# sdf2pdb (ver 1.2 delta) : an Automated, Cutting-Edge Pipeline to convert a multi-state SDFs to separate PDBs
# The present version operates with two switchable strategies activated via &use_mega trigger
# Required rdkit package installed into your python environment
# Original pipiline of rdkit-based conversion by Valentin Guillaume, 
# Code optimization, introducing &use_mega trigger by Gleb Novikov
# ####### ######### ########## ############ ######## ######## ######
## How to install rdkit on MAc in the isolate environemnt to prevent any issues:
#python3 -m venv rdkit-env
#source rdkit-env/bin/activate
##pip install --upgrade pip
#pip install rdkit
# python ./MegaSdf2pdb.py
#
# The Visual Hub - All rights reserved - 2025 (c)

import os
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem

# basic options: the folders with input and outputs
input_dir = "sdf"
output_dir = "pdb"

# advanced options:
use_mega = True #  (True) activates two-step iteration strategy (by default)

def check_output_dir(output_dir):

    # check if the output dir exists
    if not os.path.exists(output_dir):
        print(f"Creating {output_dir} repertoire for all pdbs")
        os.makedirs(output_dir)
    else:
        print(f"The directory {output_dir} actually exists.")
        print(f"Removing ...")
        shutil.rmtree(output_dir)  # Deletes old outputs and all sub-dirs
        print(f"Creating new {output_dir} repertoire for all pdbs")
        os.makedirs(output_dir)

# (i) - a single-step iteration strategy
# all sdf files from the $input_dir  will be converted to pdbs and saved within the same $output_dir
# practical when each sdf file contains only a small number of structures
def simple_converter(input_dir, output_dir):

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

# (ii) - a two-step iteration strategy
# all sdf files from the $input_dir  will be converted to pdbs and saved within the subdirs of the $output_dir
# practical when you operate with many sdf corresponded to different chemical libraries
def mega_converter(input_dir, output_dir):

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

def superb():
#def superb(use_mega=True): # ignore the global trigger
    check_output_dir(output_dir)
    # use_mega is by default
    if use_mega:
        print(f"Two-step iteration strategy is activated for {input_dir} repertoire")
        mega_converter(input_dir, output_dir)
    else:
        print(f"One-step iteration strategy is activated for {input_dir} repertoire")
        simple_converter(input_dir, output_dir)

if __name__ == "__main__":
    superb()
else:
    print("Ciao bella ciao")