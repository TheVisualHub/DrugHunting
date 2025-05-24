# da Vinci (ver 1.0) : an Automated, Cutting-Edge Pipeline to Draw 2D Images of Drugs
# Turning Big Data into Drug Discoveries
# Required rdkit package installed into your python environment
# Python coding and optimization by Valentin Guillaume, the original idea by Gleb Novikov
# The Visual Hub - all rights reserved - 2025 (c)

from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import shutil

# advanced options for visualisations
project = "my_DrugDiscovery_project"
resolution = (300,300) # resolution of images
################################################
mol_file_extension = ".pdb"
home = Path("/home/gleb/Desktop/dolce_vita/")
ligands_directory = Path(home, "ligands")
pdb_directory = Path(ligands_directory, project)
output_dir = Path(pdb_directory, "perfetto")

def draw_mol_from_dir(mol_dir: Path, ouput_dir: Path):
    if output_dir.exists():
        shutil.rmtree(output_dir)
        print("The old pictures have been removed !")
    output_dir.mkdir()
    print("The new", output_dir, "has been created !")

    for mol_file in mol_dir.iterdir():
        if mol_file.suffix == mol_file_extension:
            mol = Chem.MolFromPDBFile(str(mol_file))
            if mol:
                AllChem.Compute2DCoords(mol)
                print("The portrait of", mol_file.stem, "has been painted. Go next !")
                Draw.MolToFile(mol, Path(ouput_dir, mol_file.stem + ".png"), size=resolution)
                # the method 2: does not work yet
                #d = rdMolDraw2D.MolDraw2DCairo(resolution[0], resolution[1])
                #d.drawOptions().addStereoAnnotation = True
                #d.drawOptions().addAtomIndices = True
                #d.DrawMolecule(mol)
                #d.FinishDrawing()
                #d.WriteDrawingText(str(Path(ouput_dir, mol_file.stem + ".png")))      
            else:
                print(mol_file, " couldn't be red as molecule by rdkit.")
    print("Bravo!")

def main():
    draw_mol_from_dir(pdb_directory, output_dir)


if __name__ == "__main__":
    main()
else:
    print("Ciao bella ciao")