# Dolce Vita (ver 1.0) : an Automated, Cutting-Edge Pipeline to Identify Drugs in Docking Predictions
# Turning Big Data into Drug Discoveries
# Tested with Autodock VINA results but easily adapted for any docking software
# Python coding and optimization by Valentin Guillaume, workflow design by Gleb Novikov
# The Visual Hub - all rights reserved - 2025 (c)

import re
import statistics
from pathlib import Path
import shutil

# top data to consider: e.g. take 5% of all docking predictions
top_data = 0.05
#top_data = 1.00 # this is all big data (for the test)

#name of the directory with the results
project = "results_docking"
#######
csv_file_pattern = "ranking_"
serie_directory_name = "chemical_serie007"
protein_name = "proteinXX"
filtered_serie_directory_name = serie_directory_name + "_filteredTOP2"
write_pdbs = True # True/False : Boolean
write_for_vizu = True # True/False : Boolean

home = Path("/home/gleb/Desktop/dolce_vita/")
input_directory = Path(home, project)
pdb_directory = Path(home, "ligands")

pdb_serie_directory = Path(pdb_directory, serie_directory_name)
filtered_pdb_serie_directory = Path(pdb_directory, filtered_serie_directory_name)
filtered_pdb_vizu_directory = Path(input_directory, "top_model")

output_file_name = Path(input_directory, "top_hits.{:0.2f}.csv".format(top_data))

###################################################################
###################################################################
####################### Script ####################################
###################################################################
###################################################################


class pose_data:
    mol_id : str
    population: int = 1
    energies: list


cpnd_enum = {} # key : mol id   value : pose_data


cpnd_id_regex = re.compile("SarsCov2_(?P<cpnd_id>.+)_rep[0-9]+.*")
energy_regex = re.compile(".* (?P<energy_str>(-|\\+)?[0-9]+(.[0-9]+)?)$")

def count_poses(file_content):
    return file_content.count("\n")



def add_compound(cpnd_id: str, energy: float):
    if cpnd_id in cpnd_enum:
        cpnd_enum[cpnd_id].population += 1
        cpnd_enum[cpnd_id].energies.append(energy)
    else:
        cpnd_enum[cpnd_id] = pose_data()
        cpnd_enum[cpnd_id].mol_id = cpnd_id
        cpnd_enum[cpnd_id].energies = [energy]


def enumerate_top_cpnd(file: Path):
    file_content: str = file.read_text()
    num_lines = count_poses(file_content)

    current_num_line = 0
    for line in file_content.split("\n"):
        current_num_line += 1
        if current_num_line / num_lines <= top_data:
            # Getting compound id from line
            match_rslt = cpnd_id_regex.match(line)
            cpnd_id = match_rslt.group("cpnd_id")

            # Getting energy from line
            energy_match = energy_regex.match(line)
            energy = float(energy_match.group("energy_str"))

            # Adding line data to global dict
            add_compound(cpnd_id, energy)
        else:
            break

def compute_mean_energy(cpnd_pose_data: pose_data):
    return statistics.mean(cpnd_pose_data.energies)

def compute_stdev_energy(cpnd_pose_data: pose_data):
    return statistics.stdev(cpnd_pose_data.energies)

def rank_and_write_results(file_found):
    row_template = "{:<10} {:^5} {:^9} {:^9} {:^9} {:^9}\n"
    output_content = row_template.format(
        "Mol ID"
        , "Pop"
        , "dG mean"
        , "dG stdev"
        , "dG min"
        , "dG max"
    )
    ranking = [(k,v) for k, v in cpnd_enum.items()]
    ranking.sort(key=lambda x: x[1].population, reverse=True)
    ranked_cpnds = [] # (mol_id, enum)
    for cpnd_id, cpnd_pose_data in ranking:
        if cpnd_pose_data.population >= file_found - 1:
            output_content += row_template.format(
                cpnd_id
                , cpnd_pose_data.population
                , "{:<0.2f}".format(compute_mean_energy(cpnd_pose_data))
                , "{:<0.2f}".format(compute_stdev_energy(cpnd_pose_data))
                , "{:<0.2f}".format(min(cpnd_pose_data.energies))
                , "{:<0.2f}".format(max(cpnd_pose_data.energies))
                )
            ranked_cpnds.append((cpnd_id, cpnd_pose_data.population))
        else:
            break    

    output_file_name.write_text(output_content)
    return ranked_cpnds

def copy_top_pdbs(ranked_cpnds, file_found):
    if filtered_pdb_serie_directory.exists():
        shutil.rmtree(filtered_pdb_serie_directory)
        print("The old data has been removed.")

    filtered_pdb_serie_directory.mkdir()

    for cpnd_id, rank_enum in ranked_cpnds:
        if rank_enum >= file_found - 1 :
            mol_pdb_path = Path(pdb_serie_directory, cpnd_id+".pdb")
            if not mol_pdb_path.exists():
                print("File for molecule", cpnd_id, "is not detected.")
            else:
                shutil.copy(mol_pdb_path, Path(filtered_pdb_serie_directory, mol_pdb_path.name))

def write_pdb_vizu(ranked_cpnds, repetition_number):
    if filtered_pdb_vizu_directory.exists():
        shutil.rmtree(filtered_pdb_vizu_directory)
       
    filtered_pdb_vizu_directory.mkdir()

    complexe_file_name_pattern = "SarsCov2_{}_[^_]+_rep"
    #complexe_file_name_pattern = "SarsCov2_{}_" + protein_name + "_rep[0-9]+.pdb"

    for cpnd_id, rank_enum in ranked_cpnds:
        expected_file_name_pattern = re.compile(complexe_file_name_pattern.format(cpnd_id))
        for file in input_directory.iterdir():
            if expected_file_name_pattern.match(file.name):
                shutil.copy(file, Path(filtered_pdb_vizu_directory, file.name))




def main():
    if output_file_name.exists():
        output_file_name.unlink()
        print("the old data has been removed.")

    if not input_directory.exists():
        print("Bella ciao")
    else:
        file_found = 0
        for file in input_directory.iterdir():
            if csv_file_pattern in file.stem:
                enumerate_top_cpnd(file) 
                file_found += 1
        print(file_found, "files found")
        ranked_cpnds = rank_and_write_results(file_found)
        print("The file \n{}\n has been created.".format(output_file_name))
        if write_pdbs:
            copy_top_pdbs(ranked_cpnds, file_found)
        if write_for_vizu:
            write_pdb_vizu(ranked_cpnds, file_found)


if __name__ == "__main__":
    main()
else:
    print("Ciao bella ciao")