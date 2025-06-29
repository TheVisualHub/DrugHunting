# 🔮✨ the Royal Gamble Lite (ver 1.08 demo):
# Last update 5/06/2025: 
# The LITE version creates random three-ring linear molecules optionally joined by chiral carbons (chiral_switch=1..3)
# Conceptualized, coded and benchmarked by Gleb Novikov
# This script generates a muilti-model sdf file
# The number of prepared structures provided by the var => number_of_structures
# The current revision introduces two different strategies to select the chemical rings from the ring pool (ring_smiles)
# Finally after creation script introduce three filterins strategies applied on the crafts
# (e.g. to discard crafts contained identical rings etc)
#
# Two random strategies for ring selection:
# 🌀 Strategy 1: Pure Randomness:
# the rings are selected randomly from the pool using random.choice
# 🎲 Strategy 2: Bluff-Spin Selector (activated by the casino_switch)
# Inspired by casino bluffing, the function (bluff_spin()) adds a twist of an additional stochasticity
#
# You are free to specify which ring positions use the casino-style selector and which ones follow the straight-random approach,
# which gives you precise control over the combinatorial outcomes.
# The author does not support or promote gambling
#
## basic system modules:
import os
import sys
import glob
import shutil
## management modules
import time
#import math # ! only in full version
import random
from datetime import datetime
#import subprocess #only in full version
#import tempfile # only in full version
## rdkit modules (for sdf manips)
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
## advanced rdkit modules (for coordinate handling) #! only for full version
#from rdkit.Chem import rdMolAlign #! only in full version
#from rdkit.Chem.rdMolTransforms import ComputeCentroid #! only in full version
#from rdkit.Geometry import Point3D # ! only in full version
# advanced boolean control via dataclass
#from dataclasses import dataclass #only in full version

# System options:
script_version = '1.07 demo'

# Main options:
project_name = "CASINO" # a name of generated multi-model sdf
number_of_structures = 18 # number of generalted structures
chiral_switch = 1 # (select: 0, 1 or 2) adding one chiral linker carbon by default

# a geometrical optimization with MMFF94s
optimize_iters = 5000 # 200 is default, 5000 - for superb accuracy due to chiral centers

# ACTIVATE FILTERING:
# strategy 1 - discard triples => all crafts where all three rings the same
# strategy 2 - discard adjacent duplicate rings => all crafts where ring1==ring2 or ring2==ring3
# strategy 3 - 	allow only completely unique rings => keep only crafts with three different rings
filter_identical = True
filtering_strategy = 2 # the 2nd is more balanced

# STOCHASTIC CONTROL: activate Casino-inspired ring selector
casino_switch = True # True or False ?

# File management options:
delete_old_sdf = True

# put a creative time-stamp in the crafted sdf
timestamp = datetime.now().strftime("%d%m%Y_%H%M")

# a list of ring SMILES : 
ring_smiles = [
    # Only carbons
    ("c1ccccc1", "benzene"),
    ("C1CCCCC1", "cyclohexane"),
    ("C1CCC1", "cyclobutane"),
    ("C1CCCC1", "cyclopentane"),
    ("c1ccc2ccccc2c1", "naphthalene"),
    ("C1CC2CCC1C2", "norbornane"),

    # Nitrogen-containing rings
    ("c1ccncc1", "pyridine"),
    ("c1cncnc1", "pyrimidine"),
    ("c1ccnc2c1cccc2", "quinoline"),
    ("c1cnc[nH]1", "imidazole"),
    ("C1CNCCN1", "piperazine"),
    ("C1CCNC1", "pyrrolidine"),
    ("C1CCNCC1", "piperidine"),

    # Oxygen-containing rings
    ("C1COC1", "oxetane"),
    ("C1COCC1", "tetrahydrofuran"),         # THF
    ("c1ccoc1", "furan"),
    ("C1COC(C1)O", "1,3-dioxolane"),
    ("C1COC2OC1C2", "1,4-dioxane"),
]


# ! not implemented in the lite version !
chiral_linkers = [
    ("[C@H](F)C", "Chiral carbon bonded to fluorine and methyl"),
    ("[C@H](Cl)C", "Chiral carbon bonded to chlorine and methyl"),
    ("[C@H](Br)C", "Chiral carbon bonded to bromine and methyl"),
    ("[C@H](C)C", "Chiral carbon bonded to two methyl groups (isopropyl-like)"),
    ("[C@H](O)C", "Chiral carbon bonded to hydroxyl and methyl"),
    ("[C@H](N)C", "Chiral carbon bonded to amine and methyl"),
    ("[C@@H](C)F", "Enantiomer of [C@H](F)C (opposite stereochemistry)"),
    ("C[C@H](S)C", "Chiral carbon bonded to sulfur and methyl"),
    ("[C@H](Br)O", "Chiral carbon bonded to bromine and hydroxyl"),
    ("[C@H](CO)C", "Chiral carbon bonded to hydroxymethyl and methyl"),
]


def remove_old_sdf(delete_old_sdf=True):
    if delete_old_sdf:
        print(f"🧹Removing previusly generated SDF")
        time.sleep(0.5)
        # Get a list of all .sdf files in the current directory
        sdf_files = glob.glob("*.sdf")

        # Remove each file
        for file in sdf_files:
            try:
                os.remove(file)
                print(f"✅ Removed: {file}")
            except Exception as e:
                print(f"Error removing {file}: {e}")
    else:
        print(f"All previously generated SDF will be kept")


def print_input_info(
    chiral_info=True,
    ring_info=True,
    extra_chiral=False,
    casino_info=True,
    debug=True
):

    if chiral_info:
        """
        Print info about chiral centers added based on the current chiral setup.
        """
        if chiral_switch == 0:
            print(f"🌀 No chiral centers will be added:")
            print("🧿 Ring1 — 🧿 Ring2 — 🧿 Ring3")
        elif chiral_switch == 1:
            print(f"🌀 Just one chiral centers will be added:")
            print("🧿 Ring1 — 💠 Chiral1 — 🧿 Ring2 — 🧿 Ring3")
        elif chiral_switch == 2:
            print("🌀 Two chiral centers will be added:")
            print("🧿 Ring1 — 💠 Chiral1 — 🧿 Ring2 — 💠 Chiral2 — 🧿 Ring3")
        else:
            print("Sorry, the casion is closed today. Good bye! 🎲")
            sys.exit("💥 Script terminated due to invalid chiral input.")
        print()
        time.sleep(2)
    else:
        print("🔮 Welcome back, Master!")
    
    if ring_info:
        print("💍 The following rings will be used to craft molecules:")
        time.sleep(0.5)  # Sleep for 0.3 seconds
        for smiles, desc in ring_smiles:
            print(f" - {smiles}: {desc}")
            if debug:
                time.sleep(0.3)  # Sleep for 0.5 seconds
    else:
        print("💍 No ring info provided")
    
    if filter_identical:
        """
        Print info about activated filtering strategy.
        """
        if filtering_strategy == 1:
            print(f"🔍 THE FILTERING STRATEGY {filtering_strategy}: discard triples")
        elif filtering_strategy == 2:
            print(f"🔍 THE FILTERING STRATEGY {filtering_strategy}: discard adjacent duplicate rings")
        elif filtering_strategy == 3:
            print(f"🔍 THE FILTERING STRATEGY {filtering_strategy}: keep crafts with unique rings")
        else:
            print("Sorry, the casion is closed today. Good bye! 🎲")
            sys.exit("💥 Script terminated due to invalid filtering strategy selection.")
        print()
        time.sleep(1)
    else:
        print(f"👐 All crafts accepted without filtering")

    # optionally add an extra chiral linker to the crafted rings
    if extra_chiral:
        print("💠 The following extra chiral linkers will be used to join the rings:")
        time.sleep(0.5)  # Sleep for 0.5 seconds
        for smiles, desc in chiral_linkers:
            print(f" - {smiles}: {desc}")
            if debug:
                time.sleep(0.3)  # Sleep for 0.5 seconds
    else:
        print("👌 No extra chirality activated..")
    
    if casino_switch:
        print(f"🎰 The chemical rings will be selected using bluff-spin")
        time.sleep(1)
    else:
        print(f"🤷‍♂️ The chemical rings will be selected randomly")
        time.sleep(1)

# the first version without the fix onto the Pyrazine rings
def add_ring_to_mol(emol, ring_smiles_tuple):
    """
    Add a ring to the editable molecule emol, given a tuple (smiles, name).

    Args:
        emol (rdkit.Chem.RWMol): Editable molecule.
        ring_smiles_tuple (tuple): Tuple like (smiles, ring_name).

    Returns:
        dict: Mapping of original atom indices in ring to new atom indices in emol.
    """
    ring_smiles = ring_smiles_tuple[0]  # Extract SMILES from tuple
    ring_name = ring_smiles_tuple[1]

    ring = Chem.MolFromSmiles(ring_smiles)
    if ring is None:
        print(f"⚠️ Invalid SMILES: {ring_smiles} ({ring_name})")
        return None

    ring = Chem.AddHs(ring)
    ###ring = Chem.RemoveHs(ring)
    atom_map = {}

    for atom in ring.GetAtoms():
        new_idx = emol.AddAtom(atom)
        atom_map[atom.GetIdx()] = new_idx

    for bond in ring.GetBonds():
        emol.AddBond(atom_map[bond.GetBeginAtomIdx()], atom_map[bond.GetEndAtomIdx()], bond.GetBondType())

    return atom_map

# with chiral linker this should be modified  ..
def create_3ring_molecule(
    casino=False,
    casino_rings=[1, 2, 3], # e.g. to exclude the first ring from the bluff_spin casino_rings=[2, 3]
    chiral=0,
    ring_pool=None,
    optimize=True,
    cycles=None
):

    if ring_pool is None:
        ring_pool = ring_smiles  # calling global smile list
    if cycles is None:
        cycles = optimize_iters  # fallback to global

    emol = Chem.RWMol()
    if casino and 1 in casino_rings:
        print(f"🎲 The first ring is selected using the bluff-spin")
        time.sleep(0.1)
        ring1_smiles = bluff_spin(ring_pool)
    else:
        print(f"🌀 The first ring is selected randomly from the ring pool")
        time.sleep(0.1)
        ring1_smiles = random.choice(ring_pool)
    ring1_map = add_ring_to_mol(emol, ring1_smiles)
    if ring1_map is None:
        return None

    if chiral >= 1:
        chiral1_idx = emol.AddAtom(Chem.Atom(6))  # Carbon
        chiral1 = emol.GetAtomWithIdx(chiral1_idx)
        chiral1.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        emol.AddBond(chiral1_idx, ring1_map[0], Chem.BondType.SINGLE)

    if casino and 2 in casino_rings:
        print(f"🎲 The second ring is selected using the bluff-spin")
        time.sleep(0.1)
        ring2_smiles = bluff_spin(ring_pool)
    else:
        print(f"🌀  The second ring is selected randomly from the ring pool")
        time.sleep(0.1)
        ring2_smiles = random.choice(ring_pool)
    ring2_map = add_ring_to_mol(emol, ring2_smiles)
    if ring2_map is None:
        return None

    # add bond either between the first chiral center and the second ring
    # or just connect the both rings without chiral centers
    if chiral >= 1:
        emol.AddBond(chiral1_idx, ring2_map[0], Chem.BondType.SINGLE)
    else:
        emol.AddBond(ring1_map[0], ring2_map[0], Chem.BondType.SINGLE)

    # adding the second chiral center
    if chiral == 2:
        chiral2_idx = emol.AddAtom(Chem.Atom(6))  # Carbon
        chiral2 = emol.GetAtomWithIdx(chiral2_idx)
        chiral2.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        # Defensive: check ring2 has atom 1
        if 1 not in ring2_map:
            return None
        emol.AddBond(chiral2_idx, ring2_map[1], Chem.BondType.SINGLE)

    if casino and 3 in casino_rings:
        print(f"🎲 The third ring is selected using the bluff-spin")
        time.sleep(0.1)
        ring3_smiles = bluff_spin(ring_pool)
    else:
        print(f"🌀  The third ring is selected randomly from the ring pool")
        time.sleep(0.1)
        ring3_smiles = random.choice(ring_pool)
    ring3_map = add_ring_to_mol(emol, ring3_smiles)
    if ring3_map is None:
        return None

    # abit more comples: adding bonds
    if chiral == 2:
        emol.AddBond(chiral2_idx, ring3_map[0], Chem.BondType.SINGLE)
    elif chiral == 1:
        emol.AddBond(ring2_map[1], ring3_map[0], Chem.BondType.SINGLE)
    else:  # chiral == 0
        emol.AddBond(ring2_map[1], ring3_map[0], Chem.BondType.SINGLE)

    mol = emol.GetMol()
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)

    res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if res != 0:
        return None
    if optimize:
        #print(f"✨ Optimizing the whole structure ..")
        #AllChem.UFFOptimizeMolecule(mol, maxIters=optimize_iters)
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s', maxIters=optimize_iters)

    #return mol
    return mol, [ring1_smiles, ring2_smiles, ring3_smiles]


# ! not implemented in the lite version !
def add_chiral_linker(emol):
    linker = Chem.MolFromSmiles(chiral_linker_smiles)
    if linker is None:
        print(f"⚠️ Invalid chiral linker SMILES: {chiral_linker_smiles}")
        return None, None

    linker = Chem.AddHs(linker)
    atom_map = {}
    for atom in linker.GetAtoms():
        new_idx = emol.AddAtom(atom)
        atom_map[atom.GetIdx()] = new_idx

    for bond in linker.GetBonds():
        emol.AddBond(atom_map[bond.GetBeginAtomIdx()], atom_map[bond.GetEndAtomIdx()], bond.GetBondType())

    # Assume linker start is atom 0 and end is atom with idx = last atom
    return atom_map, emol


# a simple CASINO-inspired ring selector
# ! an more advancrd verion available in the full version !
def bluff_spin(ring_list):
    """
    A bluff-style random selector.
    70% chance: pick and keep.
    30% chance: bluff (reroll once).
    """
    first = random.choice(ring_list)
    if random.random() < 0.3:
        # Bluff: pretend we keep it, then switch
        second = random.choice(ring_list)
        return second
    else:
        return first


## FILTERING STRATEGIES:
def filter_identical_rings(
    result,
    strategy=2
):
    """
    Filtering strategy to exclude molecules with three identical rings.
    Returns True if the molecule should be kept, False if it should be discarded.
    """
    smiles_list = result[1]
    ring1, ring2, ring3 = smiles_list
    #return not (smiles_list[0] == smiles_list[1] == smiles_list[2])
    if filtering_strategy == 1:
        return not (ring1 == ring2 == ring3)         # keep unless all three are identical
    elif filtering_strategy == 2:
        return not (ring1 == ring2 or ring2 == ring3)  # discard if ring1==ring2 or ring2==ring3
    elif filtering_strategy == 3:
        return len(set([ring1, ring2, ring3])) == 3  # keep only if all three are unique !!


def main():
    print(f"👑 Welcome to the Royal Gamble ver.{script_version} 👑")
    time.sleep(0.5)
    print(f"🔮 I am going to cast {number_of_structures} structures, Master!")
    time.sleep(0.5)
    print(f"🧙‍♂️✨ But let me check the reagents first! ")
    time.sleep(0.5)
    # refresh old directories with generated sdf
    remove_old_sdf()
    # print info about used inputs
    print_input_info()

    print(f"🔆 LET'S START THE SORCERY 🔆")
    time.sleep(0.5)  # An extra  pause before the things start !
    # main workflow
    writer = SDWriter(project_name + "_" + timestamp + '.sdf')
    count = 0
    target = number_of_structures
    attempts = 0
    max_attempts = 1000

    while count < target and attempts < max_attempts:
        # Unpack the crafted rings in the result
        result = create_3ring_molecule(casino=casino_switch,chiral=chiral_switch)

        # Apply filtering if enabled
        if filter_identical and not filter_identical_rings(result, filtering_strategy):
            attempts += 1
            discarded_name = [name for _, name in result[1]]
            print(f"🃏 DISCARDED: the {attempts} attempt contains {discarded_name}") # or use directly result[1]
            continue


        #count an attempt for created molecule
        attempts += 1
        if result is not None:
            mol, ring_smiles_used = result
            writer.write(mol)
            
            count += 1
            print(f"✨ The molecule {count} is created at the {attempts} attempt")
            time.sleep(0.2)

    writer.close()
    print(f"⚗️ Work completed: totally crafted {count} molecules in {attempts} attempts.")

# call the main function
if __name__ == "__main__":
    main()
