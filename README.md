🔭🐍 The lattest development - Royal Gamble Lite - unfolding novel drug spaces through the principles of cosmic star creation. ✨ 


🎲 Royal Gamble Lite is a cheminformatics pipeline for crafting novel three-ring linear molecules using controlled stochasticity. Inspired by the randomness of the cosmos—and a dash of casino bluffing—it generates combinatorial molecules with optional chiral centers, random ring assembly, and structure filters.

This is the Lite version: streamlined, experimental, and designed for fast exploration.
🚀 Features

    🎲 Stochastic ring assembly from a curated pool of aromatic and heterocyclic rings.

    ♻️ Two ring-selection strategies:

        Pure Random: unbiased random choice

        Bluff-Spin Selector: adds controlled pseudo-randomness for combinatorial flair

    💫 Optional chiral carbon linkers

    🧹 Filter strategies to remove redundant structures (e.g., identical or adjacent rings)

    🧪 MMFF94s optimization (via RDKit)

    🗂️ SDF multi-model output with timestamped filenames

🧬 What It Does

This script generates a set of molecules with the general form:

[RING 1] — [LINKER] — [RING 2] — [LINKER] — [RING 3]

Each molecule:

    Is assembled from a ring SMILES pool

    Can include chiral centers based on settings

    Is filtered to eliminate duplicates or structurally boring outputs

🔧 Requirements

    Python 3.7+

    RDKit (tested on 2023.09.1)

Install RDKit via Conda:

conda create -c rdkit -n rgamble_env rdkit python=3.9
conda activate rgamble_env

⚙️ Configuration

Customize behavior by editing variables at the top of the script:

project_name = "CASINO"
number_of_structures = 18
chiral_switch = 1           # 0: off, 1–2: add chiral linkers
optimize_iters = 5000       # for accurate geometry
casino_switch = True        # activate Bluff-Spin selector
filter_identical = True
filtering_strategy = 2      # 1: no identical, 2: no adjacent duplicates, 3: all unique

🌀 Ring Pool

Predefined set of 18 molecular rings, including:

    Aromatic: benzene, naphthalene, pyridine, etc.

    Saturated: cyclohexane, norbornane

    Heterocycles: imidazole, pyrrolidine, piperidine, oxetane, THF

ring_smiles = [
  ("c1ccccc1", "benzene"),
  ("C1CCCCC1", "cyclohexane"),
  ...
]

📂 Output

Generates a multi-model .sdf file:

CASINO_18062025_1420.sdf

Each entry is fully 3D optimized with MMFF94s forcefield.
🧠 Behind the Scenes

    Ring selection is controlled via Python’s random and optionally modified by a bluff_spin() logic (mimicking gambling decisions).

    Molecules are connected using RDKit molecule constructors

    Filters remove boring, redundant, or chemically trivial outputs

    You can manually extend with custom ring sets, custom linkers, or property scoring


🧙 Author

Created by Gleb Novikov
Coded, conceptualized, and benchmarked in June 2025
Not affiliated with any gambling platform 😉
