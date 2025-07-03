# 🎰 Designing drug-like molecules using combinatorics inspired by casino models to unfold novel chemical spaces 💊  

🎥 Click on the image to watch the demonstration:
[![Watch the video](https://img.youtube.com/vi/gc762D6thsg/maxresdefault.jpg)](https://www.youtube.com/watch?v=gc762D6thsg)

##  🔍 Overview 
🎲 **The Royal Gamble Lite** is an innovative cheminformatics pipeline designed for crafting novel molecules using controlled stochasticity. Inspired by the univeral principles of randomness — and a dash of casino bluffing — this cutting-edge approach generates combinatorial molecules with optional chiral centers, random ring assembly, and structural-based filters. Streamlined and designed for fast exploration:    
1️⃣ Random assembly of chemical rings from **SMILES**, followed by **3D reconstruction** in **PDB** or **SDF** formats. Ready for docking!  
2️⃣ **Bluff-spin algorithms** based on casino gaming to ramp up stochastic diversity.    
3️⃣ Generation of relevant **protonation states** and **tautomers**.  
4️⃣ Addition of **chiral linkers** to expand stereochemical possibilities.  
5️⃣ Geometrical **optimization** of energetically favorable conformers.  
6️⃣ Multiple **filtering strategies** seamlessly integrated with established medicinal chemistry guidelines.  
7️⃣ **Industry-grade stability** — the code keeps crafting compounds even after a crash! 😁 

---

## 👤 Author

This script was developed and benchmarked by **Gleb Novikov**

---

## 🚀 Main Features

🎲 Stochastic ring assembly from a curated pool of aromatic and heterocyclic rings.

🛠️ The following ring-selection strategies are implemented:

-  Pure Randomness: unbiased random choice

-  Bluff-Spin Selector: adds simple bluff-function for combinatorial flair
    
-  (only in FULL version): custom Bluff-Spin functions:  
    🎰 simulating a weighted roulette spin over the pool of rings, assigning dynamic probabilities influenced by recent picks and structural diversity to avoid repetition in each sample. It then introduces a a controlled randomness layer that occasionally selects lower-probability rings, mimicking the unpredictability of a casino bluff —thereby introducing a trade-off between exploitation of high-probability candidates and exploration of under-sampled regions in the combinatorial search space. ⚖️ 

🔗 Optional chiral carbon linkers

🔍 Filter strategies to remove redundant structures (e.g., identical or adjacent rings)

🧲 MMFF94s optimization (via RDKit)

🗂️ SDF multi-model output with timestamped filenames

---


## 🧬 What It Does

➡️ The LITE version of the script generates a set of linear molecules with the three-ring pattern:

```🌀[RING 1] — 🔗[LINKER] — 🌀[RING 2] — 🔗[LINKER] — 🌀[RING 3]```

➡️ The FULL version of the script generates a set of linear molecules with the N-ring pattern:

```🌀[RING 1] — 🔗[LINKER] — 🌀[RING 2] — 🔗[LINKER] ... 🌀[RING N]```

🌟 and additionally star-shaped molecules with 4 or 5 rings:

```
         🌀[RING 4]
            │
🌀[RING 1]—🌀[RING 2]—🌀[RING 3]
            │
         🌀[RING 5]
```



In this way, each generated molecule:

⚙️ Is assembled from a ring SMILES pool

🔀 Can include chiral centers based on settings

🔍 Is filtered to eliminate duplicates or structurally boring outputs

---

## 🔧 Requirements

    Python 3.7+

    RDKit (tested on 2023.09.1)

Install RDKit via Conda:

```bash
conda create -c rdkit -n rgamble_env rdkit python=3.9
conda activate rgamble_env
```

---

## ⚙️ Configuration

### 🎛️ Main Control

Customize all variables at the top of the script. The Full version iontroduces a dedicated DATACLASS for tailored control.

```bash
project_name = "CASINO"
number_of_structures = 18
chiral_switch = 1           # 0: off, 1–2: add chiral linkers
optimize_iters = 5000       # for accurate geometry
casino_switch = True        # activate Bluff-Spin selector
filter_identical = True
filtering_strategy = 2      # 1: no identical, 2: no adjacent duplicates, 3: all unique
```

### 🌀 Ring Pool

The LITE version has predefined set of 18 chemical rings, including:

    Aromatic: benzene, naphthalene, pyridine, etc.

    Saturated: cyclohexane, norbornane

    Heterocycles: imidazole, pyrrolidine, piperidine, oxetane, THF

```bash
ring_smiles = [
  ("c1ccccc1", "benzene"),
  ("C1CCCCC1", "cyclohexane"),
  ...
]
```

### 📂 Output

**Generates a multi-model .sdf file:**

```bash
CASINO_{timestamp}.sdf
```

Each entry is fully 3D optimized with MMFF94s forcefield.
###  🧠 Behind the Scenes

By combining unmatched stochastic algorithms with physics-driven optimization and advanced medicinal chemistry filters, this tool generates more diverse, innovative molecules while ensuring reliability and seamless workflow integration — reshaping the frontiers of modern drug discovery!

- Ring selection is controlled via Python’s random and optionally modified by a bluff_spin() logic (mimicking gambling decisions).

- Molecules are connected using RDKit molecule constructors

- Filters remove redundant crafts, or chemically trivial outputs

- You can manually extend with custom ring sets, linkers, or direct output to your scoring functions for ultimate flexibility

---
👤 The Visual Hub (2025)  
Hope you enjoy it ! 🧡 ✨

