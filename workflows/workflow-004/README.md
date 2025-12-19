# Workflow-004: In-Silico Virtual Screening with improved pose prediction

## Overview

This workflow automates **virtual screening** of ligand variants against a protein target using **Gnina**, an enhanced fork of smina.  
It integrates molecular structure extraction, ligand modification, receptor preparation, and docking that includes rescoring with convoluational neural networks (CNN) to improve pose prediction.

The target system modeled here is based on **PDB ID: 4OHU**, which includes **chain A**, and the ligand of interest **2TK**.

---

## Pipeline Structure

### Step 1: Download PDB structure (`download_pdb.py`)
1. Downloads the PDB file (`4OHU.pdb`) from RCSB.

**Key Outputs:**
4OHU.pdb

**Usage Example:**
```bash
python3 download_pdb.py
```

### Step 2: Protein Preparation (`protein_preparation.py`)
Prepares the target protein for docking.
**Process:**
1. Selects and extracts chain A in 4OHU.pdb removing all non-protein heteroatoms (including ligands).
2. Fixes structural issues using PDBFixer, adds hydrogens and minimizes energy 
3. Copies final receptor to the `protein_files/` directory.

**Key Outputs:**
```
4OHU_A.pdb
4OHU_A_fixed.pdb
protein_files/4OHU_A.pdbqt
```

**Usage Example:**
```bash
python3 protein_preparation.py
```

### Step 3: Ligand Extraction
Extracts ligand coordinates from the downloaded PDB file.

**Process:**
1. Reads the `4OHU.pdb` structure.
2. Screens the protein for all co-crystallized ligands and prompts user to select one.
3. Isolates selected ligand from  PDB and saves as ligand_id variable.
4. Exports ligand_id as `{ligand_id}.sdf` (preserving its 3D coordinates).
5. Downloads an ideal ligand from RCSB (as {ligand_id}_ideal.sdf).
6. Partially fixes and aligns coordinates with {ligand_id}_ideal.sdf and rearomatize {ligand_id}.sdf and saves it as {ligand_id}_corrected_pose.sdf
7. 

**Key Outputs:**
```
{ligand_id}.sdf
{ligand_id}_ideal.sdf
{ligand_id}_corrected_pose.sdf
```

**Usage Example:**
```bash
python3 ligand_extraction_and_preparation.py
```

### Step 4: In-Silico Screening (`in_silico_screening_and_reporting.py`)
Performs docking of each ligand against the prepared receptor using Gnina.

**Process:**
1. Verifies Gnina installation.
2. Prompts user to choose docking mode and configure run parameters:

   **Docking mode**
   a) single docking
   b) batch docking
   c) flexible docking
   d) docking on unknown sites

   **Exhaustiveness**
   8, 16, 24, ..., 80
   
   **GPU usage**
   yes/no

   **Convolutional neural network (CNN) usage**
   yes/no

   If **yes** was selected, CNN score will be displayed (compute-intensive)
   
5. Runs Gnina docking with parameters selected above.
6. Writes report and saves it as csv file.

**Outputs:**
```
docking_results/
├── {ligand_id}_docked_{pdb_id}.sdf
├── multiple_ligands_docked_{pdb_id}.sdf
├── {ligand_id}_flex.sdf
├── {ligand_id}_docked_whole_{pdb_id}.sdf
├── docking_results_{pdb_id}.csv
```

**Usage Example:**
```bash
python3 in-silico_screening_and_reporting.py
```

---

## Expected Runtime

| Stage | Duration (Approx.) |
|--------|--------------------|
| Protein Chain Extraction + Preparation | ~2–4 min |
| Ligand Extraction + Preparation | ~2–4 min |
| Docking (per ligand) | ~2–5 min |
| Full Library (10–20 ligands) | ~1–2 hrs |
| Report Generation | <1 min |

---

## Output Structure

```
workflow-004/
├── 4OHU_A.pdb
├── *.sdf
├── docking_results/
│   ├── *_docked_{pdb_id}.sdf
│   ├── *_flex.sdf
└── ├── *_docked_whole_{pdb_id}.sdf
```

---

## Dependencies

Installed inside the Docker environment:
- **python 3.11**
- **RDKit**
- **Gnina**
- **BioPython**
- **OpenMM** / **PDBFixer**
- **NumPy**
- **Pandas**
- **Scipy**
- **MDAnalysis**
- **Requests**
- **Openbabel-wheel**
- **Useful-rdkit-utils**
- **Molscrub**
- **Libstdcxx-ng**

---

## References
- **PDBFixer:** https://github.com/openmm/pdbfixer
- **MDAnalysis:** https://github.com/MDAnalysis/mdanalysis
- **Molscrub:** https://github.com/forlilab/molscrub
- **Openbabel-wheel:** https://github.com/njzjz/openbabel-wheel
- **Scipy:** https://github.com/scipy/scipy
- **Useful_rdkit_utils:** https://github.com/PatWalters/useful_rdkit_utils
- **PDB101 tutorial by RCSB Protein Data Bank:** https://pdb101.rcsb.org/train/training-events/python4
- **Gnina:** McNutt, A.T., Li, Y., Meli, R. et al. GNINA 1.3: the next increment in molecular docking with deep learning. J Cheminform 17, 28 (2025). https://doi.org/10.1186/s13321-025-00973-x



