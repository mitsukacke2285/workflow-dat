# Workflow-004: In-Silico Virtual Screening with improved pose prediction

## Overview

This workflow automates **virtual screening** of ligand variants against a protein target using **Gnina**, an enhanced fork of smina.  
It integrates molecular structure extraction, ligand modification, receptor preparation, and docking that includes rescoring with convoluational neural networks (CNN) to improve pose prediction.

The target system modeled here is based on **PDB ID: 4OHU**, which includes **chain A**, and the ligand of interest **2TK**.

---

## Pipeline Structure

### Step 1: Download PDB structure (`download_pdb.py`)
1. **Downloads** the PDB file (`4OHU.pdb`) from RCSB.

**Key Outputs:**
4OHU.pdb

**Usage Example:**
```bash
python3 download_pdb.py
```

### Step 2: Protein Preparation (`protein_preparation.py`)
Prepares the **target protein** for docking.
**Process:**
1. **Selects and extracts chain A** in 4OHU.pdb removing all non-protein heteroatoms (including ligands).
2. **Fixes structural issues** using **PDBFixer**, adds hydrogens and minimizes energy 
8. Copies final receptor to the `protein_files/` directory.

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

### Step 3: Ligand Extraction and Preparation (`ligand_extraction_and_preparation.py`)
Extracts ligand coordinates from the downloaded PDB file.

**Process:**
1. Reads the `4OHU.pdb` structure.
2. Screens the protein for all co-crystallized ligands for the user to select; saves it as ligand_id.
3. Isolates selected ligand from  PDB and saves as ligand_id variable.
4. Exports ligand_id as `{ligand_id}.sdf` (preserving its 3D coordinates).
5. Downloads an ideal ligand from RCSB (as {ligand_id}_ideal.sdf).
6. Partially fixes and aligns coordinates with {ligand_id}_ideal.sdf and rearomatize {ligand_id}.sdf and saves it as {ligand_id}_corrected_pose.sdf

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

---

### Step 4: In-Silico Screening (`in_silico_screening_and_reporting.py`)
Performs docking of each ligand against the prepared receptor using **Gnina**.

**Process:**
1. Verifies Gnina installation.
2. Prompts user to select docking mode.

**Configuration:**
- Scoring: Vina
- Number of modes: 4
- Exhaustiveness: 8-80

**Outputs:**
```
docking_results/
├── variant_name_docked.sdf
├── variant_name_docking.log
```

**Usage Example:**
```bash
python3 in_silico_screening.py
```

---


## Expected Runtime

| Stage | Duration (Approx.) |
|--------|--------------------|
| Ligand Extraction + Preparation | ~2–4 min |
| Ligand Modification | ~3–5 min |
| Docking (per ligand) | ~2–5 min |
| Full Library (10–20 ligands) | ~1–2 hrs |
| Report Generation | <1 min |

---

## Dependencies

Installed inside the Docker environment:
- **Python => 3.11**
- **Numpy**
- **Requests**
- **Pandas**
- **Scipy**
- **RDKit**
- **Gnina**
- **BioPython**
- **OpenMM** / **PDBFixer**
- **MDAnalysis**

---

## References
- **PDB101 tutorial:** by RCSB Protein Data Bank: https://pdb101.rcsb.org/train/training-events/python4
- **Gnina:** McNutt, A.T., Li, Y., Meli, R. et al. GNINA 1.3: the next increment in molecular docking with deep learning. J Cheminform 17, 28 (2025). https://doi.org/10.1186/s13321-025-00973-x
- Buccheri, R.; Rescifina, A. High-Throughput, High-Quality: Benchmarking GNINA and AutoDock Vina for Precision Virtual Screening Workflow. Molecules 2025, 30, 3361. https://doi.org/10.3390/molecules30163361
- Sunseri, J., & Koes, D. R. (2021). Virtual Screening with Gnina 1.0. Molecules (Basel, Switzerland), 26(23), 7369. https://doi.org/10.3390/molecules26237369
- **PDBFixer:** https://github.com/openmm/pdbfixer  







