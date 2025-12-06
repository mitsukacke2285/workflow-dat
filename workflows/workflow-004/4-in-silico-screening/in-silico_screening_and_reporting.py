#!/usr/bin/env python
# coding: utf-8

# # Step 1 Installation of dependencies
# Gnina will run within a linux environment provided by google colab virtual machine.
# 
# 1. `useful_rdkit_utils` is a Python package written and maintained by Pat Walters that contains useful RDKit functions. We will use it for the functions `mcs_rmsd` (explained later).
# 2. `py3Dmol` is used for molecular visualization.
# 3. The RDKit is a popular cheminiformatics package we will use for processing molecules.
# 

# In[ ]:


#import subprocess
#import os
#import requests

### Check for previous gnina installation ###
#gnina_file = "gnina"

# Download URL
#gnina_url = "https://github.com/gnina/gnina/releases/download/v1.3/gnina"

# Check if file exists
#if os.path.exists(gnina_file):
    #print(f"✅ {gnina_file} already exists, skipping download.")
#else:
    #print(f"⬇️ {gnina_file} not found, downloading...")
    #subprocess.run([
        #"wget",
        #gnina_url,
        #"-O", gnina_file
    #], check=True)
    
    # Make gnina executable
    #subprocess.run(["chmod", "+x", "gnina"])
    #print("✅ Download complete.")

#pdb_id = input("Enter PDB code: ") # The Protein ID we're looking at
#ligand_id = input("Enter ligand code: ") # The ID of the co-crystallized ligand

### Working directories ###
ligand_id = "2TK"
docking_results_directory = "docking_results"

import subprocess
import os
#import useful_rdkit_utils as uru
#from rdkit import Chem
#from rdkit import RDLogger
#from rdkit.Chem import PandasTools
#from rdkit.rdBase import BlockLogs
#import pandas as pd

docking_results = ""
# ------------------------------------------------------------------------------
# 4.1 Docking modes
# ------------------------------------------------------------------------------

def single_ligand_docking():
    receptor = f"{pdb_id}_A.pdbqt"
    ligand = f"{ligand_id}.sdf"
    grid_box = f"{ligand_id}_corrected_pose.sdf"
    output_file = f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf"
    subprocess.run([
    #"./gnina",
    "gnina",
    "-r", receptor,
    "-l", ligand,
    "--autobox_ligand", grid_box,
    "-o", output_file,
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def batch_docking():
    receptor = f"{pdb_id}_A.pdbqt"
    ligand = f"ligands_to_dock.sdf"
    grid_box = f"{ligand_id}_corrected_pose.sdf"
    output_file = f"{docking_results_directory}/multiple_ligands_docked_{pdb_id}.sdf"
    subprocess.run([
    #"./gnina",
    "gnina",
    "-r", receptor,
    "-l", ligand,
    "--autobox_ligand", grid_box,
    "-o", output_file,
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def flexible_docking():
    receptor = f"{pdb_id}_A.pdbqt"
    ligand = f"ligands_to_dock.sdf"
    grid_box = f"{ligand_id}_corrected_pose.sdf"
    flexdist_ligand = f"{ligand_id}_corrected_pose.sdf"
    output_file = f"{docking_results_directory}/{ligand_id}_flex.sdf"
    subprocess.run([
    #"./gnina",
    "gnina",
    "-r", receptor,
    "-l", ligand,
    "--autobox_ligand", grid_box,
    "-o", output_file,
    "--flexdist_ligand", flexdist_ligand,
    "--flexdist", "3.59",
    "--seed", "0",
    "--exhaustiveness 64",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def unknown_site_docking():
    receptor = f"{pdb_id}_A.pdbqt"
    ligand = f"ligands_to_dock.sdf"
    grid_box = f"{pdb_id}_A.pdbqt"
    flexdist_ligand = f"{ligand_id}_corrected_pose.sdf"
    output_file = f"{docking_results_directory}/{ligand_id}_docked_whole_{pdb_id}.sdf"
    subprocess.run([
    #"./gnina",
    "gnina",
    "-r", receptor,
    "-l", ligand,
    "--autobox_ligand", grid_box,
    "-o", output_file,
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])
    
# ------------------------------------------------------------------------------
# 4.2 Main functions
# ------------------------------------------------------------------------------

def docking_main():
    # Redocking with extracted ligand
    if selection == "a":
        print("\n === Single ligand docking ===")
        single_ligand_docking()
        #rmsd_calculation()
        #report(docking_results = f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf")
    
    # Docking with multiple ligands
    elif selection == "b":
        print("\n === Batch docking ===")
        batch_docking()
        #report(docking_results = f"{docking_results_directory}/multiple_ligands_docked_{pdb_id}.sdf")
    
    # Flexible docking
    elif selection == "c":
        print("Flexible docking")
        flexible_docking()
        #report(docking_results = f"{docking_results_directory}/{ligand_id}_flex.sdf")
    
    # Docking on unknown site
    else:
        print("Docking on unknown site")
        unknown_site_docking()
        #report(docking_results = f"{docking_results_directory}/{ligand_id}_docked_whole_{pdb_id}.sdf")
        
#def report(docking_results):
    #score_columns = [
        #"minimizedAffinity",
        #"CNNscore",
        #"CNNaffinity",
        #"CNN_VS",
        #"CNNaffinity_variance",
    #]

    # Normalize input: accept either a string or a list of strings
    #if isinstance(docking_results, str):
        #sdf_paths = [docking_results]
    #else:
        #sdf_paths = docking_results
    
    #df_list = []
    #for filename in sdf_paths:
        #with BlockLogs():
            #df_list.append(PandasTools.LoadSDF(filename))
    
    #combo_df = pd.concat(df_list)
    
    # Convert score columns to float if they exist and are not entirely empty
    #for col in score_columns:
        #if col in combo_df.columns:
            # If the column is completely empty, drop it
            #if combo_df[col].isnull().all() or (combo_df[col] == "").all():
                #combo_df.drop(columns=[col], inplace=True)
            #else:
                #combo_df[col] = combo_df[col].astype(float)
    
    # Extract SMILES from the RDKit molecule column
    #combo_df["SMILES"] = combo_df["ROMol"].apply(lambda mol: Chem.MolToSmiles(mol) if mol is not None else None)
    
    # Save with SMILES included
    #output_csv = f"{docking_results_directory}/docking_results_{pdb_id}.csv"
    #combo_df.to_csv(output_csv, index=False)
    
    #print(f"Docking results with SMILES saved to {output_csv}")
    
    #combo_df.head()

# ------------------------------------------------------------------------------
# 4.3 Additional functions
# ------------------------------------------------------------------------------

#def rmsd_calculation():
    #print("\n === Running mcs (maximum common structure) rmsd calculation ===")
    #cognate = Chem.MolFromMolFile(f"{ligand_directory}/{ligand_id}_corrected_pose.sdf")
    #cognate = Chem.MolFromMolFile(f"{ligand_id}_corrected_pose.sdf")
    #poses = Chem.SDMolSupplier(f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf")
    #for i, pose in enumerate(poses):
        #RDLogger.DisableLog('rdApp.warning')
        #n_match, rmsd = uru.mcs_rmsd(cognate, pose)
        #print(f"{n_match}\t{rmsd:.2f}")


# Run gnina

os.makedirs("molecular_docking/docking_results", exist_ok=True)

#selection = input("\n === Welcome to gnina. Please select the docking mode: \
                    #\n single docking (rmsd and cnn score will be calculated): a \
                    #\n batch docking: b \
                    #\n flexible docking (maximum exhaustiveness): c \
                    #\n docking on unknown sites: d \
                #")

selection = os.getenv("PARAM_SELECTION")

inp = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80]

#ex = int(input("\n Define exhaustiveness (8, 16, 24, 32, 40, 48, 56, 64): "))
ex = os.getenv("PARAM_EX")
#ex = int(2)
if ex in inp:
    print(f"You have selected  exhaustiveness level {ex}")
    #ask_gpu = os.getenv("PARAM_ASK_GPU")
    #ask_cnn = os.getenv("PARAM_ASK_CNN")
    #ask_gpu = input("\n Run with GPU? [y/n]: ").strip().lower()
    #ask_cnn = input("\n Run with CNN? [y/n]: ").strip().lower()
    ask_gpu = 'n'
    ask_cnn = 'n'
    # Decide GPU flag
    gpu_flag = []
    if ask_gpu == 'y':
        gpu_flag = []  # run with GPU
    elif ask_gpu == 'n':
        gpu_flag = ["--no_gpu"]  # run without GPU
    else:
        print("\n Invalid input. Please try again with 'y' or 'n'.")
    
    # Decide CNN flag
    cnn_flag = []
    if ask_cnn == 'y':
        cnn_flag = []  # run with CNN
    elif ask_cnn == 'n':
        cnn_flag = ["--cnn_scoring=none"]  # run without CNN
    else:
        print("\n Invalid input. Please try again with 'y' or 'n'.")

    print(f"You have selected option {selection}.")

    if __name__ == "__main__":
        #pdb_id = "4OHU"
        pdb_id = os.getenv("PARAM_PDB_ID")
        #ligand_id = "2TK"
        ligand_id = os.getenv("PARAM_LIGAND_ID")

        docking_main()
else:
    print("\n === Invalid input. Please try again with the given options only! ===")
    print("\n === Exiting gnina ===")
    ex = ''




