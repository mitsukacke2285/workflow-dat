#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract the ligand from a PDB file and save it as ligand_name.sdf
(Preserves original coordinates from the PDB)
"""

import os
import requests
from Bio.PDB import PDBIO, PDBParser, Select
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import MDAnalysis as mda

def extract_ligand_from_pdb(pdb_path, ligand_name):
    """Extract the ligand from a PDB file and save as SDF (keeping original conformation)"""
    if not os.path.exists(pdb_path):
        print(f"File not found: {pdb_path}")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Load the original PDB
    u = mda.Universe(f"{pdb_id}_A_NAD.pdb")
    single_ligand = u.select_atoms(f"resname {ligand_name}")
    single_ligand.write(f"{ligand_name}_fromPDB.pdb")
    print(f"Ligand {ligand_name} extracted from original PDB!")

    print(f"Extracted {ligand_name} ligand and saved as: {ligand_name}.sdf")
    
def fix_and_align():
    
    print("Fixing and correcting the pose of extracted ligand...")
    # Load and remove the hydrogens that we can't map anyways yet.
    ideal_mol = Chem.MolFromMolFile(f"{ligand_name}_ideal.sdf", removeHs=True)
    pose_mol = Chem.MolFromPDBFile(f"{ligand_name}_fromPDB.pdb", removeHs=True)
    
    # Disconnect any organometal
    rdMolStandardize.DisconnectOrganometallicsInPlace(pose_mol)

    # Remove disconnected fragments
    fragmenter = rdMolStandardize.FragmentRemover()
    pose_mol_f = fragmenter.remove(pose_mol)

    # Choose largest fragment
    chooser = rdMolStandardize.LargestFragmentChooser()
    pose_mol_lf = chooser.choose(pose_mol_f)

    # Assign bond orders from the template to the pose molecule
    corrected_pose = AllChem.AssignBondOrdersFromTemplate(ideal_mol, pose_mol_lf)

    # Add hydrogens
    corrected_pose_with_H = Chem.AddHs(corrected_pose, addCoords=True)

    # Save the corrected pose to an SDF file
    ligand_corrected_pose_file = f"{ligand_name}_corrected_pose.sdf"
    writer = Chem.SDWriter(ligand_corrected_pose_file)
    writer.write(corrected_pose_with_H)
    writer.close()
    print("Extracted ligand fixed!")

def download_ideal_ligand():

    ideal_ligand_filename = f"{ligand_name}_ideal.sdf"
    print(f"Downloading ligand {ligand_name}...")
    ligand_url = f"https://files.rcsb.org/ligands/download/{ideal_ligand_filename}"
    ligand_request = requests.get(ligand_url)
    ligand_request.raise_for_status() # Check for errors

    ideal_filepath = f"{ideal_ligand_filename}"

    with open(ideal_filepath, "w") as f:
        f.write(ligand_request.text)
    print(f"Saved ligand to {ideal_filepath}")

if __name__ == "__main__":
    # pdb_id = "4OHU"
    #pdb_id = os.getenv("PARAM_PDB_ID")
    # ligand_name = "2TK"
    #ligand_name = os.getenv("PARAM_LIGAND_NAME")
    extract_ligand_from_pdb(f"{pdb_id}_A_NAD.pdb", ligand_name)
    download_ideal_ligand()
    fix_and_align()
