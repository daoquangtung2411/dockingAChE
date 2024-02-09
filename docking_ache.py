#! /usr/local/bin/env python

from rdkit import Chem
from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import subprocess
import logging
import csv

# Get information from config file

protein = '6tt0'

# Vina define

v = Vina(sf_name='vina', verbosity=1)

# Define Receptor

v.set_receptor(f'{protein}.pdbqt')

center_x = -4.556

center_y = 55.907

center_z = 299.157

box_size = 28

exhaustiveness = 8

num_modes = 10

is_macrocycles = True # state True if ligands is macrocyclyes, False for others

with open('NB_docking_list_5.csv','r') as docking_file:

    NB_docking_data = csv.reader(docking_file, delimiter=" ")

    for row in NB_docking_data:

        ligand = row[0].split(',')

        ligand_name = ligand[0]

        ligand_SMILES = ligand[1]

        # Prepare ligand with rdkit

        # Protonate smiles

        cmd = f'obabel -:"{ligand_SMILES}" -ismi -ocan -p{7.4}'

        cmd_return = subprocess.run(cmd, capture_output=True, shell=True)

        output = cmd_return.stdout.decode('utf-8')

        logging.debug(output)

        # Convert SMILES into Mol

        lig = Chem.MolFromSmiles(output.strip(), sanitize=True)

        if lig is not None:

            # Add Hydrogens

            protonated_lig = Chem.rdmolops.AddHs(lig)

            # Generate 3D coordinate using ETKDG version 3

            if is_macrocycles:

                etkdg_params = Chem.rdDistGeom.ETKDGv3()

            else:

                etkdg_params = Chem.rdDistGeom.srETKDGv3()

            Chem.rdDistGeom.EmbedMolecule(protonated_lig, etkdg_params)

            # Forcefield optimized

            Chem.rdForceFieldHelpers.MMFFSanitizeMolecule(protonated_lig)

            Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(protonated_lig, mmffVariant='MMFF94s', maxIters=10000)

            # Compute Gasteiger Charges

            Chem.rdPartialCharges.ComputeGasteigerCharges(protonated_lig)

            # Process ligand with Meeko

            # Ligand Prep

            preparator = MoleculePreparation(
                merge_these_atom_types=("H",),
                hydrate=False,
                rigid_macrocycles=True)

            mol_setups = preparator.prepare(protonated_lig)

            for setup in mol_setups:

                pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)

                if is_ok:

                    ligand_file = open(f'{ligand_name}.pdbqt', 'w')

                    content = pdbqt_string

                    n = ligand_file.write(content)

                    if n == len(content):

                        print("Success")

                    else:

                        print("Failed")

                    ligand_file.close()

            # Dock the ligand

            cmd2 = f'./mvina --receptor {protein}.pdbqt --ligand {ligand_name}.pdbqt --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {bx} --size_y {bx} --size_z {bx} --out {ligand_name}_{protein}_{bx}.pdbqt --log {ligand_name}_{protein}_{bx}.txt --num_modes {num_modes} --exhaustiveness {exhaustiveness}'

            cmd_return = subprocess.call(cmd2, shell=True)