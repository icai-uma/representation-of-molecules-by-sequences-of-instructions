"""

Demonstration script with usage of the rdkit library for the IsalChem methodology:
Lopez-Rubio, Ezequiel (2024). Instruction set and language for chemical nomenclature
https://doi.org/10.26434/chemrxiv-2024-5b4dn

Coded by Ezequiel Lopez-Rubio, February 2025.

"""

# Necessary to load Graph instances with pickle
from datastructures import Graph, Node
import rdkit
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


from isalchemmolecule import IsalChemMolecule
from isalchemutilities import random_edit, graph_to_mol, restricted_random_edit, smiles_tokenizer, \
    mol_to_graph, IsalChemCompressor
import pickle
import numpy as np
import pandas as pd


"""
# Conversion of a single SMILES string to IsalChem string

smiles_str = "C1=CC=C(C=C1)C2=CC(=O)C3=CC=CC=C3O2"
mol = Chem.MolFromSmiles(smiles_str)
molecule_graph = mol_to_graph(mol)
compressor = IsalChemCompressor(molecule_graph)
compressor.compress()
token_list = compressor.compressed_token_list
print(token_list)
"""

"""
# Batch conversion of SMILES strings to Graph objects

# Path to the CSV file
file_path = '../Zinc/250k_rndm_zinc_drugs_clean_3.csv'
# Load the CSV file
smiles_df = pd.read_csv(file_path)

# Convert SMILES strings to RDKit Mol objects and then to Graph objects
molecules_graphs = []
molecules_smiles_strings = []
molecule_counter = 0
for smiles in smiles_df.iloc[:, 0]:  # Assume SMILES strings are in the first column
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        graph = mol_to_graph(mol)
        if graph is not None:
          molecules_graphs.append(graph)
          molecules_smiles_strings.append(smiles)
          molecule_counter = molecule_counter + 1
    else:
        print(f"Invalid SMILES string: {smiles}")
    if molecule_counter >= 1000:
      break

print(f"Converted {len(molecules_graphs)} molecules to Graph objects.")

# Save the list of Graph objects to a pickle file
pickle_file_path = "molecules_graphs.pkl"
with open(pickle_file_path, "wb") as f:
    pickle.dump((molecules_graphs, molecules_smiles_strings), f)

print(f"Saved Graph objects to {pickle_file_path}")
"""



# Comparison of Levenshtein distance versus Tanimoto distance for random edits on Zinc molecules
# Also compare the number of tokens of SMILES versus IsalChem strings
# https://stackoverflow.com/questions/51681659/how-to-use-rdkit-to-calculte-molecular-fingerprint-and-similarity-of-a-list-of-s

# Load the Graph instances for some Zinc molecules
pickle_file_path = "molecules_graphs_1k.pkl"
with open(pickle_file_path, "rb") as f:
    molecules_graphs, molecules_smiles_strings = pickle.load(f)

num_edits = 21
num_samples = 500
similarity_values = {
    "tanimoto": [],
    "dice": [],
    "cosine": [],
    "sokal": [],
    "russel": [],
    "kulczynski": [],
    "mcconnaughey": []
}
similarity_functions = {
    "tanimoto": DataStructs.BulkTanimotoSimilarity,
    "dice": DataStructs.BulkDiceSimilarity,
    "cosine": DataStructs.BulkCosineSimilarity,
    "sokal": DataStructs.BulkSokalSimilarity,
    "russel": DataStructs.BulkRusselSimilarity,
    "kulczynski": DataStructs.BulkKulczynskiSimilarity,
    "mcconnaughey": DataStructs.BulkMcConnaugheySimilarity
}

for similarity in similarity_values:
    for ndx_list in range(num_edits+1):
        similarity_values[similarity].append([])

smiles_num_tokens = []
isalchem_num_tokens = []

for ndx_molecule, molecule_graph in enumerate(molecules_graphs[:num_samples]):
    # print(f"Processing molecule {ndx_molecule}: {molecules_smiles_strings[ndx_molecule]}")

    # Ions cannot be processed at this point because they are not modeled by IsalChem
    if "+" in molecules_smiles_strings[ndx_molecule] or "-" in molecules_smiles_strings[ndx_molecule]:
        continue

    # Convert Graph instances to compressed token lists
    compressor = IsalChemCompressor(molecule_graph)
    compressor.compress()
    token_list = compressor.compressed_token_list
    # print(token_list)
    smiles_num_tokens.append(len(smiles_tokenizer(molecules_smiles_strings[ndx_molecule])))
    isalchem_num_tokens.append(len(token_list))

    # Convert compressed token lists into IsalChemMolecule instances
    molecule = IsalChemMolecule(token_list)
    molecule.run()

    rdkit_molecules = []
    #rdkit_molecules.append(graph_to_mol(molecule.states[-1].graph))
    #print(MolToSmiles(rdkit_molecules[0]))

    # Perform some random edits sequentially

    #lev_distance, sequence_molecules = random_edit(token_list, num_edits)
    lev_distance, sequence_molecules = restricted_random_edit(token_list, num_edits, abundances=True)
    #print(lev_distance, sequence_molecules)
    for ndx_sequence, sequence_token_list in enumerate(sequence_molecules):
        this_molecule = IsalChemMolecule(sequence_token_list)
        for ndx_step in range(len(sequence_token_list)):
            this_molecule.step()
        rdkit_molecules.append(graph_to_mol(this_molecule.states[-1].graph))
        #print(MolToSmiles(rdkit_molecules[-1]))

    # make a list of fingerprints
    molecular_fingerprints = [FingerprintMols.FingerprintMol(x) for x in rdkit_molecules]

    # Compare all fingerprints according to the different similarities
    for similarity in similarity_functions:
        s = similarity_functions[similarity](molecular_fingerprints[0], molecular_fingerprints)
        for ndx_dist in range(len(s)):
            similarity_values[similarity][ndx_dist].append(s[ndx_dist])

    if ndx_molecule % 20 == 19:
        print(f"{ndx_molecule+1} molecules have been processed.")
        arr_smiles_num_tokens = np.array(smiles_num_tokens)
        arr_isalchem_num_tokens = np.array(isalchem_num_tokens)
        print(f"Number of tokens for SMILES: mean={np.mean(arr_smiles_num_tokens)}, std={np.std(arr_smiles_num_tokens)}")
        print(f"Number of tokens for IsalChem: mean={np.mean(arr_isalchem_num_tokens)}, std={np.std(arr_isalchem_num_tokens)}")
        for similarity in similarity_functions:
            for ndx_dist in range(num_edits+1):
                my_similarities = np.array(similarity_values[similarity][ndx_dist])
                median_similarity = np.median(my_similarities)
                mean_similarity = np.mean(my_similarities)
                std_similarity = np.std(my_similarities)
                print(f"Levenshtein distance={ndx_dist}, {similarity} similarity: median ={median_similarity}, mean ={mean_similarity}, std ={std_similarity}")

