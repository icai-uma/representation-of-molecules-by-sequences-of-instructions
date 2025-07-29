"""

Demonstration script for the IsalChem methodology:
Lopez-Rubio, Ezequiel (2024). Instruction set and language for chemical nomenclature
https://doi.org/10.26434/chemrxiv-2024-5b4dn

Coded by Ezequiel Lopez-Rubio, November 2024.

"""

# !pip install networkx python-pptx
import random

from isalchemmolecule import IsalChemMolecule
from isalchemutilities import IsalChemConverter
from isalchemutilities import IsalChemCompressor
from isalchemutilities import IsalChemComparison
import pickle

from isalchemutilities import levenshtein_distance_and_sequence
from isalchemutilities import fetch_neighbors, random_edit





# Examples for the preprint

# Propane
# token_list = ["C", "C", "C", "+", "-", ">", "<"]
# Methane
# token_list = ["C", ">", ">", ">", ">", ">"]
# Ammonia
# token_list = ["N", ">", ">", ">", ">", ">"]
# token_list = ["N", "N", "C", "O", "Cl", "Br", "C2"]
# token_list = ["C", "C2", "C", "+", "+", "+", "C2"]
# token_list = ["C", "C2", "C", "+", "+", "+", "N2", "C", "B3"]
# token_list = ["C", "C", "C", "C", "C", "J"]
# 1,4-hexadiene
# token_list = ["C", "C", "C2", "C", "C", "C2"]
# 3-propyl 4-isopropyl 1-heptene
# token_list = ["C", "C2", "C", "C", "C", "C", "C", "+", "+", "+", "+", "C", "C", "+", "+", "C", "+", "+", "C", "C", "C"]
# Cyclohexane
# token_list = ["C", "C", "C", "C", "C", "C", "J"]
# Isobutyric acid
# token_list = ["C", "C", "C", "O", "+", "O2", "+", "C"]
# Triethylamine
# token_list = ["C", "C", "N", "C", "C", "+", "+", "+", "C", "C"]
# 2,4,5-trichlorophenol
# token_list = ["C", ">", "C2", "C", "C2", "C", "C2", "J", "O2", "-", "Cl", "-", "Cl", "Cl"]
# Cyclopentene
# token_list = ["C", "C2", "-", "C", "C", "C", "J"]
# Ethanoic acid
# token_list = ["C", "C", "O2", "+", "O"]
# For debugging purposes
# token_list = ["C", "C", "O", "C", "Cl", "J2"]

# molecule = IsalChemMolecule(token_list)
# print(molecule)
# molecule.states[-1].generate_pptx("my_presentation_initial.pptx")

#for ndx_step in range(len(token_list)):
#    molecule.step()
#    print(molecule)
#    molecule.states[-1].generate_pptx(f"my_presentation_{ndx_step}.pptx")



# Example for massive generation of random strings
# random.seed(20)
# num_token_lists = 100
# max_token_list_length = 40
# instruction_set = ["A", "+", "-", ">", "<", "C", "N", "B", "O", \
#                    "C2", "N2", "B2", "O2", "C3", "N3", "B3", "J", "J2", "Cl", "Br"]
# for ndx_token_list in range(num_token_lists):
#     my_token_list_length = random.randint(0, max_token_list_length)
#     my_token_list = random.choices(instruction_set, k=my_token_list_length)
#     print(my_token_list)
#     molecule = IsalChemMolecule(my_token_list)
#     print(molecule)
#     for ndx_step in range(len(my_token_list)):
#         molecule.step()
#         print(molecule)
#     molecule.states[-1].generate_pptx(f"my_molecule_{ndx_token_list}.pptx")


# for ndx_step in range(len(token_list)+1):
#    molecule.states[ndx_step].generate_pptx(f"my_debug_{ndx_step}.pptx")


# Example for conversion of graph to string




#token_list = ["C", "C", "C", "O", "+", "O2", "+", "C"]
#token_list = ["C", "C", "N", "C", "C", "+", "+", "+", "C", "C"]
#token_list = ["C", "C", "O", "C", "Cl"]
#token_list = ["C", "C2", "C", "C", "C", "C", "C", "+", "+", "+", "+", "C", "C", "+", "+", "C", "+", "+", "C", "C", "C"]

#token_list = ["C", "C", "C", "C", "C", "C", "J"]

#token_list = ["C", "C2", "-", "C", "C", "C", "J"]

#token_list = ["C", "C", "O", "C", "Cl", "J2"]
#token_list = ['B3', '>', 'O2', 'C', '-', 'N', 'C3', 'C3', 'C3', 'J2', 'Cl', 'B2', 'Br', '<', 'Cl', 'B', 'C2', 'O2', 'C2', 'O', 'O', '>', '-', 'J2', '<', 'C3', 'B', 'O2', 'C3', 'B2', 'Br', 'N3', 'O2', '>', 'B', 'J', 'B2', 'B3', 'C', 'C3']


# Furan
#token_list = ["C", "C2", "-", "C", "C2", "O", "J"]

#Luteolin

token_list = ["C", "O", "-", "C", "C", "O", "-", "C", "C", "C", "<", "J2"]
token_list = ["C", "-", "O", "C2", "C", "O", "-", "C2", "C", "C2", "<", "-", "J2"]

token_list = ["C", "C2", "C", "C3"]
#, ">", ">", "+", "O", "C", "C2", "C", "J", "O2", "+", "+", "C", "+", "C", "C2", "O", "+", "C", "O", "-", "C2", "C", "J2"]

token_list = ["C", "C2", "-", "C", "C", "C", "N3"]





# Isovaleric acid, also known as 3-methylbutanoic acid or beta-methylbutyric acid
#token_list = ["C", "C", "C", "+", "+", "C", "C", "O2", "+", "O"]

# Isobutylamine
#token_list = ["C", "C", "C", "+", "+", "C", "N"]

# Propane
#token_list = ["C", "C", "C"]

# Propylamine
#token_list = ["C", "C", "C", "N"]

#token_list = ["C", "C", "C", "F"]

"""
# Example of neighbor generation
molecule = IsalChemMolecule(token_list)
for ndx_step in range(len(token_list)):
    molecule.step()
molecule.states[-1].generate_pptx(f"neighbors/original_molecule.pptx")

# Fetch the neighbors of the original molecule
neighbor_list = fetch_neighbors(token_list)

# Convert each token list to a tuple and use a set to filter duplicates
unique_token_tuples = list(set(tuple(my_token_list) for my_token_list in neighbor_list))

# Convert back to token lists
unique_neighbors = [list(tup) for tup in unique_token_tuples]

with open(f'neighbors/trace.txt', 'w') as file:
    for ndx_neighbor, neighbor in enumerate(unique_neighbors):
        neighbor_molecule = IsalChemMolecule(neighbor)
        for ndx_step in range(len(neighbor)):
            neighbor_molecule.step()
        neighbor_molecule.states[-1].generate_pptx(f"neighbors/neighbor_molecule{ndx_neighbor}.pptx")
        file.write(f"Neighbor {ndx_neighbor}: {neighbor}\n")
"""


"""
# Example of random edit
random.seed(20)
num_edits = 5
molecule = IsalChemMolecule(token_list)
for ndx_step in range(len(token_list)):
    molecule.step()
molecule.states[-1].generate_pptx(f"sequence/initial_molecule.pptx")
lev_distance, sequence_molecules = random_edit(token_list, num_edits)
print(lev_distance, sequence_molecules)
for ndx_sequence, sequence_token_list in enumerate(sequence_molecules):
    this_molecule = IsalChemMolecule(sequence_token_list)
    for ndx_step in range(len(sequence_token_list)):
        this_molecule.step()
    this_molecule.states[-1].generate_pptx(f"sequence/sequence_molecule{ndx_sequence}.pptx")
"""




"""
# Example of shortest path

# 4-bromopentanal
initial_token_list = ["C", "C", "Br", "C", "C", "C", "O2"]

# 2-chloro-2-hydroxyethanal
final_token_list = ["C", "Cl", "C", "O2", "O"]

lev_distance, sequence_token_lists = levenshtein_distance_and_sequence(initial_token_list, final_token_list)

print(lev_distance, sequence_token_lists)
for ndx_sequence, sequence_token_list in enumerate(sequence_token_lists):
    this_molecule = IsalChemMolecule(sequence_token_list)
    for ndx_step in range(len(sequence_token_list)):
        this_molecule.step()
    this_molecule.states[-1].generate_pptx(f"sequence/sequence_molecule{ndx_sequence}.pptx")

"""



token_list = ['J', 'C2', 'C3', 'N2', 'B2', 'N3', 'B3', 'B', 'N2', 'A', 'J2', 'Cl', 'O2', '-', 'N', '<', 'N3', 'N3', '-', 'N2', 'C3', 'J2', 'B', 'C2', 'N', 'C', '>']

token_list = ['J', 'C2', 'C3', 'N2', 'B2', 'N3', 'B3', 'B', 'N2', 'A']

# Flavone
token_list = [
  "C", "C", "C2", "C", "C2", "C", "J2",
  "+", "+", "C", "O2", "C", "C2", "O", "J",
  ">", ">", "C", "C2", "<", "<", "C", "C2", "C", "C2", "J"
]


# Example of conversion from graph to IsalChem string
molecule = IsalChemMolecule(token_list)
for ndx_step in range(len(token_list)):
    molecule.step()
print(molecule)
molecule.states[-1].generate_pptx(f"graph_to_string/input_molecule.pptx")
converter = IsalChemConverter(molecule.states[-1].graph)
converter.to_token_list()
print(converter)
converter.result.states[-1].generate_pptx(f"graph_to_string/output_molecule.pptx")
comparer = IsalChemComparison(molecule.states[-1].graph, converter.result.states[-1].graph)
print(f"Comparison result: {comparer.compare()}")

with open(f'trace.txt', 'w') as file:
    # A new file will be created
    file.write(converter.log_text)



"""
# Example for massive generation of random strings
random.seed(200)
num_token_lists = 10000
max_token_list_length = 40
instruction_set = ["A", "+", "-", ">", "<", "C", "N", "B", "O", \
                   "C2", "N2", "B2", "O2", "C3", "N3", "B3", "J", "J2", "Cl", "Br"]
for ndx_token_list in range(num_token_lists):
    my_token_list_length = random.randint(0, max_token_list_length)
    my_token_list = random.choices(instruction_set, k=my_token_list_length)
    #print(my_token_list)
    molecule = IsalChemMolecule(my_token_list)
    for ndx_step in range(len(my_token_list)):
        molecule.step()
    #if ndx_token_list==1000:
    #    molecule.states[-1].generate_pptx(f"my_input_molecule_{ndx_token_list}.pptx")
    converter = IsalChemConverter(molecule.states[-1].graph)
    try:
        converter.to_token_list()
        converter.validate()
        # converter.result.states[-1].generate_pptx(f"my_output_molecule_{ndx_token_list}.pptx")
        comparison = IsalChemComparison(molecule.states[-1].graph, converter.result.states[-1].graph)
        comparison_result = comparison.compare(True)
        if comparison_result == False:
            print(my_token_list)
            molecule.states[-1].generate_pptx(f"my_input_molecule_{ndx_token_list}.pptx")
            converter.result.states[-1].generate_pptx(f"my_output_molecule_{ndx_token_list}.pptx")
            # Open a file and use dump()
            with open(f'offending_{ndx_token_list}.pkl', 'wb') as file:
                # A new file will be created
                pickle.dump([molecule, converter, comparison], file)
            with open(f'offending_{ndx_token_list}.txt', 'w') as file:
                # A new file will be created
                file.write(converter.log_text)
            for ndx_state, state in enumerate(converter.result.states):
                state.generate_pptx(f"output_state_{ndx_state}.pptx")
            print("Conversion failed.")
            break
    except Exception as e:
        print(f"Error found: {e}, ndx_token_list={ndx_token_list}, my_token_list={my_token_list}")
    if ndx_token_list % 100 ==0:
        print(f"Converted {ndx_token_list}/{num_token_lists} molecules")
"""


"""
# Example for massive generation of pairs of random strings
random.seed(190)
num_token_list_pairs = 100
max_token_list_length = 40
instruction_set = ["A", "+", "-", ">", "<", "C", "N", "B", "O", \
                   "C2", "N2", "B2", "O2", "C3", "N3", "B3", "J", "J2", "Cl", "Br"]
for ndx_token_list in range(num_token_list_pairs):
    my_token_list_length1 = random.randint(0, max_token_list_length)
    my_token_list1 = random.choices(instruction_set, k=my_token_list_length1)
    my_token_list_length2 = random.randint(0, max_token_list_length)
    my_token_list2 = random.choices(instruction_set, k=my_token_list_length2)

    distance, transformation_sequence = levenshtein_distance_and_sequence(my_token_list1, my_token_list2)
    print(f"Inputs: {my_token_list1}\n{my_token_list2}")
    print(f"Levenshtein Distance: {distance}")
    print("Transformation Sequence:")
    for step in transformation_sequence:
        print(step)
"""


"""
token_list = ["C", "C2", "C", "C", "C", "C", "C", "+", "+", "+", "+", "C", "C", "+", "+", "C", "+", "+", "C", "C", "C"]
token_list = ['B3', '>', 'O2', 'C', '-', 'N', 'C3', 'C3', 'C3', 'J2', 'Cl', 'B2', 'Br', '<', 'Cl', 'B', 'C2', 'O2', 'C2', 'O', 'O', '>', '-', 'J2', '<', 'C3', 'B', 'O2', 'C3', 'B2', 'Br', 'N3', 'O2', '>', 'B', 'J', 'B2', 'B3', 'C', 'C3']
token_list = ["C", "C", "C", "O", "+", "O2", "+", "C"]
token_list = ["C", "C", "C", "C", "C", "C", "J"]
token_list = ["C", "C2", "-", "C", "C", "C", "J"]

# CUANDO LA LONGITUD DE LA TOKEN_LIST MAS CORTA ES 12, HAY UN FALLO EN LA CONVERSION
token_list = ["C", ">", "C2", "C", "C2", "C", "C2", "J", "O2", "-", "Cl", "-", "Cl", "Cl"]
# LA CADENA CONVERTIDA ERRONEA ES
#token_list = ['C', 'C2', 'C', 'C2', 'C', 'C2', 'Cl', 'Cl', '-', 'J2', 'Cl', 'O']



molecule = IsalChemMolecule(token_list)
for ndx_step in range(len(token_list)):
    molecule.step()
molecule.states[-1].generate_pptx(f"original_molecule.pptx")

compressor = IsalChemCompressor(molecule.states[-1].graph)
compressor.compress()
print(compressor.token_list_lengths)
print(compressor.compressed_token_list)

compressed_molecule = IsalChemMolecule(compressor.compressed_token_list)
for ndx_step in range(len(compressor.compressed_token_list)):
    compressed_molecule.step()
compressed_molecule.states[-1].generate_pptx(f"compressed_molecule.pptx")

"""




"""

molecule = IsalChemMolecule(token_list)
for ndx_step in range(len(token_list)):
    molecule.step()
    molecule.states[-1].generate_pptx(f"step_{ndx_step}_molecule.pptx")
print(molecule)
molecule.states[-1].generate_pptx(f"output_molecule.pptx")
"""