"""

Utilities for the IsalChem methodology:
Lopez-Rubio, Ezequiel (2024). Instruction set and language for chemical nomenclature
https://doi.org/10.26434/chemrxiv-2024-5b4dn

Coded by Ezequiel Lopez-Rubio, November 2024.

"""
from isalchemmolecule import IsalChemMolecule
from datastructures import Graph, Node
import numpy as np
import re

from collections import Counter

try:
    from rdkit import Chem
    from rdkit.Chem import RWMol
    from rdkit.Chem import Atom, BondType
except ImportError:
    print("The rdkit library is not available")


class IsalChemConverter:
    """ Conversion of a molecular graph into an IsalChem token list """
    def __init__(self, input_graph: Graph, start_nonhydrogen_node_id: int = 0):
        self.input_graph = input_graph      # Input molecular graph object of type Graph to be converted to a token list
        self.start_nonhydrogen_node_id = start_nonhydrogen_node_id # Node ID of the input nonhydrogen atom to start the conversion
        self.result = IsalChemMolecule([])   # The output IsalChemMolecule object to be built by the conversion
        # Dictionary to translate node IDs from self.input_graph
        # to Node objects in self.result
        self.input_to_output = dict()
        # Dictionary to translate node IDs from self.result
        # to Node objects in self.input_graph
        self.output_to_input = dict()
        self.log_text = ""

    def __repr__(self) -> str:
        input_repr = self.input_graph.__repr__()
        result_repr = self.result.__repr__()
        return f"IsalChemConverter(input_graph={input_repr}, result={result_repr}, input_to_output={self.input_to_output}, output_to_input={self.output_to_input})"

    def find_suitable_neighbors(self, x_input: Node, x_output: Node, output_to_input: dict, input_to_output: dict):
        """ Find the neighbors of x in self.input_graph
        that are suitable for adding an instruction to link them to x in self.result
        If the neighbor of x is not in self.result, then it is suitable.
        If the neighbor of x is in self.result, but the bond to x has not been added yet, then it is suitable.
        A set of Node instances is returned.

        """
        self.log_text += f"Called find_suitable_neighbors with x_input={x_input}, x_output={x_output}, input_to_output={input_to_output}\n"
        suitable_neighbors = set()
        # Consider all neighbors of x_input in self.input_graph
        for x_input_neighbor in x_input.neighbors:
            # Check if x_input_neighbor is already in self.result
            if x_input_neighbor.node_id in input_to_output:
                # x_input_neighbor is already in self.result
                x_output_neighbor = input_to_output[x_input_neighbor.node_id]
                # x_input_neighbor is a suitable neighbor only if the bond with x_input
                # is not present in self.result
                if x_output_neighbor not in x_output.neighbors:
                    suitable_neighbors.add(x_input_neighbor)
                    self.log_text += f"Added {x_input_neighbor} to suitable_neighbors, already in result, with x_output_neighbor={x_output_neighbor}, x_output.neighbors={x_output.neighbors}\n"
            else:
                # x_input_neighbor is not in self.result, so it is a suitable neighbor
                suitable_neighbors.add(x_input_neighbor)
                self.log_text += f"Added {x_input_neighbor} to suitable_neighbors, not in result, with x_input_neighbor.node_id={x_input_neighbor.node_id}\n"
        return suitable_neighbors


    def select_nonhydrogens(self, output_to_input: dict, input_to_output: dict, waiting_nonhydrogens: set):
        """
        Select the non-hydrogen atoms r in the waiting set of nonhydrogens and
        x in self.result.graph, such that r and x are
        neighbors in self.input_graph but not neighbors in self.result.graph,
        and one of the hydrogen atoms bound to x in self.result.graph
        is closest to the primary pointer self.result
        """

        # The current state
        current_state = self.result.states[-1]

        # Fetch the hydrogen atom that is pointed by the primary pointer
        primary_hydrogen = current_state.graph.get_node(current_state.primary_ptr)

        # Check if the nonhydrogen bound to the primary pointer hydrogen is suitable
        if primary_hydrogen is None:
            print(f"current_state={current_state}")
            print(f"waiting_nonhydrogens={waiting_nonhydrogens}")
            print(f"self.input_graph={self.input_graph}")
        x_output = primary_hydrogen.neighbors[0]
        x_input = output_to_input[x_output.node_id]
        suitable_neighbors = self.find_suitable_neighbors(x_input, x_output, output_to_input, input_to_output)
        self.log_text += f"Called find_suitable_neighbors with suitable_neighbors={suitable_neighbors}\n"
        r_candidates = waiting_nonhydrogens.intersection(suitable_neighbors)
        # r_candidates = waiting_nonhydrogens.intersection(set(x_input.neighbors) - set(x_output.neighbors))
        self.log_text += f"Initialized select_nonhydrogens: output_to_input={output_to_input}, primary_hydrogen={primary_hydrogen}, x_input={x_input}, x_output={x_output}, r_candidates={r_candidates}\n"
        if len(r_candidates)>0:
            r_input = r_candidates.pop()
            # Find out the edge type to be inserted (single, double, or triple bond)
            selected_edge_ndx = r_input.neighbors.index(x_input)
            selected_edge_type = r_input.edge_types[selected_edge_ndx]
            # Number of required primary pointer moves
            num_moves = 0
            if (selected_edge_type==1 and x_output==primary_hydrogen.next.neighbors[0] and
                x_output==primary_hydrogen.prev.neighbors[0]):
                # In case that this nonydrogen has a double bond yet to be created in the output,
                # we must avoid leaving two nonconsecutive hydrogens in the list that are attached to the
                # same nonhydrogen
                all_edge_types_x_input = set(x_input.edge_types)
                if 2 in all_edge_types_x_input:
                    num_moves = 1
            elif selected_edge_type>=2 and x_output!=primary_hydrogen.next.neighbors[0]:
                num_moves = -1
            elif selected_edge_type==3 and x_output!=primary_hydrogen.prev.neighbors[0]:
                num_moves = 1
            self.log_text += f"No moves required in select_nonhydrogens: r_input={r_input}, r_candidates={r_candidates}\n"
        else:
            # The nonhydrogen bound to the primary pointer hydrogen is not suitable
            # Therefore a search is launched for the first suitable nonhydrogen
            # which is carried out on the doubly linked list
            left_pointer = primary_hydrogen.prev
            right_pointer = primary_hydrogen.next
            num_moves = 1 # Number of required primary pointer movement instructions
            self.log_text += f"Prepared main loop in select_nonhydrogens: left_pointer={left_pointer}, right_pointer={right_pointer}\n"
            while True:
                # See if the left pointer is bound to a suitable nonhydrogen
                x_output = left_pointer.neighbors[0]
                x_input = output_to_input[x_output.node_id]
                suitable_neighbors = self.find_suitable_neighbors(x_input, x_output, output_to_input, input_to_output)
                self.log_text += f"Called find_suitable_neighbors with suitable_neighbors={suitable_neighbors}\n"
                r_candidates = waiting_nonhydrogens.intersection(suitable_neighbors)
                # r_candidates = waiting_nonhydrogens.intersection(set(x_input.neighbors) - set(x_output.neighbors))
                self.log_text += f"Trying left path in select_nonhydrogens: x_input={x_input}, x_output={x_output}, \
                    set(x_input.neighbors)={set(x_input.neighbors)}, \
                    set(x_output.neighbors)={set(x_output.neighbors)}, \
                    waiting_nonhydrogens={waiting_nonhydrogens}, \
                    r_candidates={r_candidates}\n"
                if len(r_candidates) > 0:
                    r_input = r_candidates.pop()
                    # Find out the edge type to be inserted (single, double, or triple bond)
                    selected_edge_ndx = r_input.neighbors.index(x_input)
                    selected_edge_type = r_input.edge_types[selected_edge_ndx]
                    # Signed number of required primary pointer moves
                    num_moves = -num_moves
                    if selected_edge_type >= 2:
                        num_moves -= 1
                    self.log_text += f"Left path chosen in select_nonhydrogens: r_input={r_input}, r_candidates={r_candidates}, num_moves={num_moves}\n"
                    break
                else:
                    left_pointer = left_pointer.prev
                # See if the right pointer is bound to a suitable nonhydrogen
                x_output = right_pointer.neighbors[0]
                x_input = output_to_input[x_output.node_id]
                suitable_neighbors = self.find_suitable_neighbors(x_input, x_output, output_to_input, input_to_output)
                self.log_text += f"Called find_suitable_neighbors with suitable_neighbors={suitable_neighbors}\n"
                r_candidates = waiting_nonhydrogens.intersection(suitable_neighbors)
                # r_candidates = waiting_nonhydrogens.intersection(set(x_input.neighbors) - set(x_output.neighbors))
                self.log_text += f"Trying right path in select_nonhydrogens: x_input={x_input}, x_output={x_output}, \
                    set(x_input.neighbors)={set(x_input.neighbors)}, \
                    set(x_output.neighbors)={set(x_output.neighbors)}, \
                    waiting_nonhydrogens={waiting_nonhydrogens}, \
                    r_candidates={r_candidates}\n"
                if len(r_candidates) > 0:
                    r_input = r_candidates.pop()
                    # Find out the edge type to be inserted (single, double, or triple bond)
                    selected_edge_ndx = r_input.neighbors.index(x_input)
                    selected_edge_type = r_input.edge_types[selected_edge_ndx]
                    # Number of required primary pointer moves
                    if selected_edge_type == 3:
                        num_moves += 1
                    self.log_text += f"Right path chosen in select_nonhydrogens: r_input={r_input}, r_candidates={r_candidates}, num_moves={num_moves}\n"
                    break
                else:
                    right_pointer = right_pointer.next
                num_moves += 1

        return num_moves, r_input, x_input, x_output, selected_edge_type


    def prepare_secondary_moves(self, r_output: Node, selected_edge_type: int):
        """
            Compute the signed number of moves required to move the secondary pointer
            to the hydrogen atom bound to r in self.result
            selected_edge_type is the type of bond: 1 for simple bond (J instruction) and 2 for double bond (J2 instruction)
        """
        # The current state
        current_state = self.result.states[-1]

        # Fetch the hydrogen atom that is pointed by the secondary pointer
        secondary_hydrogen = current_state.graph.get_node(current_state.secondary_ptr)

        self.log_text += f"Initialized prepare_secondary_moves: r_output={r_output}, secondary_hydrogen={secondary_hydrogen}, current_state={current_state}\n"
        # Check if the nonhydrogen bound to the secondary pointer hydrogen is r
        if secondary_hydrogen.neighbors[0].node_id == r_output.node_id:
            # Number of required secondary pointer moves
            num_moves = 0
            if (selected_edge_type==1 and secondary_hydrogen.neighbors[0]==secondary_hydrogen.next.neighbors[0] and
                secondary_hydrogen.neighbors[0]==secondary_hydrogen.prev.neighbors[0]):
                # In case that this nonhydrogen has a double bond yet to be created in the output,
                # this avoids leaving two nonconsecutive hydrogens in the list that are attached to the
                # same nonhydrogen
                non_hydrogen_input = self.output_to_input[secondary_hydrogen.neighbors[0].node_id]
                all_edge_types_non_hydrogen_input = set(non_hydrogen_input.edge_types)
                if 2 in all_edge_types_non_hydrogen_input:
                    num_moves = 1
            if selected_edge_type==2 and secondary_hydrogen.neighbors[0]!=secondary_hydrogen.next.neighbors[0]:
                num_moves -= 1
            self.log_text += f"No moves required in prepare_secondary_moves: secondary_hydrogen.neighbors[0].node_id={secondary_hydrogen.neighbors[0].node_id}, r_output.node_id={r_output.node_id}\n"
        else:
            # The nonhydrogen bound to the secondary pointer hydrogen is not r
            # Therefore a search is launched for r
            # which is carried out on the doubly linked list
            left_pointer = secondary_hydrogen.prev
            right_pointer = secondary_hydrogen.next
            num_moves = 1  # Number of required secondary pointer movement instructions
            self.log_text += f"Prepared main loop in prepare_secondary_moves: left_pointer={left_pointer}, right_pointer={right_pointer}\n"
            while True:
                # See if the left pointer is bound to r
                nonhydrogen = left_pointer.neighbors[0]
                self.log_text += f"Trying left path in prepare_secondary_moves: left_pointer={left_pointer}, nonhydrogen={nonhydrogen}\n"
                if nonhydrogen.node_id == r_output.node_id:
                    num_moves = -num_moves
                    if selected_edge_type == 2:
                        num_moves -= 1
                    self.log_text += f"Left path chosen in prepare_secondary_moves: nonhydrogen.node_id={nonhydrogen.node_id}, num_moves={num_moves}\n"
                    break
                else:
                    left_pointer = left_pointer.prev
                # See if the right pointer is bound to r
                nonhydrogen = right_pointer.neighbors[0]
                self.log_text += f"Trying right path in prepare_secondary_moves: right_pointer={right_pointer}, nonhydrogen={nonhydrogen}\n"
                if nonhydrogen.node_id == r_output.node_id:
                    self.log_text += f"Right path chosen in prepare_secondary_moves: nonhydrogen.node_id={nonhydrogen.node_id}, num_moves={num_moves}\n"
                    break
                else:
                    right_pointer = right_pointer.next
                num_moves += 1

        return num_moves

    def to_token_list(self, debug_output=False):
        """ Convert the molecular graph to a token list """

        non_hydrogen = None
        if self.start_nonhydrogen_node_id != 0:
            # A starting nonhydrogen atom has been specified
            if self.start_nonhydrogen_node_id in self.input_graph.nodes:
                node = self.input_graph.get_node(self.start_nonhydrogen_node_id)
                if node.data != "H":
                    non_hydrogen = node
                else:
                    raise ValueError(f"The specified node {self.start_nonhydrogen_node_id} is a hydrogen atom in {self.input_graph}")
            else:
                raise ValueError(f"The specified node {self.start_nonhydrogen_node_id} does not exist in {self.input_graph}")
        else:
            # Choose a non hydrogen atom from the input graph
            for node_id in self.input_graph.nodes:
                node = self.input_graph.get_node(node_id)
                if node.data != "H":
                    non_hydrogen = node
                    break

        # For debugging purposes
        if debug_output:
            self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")

        # The input graph does not contain non-hydrogen atoms, so conversion is finished
        if non_hydrogen is None:
            return self.result.token_list

        # Append the single-bond instruction associated to the non-hydrogen atom to the token list,
        # and execute the corresponding instruction on the current state
        self.result.token_list.append(non_hydrogen.data)
        self.result.step()
        # For debugging purposes
        # self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")

        # Initialize the dictionary to translate node IDs from self.input_graph
        # to Node objects in self.result
        self.input_to_output.update({non_hydrogen.node_id: self.result.last_nonhydrogen})

        # Initialize the dictionary to translate node IDs from self.result
        # to Node objects in self.input_graph
        self.output_to_input.update({self.result.last_nonhydrogen.node_id: non_hydrogen})

        # Initialize the set of non-hydrogen atoms waiting to be inserted into self.result
        waiting_nonhydrogens = set()
        for node in non_hydrogen.neighbors:
            if node.data != "H":
                waiting_nonhydrogens.add(node)

        self.log_text += f"Initialized to_token_list: waiting_nonhydrogens={waiting_nonhydrogens}\n"

        # Iterate while there are non-hydrogen atoms waiting to be inserted
        while len(waiting_nonhydrogens)>0:
            # Select the non-hydrogen atoms r in the waiting set of nonhydrogens and
            # x in self.result.graph

            self.log_text += f"\n\n\n\nStarted main loop iteration in to_token_list: waiting_nonhydrogens={waiting_nonhydrogens}, current_state={len(self.result.states)-1}\n"
            num_moves, r_input, x_input, x_output, selected_edge_type = \
                self.select_nonhydrogens(self.output_to_input, self.input_to_output, waiting_nonhydrogens)
            self.log_text += f"Called select_nonhydrogens in to_token_list: num_moves={num_moves}, r_input={r_input}, x_input={x_input}, x_output={x_output}\n"

            # Insert primary pointer movement instructions +/- into the token list as required
            # to move the primary pointer to the hydrogen atom bound to x in self.result
            # Execute such movement instructions on self.result
            if num_moves >= 0:
                for ndx_move in range(num_moves):
                    # Append a primary pointer instruction
                    # and execute the corresponding instruction on the current state
                    self.result.token_list.append("+")
                    self.result.step()
                    # For debugging purposes
                    if debug_output:
                        self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")
            else:
                for ndx_move in range(-num_moves):
                    # Append a primary pointer instruction
                    # and execute the corresponding instruction on the current state
                    self.result.token_list.append("-")
                    self.result.step()
                    # For debugging purposes
                    if debug_output:
                        self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")

            # Query whether r is already in self.result
            if r_input.node_id in self.input_to_output:
                # r is already in self.result so a join instruction must be executed
                r_output = self.input_to_output[r_input.node_id]

                # If x was also in the waiting set of nonydrogens, then it must be removed
                if x_input in waiting_nonhydrogens:
                    waiting_nonhydrogens.remove(x_input)

                self.log_text += f"Join instruction must be executed in to_token_list: num_moves={num_moves}, r_input={r_input}, r_output={r_output}, x_input={x_input}, x_output={x_output}\n"

                # Compute the number of steps to move the secondary pointer to the hydrogen
                # atom bound to r in self.result
                num_moves = self.prepare_secondary_moves(r_output, selected_edge_type)
                self.log_text += f"Called prepare_secondary_moves in to_token_list: num_moves={num_moves}\n"

                # Insert secondary pointer movement instructions >/< into the token list as required
                # to move the secondary pointer to the hydrogen atom bound to r in self.result
                # Execute such movement instructions on self.result
                if num_moves >= 0:
                    for ndx_move in range(num_moves):
                        # Append a secondary pointer instruction
                        # and execute the corresponding instruction on the current state
                        self.result.token_list.append(">")
                        self.result.step()
                        # For debugging purposes
                        if debug_output:
                            self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")
                else:
                    for ndx_move in range(-num_moves):
                        # Append a secondary pointer instruction
                        # and execute the corresponding instruction on the current state
                        self.result.token_list.append("<")
                        self.result.step()
                        # For debugging purposes
                        if debug_output:
                            self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")

                # Append the join instruction associated with r
                # (depending on the bond type between x and r in self.input_graph) into
                # the token list
                if selected_edge_type == 1:
                    self.result.token_list.append("J")
                elif selected_edge_type == 2:
                    self.result.token_list.append("J2")
                else:
                    raise ValueError(f"Invalid edge type when appending join instruction: {selected_edge_type}")

                # Execute the instruction on the current state. This inserts the cycle into self.result
                self.result.step()

            else:
                # r is not in self.result yet so an insertion instruction must be executed

                # Append the instruction associated with r
                # (depending on the bond type between x and r in self.input_graph) into
                # the token list
                if selected_edge_type == 1:
                    self.result.token_list.append(r_input.data)
                elif selected_edge_type == 2:
                    self.result.token_list.append(r_input.data+"2")
                elif selected_edge_type == 3:
                    self.result.token_list.append(r_input.data + "3")

                # Execute the instruction on the current state. This inserts r into self.result
                self.result.step()

                # Update the dictionary to translate node IDs from self.input_graph
                # to Node objects in self.result
                self.input_to_output.update({r_input.node_id: self.result.last_nonhydrogen})

                # Update the dictionary to translate node IDs from self.result
                # to Node objects in self.input_graph
                self.output_to_input.update({self.result.last_nonhydrogen.node_id: r_input})

                # Obtain the newly inserted node
                r_output = self.input_to_output[r_input.node_id]

            # For debugging purposes
            if debug_output:
                self.result.states[-1].generate_pptx(f"output_molecule{len(self.result.states)}.pptx")

            # Remove r from the set of nonhydrogens waiting to be inserted
            waiting_nonhydrogens.remove(r_input)


            # Refresh r_output in the new current state
            r_output = self.result.states[-1].graph.get_node(r_output.node_id)

            # Insert into waiting_nonhydrogens all the non-hydrogen neighbors
            # of r in self.input_graph that are not already neighbors of r in
            # self.result
            non_hydrogen_r_input_neighbors = set()
            for node in r_input.neighbors:
                if node.data != "H":
                    non_hydrogen_r_input_neighbors.add(node)

            non_hydrogen_r_output_neighbors = set()
            for node in r_output.neighbors:
                if node.data != "H":
                    non_hydrogen_r_output_neighbors.add(self.output_to_input[node.node_id])
            new_waiting_nonhydrogens = non_hydrogen_r_input_neighbors - non_hydrogen_r_output_neighbors
            #if len(new_waiting_nonhydrogens.intersection(waiting_nonhydrogens)) > 0:
            #    self.log_text += f"Duplication updating waiting_nonhydrogens: new_waiting_nonhydrogens={new_waiting_nonhydrogens}, waiting_nonhydrogens={waiting_nonhydrogens}\n"
            waiting_nonhydrogens = waiting_nonhydrogens.union(new_waiting_nonhydrogens)
            self.log_text += f"Updated waiting_nonhydrogens in to_token_list: waiting_nonhydrogens={waiting_nonhydrogens}, non_hydrogen_r_input_neighbors={non_hydrogen_r_input_neighbors}, r_output={r_output}, non_hydrogen_r_output_neighbors={non_hydrogen_r_output_neighbors}, new_waiting_nonhydrogens={new_waiting_nonhydrogens}\n"


    def validate(self, verbose:bool = False):
        """ Validate the result of the conversion """

        # Check that, for each input nonhydrogen, there is a matching output nonhydrogen
        for input_node_id in self.input_graph.nodes:
            input_node = self.input_graph.get_node(input_node_id)
            if input_node.data != "H":
                if input_node.node_id in self.input_to_output:
                    output_node = self.input_to_output[input_node.node_id]
                    if output_node.data != input_node.data:
                        raise ValueError(f"Node data does not match: input_node={input_node}, output_node={output_node}")
                else:
                    raise ValueError(f"Missing node: input_node={input_node}")

        # Check that, for each output nonhydrogen, there is a matching input nonhydrogen
        output_graph = self.result.states[-1].graph
        for output_node_id in output_graph.nodes:
            output_node = output_graph.get_node(output_node_id)
            if output_node.data != "H":
                if output_node.node_id in self.output_to_input:
                    input_node = self.output_to_input[output_node.node_id]
                    if input_node.data != output_node.data:
                        raise ValueError(f"Node data does not match: input_node={input_node}, output_node={output_node}")
                else:
                    raise ValueError(f"Missing node: output_node={output_node}")

        if verbose:
            print("The conversion has been validated.")


class IsalChemCompressor:
    """ Given a molecular graph, obtain a compressed IsalChem token list
     as the first one in lexicographic order among those with the minimum
      length. """
    def __init__(self, input_graph: Graph):
        self.input_graph = input_graph  # The input Graph object to be processed
        self.compressed_token_list = None # The first token list in alphabetical order, among the shortest ones
        self.token_lists = []
        self.token_list_lengths = []

    def __repr__(self) -> str:
        input_repr = self.input_graph.__repr__()
        token_lists_repr = self.token_lists.__repr__()
        token_list_lengths_repr = self.token_list_lengths.__repr__()
        compressed_token_list_repr = self.compressed_token_list.__repr__()
        return f"IsalChemCompressor(input_graph={input_repr}, token_lists_repr={token_lists_repr}, token_list_lengths_repr={token_list_lengths_repr}, compressed_token_list_repr={compressed_token_list_repr})"



    def compress(self):
        """ Obtain a compressed IsalChem token list, which is the first one in
         lexicographical order among those with minimum length """

        # For each nonhydrogen atom in the input graph, run a conversion
        # to list of tokens taking that atom as the starting nonhydrogen
        for start_node_id in self.input_graph.nodes:
            node = self.input_graph.get_node(start_node_id)
            # Only run the conversion if the node is a nonhydrogen
            if node.data != "H":
                converter = IsalChemConverter(self.input_graph, start_node_id)
                converter.to_token_list()
                this_token_list = converter.result.token_list
                # Add the obtained token list to the list of token lists
                self.token_lists.append(this_token_list)
                self.token_list_lengths.append(len(this_token_list))

        # Find the number of tokens of the shortest token lists that have been obtained
        minimum_token_list_length = min(self.token_list_lengths)
        # Obtain the shortest token lists
        short_token_lists = []
        for ndx_token_list in range(len(self.token_lists)):
            if len(self.token_lists[ndx_token_list]) == minimum_token_list_length:
                short_token_lists.append(self.token_lists[ndx_token_list])
        # Set the first shortest token list as the normalized token list
        self.compressed_token_list = min(short_token_lists)

class IsalChemComparison:
    """ Compare two IsalChem graphs for equivalence """

    def __init__(self, input_graph1: Graph, input_graph2: Graph):
        self.input_graph1 = input_graph1
        self.input_graph2 = input_graph2

    def __repr__(self) -> str:
        input_graph1_repr = self.input_graph1.__repr__()
        input_graph2_repr = self.input_graph2.__repr__()
        return f"IsalChemComparison(input_graph1={input_graph1_repr}, input_graph2={input_graph2_repr})"

    def graph_features(self, input_graph: Graph) -> (dict, dict):
        """ Compute the features of an IsalChem graph for equivalence purposes """

        node_count_graph = dict()
        edge_count_graph = dict()
        for node_id in input_graph.nodes:
            # Get the current Node object
            node = input_graph.get_node(node_id)
            # Update the node count dictionary
            node_key = node.data , len(node.neighbors)
            if node_key in node_count_graph:
                node_count_graph[node_key] += 1
            else:
                node_count_graph[node_key] = 1
            # Process all the neighbors of this node
            for ndx_neighbor, neighbor in enumerate(node.neighbors):
                edge_key = node.data, neighbor.data , node.edge_types[ndx_neighbor]
                if edge_key in edge_count_graph:
                    edge_count_graph[edge_key] += 1
                else:
                    edge_count_graph[edge_key] = 1

        return node_count_graph, edge_count_graph

    def compare(self, verbose:bool = False) -> bool:
        """ Compare two IsalChem graphs for equivalence, False indicates that they
         are not equivalent """

        node_count_graph1, edge_count_graph1 = self.graph_features(self.input_graph1)
        node_count_graph2, edge_count_graph2 = self.graph_features(self.input_graph2)

        for key, value in node_count_graph1.items():
            if key not in node_count_graph2:
                if verbose:
                    print(f"{key}, {value} not found in second graph.")
                    return False
            else:
                if node_count_graph1[key] != node_count_graph2[key]:
                    if verbose:
                        print(f"{key}, {node_count_graph1[key]} does not match {key}, {node_count_graph2[key]}")
                        return False

        for key, value in node_count_graph2.items():
            if key not in node_count_graph1:
                if verbose:
                    print(f"{key}, {value} not found in first graph.")
                    return False
            else:
                if node_count_graph1[key] != node_count_graph2[key]:
                    if verbose:
                        print(f"{key}, {node_count_graph1[key]} does not match {key}, {node_count_graph2[key]}")
                        return False

        for key, value in edge_count_graph1.items():
            if key not in edge_count_graph2:
                if verbose:
                    print(f"{key}, {value} not found in second graph.")
                    return False
            else:
                if edge_count_graph1[key] != edge_count_graph2[key]:
                    if verbose:
                        print(f"{key}, {edge_count_graph1[key]} does not match {key}, {edge_count_graph2[key]}")
                        return False

        for key, value in edge_count_graph2.items():
            if key not in edge_count_graph1:
                if verbose:
                    print(f"{key}, {value} not found in first graph.")
                    return False
            else:
                if edge_count_graph1[key] != edge_count_graph2[key]:
                    if verbose:
                        print(f"{key}, {edge_count_graph1[key]} does not match {key}, {edge_count_graph2[key]}")
                        return False

        return True


def levenshtein_distance_and_sequence(token_list1: [str], token_list2: [str]) -> (int, [[str]]):
    """
    Compute the Levenshtein distance between two token lists, and the sequence of token lists
    that leads from one token list to the other token list.

    :param token_list1: First token list
    :param token_list2: Second token list
    :return: tuple with the Levenshtein distance and the sequence to get from one token list to another
    """
    m, n = len(token_list1), len(token_list2)

    # Create dynamic programming table using NumPy
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1):
        dp[i, 0] = i
    for j in range(n + 1):
        dp[0, j] = j

    # Fill dynamic programming table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if token_list1[i - 1] == token_list2[j - 1]:
                dp[i, j] = dp[i - 1, j - 1]
            else:
                dp[i, j] = 1 + min(dp[i - 1, j], dp[i, j - 1], dp[i - 1, j - 1])

    # Reconstruct the edit sequence
    i, j = m, n
    sequence = [list(token_list1)]
    current_list = list(token_list1)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and token_list1[i - 1] == token_list2[j - 1]:
            i -= 1
            j -= 1
        elif i > 0 and dp[i, j] == dp[i - 1, j] + 1:  # Deletion
            current_list.pop(i - 1)
            sequence.append(list(current_list))
            i -= 1
        elif j > 0 and dp[i, j] == dp[i, j - 1] + 1:  # Insertion
            current_list.insert(i, token_list2[j - 1])
            sequence.append(list(current_list))
            j -= 1
        else:  # Substitution
            current_list[i - 1] = token_list2[j - 1]
            sequence.append(list(current_list))
            i -= 1
            j -= 1

    return dp[m, n], sequence


def fetch_neighbors(token_list: [str]) -> [[str]]:
    """
    Given an IsalChem string, compute all its immediate neighbors, i.e.
    the IsalChem strings that are at unit Levenshtein distance.

    :param token_list: Token list
    :return: list with the immediate neighbors of the token_list
    """
    possible_tokens = ["+", "-", ">", "<", "C", "N", "B", "O", "Cl", "Br", "F", "I", "C2", "N2", "B2", "O2", "C3", "N3", "B3", "J", "J2"]

    num_tokens = len(token_list)
    edit_types = ["insert", "replace", "delete"]

    neighbors_list = []
    for edit_type in edit_types:
        if edit_type == "insert":
            for ndx_insert in range(num_tokens+1):
                for token in possible_tokens:
                    new_neighbor = token_list.copy()
                    new_neighbor.insert(ndx_insert, token)
                    neighbors_list.append(new_neighbor)
        elif edit_type == "replace":
            for ndx_replace in range(num_tokens):
                for token in possible_tokens:
                    if token != token_list[ndx_replace]:
                        new_neighbor = token_list.copy()
                        new_neighbor[ndx_replace] = token
                        neighbors_list.append(new_neighbor)
        elif edit_type == "delete":
            for ndx_delete in range(num_tokens):
                new_neighbor = token_list.copy()
                del new_neighbor[ndx_delete]
                neighbors_list.append(new_neighbor)

    return neighbors_list


def random_edit(token_list: [str], num_edits: int) -> (int, [[str]]):
    """
    Given an IsalChem token list and a desired Levenshtein distance, draw
    a random IsalChem token list that is roughly at the specified Levenshtein distance. Please note
    that maybe a smaller Levenshtein distance is produced in the cases that
    a random edit is cancelled by another random edit.
    Please note that the output token list is almost always longer than the input one. This
    is because each token list has a majority of longer immediate neighbors over shorter ones.

    :param token_list: Token list
    :param num_edits: Number of edits (desired Levenshtein distance from the input token list
    to the final token list)
    :return: tuple with the actual Levenshtein distance and the sequence to get from the input token list
    to the final token list
    """

    current_token_list = token_list
    for ndx_edit in range(num_edits):
        possible_neighbors = fetch_neighbors(current_token_list)
        ndx_neighbor = np.random.randint(len(possible_neighbors)-1)
        current_token_list = possible_neighbors[ndx_neighbor]

    lev_distance, sequence = levenshtein_distance_and_sequence(token_list, current_token_list)

    return lev_distance, sequence


def fetch_restricted_neighbors(token_list: [str], possible_tokens:[str]=None, abundances:bool=True) -> ([[str]], [[float]]):
    """
    Given an IsalChem string, compute some of its immediate neighbors, i.e.
    some IsalChem strings that are at unit Levenshtein distance.

    :param token_list: Token list
    :param possible_tokens: Possible tokens for insert or replace edits
    :param abundances: If true, the abundances of tokens in the token list are employed to compute the generation probabilities
    :return: list with the immediate neighbors of the token_list, and list of generation probabilities for such neighbors
    """

    # Consider all tokens in case that the input is set to None
    if possible_tokens is None:
        possible_tokens = ["+", "-", ">", "<", "C", "N", "B", "O", "Cl", "Br", "F", "I", "C2", "N2", "B2", "O2", "C3", "N3", "B3", "J", "J2"]
    # possible_tokens = ["C"]

    # If token abundances must be considered, then compute the token counts in the input token list
    if abundances:
        token_counter = Counter(token_list)

    num_tokens = len(token_list)
    edit_types = ["insert", "replace", "delete"]

    neighbors_list = []
    probabilities_list = []
    for edit_type in edit_types:
        if edit_type == "insert":
            for ndx_insert in range(num_tokens+1):
                for token in possible_tokens:
                    new_neighbor = token_list.copy()
                    new_neighbor.insert(ndx_insert, token)
                    neighbors_list.append(new_neighbor)
                    if abundances:
                        probabilities_list.append(token_counter[token])
                    else:
                        probabilities_list.append(1.0)
        elif edit_type == "replace":
            for ndx_replace in range(num_tokens):
                for token in possible_tokens:
                    if token != token_list[ndx_replace]:
                        new_neighbor = token_list.copy()
                        new_neighbor[ndx_replace] = token
                        neighbors_list.append(new_neighbor)
                        if abundances:
                            probabilities_list.append(token_counter[token])
                        else:
                            probabilities_list.append(1.0)
        elif edit_type == "delete":
            for ndx_delete in range(num_tokens):
                new_neighbor = token_list.copy()
                del new_neighbor[ndx_delete]
                neighbors_list.append(new_neighbor)
                probabilities_list.append(1.0)

    # Normalize probabilities
    sum_counts = np.sum(np.array(probabilities_list))
    for ndx_token in range(len(probabilities_list)):
        probabilities_list[ndx_token]/=sum_counts

    # print(probabilities_list)

    return neighbors_list, probabilities_list

def restricted_random_edit(token_list: [str], num_edits: int, possible_tokens:[str]=None, abundances:bool=True) -> (int, [[str]]):
    """
    Given an IsalChem token list and a desired Levenshtein distance, draw
    a random IsalChem token list that is roughly at the specified Levenshtein distance.
    Only a reduced set of instructions are employed for the insertion and replacement edits.
    Please note that maybe a smaller Levenshtein distance is produced in the cases that
    a random edit is cancelled by another random edit.
    Please note that the output token list is almost always longer than the input one. This
    is because each token list has a majority of longer immediate neighbors over shorter ones.

    :param token_list: Token list
    :param num_edits: Number of edits (desired Levenshtein distance from the input token list
    to the final token list)
    :param possible_tokens: Possible tokens for insert or replace edits
    :param abundances: If true, the abundances of tokens in the token list are employed to compute the generation probabilities
    :return: tuple with the actual Levenshtein distance and the sequence to get from the input token list
    to the final token list
    """

    current_token_list = token_list
    for ndx_edit in range(num_edits):
        possible_neighbors, probabilities_neighbors = fetch_restricted_neighbors(current_token_list, possible_tokens, abundances)
        ndx_neighbor = np.random.choice(range(len(possible_neighbors)), p=probabilities_neighbors)
        current_token_list = possible_neighbors[ndx_neighbor]

    lev_distance, sequence = levenshtein_distance_and_sequence(token_list, current_token_list)

    return lev_distance, sequence

def int_bond_type_to_rdkit(bond_type: int):
    """
    Converts an integer bond type (1, 2, 3) into an RDKit BondType.
    """
    bond_map = {
        1: BondType.SINGLE,
        2: BondType.DOUBLE,
        3: BondType.TRIPLE,
    }
    return bond_map.get(bond_type, BondType.UNSPECIFIED)


def graph_to_mol(graph: Graph):
    """
    Converts an instance of the Graph class into an RDKit Mol instance.
    :param: graph is a Graph instance with the molecular graph to be converted
    :return: the computed RDKit RWMol instance
    """
    mol = RWMol()
    node_to_idx = {}  # Map node_id to RDKit atom indices

    # Add atoms to the molecule
    for node_id, node in graph.nodes.items():
        atom = Atom(node.data)  # Assume `node.data` holds the atomic symbol (e.g., 'C', 'O', 'N')
        idx = mol.AddAtom(atom)
        node_to_idx[node_id] = idx

    # Add bonds to the molecule
    for node_id, node in graph.nodes.items():
        for neighbor, bond_type in zip(node.neighbors, node.edge_types):
            # Add bonds only if the bond doesn't already exist
            if mol.GetBondBetweenAtoms(node_to_idx[node_id], node_to_idx[neighbor.node_id]) is None:
                mol.AddBond(
                    node_to_idx[node_id],
                    node_to_idx[neighbor.node_id],
                    int_bond_type_to_rdkit(bond_type)
                )

    return mol


def mol_to_graph(mol) -> Graph:
    """ Converts an RDKit Mol object into a Graph instance
    :param mol: RDKit Mol object to be converted
    :return: Graph instance
    """
    bond_mapping = {Chem.rdchem.BondType.SINGLE: 1,
                    Chem.rdchem.BondType.DOUBLE: 2,
                    Chem.rdchem.BondType.TRIPLE: 3}

    atom_symbols = {"H", "C", "N", "B", "O", "Cl", "Br", "F", "I"}

    graph = Graph()
    for atom in mol.GetAtoms():
        this_atom_symbol = atom.GetSymbol()
        # Check whether a non supported atom symbol has been found
        if this_atom_symbol not in atom_symbols:
          print(this_atom_symbol)
          return None
        graph.add_node(atom.GetIdx(), this_atom_symbol)
    for bond in mol.GetBonds():
        bond_type = bond_mapping.get(bond.GetBondType(), 1)  # Default to single bond
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type)
    return graph

def smiles_tokenizer(smiles_string:str) -> [str]:
    """
    Tokenize a SMILES molecule or reaction string
    :param smiles_string: SMILES string to tokenize
    :return: list of tokens
    """
    pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smiles_string)]
    #print(smiles_string)
    #print(''.join(tokens))
    # assert smiles_string == ''.join(tokens)
    return tokens

