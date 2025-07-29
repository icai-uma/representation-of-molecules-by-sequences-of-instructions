"""

Class to model a molecule under the IsalChem methodology:
Lopez-Rubio, Ezequiel (2024). Instruction set and language for chemical nomenclature
https://doi.org/10.26434/chemrxiv-2024-5b4dn

Coded by Ezequiel Lopez-Rubio, November 2024.

"""

from isalchemstate import IsalChemState
import copy


class IsalChemMolecule:
    """ Class to model a molecule under the IsalChem methodology """
    def __init__(self, token_list):

        self.token_list = token_list  # The list of tokens (instructions) to be executed
        self.current_token_ndx = 0  # Index of the current token in the token list
        initial_state = IsalChemState()  # The initial state of the molecule
        self.states = [initial_state] # The sequence of states that the molecule has undergone
        self.atom_counter = 3  # A counter to create new node IDs for the new atoms
        self.last_nonhydrogen = None  # The Node object of the last added non-hydrogen atom

    def __repr__(self) -> str:
        current_state_repr = self.states[-1].__repr__()
        return f"IsalChemMolecule(token_list={self.token_list}, current_state={current_state_repr})"

    def step_C(self, state: IsalChemState):
        """ Execute a C instruction: insert a carbon atom with a single bond """

        # Add the three new hydrogen atoms
        new_hydrogen_1 = state.graph.add_node(self.atom_counter, "H")
        self.atom_counter += 1
        new_hydrogen_2 = state.graph.add_node(self.atom_counter, "H")
        self.atom_counter += 1
        new_hydrogen_3 = state.graph.add_node(self.atom_counter, "H")
        self.atom_counter += 1

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Insert the three new hydrogen atoms into the list, after the old hydrogen atom
        state.dll.insert_after(old_hydrogen, new_hydrogen_1)
        state.dll.insert_after(new_hydrogen_1, new_hydrogen_2)
        state.dll.insert_after(new_hydrogen_2, new_hydrogen_3)

        # Update the primary pointer to point to the second new hydrogen atom
        state.primary_ptr = new_hydrogen_2.node_id

        # Update the secondary pointer if it was pointing to the old hydrogen atom
        if state.secondary_ptr == old_hydrogen.node_id:
            state.secondary_ptr = new_hydrogen_2.node_id

        # Remove the old hydrogen atom from the list (but not from the graph)
        state.dll.remove(old_hydrogen)

        # Reuse the old hydrogen atom as the new carbon atom
        # The old edge is now bond 0 of the new carbon atom
        old_hydrogen.data = "C"

        # The Node object of the last added non-hydrogen atom
        self.last_nonhydrogen = old_hydrogen

        # Add edges from the new carbon atom to the new hydrogen atoms into the graph
        # These are bond numbers 1, 2 and 3 of the new carbon atom
        state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_1.node_id)
        state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_2.node_id)
        state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_3.node_id)

    def step_NB(self, state: IsalChemState, current_token: str):
        """ Execute an N/B instruction: insert a nitrogen/boron atom with a single bond """

        # Add the two new hydrogen atoms
        new_hydrogen_1 = state.graph.add_node(self.atom_counter, "H")
        self.atom_counter += 1
        new_hydrogen_2 = state.graph.add_node(self.atom_counter, "H")
        self.atom_counter += 1

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Insert the two new hydrogen atoms into the list, after the old hydrogen atom
        state.dll.insert_after(old_hydrogen, new_hydrogen_1)
        state.dll.insert_after(new_hydrogen_1, new_hydrogen_2)

        # Update the primary pointer to point to the second new hydrogen atom
        state.primary_ptr = new_hydrogen_2.node_id

        # Update the secondary pointer if it was pointing to the old hydrogen atom
        if state.secondary_ptr == old_hydrogen.node_id:
            state.secondary_ptr = new_hydrogen_2.node_id

        # Remove the old hydrogen atom from the list (but not from the graph)
        state.dll.remove(old_hydrogen)

        # Reuse the old hydrogen atom as the new nitrogen/boron atom
        old_hydrogen.data = current_token

        # The Node object of the last added non-hydrogen atom
        self.last_nonhydrogen = old_hydrogen

        # Add edges from the new nitrogen/boron atom to the new hydrogen atoms into the graph
        state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_1.node_id)
        state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_2.node_id)

    def step_O(self, state: IsalChemState):
        """ Execute an O instruction: insert an oxygen atom with a single bond """

        # Add the new hydrogen atom
        new_hydrogen_1 = state.graph.add_node(self.atom_counter, "H")
        self.atom_counter += 1

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Insert the new hydrogen atom into the list, after the old hydrogen atom
        state.dll.insert_after(old_hydrogen, new_hydrogen_1)

        # Update the primary pointer to point to the new hydrogen atom
        state.primary_ptr = new_hydrogen_1.node_id

        # Update the secondary pointer if it was pointing to the old hydrogen atom
        if state.secondary_ptr == old_hydrogen.node_id:
            state.secondary_ptr = new_hydrogen_1.node_id

        # Remove the old hydrogen atom from the list (but not from the graph)
        state.dll.remove(old_hydrogen)

        # Reuse the old hydrogen atom as the new oxygen atom
        old_hydrogen.data = "O"

        # The Node object of the last added non-hydrogen atom
        self.last_nonhydrogen = old_hydrogen

        # Add edges from the new oxygen atom to the new hydrogen atom into the graph
        state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_1.node_id)

    def step_halogen(self, state: IsalChemState, current_token: str):
        """ Execute a halogen instruction: insert a halogen atom """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Update the primary pointer to point to the previous hydrogen atom in the list
        state.primary_ptr = old_hydrogen.prev.node_id

        # Update the secondary pointer if it was pointing to the old hydrogen atom
        if state.secondary_ptr == old_hydrogen.node_id:
            state.secondary_ptr = state.primary_ptr

        # Remove the old hydrogen atom from the list (but not from the graph)
        state.dll.remove(old_hydrogen)

        # Reuse the old hydrogen atom as the new halogen atom
        old_hydrogen.data = current_token

        # The Node object of the last added non-hydrogen atom
        self.last_nonhydrogen = old_hydrogen

    def step_C2(self, state: IsalChemState):
        """ Execute a C2 instruction: insert a carbon atom with a double bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the next old hydrogen atom in the list
        other_old_hydrogen = old_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_non_hydrogen = old_hydrogen.neighbors[0]

        # Check whether inserting a double bond carbon atom is feasible
        #if len(list_other_old_hydrogens) == 0:
        if (other_old_hydrogen is None or other_old_hydrogen==old_hydrogen or
                other_old_hydrogen not in old_non_hydrogen.neighbors):
            # If the insertion of the double bond is not feasible, then fall back to instruction C
            self.step_C(state)
        else:

            # Add the two new hydrogen atoms
            new_hydrogen_1 = state.graph.add_node(self.atom_counter, "H")
            self.atom_counter += 1
            new_hydrogen_2 = state.graph.add_node(self.atom_counter, "H")
            self.atom_counter += 1

            # Insert the two new hydrogen atoms into the list, after the old hydrogen atom
            state.dll.insert_after(old_hydrogen, new_hydrogen_1)
            state.dll.insert_after(new_hydrogen_1, new_hydrogen_2)

            # Update the primary pointer to point to the second new hydrogen atom
            state.primary_ptr = new_hydrogen_2.node_id

            # Update the secondary pointer if it was pointing to the old hydrogen atom
            # or the other old hydrogen atom
            old_hydrogens = {old_hydrogen.node_id, \
                             other_old_hydrogen.node_id}
            if state.secondary_ptr in old_hydrogens:
                state.secondary_ptr = new_hydrogen_2.node_id

            # Remove the old hydrogen atom from the list (but not from the graph)
            state.dll.remove(old_hydrogen)

            # Reuse the old hydrogen atom as the new carbon atom
            # The old edge is now bonds 0 and 1 of the new carbon atom
            old_hydrogen.data = "C"

            # The Node object of the last added non-hydrogen atom
            self.last_nonhydrogen = old_hydrogen

            # Change the edge type to double bond
            ndx_old_hydrogen = old_non_hydrogen.neighbors.index(old_hydrogen)
            old_non_hydrogen.edge_types[ndx_old_hydrogen] = 2
            old_hydrogen.edge_types[0] = 2

            # Remove the other old hydrogen atom from the list and the graph
            state.dll.remove(other_old_hydrogen)
            state.graph.remove_node(other_old_hydrogen.node_id)

            # Add edges from the new carbon atom to the new hydrogen atoms into the graph
            # These are bonds 2 and 3 of the new carbon atom
            state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_1.node_id)
            state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_2.node_id)

    def step_N2B2(self, state: IsalChemState, current_token: str):
        """ Execute a N2/B2 instruction: insert a nitrogen/boron atom with a double bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the next old hydrogen atom in the list
        other_old_hydrogen = old_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_non_hydrogen = old_hydrogen.neighbors[0]

        # Check whether inserting a double bond nitrogen/boron atom is feasible
        #if len(list_other_old_hydrogens) == 0:
        if (other_old_hydrogen is None or other_old_hydrogen==old_hydrogen
                or other_old_hydrogen not in old_non_hydrogen.neighbors):
            # If there are no other old hydrogen atoms, then fall back to instruction N/B
            self.step_NB(state, current_token[0])
        else:

            # Add the new hydrogen atom
            new_hydrogen_1 = state.graph.add_node(self.atom_counter, "H")
            self.atom_counter += 1

            # Insert the new hydrogen atom into the list, after the old hydrogen atom
            state.dll.insert_after(old_hydrogen, new_hydrogen_1)

            # Update the primary pointer to point to the new hydrogen atom
            state.primary_ptr = new_hydrogen_1.node_id

            # Update the secondary pointer if it was pointing to the old hydrogen atom
            # or the other old hydrogen atom
            old_hydrogens = {old_hydrogen.node_id, \
                             other_old_hydrogen.node_id}
            if state.secondary_ptr in old_hydrogens:
                state.secondary_ptr = new_hydrogen_1.node_id

            # Remove the old hydrogen atom from the list (but not from the graph)
            state.dll.remove(old_hydrogen)

            # Reuse the old hydrogen atom as the new nitrogen/boron atom
            old_hydrogen.data = current_token[0]

            # The Node object of the last added non-hydrogen atom
            self.last_nonhydrogen = old_hydrogen

            # Change the edge type to double bond
            ndx_old_hydrogen = old_non_hydrogen.neighbors.index(old_hydrogen)
            old_non_hydrogen.edge_types[ndx_old_hydrogen] = 2
            old_hydrogen.edge_types[0] = 2

            # Remove the other old hydrogen atom from the list and the graph
            state.dll.remove(other_old_hydrogen)
            state.graph.remove_node(other_old_hydrogen.node_id)

            # Add edge from the new nitrogen/boron atom to the new hydrogen atom into the graph
            state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_1.node_id)

    def step_O2(self, state: IsalChemState):
        """ Execute an O2 instruction: insert an oxygen atom with a double bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the next old hydrogen atom in the list
        other_old_hydrogen = old_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_non_hydrogen = old_hydrogen.neighbors[0]

        # Check whether inserting a double bond oxygen atom is feasible
        #if len(list_other_old_hydrogens) == 0:
        if (other_old_hydrogen is None or other_old_hydrogen==old_hydrogen
                or other_old_hydrogen not in old_non_hydrogen.neighbors):
            # If there are no other old hydrogen atoms, then fall back to instruction O
            self.step_O(state)
        else:

            # Update the primary pointer to point to the previous hydrogen atom
            old_hydrogens = {old_hydrogen.node_id, \
                             other_old_hydrogen.node_id}
            previous_hydrogens = {old_hydrogen.prev.node_id, \
                                  other_old_hydrogen.prev.node_id}

            remaining_hydrogens = list(previous_hydrogens - old_hydrogens)
            if len(remaining_hydrogens) > 0:
                remaining_hydrogen = remaining_hydrogens[0]
                state.primary_ptr = remaining_hydrogen

                # Update the secondary pointer if it was pointing to the old hydrogen atom
                # or the other old hydrogen atom
                if state.secondary_ptr in old_hydrogens:
                    state.secondary_ptr = remaining_hydrogen

            # Remove the old hydrogen atom from the list (but not from the graph)
            state.dll.remove(old_hydrogen)

            # Reuse the old hydrogen atom as the new oxygen atom
            old_hydrogen.data = "O"

            # The Node object of the last added non-hydrogen atom
            self.last_nonhydrogen = old_hydrogen

            # Change the edge type to double bond
            ndx_old_hydrogen = old_non_hydrogen.neighbors.index(old_hydrogen)
            old_non_hydrogen.edge_types[ndx_old_hydrogen] = 2
            old_hydrogen.edge_types[0] = 2

            # Remove the other old hydrogen atom from the list and the graph
            state.dll.remove(other_old_hydrogen)
            state.graph.remove_node(other_old_hydrogen.node_id)

    def step_C3(self, state: IsalChemState):
        """ Execute a C3 instruction: insert a carbon atom with a triple bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the previous and next old hydrogen atoms in the list
        first_other_old_hydrogen = old_hydrogen.prev
        second_other_old_hydrogen = old_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_non_hydrogen = old_hydrogen.neighbors[0]

        # Check whether inserting a triple bond carbon atom is feasible
        if (first_other_old_hydrogen is None or second_other_old_hydrogen is None or
                first_other_old_hydrogen==old_hydrogen or second_other_old_hydrogen==old_hydrogen
                or first_other_old_hydrogen==second_other_old_hydrogen or
                first_other_old_hydrogen not in old_non_hydrogen.neighbors or
                second_other_old_hydrogen not in old_non_hydrogen.neighbors):
            # If there are no other old hydrogen atoms, then fall back to instruction C2
            self.step_C2(state)
        else:

            # Add the new hydrogen atom
            new_hydrogen_1 = state.graph.add_node(self.atom_counter, "H")
            self.atom_counter += 1

            # Insert the new hydrogen atom into the list, after the old hydrogen atom
            state.dll.insert_after(old_hydrogen, new_hydrogen_1)

            # Update the primary pointer to point to the new hydrogen atom
            state.primary_ptr = new_hydrogen_1.node_id

            # Update the secondary pointer if it was pointing to the old hydrogen atom
            # or the two other old hydrogen atoms
            old_hydrogens = {old_hydrogen.node_id, \
                             first_other_old_hydrogen.node_id, \
                             second_other_old_hydrogen.node_id}
            if state.secondary_ptr in old_hydrogens:
                state.secondary_ptr = new_hydrogen_1.node_id

            # Remove the old hydrogen atom from the list (but not from the graph)
            state.dll.remove(old_hydrogen)

            # Reuse the old hydrogen atom as the new carbon atom
            # The old edge is now bonds 0, 1 and 2 of the new carbon atom
            old_hydrogen.data = "C"

            # The Node object of the last added non-hydrogen atom
            self.last_nonhydrogen = old_hydrogen

            # Change the edge type to triple bond
            ndx_old_hydrogen = old_non_hydrogen.neighbors.index(old_hydrogen)
            old_non_hydrogen.edge_types[ndx_old_hydrogen] = 3
            old_hydrogen.edge_types[0] = 3

            # Remove the other old hydrogen atoms from the list and the graph
            state.dll.remove(first_other_old_hydrogen)
            state.graph.remove_node(first_other_old_hydrogen.node_id)
            state.dll.remove(second_other_old_hydrogen)
            state.graph.remove_node(second_other_old_hydrogen.node_id)

            # Add an edge from the new carbon atom to the new hydrogen atom into the graph
            # This is bond 3 of the new carbon atom
            state.graph.add_edge(old_hydrogen.node_id, new_hydrogen_1.node_id)


    def step_N3B3(self, state: IsalChemState, current_token: str):
        """ Execute an N3/B3 instruction: insert a nitrogen/boron atom with a triple bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the previous and next old hydrogen atoms in the list
        first_other_old_hydrogen = old_hydrogen.prev
        second_other_old_hydrogen = old_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_non_hydrogen = old_hydrogen.neighbors[0]

        # Check whether inserting a triple bond nitrogen/boron atom is feasible
        if (first_other_old_hydrogen is None or second_other_old_hydrogen is None or
                first_other_old_hydrogen == old_hydrogen or second_other_old_hydrogen == old_hydrogen
                or first_other_old_hydrogen == second_other_old_hydrogen or
                first_other_old_hydrogen not in old_non_hydrogen.neighbors or
                second_other_old_hydrogen not in old_non_hydrogen.neighbors):
            # If there are no other old hydrogen atoms, then fall back to instruction N2/B2
            self.step_N2B2(state, current_token)
        else:

            # Update the primary pointer to point to the previous hydrogen atom
            old_hydrogens = {old_hydrogen.node_id, \
                             first_other_old_hydrogen.node_id, \
                             second_other_old_hydrogen.node_id}
            previous_hydrogens = {old_hydrogen.prev.node_id, \
                                  first_other_old_hydrogen.prev.node_id, \
                                  second_other_old_hydrogen.prev.node_id}

            remaining_hydrogens = list(previous_hydrogens - old_hydrogens)
            if len(remaining_hydrogens) > 0:
                remaining_hydrogen = remaining_hydrogens[0]
                state.primary_ptr = remaining_hydrogen

                # Update the secondary pointer if it was pointing to the old hydrogen atom
                # or the two other old hydrogen atoms
                if state.secondary_ptr in old_hydrogens:
                    state.secondary_ptr = remaining_hydrogen

            # Remove the old hydrogen atom from the list (but not from the graph)
            state.dll.remove(old_hydrogen)

            # Reuse the old hydrogen atom as the new carbon atom
            old_hydrogen.data = current_token[0]

            # The Node object of the last added non-hydrogen atom
            self.last_nonhydrogen = old_hydrogen

            # Change the edge type to triple bond
            ndx_old_hydrogen = old_non_hydrogen.neighbors.index(old_hydrogen)
            old_non_hydrogen.edge_types[ndx_old_hydrogen] = 3
            old_hydrogen.edge_types[0] = 3

            # Remove the other old hydrogen atoms from the list and the graph
            state.dll.remove(first_other_old_hydrogen)
            state.graph.remove_node(first_other_old_hydrogen.node_id)
            state.dll.remove(second_other_old_hydrogen)
            state.graph.remove_node(second_other_old_hydrogen.node_id)

    def step_J(self, state: IsalChemState):
        """ Execute a J instruction: join two atoms with a simple bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_primary_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_primary_non_hydrogen = old_primary_hydrogen.neighbors[0]

        # Check for the case that there are just two hydrogen atoms in the molecule
        if old_primary_non_hydrogen.data == "H":
            return

        # Fetch the old hydrogen atom that is pointed by the secondary pointer
        old_secondary_hydrogen = state.graph.get_node(state.secondary_ptr)

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the secondary pointer
        old_secondary_non_hydrogen = old_secondary_hydrogen.neighbors[0]

        # If both pointers are associated with the same non hydrogen atom,
        # then skip the instruction
        if old_primary_non_hydrogen == old_secondary_non_hydrogen:
            return

        # If both nonhydrogens are already bound, then skip the instruction
        if old_primary_non_hydrogen in old_secondary_non_hydrogen.neighbors or \
            old_secondary_non_hydrogen in old_primary_non_hydrogen.neighbors:
            return

        # Join the non hydrogen atoms
        ndx_old_primary_hydrogen = old_primary_non_hydrogen.neighbors.index(old_primary_hydrogen)
        ndx_old_secondary_hydrogen = old_secondary_non_hydrogen.neighbors.index(old_secondary_hydrogen)
        old_primary_non_hydrogen.neighbors[ndx_old_primary_hydrogen] = old_secondary_non_hydrogen
        old_primary_non_hydrogen.edge_types[ndx_old_primary_hydrogen] = 1
        old_secondary_non_hydrogen.neighbors[ndx_old_secondary_hydrogen] = old_primary_non_hydrogen
        old_secondary_non_hydrogen.edge_types[ndx_old_secondary_hydrogen] = 1

        # Update the pointers to their respective previous hydrogen atoms in the list
        if old_primary_hydrogen.prev.node_id != old_secondary_hydrogen.node_id:
            state.primary_ptr = old_primary_hydrogen.prev.node_id
        else:
            state.primary_ptr = old_primary_hydrogen.prev.prev.node_id
        if old_secondary_hydrogen.prev.node_id != old_primary_hydrogen.node_id:
            state.secondary_ptr = old_secondary_hydrogen.prev.node_id
        else:
            state.secondary_ptr = old_secondary_hydrogen.prev.prev.node_id

        # Remove the old hydrogen atoms from the list and the graph
        state.dll.remove(old_primary_hydrogen)
        del state.graph.nodes[old_primary_hydrogen.node_id]
        state.dll.remove(old_secondary_hydrogen)
        del state.graph.nodes[old_secondary_hydrogen.node_id]

    def step_J2(self, state: IsalChemState):
        """ Execute a J2 instruction: join two atoms with a double bond """

        # Fetch the old hydrogen atom that is pointed by the primary pointer
        old_primary_hydrogen = state.graph.get_node(state.primary_ptr)

        # Fetch the next hydrogen atom in the list
        other_old_primary_hydrogen = old_primary_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the primary pointer
        old_primary_non_hydrogen = old_primary_hydrogen.neighbors[0]

        # Check for the case that there are just two hydrogen atoms in the molecule
        if old_primary_non_hydrogen.data == "H":
            return

        # Fetch the old hydrogen atom that is pointed by the secondary pointer
        old_secondary_hydrogen = state.graph.get_node(state.secondary_ptr)

        # Fetch the next hydrogen atom in the list
        other_old_secondary_hydrogen = old_secondary_hydrogen.next

        # Fetch the old non hydrogen atom that is bound to the hydrogen atom
        # pointed by the secondary pointer
        old_secondary_non_hydrogen = old_secondary_hydrogen.neighbors[0]

        # If both pointers are associated with the same non hydrogen atom,
        # then skip the instruction
        if old_primary_non_hydrogen == old_secondary_non_hydrogen:
            return

        # If both nonhydrogens are already bound, then skip the instruction
        if old_primary_non_hydrogen in old_secondary_non_hydrogen.neighbors or \
            old_secondary_non_hydrogen in old_primary_non_hydrogen.neighbors:
            return

        # If the hydrogen atom pointed by the primary pointer is bound to an atom
        # with no other hydrogen atoms bound or the secondary pointer is bound to
        # an atom with no other hydrogen atoms bound, then instruction J is executed
        if (other_old_primary_hydrogen is None or other_old_secondary_hydrogen is None or
            other_old_primary_hydrogen not in old_primary_non_hydrogen.neighbors or
            other_old_secondary_hydrogen not in old_secondary_non_hydrogen.neighbors):
            self.step_J(state)
        else:

            # Join the non hydrogen atoms
            ndx_old_primary_hydrogen = old_primary_non_hydrogen.neighbors.index(old_primary_hydrogen)
            ndx_old_secondary_hydrogen = old_secondary_non_hydrogen.neighbors.index(old_secondary_hydrogen)
            old_primary_non_hydrogen.neighbors[ndx_old_primary_hydrogen] = old_secondary_non_hydrogen
            old_primary_non_hydrogen.edge_types[ndx_old_primary_hydrogen] = 2
            old_secondary_non_hydrogen.neighbors[ndx_old_secondary_hydrogen] = old_primary_non_hydrogen
            old_secondary_non_hydrogen.edge_types[ndx_old_secondary_hydrogen] = 2

            # The set of Node IDs of the old hydrogen atoms that will be removed from the list and the graph
            old_hydrogens = {old_primary_hydrogen.node_id, \
                                     other_old_primary_hydrogen.node_id, \
                                     old_secondary_hydrogen.node_id, \
                                     other_old_secondary_hydrogen.node_id
                                     }

            # Update the primary pointer to point to the previous hydrogen atom
            remaining_primary_hydrogen = old_primary_hydrogen.prev.node_id
            # Go back in the list until a hydrogen atom is found which is not going to be removed
            while remaining_primary_hydrogen in old_hydrogens and remaining_primary_hydrogen != state.primary_ptr:
                remaining_primary_hydrogen = state.graph.get_node(remaining_primary_hydrogen).prev.node_id
            # Check whether a hydrogen atom was found that is not going to be removed
            if remaining_primary_hydrogen in old_hydrogens:
                # No remaining hydrogen exists in the state
                state.primary_ptr = None
            else:
                state.primary_ptr = remaining_primary_hydrogen

            # Update the secondary pointer to point to the previous hydrogen atom
            remaining_secondary_hydrogen = old_secondary_hydrogen.prev.node_id
            # Go back in the list until a hydrogen atom is found which is not going to be removed
            while remaining_secondary_hydrogen in old_hydrogens and remaining_secondary_hydrogen != state.secondary_ptr:
                remaining_secondary_hydrogen = state.graph.get_node(remaining_secondary_hydrogen).prev.node_id
            # Check whether a hydrogen atom was found that is not going to be removed
            if remaining_secondary_hydrogen in old_hydrogens:
                # No remaining hydrogen exists in the state
                state.secondary_ptr = None
            else:
                state.secondary_ptr = remaining_secondary_hydrogen

            # Remove the old hydrogen atoms from the list and the graph
            state.dll.remove(old_primary_hydrogen)
            del state.graph.nodes[old_primary_hydrogen.node_id]
            state.dll.remove(other_old_primary_hydrogen)
            ndx_other_old_primary_hydrogen = old_primary_non_hydrogen.neighbors.index(other_old_primary_hydrogen)
            del old_primary_non_hydrogen.neighbors[ndx_other_old_primary_hydrogen]
            del old_primary_non_hydrogen.edge_types[ndx_other_old_primary_hydrogen]
            del state.graph.nodes[other_old_primary_hydrogen.node_id]
            state.dll.remove(old_secondary_hydrogen)
            del state.graph.nodes[old_secondary_hydrogen.node_id]
            state.dll.remove(other_old_secondary_hydrogen)
            ndx_other_old_secondary_hydrogen = old_secondary_non_hydrogen.neighbors.index(other_old_secondary_hydrogen)
            del old_secondary_non_hydrogen.neighbors[ndx_other_old_secondary_hydrogen]
            del old_secondary_non_hydrogen.edge_types[ndx_other_old_secondary_hydrogen]
            del state.graph.nodes[other_old_secondary_hydrogen.node_id]

    def step(self):
        """ Execute the next instruction in the token list """

        if self.current_token_ndx >= len(self.token_list):
            raise RuntimeError("No more tokens to execute in the token list")

        # The previous state
        previous_state = self.states[-1]
        # The current state is initialized as a deep copy of the previous state
        state = copy.deepcopy(previous_state)

        # If the hydrogen atom list is empty, then no operation is performed
        if state.dll.head is not None:
            # Analyze the current token
            current_token = self.token_list[self.current_token_ndx]
            if current_token == "A":
                pass
            elif current_token == "+":
                state.primary_ptr = state.graph.get_node(state.primary_ptr).next.node_id
            elif current_token == "-":
                state.primary_ptr = state.graph.get_node(state.primary_ptr).prev.node_id
            elif current_token == ">":
                state.secondary_ptr = state.graph.get_node(state.secondary_ptr).next.node_id
            elif current_token == "<":
                state.secondary_ptr = state.graph.get_node(state.secondary_ptr).prev.node_id
            elif current_token == "C":
                self.step_C(state)
            elif current_token == "N" or current_token == "B":
                self.step_NB(state, current_token)
            elif current_token == "O":
                self.step_O(state)
            elif current_token == "Cl" or current_token == "Br" or current_token == "F" or current_token == "I":
                self.step_halogen(state, current_token)
            elif current_token == "C2":
                self.step_C2(state)
            elif current_token == "N2" or current_token == "B2":
                self.step_N2B2(state, current_token)
            elif current_token == "O2":
                self.step_O2(state)
            elif current_token == "C3":
                self.step_C3(state)
            elif current_token == "N3" or current_token == "B3":
                self.step_N3B3(state, current_token)
            elif current_token == "J":
                self.step_J(state)
            elif current_token == "J2":
                self.step_J2(state)
            else:
                raise ValueError(f"Invalid token in the token list: {current_token}")

        primary_ptr_node = state.graph.get_node(state.primary_ptr)
        if primary_ptr_node is None and state.dll.head is not None:
            print(f"Invalid primary ptr: {state}")
            print(f"Current token: {current_token}")
            print(f"Previous state: {self.states[-1]}")
            state.generate_pptx("error_current.pptx")
            self.states[-1].generate_pptx("error_previous.pptx")
            raise RuntimeError(f"Empty list of non-hydrogens for token list: {current_token}")

        # Append the new state to the state list
        self.states.append(state)
        # Update the current token pointer into the token list
        self.current_token_ndx += 1

    def run(self):
        """
        Execute the entire token list.
        """
        for ndx_step in range(len(self.token_list)):
            self.step()
