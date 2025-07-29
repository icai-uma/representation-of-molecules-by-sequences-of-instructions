"""

Basic data structures for the IsalChem methodology:
Lopez-Rubio, Ezequiel (2024). Instruction set and language for chemical nomenclature
https://doi.org/10.26434/chemrxiv-2024-5b4dn

Coded by Ezequiel Lopez-Rubio, November 2024.

"""

class Node:
    """ A node for the IsalChem lists and graphs """
    def __init__(self, node_id: int, data: str = ""):
        self.node_id = node_id
        self.data = data  # A string to identify the element of the atom
        self.neighbors = []  # List of adjacent Node objects in the graph
        self.edge_types = []  # List of bond types for the adjacent Node objects in the graph
        self.prev = None  # Previous Node object in the doubly linked circular list
        self.next = None  # Next Node object in the doubly linked circular list
        self.in_list = False  # Flag to check if the node is currently in the list

    def __repr__(self) -> str:
        return f"Node(node_id={self.node_id}, data='{self.data}')"

    # Nodes must be compared considering their Node IDs only
    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
            getattr(other, 'node_id', None) == self.node_id)

    # Nodes must be hashed considering their Node IDs only
    def __hash__(self):
        return hash(self.node_id)


class Graph:
    """ An IsalChem graph that contains a Node object for each atom of the molecule """
    def __init__(self):
        # Dictionary to store the nodes.
        # Each key is a node ID (an integer)
        # Each value is the associated Node object
        self.nodes = {}

    def add_node(self, node_id: int, data: str = "") -> Node:
        """  Given an unused node ID, create a node with that node ID, and add the node to the graph.
        The new node is returned. """
        if node_id not in self.nodes:
            new_node = Node(node_id, data)
            self.nodes[node_id] = new_node
            return new_node
        else:
            raise RuntimeError(f"Node {node_id} already exists.")
            return None

    def add_edge(self, src_node_id: int, dest_node_id: int, edge_type: int = 1):
        """ Add an edge to the graph """
        # Check if the nodes already exist
        if src_node_id in self.nodes and dest_node_id in self.nodes:
            src_node = self.nodes[src_node_id]
            dest_node = self.nodes[dest_node_id]

            # Check if the edge already exists before inserting the edge
            if dest_node not in src_node.neighbors:
                src_node.neighbors.append(dest_node)
                dest_node.neighbors.append(src_node)
                src_node.edge_types.append(edge_type)
                dest_node.edge_types.append(edge_type)
            else:
                raise RuntimeError(f"Edge between {src_node_id} and {dest_node_id} already exists.")

    def remove_edge(self, src_node_id: int, dest_node_id: int):
        """ Remove an edge from the graph """
        if src_node_id in self.nodes and dest_node_id in self.nodes:
            src_node = self.nodes[src_node_id]
            dest_node = self.nodes[dest_node_id]

            # Remove dest_node from src_node's neighbors, if it exists
            if dest_node in src_node.neighbors:
                ndx_dest = src_node.neighbors.index(dest_node)
                src_node.neighbors.pop(ndx_dest)
                src_node.edge_types.pop(ndx_dest)
                ndx_src = dest_node.neighbors.index(src_node)
                dest_node.neighbors.pop(ndx_src)
                dest_node.edge_types.pop(ndx_src)
            else:
                raise RuntimeError(f"No edge exists between {src_node_id} and {dest_node_id}.")
        else:
            raise RuntimeError(f"One or both nodes {src_node_id}, {dest_node_id} do not exist.")

    def remove_node(self, node_id: int):
        """  Remove a node from the graph. All edges connected to the node are also removed
            from the graph"""
        if node_id in self.nodes:
            # Remove all edges connected to this node
            node_to_remove = self.nodes[node_id]
            for neighbor in node_to_remove.neighbors:
                ndx = neighbor.neighbors.index(node_to_remove)
                neighbor.neighbors.pop(ndx)
                neighbor.edge_types.pop(ndx)

            # Finally, delete the node from the graph
            del self.nodes[node_id]
        else:
            raise RuntimeError(f"Node {node_id} does not exist in the graph.")

    # Get a Node object from the node ID
    def get_node(self, node_id: int) -> Node:
        return self.nodes.get(node_id)

    def __repr__(self) -> str:
        edges = []
        for node_id, node in self.nodes.items():
            for neighbor in node.neighbors:
                edges.append((node_id, neighbor.node_id))

        return f"Graph({list(self.nodes.keys())}, {edges})"


class DoublyLinkedCircularList:
    """ An IsalChem doubly linked circular list of hydrogen atoms represented by Node objects """
    def __init__(self):
        self.head = None

    def insert(self, node: Node):
        """ Insert a new node into the list, before the head """
        if node.in_list:
            raise RuntimeError(f"Node {node.node_id} is already in the list.")
            return

        node.in_list = True
        if self.head is None:
            # Initialize with the first node in the list
            node.next = node
            node.prev = node
            self.head = node
        else:
            # Insert node before the head
            tail = self.head.prev
            tail.next = node
            node.prev = tail
            node.next = self.head
            self.head.prev = node

    def insert_after(self, existing_node: Node, new_node: Node):
        """  Insert a new node into the list, exactly after a node that is already in
             list """
        if not self.head:
            raise RuntimeError("The list is empty. Cannot insert after a non-existent node.")

        # Insert the new node after the existing node
        new_node.prev = existing_node
        new_node.next = existing_node.next
        existing_node.next.prev = new_node
        existing_node.next = new_node
        new_node.in_list = True

    def remove(self, node: Node):
        """ Remove a node from the list """
        if not node.in_list:
            raise RuntimeError(f"Node {node.node_id} is not in the list.")

        node.in_list = False
        if node.next == node:  # Only one node in the list
            self.head = None
        else:
            node.prev.next = node.next
            node.next.prev = node.prev
            if self.head == node:
                self.head = node.next

        node.next = None
        node.prev = None

    def __repr__(self) -> str:
        nodes = []
        current = self.head
        if current is not None:
            while True:
                nodes.append(current.node_id)
                current = current.next
                if current == self.head:
                    break
        return f"DoublyLinkedCircularList({nodes})"

