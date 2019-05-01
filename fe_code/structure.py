from .node import Node

class Structure:
    def __init__(self):
        self._nodes = dict()
        self._elements = dict()

    def add_node(self, node_id, x_pos, y_pos, z_pos):
        self._nodes[node_id] = Node(node_id, x_pos, y_pos, z_pos)

    def add_fiber_beam_element(self, element_id, element):
        self._elements[element_id] = element
