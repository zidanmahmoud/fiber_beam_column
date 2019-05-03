import numpy as np

from .node import Node
from .fiber_beam import FiberBeam

DOF_INDEX_MAP = {
    "u": 0,
    "v": 1,
    "w": 2,
    "x": 3,
    "y": 4,
    "z": 5,
}

class Structure:
    def __init__(self):
        self._nodes = dict()
        self._elements = dict()
        self._dirichlet_conditions = dict()
        self._newmann_conditions = dict()
        self._tolerance = 1e-7

    @property
    def tolerance(self):
        return self._tolerance
    @tolerance.setter
    def tolerance(self, value):
        self._tolerance = value

    @property
    def nodes(self):
        return self._nodes.values()

    def get_node(self, node_id):
        return self._nodes[node_id]

    @property
    def elements(self):
        return self._elements.values()

    def get_element(self, element_id):
        return self._elements[element_id]


    def index_from_dof(self, dof):
        node_id, dof_type = dof
        return 3*(node_id - 1) + DOF_INDEX_MAP[dof_type]

    def dof_from_index(self, index):
        raise NotImplementedError


    def add_node(self, node_id, x_pos, y_pos, z_pos):
        self._nodes[node_id] = Node(node_id, x_pos, y_pos, z_pos)

    def add_fiber_beam_element(self, element_id, node1_id, node2_id):
        self._elements[element_id] = FiberBeam(self.get_node(node1_id), self.get_node(node2_id))

    def add_dirichlet_condition(self, node_id, dof_types, value):
        for dof_type in dof_types:
            dof = (node_id, dof_type)
            self._dirichlet_conditions[dof] = value

    def add_neumann_condition(self, node_id, dof_types, value):
        for dof_type in dof_types:
            dof = (node_id, dof_type)
            self._newmann_conditions[dof] = value

    def initialize(self):
        for element in self.elements:
            element.initialize()
        # self.calculate_stiffness_matrix()

    def calculate_stiffness_matrix(self):
        dof_per_node = 6
        no_dofs = len(self.nodes)*dof_per_node
        stiffness_matrix = np.zeros((no_dofs, no_dofs))

        for element in self.elements:
            k_e = element.calculate_global_stiffness_matrix()
            i = np.array([self.index_from_dof(dof) for dof in element.dofs], dtype=int)
            stiffness_matrix[np.ix_(i, i)] = k_e

        return k_e

    def solve_displacement_control(self):
        pass

