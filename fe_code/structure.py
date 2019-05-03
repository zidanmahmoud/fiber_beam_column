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

        self._stiffness = None

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
        return 6*(node_id - 1) + DOF_INDEX_MAP[dof_type]

    def dof_from_index(self, index):
        raise NotImplementedError

    @property
    def no_dofs(self):
        return len(self._nodes) * 6


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

    @property
    def tangent_stiffness(self):
        """ Last NR iterations """
        return self._stiffness

    @tangent_stiffness.setter
    def tangent_stiffness(self, value):
        if isinstance(value, np.ndarray):
            if value.shape == (self.no_dofs, self.no_dofs):
                self._stiffness = value
            else:
                raise ValueError("Structure stiffness is of wrong shape.")
        else:
            raise ValueError("Structure stiffness must be a numpy.ndarray")

    def calculate_stiffness_matrix(self):
        stiffness_matrix = np.zeros((self.no_dofs, self.no_dofs))

        for element in self.elements:
            k_e = element.calculate_global_stiffness_matrix()
            i = [self.index_from_dof(dof) for dof in element.dofs]
            stiffness_matrix[np.ix_(i, i)] = k_e

        self.tangent_stiffness = stiffness_matrix

    def solve_displacement_control(self):
        pass

