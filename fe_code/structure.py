import numpy as np

from .node import Node
from .fiber_beam import FiberBeam

def debug(message, *args, **kwargs):
    print(message, *args, **kwargs)

def warning(message, *args, **kwargs):
    print('\33[93m'+"WARNING: "+message.upper()+'\33[0m', *args, **kwargs)

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

        self.load_factor = 0

        self._stiffness = None
        self._unbalanced_forces = None
        self._displacement_increment = None

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

    @property
    def displacement_increment(self):
        if self._displacement_increment is None:
            return np.zeros(self.no_dofs)
        return self._displacement_increment

    @displacement_increment.setter
    def displacement_increment(self, value):
        self._displacement_increment = value

    def _calculate_stiffness_matrix(self):
        stiffness_matrix = np.zeros((self.no_dofs, self.no_dofs))

        for element in self.elements:
            k_e = element.calculate_global_stiffness_matrix()
            i = [self.index_from_dof(dof) for dof in element.dofs]
            stiffness_matrix[np.ix_(i, i)] = k_e

        self.tangent_stiffness = stiffness_matrix

    def _calculate_force_vector(self):
        forces = np.zeros(self.no_dofs)
        for condition in self._newmann_conditions.items():
            dof, value = condition
            forces[self.index_from_dof(dof)] += self.load_factor * value
        return forces

    def _construct_unbalance_forces_first_iteration(self):
        self._unbalanced_forces = np.zeros(self.no_dofs)
        for condition in self._newmann_conditions.items():
            dof, value = condition
            self._unbalanced_forces[self.index_from_dof(dof)] += value

    def solve(self, max_ele_iterations):
        """
        main solution loop until element convergence
        """
        if self._unbalanced_forces is None: #First NR iteration
            self._construct_unbalance_forces_first_iteration()

        self._calculate_stiffness_matrix()

        for condition in self._dirichlet_conditions.items():
            dof, value = condition
            i = self.index_from_dof(dof)
            self.tangent_stiffness[:, i] = 0
            self.tangent_stiffness[i, :] = 0
            self.tangent_stiffness[i, i] = 1
            self._unbalanced_forces[i] = value

        change_in_displacement_increment = np.linalg.solve(
            self.tangent_stiffness, self._unbalanced_forces
        )
        self.displacement_increment += change_in_displacement_increment

        for element in self.elements:
            indices = [self.index_from_dof(dof) for dof in element.dofs]
            element.calculate_displacement_increment_from_structure(
                change_in_displacement_increment[indices]
            )

        conv = False
        for j in range(1, max_ele_iterations+1):
            for element in self.elements:
                element.calculate_force_increment()
                element.increment_resisting_forces()
                element.update_stiffness()

            for element in self.elements:
                conv += element.check_convergence()

            if conv:
                debug(f"Elements have converged with {j} iteration(s).")
                break

            for element in self.elements:
                element.update_chng_displacement_increment()

            if j == max_ele_iterations:
                warning(f"ELEMENTS DID NOT CONVERGE WITH {max_ele_iterations} ITERATIONS")

    def check_nr_convergence(self):
        resisting_forces = np.zeros(self.no_dofs)
        for element in self.elements:
            f_e = element.l_e @ element.resisting_forces
            i = [self.index_from_dof(dof) for dof in element.dofs]
            # stiffness_matrix[np.ix_(i, i)] = k_e
            resisting_forces[i] += f_e

        external_forces = self._calculate_force_vector()
        self._unbalanced_forces = external_forces - resisting_forces
        return abs(np.linalg.norm(self._unbalanced_forces)) < self._tolerance
