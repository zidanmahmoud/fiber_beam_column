"""
Module contains the structure class
"""
import numpy as np

from .node import Node
from .dof import DoF
from .fiber_beam import FiberBeam
from .io import warning


DOF_INDEX_MAP = {"u": 0, "v": 1, "w": 2, "x": 3, "y": 4, "z": 5}


def index_from_dof(dof):
    """get the index from dof"""
    return 6 * (dof.node_id - 1) + DOF_INDEX_MAP[dof.type]


def dof_from_index(index):
    """get the dof from index"""
    node_id = 1 + index // 6
    for key, value in DOF_INDEX_MAP.items():
        if value == index % 6:
            dof_type = key
    return DoF(node_id, dof_type)


class Structure:
    """
    Structure class

    Attributes
    ----------
    nodes : dict_values

    elements : dict_values

    controled_dof : Dof object

    no_dofs : flaot
        number of dofs

    converged_load_factor : float

    converged_displacement : ndarray

    controled_dof_increment : float
        used in the displacement-control solver

    tolerance : float
        used to check convergence
        default is 1e-7
    """

    def __init__(self):
        self._nodes = dict()
        self._elements = dict()
        self._dirichlet_conditions = dict()
        self._newmann_conditions = dict()
        self._tolerance = 1e-7

        self.controled_dof = None

        self._load_factor_increment = 0.0
        self._load_factor = 0.0
        self.converged_load_factor = 0.0
        self.controled_dof_increment = 0.0

        # initialized as None because the number of dofs is not yet determined
        self._stiffness = None
        self._unbalanced_forces = None
        self._displacement_increment = None
        self._displacement = None
        self.converged_displacement = None

    @property
    def tolerance(self):
        """ tolerance to check convergence """
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value):
        self._tolerance = value

    def set_section_tolerance(self, value):
        """ set convergence tolerance for sections """
        for element in self.elements:
            for section in element.sections:
                section.tolerance = value

    @property
    def nodes(self):
        """ nodes """
        return self._nodes.values()

    def get_node(self, node_id):
        """ get node with id """
        return self._nodes[node_id]

    @property
    def elements(self):
        """ elements """
        return self._elements.values()

    def get_element(self, element_id):
        """ get element with id """
        return self._elements[element_id]

    @property
    def no_dofs(self):
        """ number of dofs """
        return len(self._nodes) * 6

    def add_node(self, node_id, x_pos, y_pos, z_pos):
        """ add a node """
        self._nodes[node_id] = Node(node_id, x_pos, y_pos, z_pos)

    def add_fiber_beam_element(self, element_id, node1_id, node2_id):
        """ add an element """
        node1 = self.get_node(node1_id)
        node2 = self.get_node(node2_id)
        self._elements[element_id] = FiberBeam(element_id, node1, node2)

    def add_dirichlet_condition(self, node_id, dof_types, value):
        """ add a dirichlet boundary condition """
        for dof_type in dof_types:
            dof = DoF(node_id, dof_type)
            self._dirichlet_conditions[dof] = value

    def add_neumann_condition(self, node_id, dof_types, value):
        """ add a Neumann boundary condition """
        for dof_type in dof_types:
            dof = DoF(node_id, dof_type)
            self._newmann_conditions[dof] = value

    def set_controled_dof(self, node_id, dof_type):
        """ sets the controled dof """
        self.controled_dof = DoF(node_id, dof_type)

    @property
    def converged_controled_dof(self):
        """ converged_controled_dof """
        return self.converged_displacement[index_from_dof(self.controled_dof)]

    def initialize(self):
        """ initialize all arrays and stuff """
        self._stiffness = np.zeros((self.no_dofs, self.no_dofs))
        self._displacement_increment = np.zeros(self.no_dofs)
        self._unbalanced_forces = np.zeros(self.no_dofs)
        self._displacement_increment = np.zeros(self.no_dofs)
        self._displacement = np.zeros(self.no_dofs)
        self.converged_displacement = np.zeros(self.no_dofs)

        for element in self.elements:
            element.initialize()
        self._calculate_stiffness_matrix()

    def _calculate_stiffness_matrix(self):
        stiffness_matrix = np.zeros((self.no_dofs, self.no_dofs))

        for element in self.elements:
            k_e = element.get_global_stiffness_matrix()
            i = [index_from_dof(dof) for dof in element.dofs]
            stiffness_matrix[np.ix_(i, i)] = k_e  # stiffness_matrix[i][:i] = k_e

        self._stiffness = stiffness_matrix

    def _calculate_force_vector(self):
        forces = np.zeros(self.no_dofs)
        dof = self.controled_dof
        forces[index_from_dof(dof)] += self._load_factor
        return forces

    def solve(self, max_ele_iterations):
        """
        main solution loop until element convergence

        steps 3-17
        """
        dofs = self.no_dofs
        lhs = np.zeros((dofs + 1, dofs + 1))
        lhs[:dofs, :dofs] = self._stiffness
        for dof, value in self._dirichlet_conditions.items():
            if value == 0:
                i = index_from_dof(dof)
                lhs[:, i] = 0
                lhs[i, :] = 0
                lhs[i, i] = 1

        vector = np.zeros(dofs)
        vector[index_from_dof(self.controled_dof)] = 1.0
        lhs[:dofs, -1] = -vector
        lhs[-1, :dofs] = -vector
        rhs = np.zeros(dofs + 1)
        rhs[:dofs] = self._unbalanced_forces
        rhs[-1] = vector @ self._displacement_increment - self.controled_dof_increment

        change_in_increments = np.linalg.solve(lhs, rhs)
        self._displacement_increment += change_in_increments[:dofs]
        self._load_factor_increment += change_in_increments[-1]
        self._displacement = self.converged_displacement + self._displacement_increment
        self._load_factor = self.converged_load_factor + self._load_factor_increment

        # STEP 4
        for element in self.elements:
            indices = [index_from_dof(dof) for dof in element.dofs]
            element.calculate_displacement_increment_from_structure(
                change_in_increments[:dofs][indices]
            )

        # STEP 5
        for j in range(1, max_ele_iterations + 1):
            conv, rev = self.element_loop()
            if rev > 0:
                print(f" --- reversed {rev} fibers --- ")

            # STEP 17
            if conv:  # all elements converged
                self._calculate_stiffness_matrix()
                for element in self.elements:
                    element.displacement_residual.fill(0.0)
                    for section in element.sections:
                        section.residual.fill(0.0)
                print(f"Elements converged with {j} iteration(s).")
                break

            if j == max_ele_iterations:
                warning(f"ELEMENTS DID NOT CONVERGE WITH {max_ele_iterations} ITERATIONS")

    def element_loop(self):
        """
        FIXME
        """
        conv = True
        rev = 0
        for element in self.elements:
            element.calculate_forces()
            rev += element.state_determination()
            element.calculate_displacement_residuals()
            # print(np.linalg.norm(element.displacement_residual))
            conv *= element.check_convergence()
        return conv, rev

    def check_nr_convergence(self):
        """ steps 18-20 """
        resisting_forces = np.zeros(self.no_dofs)
        for element in self.elements:
            f_e = element.get_global_resisting_forces()
            i = [index_from_dof(dof) for dof in element.dofs]
            resisting_forces[i] += f_e
        for dof, value in self._dirichlet_conditions.items():
            if value == 0:
                resisting_forces[index_from_dof(dof)] = value

        external_forces = self._calculate_force_vector()

        self._unbalanced_forces = external_forces - resisting_forces
        res = abs(np.linalg.norm(self._unbalanced_forces))
        return res < self._tolerance, res

    def finalize_load_step(self):
        """ step 21 """
        self._displacement_increment.fill(0.0)
        self._load_factor_increment = 0.0
        self.converged_displacement = self._displacement
        self.converged_load_factor = self._load_factor
        for node in self.nodes:
            indices = [index_from_dof(DoF(node.id, dof_type)) for dof_type in "uvw"]
            node.u, node.v, node.w = self.converged_displacement[indices]
        for element in self.elements:
            element.finalize_load_step()
