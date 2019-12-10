"""
Module contains the structure class
"""
import numpy as np
from copy import copy

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

    controlled_dof : Dof object

    no_dofs : flaot
        number of dofs

    controlled_dof_increment : float
        used in the displacement-control solver
    """

    def __init__(self):
        self._nodes = dict()
        self._elements = dict()
        self._dirichlet_conditions = dict()
        self._neumann_conditions = dict()
        self._tolerance = 1e-7

        self._controlled_dof = None

        self._load_factor_increment = 0.0
        self._load_factor = 0.0
        self._converged_load_factor = 0.0
        self.controlled_dof_increment = 0.0

        # initialized as None because the number of dofs is not yet determined
        self._stiffness_matrix = None
        self._resisting_forces = None
        self._unbalanced_forces = None
        self._previous_dx = None
        self._displacement_increment = None
        self._displacement = None
        self._converged_displacement = None

    def set_tolerance(self, value):
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
            self._neumann_conditions[dof] = value

    def set_controlled_dof(self, node_id, dof_type):
        """ sets the controlled dof """
        self._controlled_dof = DoF(node_id, dof_type)

    def get_force(self, dof):
        node_id, dof_type = dof
        i = index_from_dof(DoF(node_id, dof_type))
        return self._resisting_forces[i]

    def get_forces(self):
        return self._resisting_forces

    def get_displacement(self, dof):
        node_id, dof_type = dof
        i = index_from_dof(DoF(node_id, dof_type))
        return self._converged_displacement[i]

    def get_displacements(self):
        return self._converged_displacement

    def get_dof_value(self, dof):
        if isinstance(dof, DoF):
            node_id = dof.node_id
            dof_type = dof.type
        elif isinstance(dof, tuple):
            node_id, dof_type = dof
        else:
            raise
        i = index_from_dof(DoF(node_id, dof_type))
        return self._converged_displacement[i]

    ####################################################################################

    def initialize(self):
        """ initialize all arrays and stuff """
        #== step 1 ==#
        self._stiffness_matrix = np.zeros((self.no_dofs, self.no_dofs))
        self._displacement_increment = np.zeros(self.no_dofs)
        self._unbalanced_forces = np.zeros(self.no_dofs)
        self._previous_dx = np.zeros(self.no_dofs + 1)
        self._displacement_increment = np.zeros(self.no_dofs)
        self._displacement = np.zeros(self.no_dofs)
        self._converged_displacement = np.zeros(self.no_dofs)
        self._resisting_forces = np.zeros(self.no_dofs)
        for element in self.elements:
            element.initialize()
        self._update_stiffness_matrix()

    def solve_NR_iteration(self, max_ele_iterations):
        """
        main solution loop until element convergence
        """
        #== step 4 ==#
        lhs, rhs = self._build_system_NR_displacement_control()
        self._apply_homogenuous_dirichlet_BCs(lhs, rhs)

        dx = np.linalg.solve(lhs, rhs)
        change_in_increments = dx - self._previous_dx
        self._previous_dx = dx
        self._displacement_increment = dx[:self.no_dofs]
        self._load_factor_increment = dx[-1]
        self._displacement = self._converged_displacement + self._displacement_increment
        self._load_factor = self._converged_load_factor + self._load_factor_increment

        #== steps 5-14 ==#
        for element in self.elements:
            indices = [index_from_dof(dof) for dof in element.dofs]
            element.state_determination(
                change_in_increments[:self.no_dofs][indices], max_ele_iterations)

        #== step 15 ==#
        self._update_stiffness_matrix()
        for element in self.elements:
            element.reset_section_residuals()

        self._resisting_forces.fill(0.0)
        for element in self.elements:
            f_e = element.get_global_resisting_forces()
            i = [index_from_dof(dof) for dof in element.dofs]
            self._resisting_forces[i] += f_e

        external_forces = self._get_external_force_vector() * self._load_factor

        self._unbalanced_forces = external_forces - self._resisting_forces
        unbalance = copy(self._unbalanced_forces)
        self._apply_homogenuous_dirichlet_BCs(vector=unbalance)

        #== step 16 ==#
        res = abs(np.linalg.norm(unbalance))
        return res < self._tolerance, res


    def finalize_load_step(self):
        self._displacement_increment.fill(0.0)
        self._load_factor_increment = 0.0
        self._converged_displacement = self._displacement
        self._converged_load_factor = self._load_factor
        for node in self.nodes:
            indices = [index_from_dof(DoF(node.id, dof_type)) for dof_type in "uvw"]
            node.u, node.v, node.w = self._converged_displacement[indices]
        for element in self.elements:
            element.finalize_load_step()


    ####################################################################################


    def _update_stiffness_matrix(self):
        self._stiffness_matrix.fill(0.0)
        for element in self.elements:
            k_e = element.get_global_stiffness_matrix()
            i = [index_from_dof(dof) for dof in element.dofs]
            self._stiffness_matrix[np.ix_(i, i)] = k_e  # stiffness_matrix[i][:i] = k_e

    def _apply_homogenuous_dirichlet_BCs(self, matrix=None, vector=None):
        if matrix is not None:
            for dof, value in self._dirichlet_conditions.items():
                if value == 0:
                    i = index_from_dof(dof)
                    matrix[:, i] = 0
                    matrix[i, :] = 0
                    matrix[i, i] = 1
        if vector is not None:
            for dof, value in self._dirichlet_conditions.items():
                if value == 0:
                    i = index_from_dof(dof)
                    vector[i] = 0

    def _get_external_force_vector(self):
        external_forces = np.zeros(self.no_dofs)
        for dof, value in self._neumann_conditions.items():
            external_forces[index_from_dof(self._controlled_dof)] += value
        return external_forces

    def _build_system_NR_displacement_control(self):
        i = index_from_dof(self._controlled_dof)
        dofs = self.no_dofs

        lhs = np.zeros((dofs + 1, dofs + 1))
        lhs[:dofs, :dofs] = self._stiffness_matrix # dr/du
        lhs[:dofs, -1] = -self._get_external_force_vector() # dr/d lam
        lhs[-1, i] = 1.0 # dr/dC

        rhs = np.zeros(dofs + 1)
        rhs[:dofs] = self._unbalanced_forces # r
        rhs[-1] = self.controlled_dof_increment - self._displacement[i]# C
        return lhs, rhs
