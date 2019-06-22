"""
fiber_beam
==========

Module contains the fiber beam element class
"""
import numpy as np

from .section import Section
from .dof import DoF
from .gauss_lobatto import GaussLobatto


class FiberBeam:
    """
    fiber-beam element class

    Attributes
    ----------
    nodes : dict_values

    sections : dict_values

    transform_matrix : ndarray

    converged_resisting_forces : ndarray
        converged in last load step
    resisting_forces : ndarray
        current
    displacement_residual : ndarray
    """

    def __init__(self, node1, node2):
        self._nodes = [node1, node2]
        self._sections = dict()

        self._local_stiffness_matrix = np.zeros((5, 5))
        self._transform_matrix = np.zeros((12, 5))

        self._chng_disp_incr = np.zeros(5)
        self._chng_force_increment = np.zeros(5)
        self._displacement_increment = np.zeros(5)
        self._force_increment = np.zeros(5)
        self.converged_resisting_forces = np.zeros(5)
        self.resisting_forces = np.zeros(5)
        self.displacement_residual = np.zeros(5)

        dof_types = "uvwxyz"
        self.dofs = [DoF(node.id, dof_type) for node in self._nodes for dof_type in dof_types]

    @property
    def nodes(self):
        """
        nodes as a dict_values
        """
        return self._nodes

    @property
    def sections(self):
        """
        sections as a dict_values
        """
        return self._sections.values()

    def add_section(self, section_id):
        """
        add a section
        """
        if section_id in self._sections:
            raise RuntimeError(f"Structure has already a section with id {section_id}")
        self._sections[section_id] = Section()

    def initialize(self):
        """
        initialize matrices
        """
        points, weights = GaussLobatto(len(self.sections))
        for i, section in enumerate(self.sections):
            section.initialize()
            section.position = points[i]
            section.weight = weights[i]
        self._calculate_transform_matrix()
        self.update_local_stiffness_matrix()

    def get_global_stiffness_matrix(self):
        """
        l_e^T . K_e . l_e
        """
        return self._transform_matrix @ self._local_stiffness_matrix @ self._transform_matrix.T

    def get_global_resisting_forces(self):
        """
        l_e^T . D_R
        """
        return self._transform_matrix @ self.resisting_forces

    def update_local_stiffness_matrix(self):
        """
        update_local_stiffness_matrix based on the section iterations
        """
        local_flexibility_matrix = np.zeros((5, 5))

        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)

        for section in self.sections:
            section_flexibility_matrix = np.linalg.inv(section.stiffness_matrix)
            b_matrix = _calculate_b_matrix(section.position)
            local_flexibility_matrix += (
                reference_length
                / 2
                * section.weight
                * (b_matrix.T @ section_flexibility_matrix @ b_matrix)
            )

        self._local_stiffness_matrix = np.linalg.inv(local_flexibility_matrix)

    def _calculate_transform_matrix(self):
        """
        called only in initialize
        """
        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)
        triad = _get_triad(reference_local_vector)
        e_1, e_2, e_3 = triad
        # Forces of first node
        self._transform_matrix[0:3, 0] = e_2 / reference_length
        self._transform_matrix[0:3, 1] = e_2 / reference_length
        self._transform_matrix[0:3, 2] = -e_3 / reference_length
        self._transform_matrix[0:3, 3] = -e_3 / reference_length
        self._transform_matrix[0:3, 4] = -e_1
        # Moments of first node
        self._transform_matrix[3:6, 0] = e_3
        self._transform_matrix[3:6, 2] = e_2
        # Forces of second node
        self._transform_matrix[6:9, :] = -self._transform_matrix[0:3, :]
        # Moments of second node
        self._transform_matrix[9:12, 1] = e_3
        self._transform_matrix[9:12, 3] = e_2

    def calculate_displacement_increment_from_structure(self, structure_chng_disp_incr):
        """step 4"""
        self._calculate_transform_matrix()
        self._chng_disp_incr = self._transform_matrix.T @ structure_chng_disp_incr
        self._displacement_increment += self._chng_disp_incr

    def calculate_forces(self):
        """ steps 6 & 7 """
        if not np.any(self.displacement_residual):
            self._chng_force_increment = self._local_stiffness_matrix @ self._chng_disp_incr
        else:
            self._chng_force_increment = -self._local_stiffness_matrix @ self.displacement_residual
        self._force_increment += self._chng_force_increment
        self.resisting_forces = self.converged_resisting_forces + self._force_increment

    def state_determination(self):
        """ steps 8-12 """
        for section in self.sections:
            section.calculate_force_increment_from_element(self._chng_force_increment)
            section.increment_section_forces()
            section.calculate_deformation_increment()
            section.calculate_fiber_deformation_increment()
            section.update_stiffness_matrix()

        self.update_local_stiffness_matrix()

    def check_convergence(self):
        """
        check all sections convergence
        """
        for section in self.sections:
            if not section.check_convergence():
                return False
        return True

    def calculate_displacement_residuals(self):
        """
        residuals if convergence is not satisfied
        """
        for section in self.sections:
            section.calculate_residuals()

        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)

        self.displacement_residual.fill(0.0)
        for section in self.sections:
            self.displacement_residual += (
                reference_length
                / 2.0
                * section.weight
                * _calculate_b_matrix(section.position).T
                @ section.residual
            )

    def finalize_load_step(self):
        """
        finalize for next load step
        """
        self._displacement_increment = np.zeros(5)
        self._force_increment = np.zeros(5)
        for section in self.sections:
            section.finalize_load_step()
        self.converged_resisting_forces = self.resisting_forces


def _calculate_b_matrix(gauss_point):
    b_matrix = np.zeros([3, 5])
    b_matrix[0, 0] = gauss_point / 2 - 1 / 2
    b_matrix[0, 1] = gauss_point / 2 + 1 / 2
    b_matrix[1, 2] = gauss_point / 2 - 1 / 2
    b_matrix[1, 3] = gauss_point / 2 + 1 / 2
    b_matrix[2, 4] = 1
    return b_matrix


def _get_triad(reference_local_vector):
    e_1 = reference_local_vector / np.linalg.norm(reference_local_vector)
    e_2 = np.cross([0, 0, 1], e_1)
    e_2 /= np.linalg.norm(e_2)
    e_3 = np.cross(e_1, e_2)
    return e_1, e_2, e_3
