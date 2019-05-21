"""
fiber_beam
==========

Module contains the fiber beam element class
"""
import numpy as np

from .element import Element
from .section import Section
from .gauss_lobatto import GaussLobatto


class FiberBeam(Element):
    def __init__(self, node1, node2):
        self._nodes = [node1, node2]
        self._sections = dict()

        self._local_stiffness_matrix = np.zeros((5, 5))

        self.chng_disp_incr = None
        self.chng_force_increment = None
        self._displacement_increment = None
        self._force_increment = None
        self.converged_resisting_forces = np.zeros(5)
        self.resisting_forces = np.zeros(5)
        # self._displacement_residual = None

    @property
    def nodes(self):
        return self._nodes

    @property
    def sections(self):
        return self._sections.values()

    @property
    def dofs(self):
        dof_types = "uvwxyz"
        return [(node.id, dof) for node in self.nodes for dof in dof_types]

    def add_section(self, section_id):
        if section_id in self._sections:
            raise RuntimeError(f"Structure has already a section with id {section_id}")
        self._sections[section_id] = Section()

    def initialize(self):
        points, weights = GaussLobatto(len(self.sections))
        for i, section in enumerate(self.sections):
            section.initialize()
            section.position = points[i]
            section.weight = weights[i]
        self.update_local_stiffness_matrix()
        # save NR

    def calculate_global_stiffness_matrix(self):
        l_e = self._calculate_transform_matrix()
        return l_e @ self._local_stiffness_matrix @ l_e.T

    def update_local_stiffness_matrix(self):
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
        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)
        triad = _get_triad(reference_local_vector)
        e_1, e_2, e_3 = triad
        reference_transform_matrix = np.zeros([12, 5])
        # Forces of first node
        reference_transform_matrix[0:3, 0] = e_2 / reference_length
        reference_transform_matrix[0:3, 1] = e_2 / reference_length
        reference_transform_matrix[0:3, 2] = -e_3 / reference_length
        reference_transform_matrix[0:3, 3] = -e_3 / reference_length
        reference_transform_matrix[0:3, 4] = -e_1
        # Moments of first node
        reference_transform_matrix[3:6, 0] = e_3
        reference_transform_matrix[3:6, 2] = e_2
        # Forces of second node
        reference_transform_matrix[6:9, :] = -reference_transform_matrix[0:3, :]
        # Moments of second node
        reference_transform_matrix[9:12, 1] = e_3
        reference_transform_matrix[9:12, 3] = e_2

        return reference_transform_matrix

    @property
    def displacement_increment(self):
        if self._displacement_increment is None:
            return np.zeros(5)
        return self._displacement_increment

    @displacement_increment.setter
    def displacement_increment(self, value):
        self._displacement_increment = value

    @property
    def force_increment(self):
        if self._force_increment is None:
            return np.zeros(5)
        return self._force_increment

    @force_increment.setter
    def force_increment(self, value):
        self._force_increment = value

    @property
    def l_e(self):
        return self._calculate_transform_matrix()

    # @property
    # def displacement_residual(self):
    #     if self._displacement_residual is None:
    #         return np.zeros(5)
    #     return self._displacement_residual
    # @displacement_residual.setter
    # def displacement_residual(self, value):
    #     self._displacement_residual = value

    def calculate_displacement_increment_from_structure(self, structure_chng_disp_incr):
        """step 4"""
        l_e = self._calculate_transform_matrix()
        self.chng_disp_incr = l_e.T @ structure_chng_disp_incr
        self.displacement_increment += self.chng_disp_incr

    def calculate_force_increment(self):
        """ steps 6 & 7 """
        self.chng_force_increment = self._local_stiffness_matrix @ self.chng_disp_incr
        # print(self._local_stiffness_matrix)
        # print("\t@")
        # print(self.chng_disp_incr)
        # print("\t=")
        # print(self.chng_force_increment)
        # input()
        self.force_increment += self.chng_force_increment

    def increment_resisting_forces(self):
        """ step 7 """
        self.resisting_forces = self.converged_resisting_forces + self.force_increment

    def state_determination(self):
        """ steps 8-12 """
        for section in self.sections:
            section.calculate_force_increment_from_element(self.chng_force_increment)
            section.increment_section_forces()
            section.calculate_deformation_increment()
            section.calculate_fiber_deformation_increment()
            section.update_stiffness_matrix()

        self.update_local_stiffness_matrix()

    def check_convergence(self):
        # for section in self.sections:
        #     if not section.check_convergence():
        #         return False
        # return True
        conv = True
        for section in self.sections:
            conv *= section.check_convergence()
            # print(section.forces)
        # input()
        return conv

    def calculate_displacement_residuals(self):

        for section in self.sections:
            section.calculate_displacement_residuals()

        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)

        residual = np.zeros(5)
        for section in self.sections:
            residual += (
                reference_length
                / 2.0
                * section.weight
                * _calculate_b_matrix(section.position).T
                @ section.residual
            )
        self.chng_disp_incr = -1 * residual

    def finalize_load_step(self):
        self._force_increment = None
        self._displacement_increment = None
        self.converged_resisting_forces = self.resisting_forces
        for section in self.sections:
            section.finalize_load_step()


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
    # return np.column_stack((e_1, e_2, e_3))
    return e_1, e_2, e_3
