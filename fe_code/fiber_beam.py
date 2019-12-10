"""
fiber_beam
==========

Module contains the fiber beam element class
"""
import numpy as np

from .section import Section
from .dof import DoF
from .gauss_lobatto import gauss_lobatto
from .io import warning


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

    def __init__(self, element_id, node1, node2):
        self._id = element_id
        self._nodes = [node1, node2]
        self._sections = dict()

        self._force_increment = np.zeros(5)
        self.resisting_forces = np.zeros(5)
        self.converged_resisting_forces = np.zeros(5)
        self._displacement_residual = np.zeros(5)

        self._local_stiffness_matrix = np.zeros((5, 5))
        self._transform_matrix = np.zeros((12, 5))

        dof_types = "uvwxyz"
        self.dofs = [DoF(node.id, dof_type) for node in self._nodes for dof_type in dof_types]

    @property
    def id(self):
        return self._id

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
        self._sections[section_id] = Section(section_id)

    def get_section(self, section_id):
        return self._sections[section_id]


    ####################################################################################


    def initialize(self):
        """
        initialize matrices
        """
        points, weights = gauss_lobatto(len(self.sections))
        for i, section in enumerate(self.sections):
            section.position = points[i]
            section.weight = weights[i]
            section.initialize()
        self._calculate_transform_matrix()
        self._update_local_stiffness_matrix()

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

    def state_determination(self, structure_chng_disp_incr, max_ele_iterations):
        #== step 6 ==#
        chng_disp_incr = self._transform_matrix.T @ structure_chng_disp_incr
        a=5
        for j in range(1, max_ele_iterations + 1):
            #== step 7 ==#
            if j==1:
                chng_force_increment = self._local_stiffness_matrix @ chng_disp_incr
            else:
                chng_force_increment = -self._local_stiffness_matrix @ self._displacement_residual
            self._force_increment += chng_force_increment
            self.resisting_forces = self.converged_resisting_forces + self._force_increment
            #== steps 8-12 ==#
            conv = True
            for section in self.sections:
                conv *= section.state_determination(chng_force_increment)
                # print(f"section {section.id}:")
                # for fiber in section.fibers:
                #     print(fiber.strain)
                # print()
            #== step 13 ==#
            self._update_local_stiffness_matrix()
            #== step 14 ==#
            if conv:
                print(f"Element {self._id} converged with {j} iteration(s).")
                return # FIXME:break piece of shite
            else:
                J = self._get_jacobian_determinant()
                self._displacement_residual.fill(0.0)
                for section in self.sections:
                    self._displacement_residual += J * section.weight * section.get_global_residuals()
        warning(f"ELEMENTS DID NOT CONVERGE WITH {max_ele_iterations} ITERATIONS")

    def reset_section_residuals(self):
        for section in self.sections:
            section.reset_residual()

    def finalize_load_step(self):
        """
        finalize for next load step
        """
        self.converged_resisting_forces = self.resisting_forces
        self._force_increment.fill(0.0)
        for section in self.sections:
            section.finalize_load_step()


    ####################################################################################


    def _update_local_stiffness_matrix(self):
        """
        update_local_stiffness_matrix based on the section iterations
        """
        local_flexibility_matrix = np.zeros((5, 5))
        J = self._get_jacobian_determinant()
        for section in self.sections:
            local_flexibility_matrix += J * section.weight * section.get_global_flexibility_matrix()
        self._local_stiffness_matrix = np.linalg.inv(local_flexibility_matrix)

    def _get_jacobian_determinant(self):
        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)
        return reference_length / 2.0

    def _calculate_transform_matrix(self):
        """
        called only in initialize
        """
        reference_local_vector = (
            self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        )
        reference_length = np.linalg.norm(reference_local_vector)
        triad = self._get_triad()
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

    def _get_triad(self):
        node1_coords = self._nodes[0].get_reference_location()
        node2_coords = self._nodes[1].get_reference_location()
        v1 = node2_coords - node1_coords
        e1 = v1 / np.linalg.norm(v1)
        # TODO: Generalize e2 ...
        e2 = np.array([0, 1, 0])
        e3 = np.cross(e1, e2)
        return e1, e2, e3
