import numpy as np

from .element import Element
from .section import Section
from .gauss_lobatto import GaussLobatto

class FiberBeam(Element):
    def __init__(self, node1, node2):
        self._nodes = [node1, node2]
        self._sections = dict()

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
        for section in self.sections:
            section.initialize()
        # stiffness matrix
        # save NR

    def calculate_global_stiffness_matrix(self):
        k_e_local = self._calculate_local_stiffness_matrix()
        l_e = self._calculate_transform_matrix()
        return l_e @ k_e_local @ l_e.T

    def _calculate_local_stiffness_matrix(self):
        points, weights = GaussLobatto(len(self.sections))
        local_flexibility_matrix = np.zeros((5, 5))

        reference_local_vector = self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        reference_length = np.linalg.norm(reference_local_vector)

        for i, section in enumerate(self.sections):
            section_flexibility_matrix = np.linalg.inv(section.calculate_stiffness_matrix())
            b_matrix = _calculate_b_matrix(points[i])
            local_flexibility_matrix += reference_length / 2 * weights[i] * (b_matrix.T @ section_flexibility_matrix @ b_matrix)

        return np.linalg.inv(local_flexibility_matrix)

    def _calculate_transform_matrix(self):
        reference_local_vector = self.nodes[1].get_reference_location() - self.nodes[0].get_reference_location()
        reference_length = np.linalg.norm(reference_local_vector)
        triad = _get_triad(reference_local_vector)
        e_1, e_2, e_3 = triad
        reference_transform_matrix = np.zeros([12, 5])
        # Forces of first node
        reference_transform_matrix[0:3, 0] = e_2 /reference_length
        reference_transform_matrix[0:3, 1] = e_2 /reference_length
        reference_transform_matrix[0:3, 2] = -e_3 /reference_length
        reference_transform_matrix[0:3, 3] = -e_3 /reference_length
        reference_transform_matrix[0:3, 4] = -e_1
        # Moments of first node
        reference_transform_matrix[3:6, 0] = e_3
        reference_transform_matrix[3:6, 2] = e_2
        # Forces of second node
        reference_transform_matrix[6:9, :] = - reference_transform_matrix[0:3, :]
        # Moments of second node
        reference_transform_matrix[9:12, 1] = e_3
        reference_transform_matrix[9:12, 3] = e_2

        return reference_transform_matrix

def _calculate_b_matrix(gauss_point):
    b_matrix = np.zeros([3, 5])
    b_matrix[0, 0] = gauss_point / 2 - 1/2
    b_matrix[0, 1] = gauss_point / 2 + 1/2
    b_matrix[1, 2] = gauss_point / 2 - 1/2
    b_matrix[1, 3] = gauss_point / 2 + 1/2
    b_matrix[2, 4] = 1
    return b_matrix

def _get_triad(reference_local_vector):
    e_1 = reference_local_vector / np.linalg.norm(reference_local_vector)
    e_2 = np.cross([0, 0, 1], e_1)
    e_2 /= np.linalg.norm(e_2)
    e_3 = np.cross(e_1, e_2)
    # return np.column_stack((e_1, e_2, e_3))
    return e_1, e_2, e_3
