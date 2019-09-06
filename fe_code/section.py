"""
Section
"""
import numpy as np

from .fiber import Fiber
from .material_laws import Material


class Section:
    """ Section class

    Attributes
    ----------
    fibers : dict_values

    position : float
        position based on Gauss-Lobatto rule
    weight : float
        weight based on Gauss-Lobatto rule
    residual : ndarray
        norm of the unbalance forces
    """

    def __init__(self, section_id):
        self._id = section_id
        self._fibers = dict()
        self._tolerance = 1e-7

        self.position = None
        self.weight = None

        self.stiffness_matrix = np.zeros((3, 3))

        self._chng_def_increment = np.zeros(3)
        self._chng_force_increment = np.zeros(3)
        self._deformation_increment = np.zeros(3)
        self._force_increment = np.zeros(3)
        self._converged_section_forces = np.zeros(3)
        self._forces = np.zeros(3)
        self._unbalance_forces = np.zeros(3)
        self.residual = np.zeros(3)
        self._b_matrix = None

    @property
    def id(self):
        return self._id

    @property
    def tolerance(self):
        """N-R iterations tolerance"""
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value):
        self._tolerance = value

    @property
    def fibers(self):
        """fibers list"""
        return self._fibers.values()

    @property
    def b_matrix(self):
        """ b_matrix """
        if self._b_matrix is None:
            self._b_matrix = _calculate_b_matrix(self.position)
        return self._b_matrix

    def add_fiber(self, fiber_id, y, z, area, material_class, w, h):
        """add a fiber to the section

        Parameters
        ----------
        fiber_id : int
            id of the fiber
        y : float
            y coordinate
        z : float
            z coordinate
        area : float
            area of the fiber
        material_class : object of type Material
            material model
        """
        # if not isinstance(material_class, Material):
        #     raise ValueError("material_class is not of type : Material")
        self._fibers[fiber_id] = Fiber(fiber_id, y, z, area, material_class, w, h)

    def initialize(self):
        """initialize stiffness matrix"""
        self.update_stiffness_matrix()

    def update_stiffness_matrix(self):
        """ section stiffness matrix """
        self.stiffness_matrix.fill(0.0)
        for fiber in self.fibers:
            EA = fiber.tangent_stiffness * fiber.area
            self.stiffness_matrix += EA * fiber.direction_matrix

    def calculate_force_increment_from_element(self, ele_chng_force_increment):
        """ step 8 """
        b_matrix = self.b_matrix
        self._chng_force_increment = b_matrix @ ele_chng_force_increment
        self._force_increment += self._chng_force_increment

    def increment_section_forces(self):
        """ step 8 """
        self._forces = self._converged_section_forces + self._force_increment

    def calculate_deformation_increment(self):
        """ step 9 """
        f_e = np.linalg.inv(self.stiffness_matrix)
        self._chng_def_increment = self.residual + f_e @ self._chng_force_increment
        self._deformation_increment += self._chng_def_increment

    def calculate_fiber_deformation_increment(self):
        """ step 10 """
        rev = 0
        for fiber in self.fibers:
            reversal = fiber.calculate_strain_increment_from_section(self._chng_def_increment)
            rev += reversal
            fiber.increment_strain()
            fiber.calculate_stress()
        return rev

    def check_convergence(self):
        """ steps 13-15 """
        return abs(np.linalg.norm(self._unbalance_forces)) < self._tolerance

    def calculate_residuals(self):
        """
        residuals (norm of the unbalance forces)
        """
        resisting_forces = np.zeros(3)
        for fiber in self.fibers:
            resisting_forces += fiber.stress * fiber.area * fiber.direction
        self._unbalance_forces = self._forces - resisting_forces
        self.residual = np.linalg.inv(self.stiffness_matrix) @ self._unbalance_forces

    def finalize_load_step(self):
        """
        finalize for next load step
        """
        self._converged_section_forces = self._forces
        self._deformation_increment = np.zeros(3)
        self._force_increment = np.zeros(3)
        for fiber in self.fibers:
            fiber.finalize_load_step()


def _calculate_b_matrix(gauss_point):
    b_matrix = np.zeros([3, 5])
    b_matrix[0, 0] = gauss_point / 2 - 1 / 2
    b_matrix[0, 1] = gauss_point / 2 + 1 / 2
    b_matrix[1, 2] = gauss_point / 2 - 1 / 2
    b_matrix[1, 3] = gauss_point / 2 + 1 / 2
    b_matrix[2, 4] = 1
    return b_matrix
