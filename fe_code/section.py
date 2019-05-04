"""
Section
"""
import numpy as np

from .fiber import Fiber

class Section:
    """ Section class
    Attributes
    ----------
    """
    def __init__(self):
        self._fibers = dict()
        self._tolerance = 1e-7

        self.position = None
        self.weight = None

        self._force_increment = None
        self._deformation_increment = None
        self.converged_section_forces = np.zeros(3)
        self.forces = np.zeros(3)
        self.residual = np.zeros(3)

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

    def add_fiber(self, fiber_id, y, z, area, material_class):
        """add a fiber to the section

        Parameters
        ----------
        """
        self._fibers[fiber_id] = Fiber(y, z, area, material_class)

    def initialize(self):
        """probably unnecessary"""
        for fiber in self.fibers:
            fiber.initialize()
        # stiffness matrix

    @property
    def force_increment(self):
        """force increment vector"""
        if self._force_increment is None:
            return np.zeros(3)
        return self._force_increment
    @force_increment.setter
    def force_increment(self, value):
        self._force_increment = value

    @property
    def deformation_increment(self):
        """deformation_increment vector"""
        if self._deformation_increment is None:
            return np.zeros(3)
        return self._deformation_increment
    @deformation_increment.setter
    def deformation_increment(self, value):
        self._deformation_increment = value


    def calculate_stiffness_matrix(self):
        """section stiffness matrix

        Returns
        -------
        stiffness_matrix : ndarray(3x3)
        """
        stiffness_matrix = np.zeros((3, 3))
        for fiber in self.fibers:
            EA = fiber.tangent_stiffness * fiber.area
            stiffness_matrix += EA * np.outer(fiber.direction, fiber.direction)
        return stiffness_matrix

    def calculate_force_increment_from_element(self, ele_chng_force_increment):
        """ step 8 """
        b_matrix = _calculate_b_matrix(self.position)
        self.chng_force_increment = b_matrix @ ele_chng_force_increment
        self.force_increment += self.chng_force_increment

    def increment_section_forces(self):
        self.forces = self.converged_section_forces + self.force_increment

    def calculate_deformation_increment(self):
        """ step 9 """
        f_e = np.linalg.inv(self.calculate_stiffness_matrix())
        self.chng_def_increment = self.residual + f_e @ self.chng_force_increment
        self.deformation_increment += self.chng_def_increment

    def calculate_fiber_deformation_increment(self):
        """ step 10 """
        for fiber in self.fibers:
            fiber.calculate_strain_increment_from_section(
                self.chng_def_increment
            )
            fiber.increment_strain()
            fiber.calculate_stress()

    def check_convergence(self):
        """ steps 11 & 12 """
        resisting_forces = np.zeros(3)
        for fiber in self.fibers:
            resisting_forces += fiber.stress * fiber.area * fiber.direction
        self.residual = np.linalg.inv(self.calculate_stiffness_matrix()) @ (self.forces - resisting_forces)
        return abs(np.linalg.norm(self.forces - resisting_forces)) < self._tolerance



def _calculate_b_matrix(gauss_point):
    b_matrix = np.zeros([3, 5])
    b_matrix[0, 0] = gauss_point / 2 - 1/2
    b_matrix[0, 1] = gauss_point / 2 + 1/2
    b_matrix[1, 2] = gauss_point / 2 - 1/2
    b_matrix[1, 3] = gauss_point / 2 + 1/2
    b_matrix[2, 4] = 1
    return b_matrix
