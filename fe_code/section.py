"""
Section
"""
import numpy as np

from .fiber import Fiber
from .material_laws import UniaxialIncrementalMaterial


class Section:
    """ Section class

    Attributes
    ----------
    fibers : dict_values

    position : float
        position based on Gauss-Lobatto rule
    weight : float
        weight based on Gauss-Lobatto rule
    """

    def __init__(self, section_id):
        self._id = section_id
        self._fibers = dict()
        self._tolerance = 1e-7

        self.position = None
        self.weight = None

        self._residual = np.zeros(3)
        self._force_increment = np.zeros(3)
        self._forces = np.zeros(3)
        self._converged_section_forces = np.zeros(3)
        self._unbalance_forces = np.zeros(3)

        self._flexibility_matrix = np.zeros((3, 3))
        self._b_matrix = np.zeros((3, 5))

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
        if not isinstance(material_class, UniaxialIncrementalMaterial):
            raise ValueError("material_class is not of type : UniaxialIncrementalMaterial")
        if fiber_id in self._fibers:
            raise RuntimeError(f"Section already contains fiber with id {fiber_id}")
        self._fibers[fiber_id] = Fiber(fiber_id, y, z, area, material_class, w, h)

    def get_fiber(self, fiber_id):
        return self._fibers[fiber_id]


    ####################################################################################


    def initialize(self):
        """initialize stiffness matrix"""
        self._calculate_b_matrix()
        self._update_flexibility_matrix()

    def get_global_flexibility_matrix(self):
        return self._b_matrix.T @ self._flexibility_matrix @ self._b_matrix

    def get_global_residuals(self):
        return self._b_matrix.T @ self._residual

    def state_determination(self, ele_chng_force_increment):
        #== step 8 ==#
        chng_force_increment = self._b_matrix @ ele_chng_force_increment
        self._force_increment += chng_force_increment
        self._forces = self._converged_section_forces + self._force_increment
        #== step 9 ==#
        chng_def_increment = self._residual + self._flexibility_matrix @ chng_force_increment
        #== step 10 ==#
        for fiber in self.fibers:
            fiber.state_determination(chng_def_increment)
        #== step 11 ==#
        self._update_flexibility_matrix()
        #== step 12 ==#
        resisting_forces = np.zeros(3)
        for fiber in self.fibers:
            resisting_forces += fiber.stress * fiber.area * fiber.direction
        self._unbalance_forces = self._forces - resisting_forces
        self._residual = self._flexibility_matrix @ self._unbalance_forces
        return abs(np.linalg.norm(self._unbalance_forces)) < self._tolerance

    def reset_residual(self):
        self._residual.fill(0.0)

    def finalize_load_step(self):
        """
        finalize for next load step
        """
        self._converged_section_forces = self._forces
        self._force_increment.fill(0.0)
        for fiber in self.fibers:
            fiber.finalize_load_step()


    ####################################################################################


    def _update_flexibility_matrix(self):
        """ section stiffness matrix """
        stiffness_matrix = np.zeros((3, 3))
        for fiber in self.fibers:
            EA = fiber.tangent_stiffness * fiber.area
            stiffness_matrix += EA * fiber.direction_matrix
        self._flexibility_matrix = np.linalg.inv(stiffness_matrix)

    def _calculate_b_matrix(self):
        self._b_matrix[0, 0] = self.position / 2 - 1 / 2
        self._b_matrix[0, 1] = self.position / 2 + 1 / 2
        self._b_matrix[1, 2] = self.position / 2 - 1 / 2
        self._b_matrix[1, 3] = self.position / 2 + 1 / 2
        self._b_matrix[2, 4] = 1
