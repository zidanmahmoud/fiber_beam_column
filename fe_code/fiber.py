"""
Module contains only the fiber class
"""
import numpy as np


class Fiber:
    """
    Fiber class

    Parameters
    ----------
    y : float
        y location
    z : float
        z location
    area : float
        section area of the fiber
    material_class : Material object

    Attributes
    ----------
    converged_strain : float
        converged from last load step
    strain : float
        current
    direction : ndarray
        fiber to section variables
    area : float
        section area
    tangent_stiffness : float
        material stiffness
    """

    def __init__(self, fiber_id, y, z, area, material_class, w, h):
        self._id = fiber_id
        self.direction = np.array([-y, z, 1.0])
        self.direction_matrix = np.outer(self.direction, self.direction)
        self.area = area
        self._material = material_class

        self.w = w
        self.h = h

        self.strain = 0.0
        self.stress = 0.0

        self._rev = True

    @property
    def id(self):
        return self._id

    @property
    def tangent_stiffness(self):
        """
        material stiffness
        """
        return self._material.tangent_modulus

    def state_determination(self, sec_def_increment):
        """ step 10 + 11 """
        strain_increment = self.direction @ sec_def_increment
        self.strain += strain_increment
        self._material.update_strain(self.strain)
        self.stress = self._material.stress

    def finalize_load_step(self):
        """
        finalize for next load step
        """
        self._material.finalize_load_step()
