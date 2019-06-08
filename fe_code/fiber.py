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
    ny : int
        fiber number in y direction
    nz : int
        fiber number in z direction
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

    def __init__(self, y, z, ny, nz, area, material_class):
        self.direction = np.array([-y, z, 1.0])
        self.direction_matrix = np.outer(self.direction, self.direction)
        # TODO: get rid of ny, nz
        self._ny = ny
        self._nz = nz
        self.area = area
        self._material_class = material_class

        self._first_iteration = True  # i=1
        self._chng_strain_increment = 0.0
        self._strain_increment = 0.0
        self.converged_strain = 0.0
        self.strain = 0.0
        self.stress = 0.0

    @property
    def tangent_stiffness(self):
        """
        material stiffness
        """
        return self._material_class.get_material_tangent_modulus()

    def initialize(self):
        """
        TODO: get rid of this
        """
        self._material_class.determin_direction(self._nz)

    def calculate_strain_increment_from_section(self, sec_chng_def_increment):
        """ step 10 + 11 """
        self._chng_strain_increment = self.direction @ sec_chng_def_increment
        self._strain_increment += self._chng_strain_increment

        self._material_class.update_change_in_material_strain_incr(self._chng_strain_increment)
        self._material_class.update_material_strain_incr()
        self._material_class.update_material_strain()
        if self._first_iteration:  # FIXME: this should be a check reversal
            self._material_class.update_model_parameters(self._nz)
        self._first_iteration = False
        self._material_class.update_material_stress()
        self._material_class.update_material_tangent_modulus()

    def increment_strain(self):
        """ step 10 """
        self.strain = self.converged_strain + self._strain_increment

    def calculate_stress(self):
        """ step 11 """
        self.stress = self._material_class.get_material_stress()

    def finalize_load_step(self):
        """
        finalize for next load step
        """
        self._material_class.save_to_last_loadstep_material_strain_incr(True)
        self._material_class.save_to_material_strain_last_loadstep(True)
        self.converged_strain = self.strain
        self._material_class.initialize_material_strain_incr()
        self._strain_increment = 0.0

    def reverse_loading(self):
        """
        to reverse material stuff ONLY in the next first iteration
        """
        self._first_iteration = True
