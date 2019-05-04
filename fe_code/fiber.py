import numpy as np

class Fiber:
    def __init__(self, y, z, area, material_class):
        self._y = y
        self._z = z
        self._area = area
        self._material_class = material_class

        self._first_iteration = True #i=1 & j=1
        self._strain_increment = None
        self.converged_strain = 0.0
        self.strain = 0.0

    @property
    def direction(self):
        return np.array([self._y, self._z, 1.0])

    @property
    def area(self):
        return self._area

    @property
    def tangent_stiffness(self):
        return self._material_class.get_material_tangent_modulus()

    def initialize(self):
        nz = 0
        self._material_class.determin_direction(nz)


    @property
    def strain_increment(self):
        if self._strain_increment is None:
            return 0
        return self._strain_increment
    @strain_increment.setter
    def strain_increment(self, value):
        self._strain_increment = value

    def calculate_strain_increment_from_section(self, sec_chng_def_increment):
        """ step 10 """
        self.chng_strain_increment = self.direction @ sec_chng_def_increment
        self.strain_increment += self.chng_strain_increment

        self._material_class.update_change_in_material_strain_incr(self.chng_strain_increment)
        self._material_class.update_material_strain_incr()
        self._material_class.update_material_strain()
        if self._first_iteration: #FIXME: make it somethin from inside the material class.. maybe
            self._material_class.update_model_parameters(nz=0)
        self._material_class.update_material_stress()
        self._material_class.update_material_tangent_modulus()

    def increment_strain(self):
        """ step 10 """
        self.strain = self.converged_strain + self.strain_increment

    def calculate_stress(self):
        self.stress = self._material_class.get_material_stress()
