class Fiber:
    def __init__(self, y, z, area, material_class):
        self._y = y
        self._z = z
        self._area = area
        self._material_class = material_class

    @property
    def direction(self):
        return [self._y, self._z, 1.0]

    @property
    def area(self):
        return self._area

    @property
    def tangent_stiffness(self):
        return self._material_class.get_material_tangent_modulus()

    def initialize(self):
        nz = 0
        self._material_class.determin_direction(nz)
