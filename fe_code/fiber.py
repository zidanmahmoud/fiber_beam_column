class Fiber:
    def __init__(self, y, z, area, material_class):
        self._y = y
        self._z = z
        self._area = area
        self._material_class = material_class

    def initialize(self):
        nz = 0
        self._material_class.determin_direction(nz)
