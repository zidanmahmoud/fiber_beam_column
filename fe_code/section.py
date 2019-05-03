import numpy as np

from .fiber import Fiber

class Section:
    def __init__(self):
        self._fibers = dict()
        self._tolerance = 1e-7

    @property
    def tolerance(self):
        return self._tolerance
    @tolerance.setter
    def tolerance(self, value):
        self._tolerance = value

    @property
    def fibers(self):
        return self._fibers.values()

    def add_fiber(self, fiber_id, y, z, area, material_class):
        self._fibers[fiber_id] = Fiber(y, z, area, material_class)

    def initialize(self):
        for fiber in self.fibers:
            fiber.initialize()
        # stiffness matrix

    def calculate_stiffness_matrix(self):
        stiffness_matrix = np.zeros((3, 3))
        for fiber in self.fibers:
            EA = fiber.tangent_stiffness * fiber.area
            stiffness_matrix += EA * np.outer(fiber.direction, fiber.direction)
        return stiffness_matrix
