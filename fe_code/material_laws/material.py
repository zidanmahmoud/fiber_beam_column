"""
Module contains the class Material
"""

from abc import ABC, abstractclassmethod

class Material(ABC):
    """
    Material abstract class to be used in fiber-beam-column element
    """
    @abstractclassmethod
    def calculate_strain_from_fiber(self, fiber_chng_strain_increment):
        pass

    @abstractclassmethod
    def reverse(self, nz):
        pass

    @abstractclassmethod
    def calculate_stress_and_tangent_modulus(self):
        pass

    @abstractclassmethod
    def finalize_load_step(self):
        pass
