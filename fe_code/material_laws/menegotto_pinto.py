"""
Module contains only MenegottoPinto class
"""

import numpy as np
from .material import Material

class MenegottoPinto(Material):
    """
    Steel uniaxial material law

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, youngs_modulus, asymptotic_modulus, yield_stress, transition, a1, a2):
        self.youngs_modulus = youngs_modulus
        self.asymptotic_modulus = asymptotic_modulus
        self.transition = transition
        self.yield_stress = yield_stress
        self.a1 = a1
        self.a2 = a2

        self.yield_strain = yield_stress / youngs_modulus
        self.b = asymptotic_modulus / youngs_modulus
        self.cosi = 0.0 #TODO: understand this!
        self.transition_0 = transition

        self.strain_r = 0.0
        self.converged_strain_r = 0.0
        self.stress_r = 0.0
        self.converged_stress_r = 0.0
        self.strain_0 = self.yield_strain
        self.stress_0 = self.yield_stress

        self.chng_strain_increment = 0.0
        self.strain_increment = 0.0
        self.strain = 0.0
        self.converged_strain = 0.0

        self.chng_stress_increment = 0.0
        self.stress_increment = 0.0
        self.stress = 0.0
        self.converged_stress = 0.0

        self.direction = 0

    def calculate_strain_increment_from_fiber(self, fiber_chng_strain_increment):
        """
        FIXME
        """
        self.chng_strain_increment = fiber_chng_strain_increment
        self.strain_increment += self.chng_strain_increment
