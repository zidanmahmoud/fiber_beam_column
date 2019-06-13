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
        self.tangnt_modulus = youngs_modulus  # TODO: Do we actually need dis?
        self.asymptotic_modulus = asymptotic_modulus
        self.transition = transition
        self.yield_stress = yield_stress
        self.a1 = a1
        self.a2 = a2

        self.yield_strain = yield_stress / youngs_modulus
        self.b = asymptotic_modulus / youngs_modulus
        self.cosi = 0.0  # TODO: understand this!
        self.transition_0 = transition

        self.strain_r = 0.0
        self.converged_strain_r = 0.0
        self.stress_r = 0.0
        self.converged_stress_r = 0.0
        self.strain_0 = self.yield_strain
        self.stress_0 = self.yield_stress

        self.chng_strain_increment = 0.0
        self.strain_increment = 0.0
        self.converged_strain_increment = 0.0  # TODO: understand why do we need dis
        self.strain = 0.0
        self.converged_strain = 0.0

        self.chng_stress_increment = 0.0
        self.stress_increment = 0.0
        self.stress = 0.0
        self.converged_stress = 0.0

        self.direction = 0

    def calculate_strain_from_fiber(self, fiber_chng_strain_increment):
        """
        FIXME
        """
        self.chng_strain_increment = fiber_chng_strain_increment
        self.strain_increment += self.chng_strain_increment
        self.strain = self.converged_strain + self.strain_increment

    def reverse(self, nz):  # TODO: FIX THIS WHOLE FFUUNNCCTTIIOONN
        """
        FIXME
        """

        E = self.youngs_modulus
        b = self.b

        if self.converged_strain_increment == 0:
            if nz > 4:
                self.strain_0 *= -1
                self.stress_0 *= -1
        else:
            self.converged_strain_r = self.strain_r
            self.converged_stress_r = self.stress_r
            self.strain_r = self.strain
            self.stress_r = self.stress
            self.direction *= -1

            epr = self.strain_r
            sgr = self.stress_r
            csgr = self.converged_stress_r
            cepr = self.converged_strain_r
            di = self.direction
            sgy = self.yield_stress

            self.strain_0 = ((E * epr) - sgr + (di * sgy * (1 - b))) / (E * (1 - b))
            eps = ((sgr - csgr) + E * b * cepr - E * epr) / (E * (1 - b))
            self.cosi = abs(eps - cepr)
            self.transition = self.transition_0 - self.a1 * self.cosi / (self.a2 + self.cosi)

    def calculate_stress_and_tangnt_modulus(self):
        """
        FIXME
        """
        b = self.b
        eps = self.strain
        epr = self.strain_r
        ep0 = self.strain_0
        sgr = self.stress_r
        sg0 = self.stress_0
        eps_star = (eps - epr) / (ep0 - epr)
        if eps_star < 0:  # TODO: fix dis
            ep0 = 1e-10
        dg = (eps_star ** self.transition).real
        sg_star = ((b * eps_star) + (1 - b) * eps_star) / (1 + dg) ** self.transition
        self.stress = sg_star * (sg0 - sgr) + sgr
        self.tangnt_modulus = b + (1 - b) * (
            (1 + dg) ** (-1 / self.transition) - dg * (1 + dg) ** (-1 - 1 / self.transition)
        ) * (sg0 - sgr) / (ep0 - epr)

    def finalize_load_step(self):
        """
        FIXME
        """
        self.converged_strain_increment = self.strain_increment
        self.strain_increment = 0.0
        self.converged_strain = self.strain
        self.converged_stress = self.stress
