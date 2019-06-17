"""
Module contains only KentPark class
"""

import numpy as np
from .material import Material


class KentPark(Material):
    """
    Steel uniaxial material law

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, compressive_strength, confinement_factor, epu, Z):
        self._compressive_strength = compressive_strength
        self._confinement_factor = confinement_factor
        self._strain_0 = -0.0027 * confinement_factor
        self._strain_ultimate = epu
        self._Z = Z

        self._strain_r = 0.0
        self._stress_r = 0.0
        self._strain_p = 0.0

        self._chng_strain_increment = 0.0
        self._converged_strain_increment = 0.0
        self._converged_strain = 0.0
        self._converged_stress = 0.0
        self._strain_increment = 0.0
        self._strain = 0.0
        self.stress = 0.0
        self.tangent_modulus = (
            -2 * self._compressive_strength * self._confinement_factor / self._strain_0
        )
        self._stress = 0.0

    def calculate_strain_from_fiber(self, fiber_chng_strain_increment):
        """
        FIXME
        """
        self._chng_strain_increment = fiber_chng_strain_increment
        self._strain_increment += self._chng_strain_increment
        self._strain = self._converged_strain + self._strain_increment

    def reverse(self, nz):
        """
        FIXME
        """
        if self._strain < 0:
            self._strain_r = self._converged_strain
            self._stress_r = self._converged_stress

            critical_point = self._strain_r / self._strain_0
            ep0 = self._strain_0
            if critical_point >= 2:
                self._strain_p = ep0 * (0.707 * (critical_point - 2) + 0.834)
            elif critical_point > 0:
                self._strain_p = ep0 * (0.145 * critical_point ** 2 + 0.13 * critical_point)

    def calculate_stress_and_tangent_modulus(self):
        """
        FIXME
        """
        eps = self._strain
        ep0 = self._strain_0
        epp = self._strain_p
        epr = self._strain_r
        epu = self._strain_ultimate
        sgr = self._stress_r
        K = self._confinement_factor
        Z = self._Z
        fc = self._compressive_strength

        if self._strain_increment < 0:
            if eps >= epp:
                sg = 0.0
                Et = 0.0
            elif epr <= eps < epp:
                sg = sgr / (epr - epp) * (eps - epr) + sgr
                Et = sgr / (epr - epp)
            elif eps < epr:
                if eps > 0:
                    sg = 0.0
                    Et = 0.0
                elif eps > ep0:
                    sg = K * fc * (-2 * (eps / ep0) + (eps / ep0) ** 2)
                    Et = K * fc * (-(2 / ep0) + (eps / ep0) * 2 / ep0)
                elif eps > epu:
                    # sg = max(K * fc * (-1 - Z * (eps - ep0)), -0.2 * K * fc)
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    if sg < -0.2 * K * fc:
                        Et = -K * fc * Z
                    else:
                        sg = -0.2 * K * fc
                        Et = 0.0
                else:
                    sg = 0.0
                    Et = 0.0
        else:
            if eps <= epr:
                if eps > 0:
                    sg = 0.0
                    Et = 0.0
                elif eps > ep0:
                    sg = K * fc * (-2 * (eps / ep0) + (eps / ep0) ** 2)
                    Et = K * fc * (-(2 / ep0) + (eps / ep0) * 2 / ep0)
                elif eps > epu:
                    # sg = max(K * fc * (-1 - Z * (eps - ep0)), -0.2 * K * fc)
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    if sg < -0.2 * K * fc:
                        Et = -K * fc * Z
                    else:
                        sg = -0.2 * K * fc
                        Et = 0.0
                else:
                    sg = 0.0
                    Et = 0.0

            elif eps < epp:
                sg = sgr / (epr - epp) * (eps - epr) + sgr
                Et = sgr / (epr - epp)
            else:
                sg = 0.0
                Et = 0.0

        self.stress = sg
        self.tangent_modulus = Et

    def finalize_load_step(self):
        """
        FIXME
        """
        self._converged_strain_increment = self._strain_increment
        self._converged_strain = self._strain
        self._converged_stress = self.stress
        self._strain_increment = 0.0

    def determin_direction(self, nz):
        """
        FIXME
        """
        if nz > 4:
            self._stress = -1
        elif nz <= 4:
            self._stress = 1
