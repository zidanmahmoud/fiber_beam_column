"""
Module contains only MenegottoPinto class
"""

from .material import Material

TOLERANCE = 1e-4


class MenegottoPinto(Material):
    """
    Steel uniaxial material law

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, youngs_modulus, asymptotic_modulus, yield_stress, transition, a1, a2):
        self._youngs_modulus = youngs_modulus
        self.tangent_modulus = youngs_modulus
        self._asymptotic_modulus = asymptotic_modulus
        self._transition = transition
        self._yield_stress = yield_stress
        self._a1 = a1
        self._a2 = a2

        self._yield_strain = yield_stress / youngs_modulus
        self._b = asymptotic_modulus / youngs_modulus
        self._xi = 0.0
        self._transition_0 = transition

        self._strain_r = 0.0
        self._last_strain_r = 0.0
        self._stress_r = 0.0
        self._last_stress_r = 0.0
        self._strain_0 = self._yield_strain
        self._stress_0 = self._yield_stress

        self._chng_strain_increment = 0.0
        self._strain_increment = 0.0
        self._strain = 0.0
        self._converged_strain = 0.0
        self._last_converged_strain = 0.0

        self.stress = 0.0

        self._direction = 0

    def calculate_strain_from_fiber(self, fiber_chng_strain_increment):
        """
        FIXME
        """
        self._chng_strain_increment = fiber_chng_strain_increment
        self._strain_increment += self._chng_strain_increment
        self._strain = self._converged_strain + self._strain_increment

        if self._check_reversal():
            self._reverse()


    def _check_reversal(self):
        # if strain was increasing
        if self._converged_strain - self._last_converged_strain > 0:
            # if the new strain is decreasing
            if self._strain < self._converged_strain:
                return True

        # else if strain was decreasing
        if self._converged_strain - self._last_converged_strain < 0:
            # if the new strain is increasing
            if self._strain > self._converged_strain:
                return True

        return False

    def _reverse(self):
        """
        FIXME
        """

        E = self._youngs_modulus
        b = self._b
        self._last_strain_r = self._strain_r
        self._last_stress_r = self._stress_r
        self._strain_r = self._strain
        self._stress_r = self.stress
        self._direction *= -1

        epr = self._strain_r
        sgr = self._stress_r
        sgy = self._direction * self._yield_stress
        lepr = self._last_strain_r

        self._strain_0 = (E * epr - sgr + sgy * (1 - b)) / (E * (1 - b))
        self._stress_0 = b * E * self._strain_0 + sgy * (1 - b)
        eps_intersect = ((sgr - lepr) + E * b * lepr - E * epr) / (E * (b - 1))
        self._xi = abs(eps_intersect - lepr)
        self._transition = self._transition_0 - self._a1 * self._xi / (self._a2 + self._xi)
        self._xi = 0.0

    def calculate_stress_and_tangent_modulus(self):
        """
        FIXME
        """
        b = self._b
        eps = self._strain
        epr = self._strain_r
        ep0 = self._strain_0
        sgr = self._stress_r
        sg0 = self._stress_0
        R = self._transition

        eps_star = (eps - epr) / (ep0 - epr)
        if eps_star < 0:  # TODO: fix dis
            eps_star = 1e-10
        dg = (eps_star ** R).real
        sg_star = b * eps_star + (1 - b) * eps_star / (1 + dg) ** (1 / R)
        self.stress = sg_star * (sg0 - sgr) + sgr
        self.tangent_modulus = b + (1 - b) * ((1 + dg) ** (-1 / R) - dg * (1 + dg) ** (-1 - 1 / R))
        self.tangent_modulus *= (sg0 - sgr) / (ep0 - epr)

    def finalize_load_step(self):
        """
        FIXME
        """
        self._strain_increment = 0.0
        self._last_converged_strain = self._converged_strain
        self._converged_strain = self._strain

    def determin_direction(self, nz):
        """
        FIXME
        """
        self._direction = 1
