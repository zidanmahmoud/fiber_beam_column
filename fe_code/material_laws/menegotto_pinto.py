"""
Module contains only MenegottoPinto class
"""

from .material import Material

TOLERANCE = 0


class MenegottoPinto(Material):
    """
    Steel uniaxial material law

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, E, b, fy, R, a1, a2):
        self._E = E
        self._Et = E
        self._b = b
        self._R0 = R
        self._R = R
        self._fy = fy
        self._a1 = a1
        self._a2 = a2

        # Loading index:
        #   0 => initial state
        #   1 => increasing strain
        #   2 => decreasing strain
        #   3 => strain not changing
        self._loading_index = 0
        self._strain_0 = fy / E
        self._stress_0 = fy
        self._strain_r = 0.0
        self._stress_r = 0.0
        self._last_strain_r = 0.0
        self._last_stress_r = 0.0
        self._strain = 0.0
        self._stress = 0.0
        self._xi = 0.0

        # Converged Variables
        self._c_loading_index = 0
        self._c_Et = E
        self._c_strain_0 = self._strain_0
        self._c_stress_0 = self._stress_0
        self._c_strain_r = 0.0
        self._c_stress_r = 0.0
        self._c_strain = 0.0
        self._c_stress = 0.0
        self._c_xi = 0.0

    @property
    def tangent_modulus(self):
        return self._Et

    @property
    def stress(self):
        return self._stress

    @property
    def strain(self):
        return self._strain

    def update_strain(self, fiber_strain):
        """
        FIXME
        """
        self._strain = fiber_strain
        self._set_trial_state()

    def _set_trial_state(self):
        deps = self._strain - self._c_strain
        if self._loading_index == 0 or self._loading_index == 3:
            if abs(deps) < 1e-15:  # nearly zero
                self._Et = self._E
                self._stress = 0
                self._loading_index = 3
            else:
                if deps < 0:
                    self._loading_index = 2
                    self._strain_0 = -self._fy / self._E
                    self._stress_0 = -self._fy
                else:
                    self._loading_index = 1
                    self._strain_0 = self._fy / self._E
                    self._stress_0 = self._fy
        reversal = self.check_reversal()
        if reversal:
            self.reverse()

    def check_reversal(self):
        deps = self._strain - self._c_strain
        if self._loading_index == 2 and deps > 0:
            self._loading_index = 1
            return True
        if self._loading_index == 1 and deps < 0:
            self._loading_index = 2
            return True
        return False

    def reverse(self):
        """
        FIXME
        """
        self._last_strain_r = self._strain_r
        self._last_stress_r = self._stress_r
        self._strain_r = self._c_strain
        self._stress_r = self._c_stress
        E = self._E
        b = self._b
        epr = self._strain_r
        sgr = self._stress_r
        if self._loading_index == 1:
            sgy = self._fy
        else:
            sgy = -self._fy
        lepr = self._last_strain_r
        self._strain_0 = (E * epr - sgr + sgy * (1 - b)) / (E * (1 - b))
        self._stress_0 = b * E * self._strain_0 + sgy * (1 - b)
        eps_intersect = ((sgr - lepr) + E * b * lepr - E * epr) / (E * (b - 1))
        self._xi = abs(eps_intersect - lepr)
        self._R = self._R0 - self._a1 * self._xi / (self._a2 + self._xi)

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
        R = self._R

        eps_star = (eps - epr) / (ep0 - epr)
        dum1 = 1.0 + (abs(eps_star))**R
        dum2 = (dum1)**(1.0/R)
        sg_star = b * eps_star + (1.0 - b) * eps_star / dum2
        self._stress = sg_star * (sg0 - sgr) + sgr
        self._Et = b + (1 - b) / (dum1 * dum2)
        self._Et *= (sg0 - sgr) / (ep0 - epr)

    def finalize_load_step(self):
        """
        FIXME
        """
        self._c_loading_index = self._loading_index
        self._c_Et = self._Et
        self._c_strain_0 = self._strain_0
        self._c_stress_0 = self._stress_0
        self._c_strain_r = self._strain_r
        self._c_stress_r = self._stress_r
        self._c_strain = self._strain
        self._c_stress = self._stress
        self._c_xi = self._xi

    def determin_direction(self, nz):
        """
        FIXME
        """
        pass
