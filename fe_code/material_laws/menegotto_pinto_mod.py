"""
Module contains only MenegottoPinto class
"""

from .material import Material


class MenegottoPintoMod:
    """
    Steel uniaxial material law

    Parameters
    ----------
    E : float
        Young's Modulus
    b : flaot
        hardening ratio
    fy : float
        yield strength
    R : float
        initial transition variable
    a1 : float
        imperical parameter
    a2 : float
        imperical parameter


    Attributes
    ----------
    tangent_modulus : float
    stress : float
    strain : float
    """

    def __init__(self, E, b, fy, R0, a1, a2):
        self._E = E
        self._b = b
        self._R0 = R0
        self._fy = fy
        self._a1 = a1
        self._a2 = a2
        self._stress_initial = 0.0

        # Loading index:
        #   0 => initial state
        #   1 => increasing strain
        #   2 => decreasing strain
        #   3 => strain not changing
        self._loading_index = 0
        self._Et = E
        self._strain_0 = 0.0
        self._stress_0 = 0.0
        self._strain_r = 0.0
        self._stress_r = 0.0
        self._strain_plastic = 0.0
        self._strain_max = fy / E
        self._strain_min = -self._strain_max
        self._strain = 0.0
        self._stress = 0.0

        # Converged Variables
        self._c_loading_index = 0
        self._c_Et = E
        self._c_strain_0 = 0.0
        self._c_stress_0 = 0.0
        self._c_strain_r = 0.0
        self._c_stress_r = 0.0
        self._c_strain_plastic = 0.0
        self._c_strain_max = fy / E
        self._c_strain_min = -self._c_strain_max
        self._c_strain = 0.0
        self._c_stress = 0.0

    @property
    def tangent_modulus(self):
        """ current tangent modulus """
        return self._Et

    @property
    def stress(self):
        """ current stress """
        return self._stress

    @property
    def strain(self):
        """ current strain """
        return self._strain

    def prestress(self, value):
        self._stress_initial = value
        self._c_strain = value / self._E
        self._c_stress = value

    def update_strain(self, fiber_strain):
        """
        set new strain and test for reversal

        Parameters
        ----------
        fiber_strain : float
            fiber strain

        Returns
        -------
        reversal : bool
            flag True if reversed
        """
        reversal = self._set_trial_state(fiber_strain)
        return reversal

    def _set_trial_state(self, new_strain):
        fy = self._fy
        E = self._E
        b = self._b
        E_inf = b * E
        strain_y = fy / E

        if self._stress_initial != 0:
            self._strain = self._stress_initial / E + new_strain
        else:
            self._strain = new_strain

        deps = self._strain - self._c_strain

        self._strain_max = self._c_strain_max
        self._strain_min = self._c_strain_min
        self._strain_plastic = self._c_strain_plastic
        self._strain_0 = self._c_strain_0
        self._stress_0 = self._c_stress_0
        self._strain_r = self._c_strain_r
        self._stress_r = self._c_stress_r
        self._loading_index = self._c_loading_index

        if self._loading_index == 0 or self._loading_index == 3:

            if abs(deps) < 1e-15:  # nearly zero
                self._Et = E
                self._stress = 0.0
                self._loading_index = 3
                return False

            else:
                self._strain_max = strain_y
                self._strain_min = -strain_y

                if deps < 0:
                    self._loading_index = 2
                    self._strain_0 = self._strain_min
                    self._stress_0 = -fy
                    self._strain_plastic = self._strain_min

                else:
                    self._loading_index = 1
                    self._strain_0 = self._strain_max
                    self._stress_0 = fy
                    self._strain_plastic = self._strain_max

        reversal = False

        if self._loading_index == 2 and deps > 0:
            reversal = True
            self._loading_index = 1
            self._strain_r = self._c_strain
            self._stress_r = self._c_stress
            if self._c_strain < self._strain_min:
                self._strain_min = self._c_strain
            self._strain_0 = (fy - E_inf * strain_y - self._stress_r + E * self._strain_r) / (
                E - E_inf
            )
            self._stress_0 = fy + E_inf * (self._strain_0 - strain_y)
            self._strain_plastic = self._strain_max

        elif self._loading_index == 1 and deps < 0:
            reversal = True
            self._loading_index = 2
            self._strain_r = self._c_strain
            self._stress_r = self._c_stress
            if self._c_strain > self._strain_max:
                self._strain_max = self._c_strain
            self._strain_0 = (-fy + E_inf * strain_y - self._stress_r + E * self._strain_r) / (
                E - E_inf
            )
            self._stress_0 = -fy + E_inf * (self._strain_0 + strain_y)
            self._strain_plastic = self._strain_min

        xi = abs((self._strain_plastic - self._strain_0) / strain_y)
        R = self._R0 - self._a1 * xi / (self._a2 + xi)
        eps_star = (self._strain - self._strain_r) / (self._strain_0 - self._strain_r)
        dum1 = 1.0 + (abs(eps_star)) ** R
        dum2 = (dum1) ** (1.0 / R)
        sg_star = b * eps_star + (1.0 - b) * eps_star / dum2
        self._stress = sg_star * (self._stress_0 - self._stress_r) + self._stress_r
        self._Et = b + (1.0 - b) / (dum1 * dum2)
        self._Et *= (self._stress_0 - self._stress_r) / (self._strain_0 - self._strain_r)

        return reversal

    def finalize_load_step(self):
        """
        update the converged variables. Called when the
        whole structure is converged at the load step
        """
        self._c_loading_index = self._loading_index
        self._c_Et = self._Et
        self._c_strain_0 = self._strain_0
        self._c_stress_0 = self._stress_0
        self._c_strain_r = self._strain_r
        self._c_stress_r = self._stress_r
        self._c_strain_plastic = self._strain_plastic
        self._c_strain_max = self._strain_max
        self._c_strain_min = self._strain_min
        self._c_strain = self._strain
        self._c_stress = self._stress
