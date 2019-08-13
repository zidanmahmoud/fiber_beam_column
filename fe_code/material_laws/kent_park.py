"""
Module contains only KentPark class
"""

from .material import Material


class KentPark(Material):
    """
    Concrete uniaxial material law

    Parameters
    ----------
    fc : float
        yield compressive strength
    K : float
        confinement factor
    Z : float
        softening slope
    e0 : float
        strain at maximum stress. default 0.002
    eu : float, optional
        ultimate strain (crushing)

    Note :: can be construced either using
    KentPark(fc, K, Z, e0) or
    KentPark(fc, K, eu, e0)

    Attributes
    ----------
    tangent_modulus : float
    stress : float
    strain : float
    """

    def __init__(self, fc, K, Z, e0=0.002):
        self._fc = fc
        self._K = K
        self._strain_0 = -e0 * K
        self._Z = Z
        self._Et = 2 * fc / self._strain_0

        # Loading index:
        #   0 => initial state
        #   1 => increasing strain
        #   2 => decreasing strain
        #   3 => strain not changing
        self._loading_index = 0
        self._strain_p = 0.0
        self._strain_r = 0.0
        self._stress_r = 0.0
        self._strain = 0.0
        self._stress = 0.0

        # == converged variables
        self._c_loading_index = 0
        self._c_strain_p = 0.0
        self._c_strain_r = 0.0
        self._c_stress_r = 0.0
        self._c_strain = 0.0
        self._c_stress = 0.0

    @classmethod
    def eu(cls, fc, K, eu, e0=0.002):
        """
        Overloading the default constructor with eu
        instead of Z
        """
        Z = 0.8 / (eu - e0)
        return cls(fc, K, Z, e0)

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
        self._strain = fiber_strain
        return self._set_trial_state()

    def _set_trial_state(self):
        deps = self._strain - self._c_strain
        if abs(deps) < 1e-15:  # nearly zero
            self._loading_index = 3
        else:
            if deps < 0:
                self._loading_index = 2
            else:
                self._loading_index = 1
        reversal = self.check_reversal()
        if reversal:
            self.reverse()
        return reversal

    def check_reversal(self):
        """
        check for reversal
        """
        if abs(self._strain) > 1e-15:
            if self._strain < 0:
                if self._c_loading_index in (2, 3):
                    if self._loading_index == 1:
                        return True
        return False

    def reverse(self):
        """
        reverse the material parameters
        """
        self._strain_r = self._c_strain
        self._stress_r = self._c_stress

        crit = self._strain_r / self._strain_0
        ep0 = self._strain_0
        if crit >= 2:
            self._strain_p = ep0 * (0.707 * (crit - 2) + 0.834)
        elif crit < 2:
            self._strain_p = ep0 * (0.145 * crit ** 2 + 0.13 * crit)

    def calculate_stress_and_tangent_modulus(self):
        """
        update the current stress and tangent modulus
        based on the current strain
        """
        eps = self._strain
        ep0 = self._strain_0
        epp = self._strain_p
        epr = self._strain_r
        sgr = self._stress_r
        K = self._K
        Z = self._Z
        fc = self._fc

        # == inequality signs are reversed compared to theory becuase of the negative signs

        # positive strain
        if eps >= 0:
            self._stress = 0.0
            self._Et = 0.0
            return

        # loading path
        if eps <= epr:
            if eps >= ep0:
                stress = K * fc * (2 * eps / ep0 - (eps / ep0) ** 2)
                tangen = K * fc * (2 / ep0 - 2 * (eps / ep0 ** 2))
            else:
                stress = K * fc * (1 + Z * (eps - ep0))
                if stress < 0.2 * K * fc:
                    stress = 0.2 * K * fc
                    tangen = 0
                else:
                    tangen = K * fc * Z

        # unloading path
        else:
            if eps >= epp:
                self._stress = 0.0
                self._Et = 0.0
                return
            stress = -(sgr * eps - epp * sgr) / (epr - epp)
            tangen = -sgr / (epr - epp)

        self._stress = -1 * stress
        self._Et = -1 * tangen

    def finalize_load_step(self):
        """
        update the converged variables. Called when the
        whole structure is converged at the load step
        """
        self._c_loading_index = self._loading_index
        self._c_strain_p = self._strain_p
        self._c_strain_r = self._strain_r
        self._c_stress_r = self._stress_r
        self._c_strain = self._strain
        self._c_stress = self._stress
