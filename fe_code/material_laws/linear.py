"""
Module contains only MenegottoPinto class
"""

from .material import UniaxialIncrementalMaterial


class LinearElastic(UniaxialIncrementalMaterial):
    """
    Linear Elastic uniaxial material law
    """

    def __init__(self, E):
        self._E = E

    @property
    def tangent_modulus(self):
        """ current tangent modulus """
        return self._E

    @property
    def stress(self):
        """ current stress """
        return self._stress

    @property
    def strain(self):
        """ current strain """
        return self._strain

    def update_strain(self, fiber_strain):
        """
        set new strain and test for reversal

        Parameters
        ----------
        fiber_strain : float
            fiber strain
        """
        self._set_trial_state(fiber_strain)

    def _set_trial_state(self, new_strain):
        self._strain = new_strain
        self._stress = self._E * self._strain

    def finalize_load_step(self):
        """
        update the converged variables. Called when the
        whole structure is converged at the load step
        """
        self._c_strain = self._strain
        self._c_stress = self._stress
