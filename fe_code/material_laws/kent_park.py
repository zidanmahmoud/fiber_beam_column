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
        self._chng_strain_increment = fiber_chng_strain_increment
        self._strain_increment += self._chng_strain_increment
        self._strain = self._converged_strain + self._strain_increment

    def reverse(self, nz):
        if self._strain < 0:
            self._strain_r = self._converged_strain
            self._stress_r = self._converged_stress

            critical_point = self._strain_r / self._strain_0
            ep0 = self._strain_0
            if critical_point >= 2:
                self._strain_p = ep0 * (0.707 * (critical_point - 2) + 0.834)
            elif critical_point > 0:
                self._strain_p = ep0 * (0.145 * critical_point**2 + 0.13 * critical_point)

    def calculate_stress_and_tangent_modulus(self):
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
        self.save_to_last_loadstep_material_strain_incr(True)
        self.save_to_material_strain_last_loadstep(True)
        self.initialize_material_strain_incr()

    def update_change_in_material_strain_incr(self, change_in_fiber_strain_incr):

        self._chng_strain_increment = change_in_fiber_strain_incr

    def save_to_last_loadstep_material_strain_incr(self, load_step_convergence):

        if load_step_convergence:
            self._converged_strain_increment = self._strain_increment

    def get_last_loadstep_material_strain_incr(self):

        return self._converged_strain_increment

    def determin_direction(self, nz):

        if nz > 4:
            self._stress = -1
        elif nz <= 4:
            self._stress = 1

    """
    material_strain_incr :
         If i == 1 and j == 1: it is initialized to be zero
         Else: it is updated from change_in_fiber_strain_incr

         Function needed:
         change_in_fiber_strain_incr
    """

    def get_material_strain_incr(self):

        return self._strain_increment

    def initialize_material_strain_incr(self):

        self._strain_increment = 0.0

    def update_material_strain_incr(self):

        self._strain_increment += self._chng_strain_increment

    """
    material_strain :

         It is updated from material_strain_incr
    """

    def get_material_strain(self):

        return self._strain

    def initialize_material_strain(self):

        self._strain = 0.0

    def update_material_strain(self):

        self._strain = self._converged_strain + self._strain_increment

    def save_to_material_strain_last_loadstep(self, load_step_convergence):

        if load_step_convergence:
            self._converged_strain = self._strain
            self._converged_stress = self.stress

    """
    material_stress :

         It is calculated from
         material_strain

         Function needed:
         calculate_material_stress()
    """

    def get_material_stress(self):

        return self.stress

    def initialize_material_stress(self):

        self.stress = 0.0

    def update_material_stress(self):

        self.stress = self.calculate_material_stress()

    """
    material_tangent_modulus :

         Function needed:
         calculate_material_tangent_modulus()
    """

    def get_material_tangent_modulus(self):

        return self.tangent_modulus

    def initialize_material_tangent_modulus(self):

        self.tangent_modulus = (
            -2 * self._compressive_strength * self._confinement_factor / self._strain_0
        )

    def update_material_tangent_modulus(self):

        self.tangent_modulus = self.calculate_material_tangent_modulus()

    """
    self._strain_p :

         Function needed:
         calculate_material_cosi()
    """

    def get_strain_p(self):

        return self._strain_p

    def initialize_strain_p(self):

        self._strain_p = 0.0

    def update_strain_p(self):

        critical_point = self._strain_r / self._strain_0
        if critical_point < 2 and critical_point > 0:
            self._strain_p = self._strain_0 * (0.145 * critical_point ** 2 + 0.13 * critical_point)
        elif critical_point >= 2:
            self._strain_p = self._strain_0 * (0.707 * (critical_point - 2) + 0.834)

    """
    Update model parameters:

         ep0
         sg0
         self._strain_r
         sgr
         R
         cosi

    """

    def update_model_parameters(self, nz):

        if self._strain < 0:
            self._strain_r = self._converged_strain
            self._stress_r = self._converged_stress
            self.update_strain_p()
        print_ = False
        if print_:
            print("self._strain_r = ", self._strain_r)
            print("sgr = ", self._stress_r)
            print("self._strain_p = ", self._strain_p)
            print("last_loadstep_material_strain_incr = ", self._converged_strain_increment)
            print("material_strain_incr = ", self._strain_increment)
            print("material_strain = ", self._strain)

    """
    other:
    """

    def get_material_PEEQ(self):

        return self._strain_p

    ################################################################################
    ################################################################################
    """
    The second block provides :
    functions which are needed in the first block.
    """
    ################################################################################
    def check_reversal(self, do_reversal, j, nx, ny, nz):
        print_ = False

        if do_reversal and (j == 1):
            reversal = True
        else:
            reversal = False

        return reversal

    def check_load(self):

        loading = True
        loading *= self._strain_increment < 0
        return loading

    def calculate_material_stress(self):

        eps = self._strain
        K = self._confinement_factor
        fc = self._compressive_strength
        ep0 = self._strain_0
        Z = self._Z
        self._strain_r = self._strain_r
        sgr = self._stress_r
        self._strain_p = self._strain_p
        epu = self._strain_ultimate
        if self._strain_increment < 0:
            if eps > self._strain_p:
                sg = 0
            elif eps > self._strain_r and eps < self._strain_p:
                sg = sgr / (self._strain_r - self._strain_p) * (eps - self._strain_r) + sgr
            elif eps <= self._strain_r:
                if eps > 0:
                    sg = 0
                elif eps > ep0:
                    sg = K * fc * (-2 * (eps / ep0) + (eps / ep0) ** 2)
                elif eps > epu:
                    sg = max(K * fc * (-1 - Z * (eps - ep0)), -0.2 * K * fc)
                else:
                    sg = 0
        else:
            if eps <= self._strain_r:
                if eps > 0:
                    sg = 0
                elif eps > ep0:
                    sg = K * fc * (-2 * (eps / ep0) + (eps / ep0) ** 2)
                elif eps > epu:
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    sg = sg * (sg < -0.2 * K * fc) + -0.2 * K * fc * (sg >= -0.2 * K * fc)
                else:
                    sg = 0

            elif eps < self._strain_p:
                sg = sgr / (self._strain_r - self._strain_p) * (eps - self._strain_r) + sgr
            else:
                sg = 0

        return sg

    def calculate_material_tangent_modulus(self):

        eps = self._strain
        K = self._confinement_factor
        fc = self._compressive_strength
        ep0 = self._strain_0
        Z = self._Z
        self._strain_r = self._strain_r
        sgr = self._stress_r
        self._strain_p = self._strain_p
        epu = self._strain_ultimate
        if self.check_load():
            if eps > self._strain_p:
                Et = 0
            elif eps > self._strain_r and eps < self._strain_p:
                Et = sgr / (self._strain_r - self._strain_p)
            elif eps <= self._strain_r:
                if eps > 0:
                    Et = 0
                elif eps > ep0:
                    Et = K * fc * (-(2 / ep0) + (eps / ep0) * 2 / ep0)
                elif eps > epu:
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    Et = -K * fc * Z * (sg < -0.2 * K * fc)
                else:
                    Et = 0
        else:
            if eps <= self._strain_r:
                # print('warning: strain exceeds reverse point.')
                if eps > 0:
                    Et = 0
                elif eps > ep0:
                    Et = K * fc * (-(2 / ep0) + (eps / ep0) * 2 / ep0)
                elif eps > epu:
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    Et = -K * fc * Z * (sg < -0.2 * K * fc)
                else:
                    Et = 0
            elif eps < self._strain_p:
                Et = sgr / (self._strain_r - self._strain_p)
            else:
                Et = 0

        return Et
