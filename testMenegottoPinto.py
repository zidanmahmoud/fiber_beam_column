"""
test file
"""

import numpy as np
import matplotlib.pyplot as plt


class MenegottoPinto:
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

    def update_strain(self, value):
        """
        FIXME
        """
        self._strain = value
        reversal = self._set_trial_state()
        return reversal

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
        reversal = self._check_reversal()
        if reversal:
            self._reverse()
        return reversal

    def _check_reversal(self):
        deps = self._strain - self._c_strain
        if self._loading_index == 2 and deps > 0:
            self._loading_index = 1
            return True
        if self._loading_index == 1 and deps < 0:
            self._loading_index = 2
            return True
        return False

    def _reverse(self):
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

        # global ax
        # ax.plot(self._strain_0, self._stress_0, "-o", color="black")
        # ax.plot(self._strain_r, self._stress_r, "-o", color="black", markerfacecolor="none")

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
        dum1 = 1.0 + (abs(eps_star)) ** R
        dum2 = (dum1) ** (1.0 / R)
        sg_star = b * eps_star + (1.0 - b) * eps_star / dum2
        self._stress = sg_star * (sg0 - sgr) + sgr
        self._Et = b + (1 - b) / (dum1 * dum2)
        self._Et *= (sg0 - sgr) / (ep0 - epr)

    def finalize(self):
        self._c_loading_index = self._loading_index
        self._c_Et = self._Et
        self._c_strain_0 = self._strain_0
        self._c_stress_0 = self._stress_0
        self._c_strain_r = self._strain_r
        self._c_stress_r = self._stress_r
        self._c_strain = self._strain
        self._c_stress = self._stress
        self._c_xi = self._xi


fiber = MenegottoPinto(E=29000, b=0.08, fy=60, R=15, a1=8.5, a2=0.0002)  # 0.0042  # 20


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(fiber._strain_0, fiber._stress_0, "-o", color="black")

# strains = np.concatenate(
#     (
#         np.linspace(0.000, 0.005),
#         np.linspace(0.005, -0.01),
#         np.linspace(-0.01, -0.0005),
#     )
# )


strains = np.array([0, 0.001, 0.002, 0.0035, 0.0015, 0.0031, 0.003, 0.004, 0.005])
f = [0, 1, 2, 6, 7, 8]
nf = [3, 4, 5]
stresses = []
for i, strain in enumerate(strains):
    reversal = fiber.update_strain(strain)
    fiber.calculate_stress_and_tangent_modulus()
    if i in f:
        fiber.finalize()
    stresses.append(fiber.stress)
    if reversal:
        print(fiber._stress_0)
stresses = np.array(stresses)
ax.plot(strains[f], stresses[f], "-o", color="black", markerfacecolor="none")
ax.plot(strains[nf], stresses[nf], "o", color="orange")
ax.grid()
ax.axhline(linewidth=3, color="black")
ax.axvline(linewidth=3, color="black")
ax.set(xlabel="STEEL STRAIN", ylabel="STEEL STRESS")
plt.show()
