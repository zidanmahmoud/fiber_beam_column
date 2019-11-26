"""
test file
"""

import numpy as np
import matplotlib.pyplot as plt


class MenegottoPinto:
    def __init__(self, ax, E, b, fy, R0, a1, a2):
        self.ax=ax

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
        return self._Et

    @property
    def stress(self):
        return self._stress

    @property
    def strain(self):
        return self._strain

    def prestress(self, value):
        self._stress_initial = value
        self._c_strain = value / self._E
        self._c_stress = value

    def update_strain(self, value):
        """
        FIXME
        """
        reversal, eps_star, sg_star = self._set_trial_state(value)
        return reversal, eps_star, sg_star

    def _set_trial_state(self, new_strain):
        E_inf = self._b * self._E
        strain_y = self._fy / self._E

        if self._stress_initial != 0:
            self._strain = self._stress_initial / self._E + new_strain
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
                self._Et = self._E
                self._stress = 0.0
                self._loading_index = 3
                return False, 0, 0

            else:
                self._strain_max = strain_y
                self._strain_min = -strain_y

                if deps < 0:
                    self._loading_index = 2
                    self._strain_0 = self._strain_min
                    self._stress_0 = -self._fy
                    self._strain_plastic = self._strain_min
                    # self.update_plot()

                else:
                    self._loading_index = 1
                    self._strain_0 = self._strain_max
                    self._stress_0 = self._fy
                    self._strain_plastic = self._strain_max
                    # self.update_plot()

        reversal = False

        if self._loading_index == 2 and deps > 0:
            reversal = True
            self._loading_index = 1
            self._strain_r = self._c_strain
            self._stress_r = self._c_stress
            if self._c_strain < self._strain_min:
                self._strain_min = self._c_strain
            self._strain_0 = (
                self._E * self._strain_r - self._stress_r + self._fy * (1 - self._b)
            ) / (self._E * (1 - self._b))
            self._stress_0 = self._b * self._E * self._strain_0 + self._fy * (1 - self._b)
            self._strain_plastic = self._strain_max

        elif self._loading_index == 1 and deps < 0:
            reversal = True
            self._loading_index = 2
            self._strain_r = self._c_strain
            self._stress_r = self._c_stress
            if self._c_strain > self._strain_max:
                self._strain_max = self._c_strain
            self._strain_0 = (
                self._E * self._strain_r - self._stress_r - self._fy * (1 - self._b)
            ) / (self._E * (1 - self._b))
            self._stress_0 = self._b * self._E * self._strain_0 - self._fy * (1 - self._b)
            self._strain_plastic = self._strain_min

        xi = abs((self._strain_plastic - self._strain_0) / strain_y)
        R = self._R0 - self._a1 * xi / (self._a2 + xi)
        eps_star = (self._strain - self._strain_r) / (self._strain_0 - self._strain_r)
        if reversal:
            print(f"xi : {xi}")
            print(f"eps* : {eps_star}")
            print(f"R : {R}")
            print()
        # print(f"pl* : {self._strain_plastic/strain_y}")
        dum1 = 1.0 + (abs(eps_star)) ** R
        dum2 = (dum1) ** (1.0 / R)
        sg_star = self._b * eps_star + (1.0 - self._b) * eps_star / dum2
        self._stress = sg_star * (self._stress_0 - self._stress_r) + self._stress_r
        self._Et = self._b + (1 - self._b) / (dum1 * dum2)
        self._Et *= (self._stress_0 - self._stress_r) / (self._strain_0 - self._strain_r)

        # if reversal:
        #     self.update_plot()

        return reversal, self._strain/strain_y, self._stress/self._fy

    def update_plot(self):
        self.ax.plot(self._strain_r, self._stress_r, "o", color="blue")
        self.ax.plot(self._strain_0, self._stress_0, "o", color="red")

    def finalize(self):
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


fig = plt.figure()
ax = fig.add_subplot(111)

fiber = MenegottoPinto(ax, E=29000, b=0.08, fy=60, R0=15, a1=8.5, a2=0.0002)  # 0.0042  # 20
# fiber = MenegottoPinto(ax, E=29000, b=0.0042, fy=60, R0=20, a1=18.5, a2=0.0002)

# ax.plot(fiber._strain_0, fiber._stress_0, "-o", color="black")
# fiber.update_plot()

strains = np.concatenate(
    (
        np.linspace(0.000, 0.005),
        np.linspace(0.005, -0.01),
        np.linspace(-0.01, -0.0005)
    )
)
# f = np.sort(np.random.randint(0, strains.size, strains.size // 2))
f = np.linspace(0, strains.size - 1, strains.size, dtype=int)
nf = []


# strains = np.array([0, 0.001, 0.002, 0.0035, 0.0015, 0.0031, 0.003, 0.004, 0.005])
# f = [0, 1, 2, 6, 7, 8]
# nf = [3, 4, 5]
stresses = []
sgrs = [0]
epss = [0]
for i, strain in enumerate(strains):
    reversal, eps_star, sg_star = fiber.update_strain(strain)
    if i in f:
        fiber.finalize()
    stresses.append(fiber.stress)
    epss.append(eps_star)
    sgrs.append(sg_star)
stresses = np.array(stresses)

ax.plot(strains[f], stresses[f], "-", color="black", markerfacecolor="none")
# ax.plot(epss, sgrs, "-", color="black", markerfacecolor="none")
# ax.plot(strains[nf], stresses[nf], "o", color="orange")
ax.grid()
ax.axhline(linewidth=3, color="black")
ax.axvline(linewidth=3, color="black")
# ax.set(xlabel="STEEL STRAIN", ylabel="STEEL STRESS")
plt.show()
