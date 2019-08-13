"""
test file
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# matplotlib.use("Qt5Agg", warn=True)

class KentPark:
    """
    Steel uniaxial material law

    Parameters
    ----------

    Attributes
    ----------

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

        self._c_loading_index = 0
        self._c_strain_p = 0.0
        self._c_strain_r = 0.0
        self._c_stress_r = 0.0
        self._c_strain = 0.0
        self._c_stress = 0.0


    @classmethod
    def eu(cls, K, fc, eu, e0=0.002):
        Z = 0.8 / (eu - e0)
        return cls(fc, K, Z, e0)


    def __str__(self):
        if self._K == 1:
            string = "Unconfined"
        else:
            string = "Confined"
        string += " Concrete Material:\n"
        string += "--------------------------------\n"
        string += f"K\t:\t{self._K}\n"
        string += f"fc\t:\t{self._fc}\n"
        string += f"e_0\t:\t{self._strain_0}\n"
        string += f"Z\t:\t{self._Z}\n"
        return string


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
        self._set_trial_state()


    def _set_trial_state(self):
        deps = self._strain - self._c_strain
        if abs(deps) < 1e-15:  # nearly zero
            self._loading_index = 3
        else:
            if deps < 0:
                self._loading_index = 2
            else:
                self._loading_index = 1
        reversal = self._check_reversal()
        if reversal:
            self._reverse()


    def _check_reversal(self):
        if abs(self._strain) > 1e-15:
            if self._strain < 0:
                if self._c_loading_index in (2, 3):
                    if self._loading_index == 1:
                        return True
        return False

    def _reverse(self):
        """
        FIXME
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
                stress = K * fc * (2 * eps / ep0 - (eps / ep0)**2)
                tangen = K * fc * (2 / ep0 - 2 * (eps / ep0**2))
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
            stress = - (sgr * eps - epp * sgr) / (epr - epp)
            tangen = - sgr / (epr - epp)

        self._stress = -1 * stress
        self._Et = -1 * tangen


    def finalize(self):
        self._c_loading_index = self._loading_index
        self._c_strain_p = self._strain_p
        self._c_strain_r = self._strain_r
        self._c_stress_r = self._stress_r
        self._c_strain = self._strain
        self._c_stress = self._stress


fiber = KentPark(
    fc=6.95,
    K=1,
    Z=770,
    # eu=0.0037,
    e0=0.0027,
)

fig = plt.figure()
ax = fig.add_subplot(111)
strains = np.concatenate((
    np.linspace(0, -0.0025, num=100),
    np.linspace(-0.0025, -0.001, num=100),
    np.linspace(-0.001, -0.003, num=100),
    np.linspace(-0.003, 0, num=100),
    np.linspace(0, -0.009, num=100),
    np.linspace(-0.009, 0, num=100),
    np.linspace(0, -0.01, num=100),
))

stresses = list()
for i, strain in enumerate(strains):
    fiber.update_strain(strain)
    fiber.calculate_stress_and_tangent_modulus()
    fiber.finalize()
    stresses.append(fiber.stress)

line, = ax.plot(strains, stresses, "o", color="black", markerfacecolor="none", label="unconfined")
ax.invert_yaxis()
ax.invert_xaxis()
ax.grid()
ax.axhline(linewidth=3, color="black")
ax.axvline(linewidth=3, color="black")
ax.set(
    xlabel="CONCRETE STRAIN",
    ylabel="CONCRETE STRESS"
)
def update(frame):
    line.set_data(strains[:frame], stresses[:frame])
    return line,
ani = ani.FuncAnimation(fig, update, len(strains), interval=25, blit=True)
plt.show()
