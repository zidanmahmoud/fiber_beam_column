"""
test file
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Qt5Agg", warn=True)

class KentPark:
    """
    Steel uniaxial material law

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, compressive_strength, confinement_factor, Z):
        self._compressive_strength = compressive_strength
        self._confinement_factor = confinement_factor
        self._strain_0 = -0.002 * confinement_factor
        self._Z = Z

        self._strain_r = 0.0
        self._stress_r = 0.0
        self._strain_p = 0.0

        self._last_strain = 0.0
        self._strain = 0.0
        self.tangent_modulus = 2 * compressive_strength / self._strain_0
        self.stress = 0.0

    def update_strain(self, value):
        """
        FIXME
        """

        if self._check_reversal(value):
            self.reverse()

            global ax
            ax.plot(self._strain_r, self._stress_r, "-o", color="black", markerfacecolor="none")
            ax.plot(self._strain_p, 0, "-o", color="black")

        self._last_strain = self._strain
        self._strain = value

    def _check_reversal(self, new_strain):

        # workaround
        if self._last_strain == 0:
            return False
        # if strain was decreasing
        if self._strain - self._last_strain < 0:
            # if the new strain is increasing
            if new_strain > self._strain:
                return True
        return False

    def reverse(self):
        """
        FIXME
        """
        self._strain_r = self._strain
        self._stress_r = self.stress

        critical_point = self._strain_r / self._strain_0
        ep0 = self._strain_0
        if critical_point >= 2:
            self._strain_p = ep0 * (0.707 * (critical_point - 2) + 0.834)
        elif critical_point < 2:
            self._strain_p = ep0 * (0.145 * critical_point ** 2 + 0.13 * critical_point)

    def calculate_stress_and_tangent_modulus(self):
        eps = self._strain
        ep0 = self._strain_0
        epp = self._strain_p
        epr = self._strain_r
        sgr = self._stress_r
        K = self._confinement_factor
        Z = self._Z
        fc = self._compressive_strength

        # == inequality signs are reversed compared to theory becuase of the negative signs

        # positive strain
        if eps >= 0:
            self.stress = 0.0
            self.tangent_modulus = 0.0
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
                self.stress = 0.0
                self.tangent_modulus = 0.0
                return
            stress = - (sgr * eps - epp * sgr) / (epr - epp)
            tangen = - sgr / (epr - epp)

        self.stress = -1 * stress
        self.tangent_modulus = -1 * tangen


fiber = KentPark(
    compressive_strength=6.95,
    confinement_factor=1,
    Z=770
)

fig = plt.figure()
ax = fig.add_subplot(111)

strains = np.concatenate((
    # np.linspace(-0.000001, -0.004, num=200),
    # np.linspace(-0.0039, -0.000001, num=200),
    # np.linspace(-0.0000011, -0.007, num=200),
    # np.linspace(-0.0069, -0.000001, num=200),
    np.linspace(-0.0000011, -0.01, num=200),
))
# strains.reshape(strains.size)
stresses = list()
for i, strain in enumerate(strains):
    fiber.update_strain(strain)
    fiber.calculate_stress_and_tangent_modulus()
    stresses.append(fiber.stress)

ax.plot(strains, stresses, "-", color="black")

ax.invert_yaxis()
ax.invert_xaxis()
ax.grid()
ax.axhline(linewidth=3, color="black")
ax.axvline(linewidth=3, color="black")
ax.set(
    xlabel="CONCRETE STRAIN",
    ylabel="CONCRETE STRESS"
)
plt.show()
