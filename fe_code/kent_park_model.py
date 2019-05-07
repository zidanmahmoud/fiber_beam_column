from .material import Material


class KentParkModel(Material):
    """
    Material KentParkModel provides the constitutive model for concretes.

    Attributes
    ----------
    id : int

    compressive_strength : float
        fc'

    confinement_factor : float
        K

    epu : float
        Ultimate strain

    Z : float
        Softening slope

    ep0 : float
        Beyond which the material starts softening

    epr : float
        Unloading strain

    sgr : float
        Unloading stress

    epp : float
        Complete unloading strain, where open crack forms.

    last_NR_change_in_material_strain_incr : float
    change_in_material_strain_incr : float
    material_strain_incr : float
    material_strain_last_loadstep : float
    material_strain : float

    material_stress : float

    material_tangent_modulus : float
    """

    def __init__(self, compressive_strength, confinement_factor, epu, Z):
        """
        Create a new material.

        Parameters
        ----------
        id : int

        compressive_strength : float
            fc'

        confinement_factor : float
            K

        epu : float
            Ultimate strain

        Z : float
            Softening slope
        """

        self.compressive_strength = compressive_strength
        self.confinement_factor = confinement_factor
        self.ep0 = -0.0027 * confinement_factor
        self.epu = epu
        self.Z = Z

        self.epr = 0
        self.sgr = 0
        self.epp = 0

        self.change_in_material_strain_incr = 0.0
        self.last_loadstep_material_strain_incr = 0.0
        self.material_strain_last_loadstep = 0.0
        self.material_stress_last_loadstep = 0.0
        self.material_strain_incr = 0.0
        self.material_strain = 0.0
        self.material_PEEQ = 0.0
        self.material_stress = 0.0
        self.material_tangent_modulus = (
            -2 * self.compressive_strength * self.confinement_factor / self.ep0
        )
        self.direction = 0.0

    ################################################################################
    ################################################################################
    """
    The first block implements functions for attributes
    which are referred in the steps :
    To get value
    To initialize value
    To update value
    """
    ################################################################################
    def get_change_in_material_strain_incr(self):

        return self.change_in_material_strain_incr

    def update_change_in_material_strain_incr(self, change_in_fiber_strain_incr):

        self.change_in_material_strain_incr = change_in_fiber_strain_incr

    def save_to_last_loadstep_material_strain_incr(self, load_step_convergence):

        if load_step_convergence:
            self.last_loadstep_material_strain_incr = self.material_strain_incr

    def get_last_loadstep_material_strain_incr(self):

        return self.last_loadstep_material_strain_incr

    def determin_direction(self, nz):

        if nz > 4:
            self.direction = -1
        elif nz <= 4:
            self.direction = 1

    """
    material_strain_incr :
         If i == 1 and j == 1: it is initialized to be zero
         Else: it is updated from change_in_fiber_strain_incr

         Function needed:
         change_in_fiber_strain_incr
    """

    def get_material_strain_incr(self):

        return self.material_strain_incr

    def initialize_material_strain_incr(self):

        self.material_strain_incr = 0.0

    def update_material_strain_incr(self):

        self.material_strain_incr += self.change_in_material_strain_incr

    """
    material_strain :

         It is updated from material_strain_incr
    """

    def get_material_strain(self):

        return self.material_strain

    def initialize_material_strain(self):

        self.material_strain = 0.0

    def update_material_strain(self):

        self.material_strain = self.material_strain_last_loadstep + self.material_strain_incr

    def save_to_material_strain_last_loadstep(self, load_step_convergence):

        if load_step_convergence:
            self.material_strain_last_loadstep = self.material_strain
            self.material_stress_last_loadstep = self.material_stress

    """
    material_stress :

         It is calculated from
         material_strain

         Function needed:
         calculate_material_stress()
    """

    def get_material_stress(self):

        return self.material_stress

    def initialize_material_stress(self):

        self.material_stress = 0.0

    def update_material_stress(self):

        self.material_stress = self.calculate_material_stress()

    """
    material_tangent_modulus :

         Function needed:
         calculate_material_tangent_modulus()
    """

    def get_material_tangent_modulus(self):

        return self.material_tangent_modulus

    def initialize_material_tangent_modulus(self):

        self.material_tangent_modulus = (
            -2 * self.compressive_strength * self.confinement_factor / self.ep0
        )

    def update_material_tangent_modulus(self):

        self.material_tangent_modulus = self.calculate_material_tangent_modulus()

    """
    epp :

         Function needed:
         calculate_material_cosi()
    """

    def get_epp(self):

        return self.epp

    def initialize_epp(self):

        self.epp = 0.0

    def update_epp(self):

        critical_point = self.epr / self.ep0
        if critical_point < 2 and critical_point > 0:
            self.epp = self.ep0 * (0.145 * critical_point ** 2 + 0.13 * critical_point)
        elif critical_point >= 2:
            self.epp = self.ep0 * (0.707 * (critical_point - 2) + 0.834)

    """
    Update model parameters:

         ep0
         sg0
         epr
         sgr
         R
         cosi

    """

    def update_model_parameters(self, nz):

        if self.material_strain < 0:
            self.epr = self.material_strain_last_loadstep
            self.sgr = self.material_stress_last_loadstep
            self.update_epp()
        print_ = False
        if print_:
            print("epr = ", self.epr)
            print("sgr = ", self.sgr)
            print("epp = ", self.epp)
            print("last_loadstep_material_strain_incr = ", self.last_loadstep_material_strain_incr)
            print("material_strain_incr = ", self.material_strain_incr)
            print("material_strain = ", self.material_strain)

    """
    other:
    """

    def get_material_PEEQ(self):

        return self.epp

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
            reversal = Fals

        return reversal

    def check_load(self):

        loading = True
        loading *= self.material_strain_incr < 0
        return loading

    def print_warning(self, nx, ny, nz):

        eps = self.material_strain
        epr = self.epr

    def calculate_material_stress(self):

        eps = self.material_strain
        K = self.confinement_factor
        fc = self.compressive_strength
        ep0 = self.ep0
        Z = self.Z
        epr = self.epr
        sgr = self.sgr
        epp = self.epp
        epu = self.epu
        if self.check_load():
            if eps > epp:
                sg = 0
            elif eps > epr and eps < epp:
                sg = sgr / (epr - epp) * (eps - epr) + sgr
            elif eps <= epr:
                if eps > 0:
                    sg = 0
                elif eps > ep0:
                    sg = K * fc * (-2 * (eps / ep0) + (eps / ep0) ** 2)
                elif eps > epu:
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    sg = sg * (sg < -0.2 * K * fc) + -0.2 * K * fc * (sg >= -0.2 * K * fc)
                else:
                    sg = 0
        else:
            if eps <= epr:
                if eps > 0:
                    sg = 0
                elif eps > ep0:
                    sg = K * fc * (-2 * (eps / ep0) + (eps / ep0) ** 2)
                elif eps > epu:
                    sg = K * fc * (-1 - Z * (eps - ep0))
                    sg = sg * (sg < -0.2 * K * fc) + -0.2 * K * fc * (sg >= -0.2 * K * fc)
                else:
                    sg = 0

            elif eps < epp:
                sg = sgr / (epr - epp) * (eps - epr) + sgr
            else:
                sg = 0

        return sg

    def calculate_material_tangent_modulus(self):

        eps = self.material_strain
        K = self.confinement_factor
        fc = self.compressive_strength
        ep0 = self.ep0
        Z = self.Z
        epr = self.epr
        sgr = self.sgr
        epp = self.epp
        epu = self.epu
        if self.check_load():
            if eps > epp:
                Et = 0
            elif eps > epr and eps < epp:
                Et = sgr / (epr - epp)
            elif eps <= epr:
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
            if eps <= epr:
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
            elif eps < epp:
                Et = sgr / (epr - epp)
            else:
                Et = 0

        return Et
