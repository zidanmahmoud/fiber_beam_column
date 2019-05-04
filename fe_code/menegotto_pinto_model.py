from .material import Material

class MenegottoPintoModel(Material):
    """
    Material MenegottoPintoModel provides the constitutive model for steel.

    Attributes
    ----------
    id : int

    youngs_modulus : float

    hardening_modulus : float

    yield_stress : float

    change_in_material_strain_incr : float
    material_strain_incr : float
    material_strain_last_loadstep : float
    material_strain : float

    material_PEEQ_incr : float
    material_PEEQ_last_loadstep :float
    material_PEEQ : float

    material_stress_lasttime_j: float
    material_stress : float

    material_bstress : float

    material_tangent_modulus : float

    To compute
    -----------
    ep_state: bool
        To determine whether the stress state is beyond elastic region.
    """

    def __init__(self, youngs_modulus, asymtopic_modulus, yield_stress, R0, a1, a2):
        """
        Create a new material.

        Parameters
        ----------
        id : int

        youngs_modulus : float

        hardening_modulus : float

        yield_stress : float
        """
        self.youngs_modulus = youngs_modulus
        self.asymtopic_modulus = asymtopic_modulus
        self.yield_stress = yield_stress
        self.yield_strain = yield_stress / youngs_modulus
        self.R0 = R0
        self.a1 = a1
        self.a2 = a2

        self.b = asymtopic_modulus / youngs_modulus
        self.cosi = 0
        self.R = R0
        self.epr = 0
        self.sgr = 0
        self.last_epr = 0
        self.last_sgr = 0
        self.ep0 = self.yield_strain
        self.sg0 = self.yield_stress

        self.change_in_material_strain_incr = 0.
        self.last_loadstep_material_strain_incr = 0.
        self.material_strain_incr = 0.
        self.material_strain_last_loadstep = 0.
        self.material_stress_last_loadstep = 0.
        self.material_strain = 0.
        self.material_PEEQ = 0.
        self.material_stress = 0.
        self.material_tangent_modulus = youngs_modulus

        self.direction = 0.

################################################################################
################################################################################
    '''
    The first block implements functions for attributes
    which are referred in the steps :
    To get value
    To initialize value
    To update value
    '''
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
        # Question: Why is it 1 & -1? It is called in Fiber.initialize_at_beginning(nz)
        if nz > 4:
            self.direction = -1
        elif nz <= 4:
            self.direction = 1

    '''
    material_strain_incr :
         If i == 1 and j == 1: it is initialized to be zero
         Else: it is updated from change_in_fiber_strain_incr

         Function needed:
         change_in_fiber_strain_incr
    '''
    def get_material_strain_incr(self):

        return self.material_strain_incr

    def initialize_material_strain_incr(self):

        self.material_strain_incr = 0.

    def update_material_strain_incr(self):

        self.material_strain_incr += self.change_in_material_strain_incr

    '''
    material_strain :

         It is updated from material_strain_incr
    '''
    def get_material_strain(self):

        return self.material_strain

    def initialize_material_strain(self):

        self.material_strain = 0.

    def update_material_strain(self):

        # if print_:
        #     print('update_material_strain called')
        #     print('material_strain_last_loadstep', self.material_strain_last_loadstep)
        #     print('material_strain_incr', self.material_strain_incr)
        #     print('material_strain', self.material_strain)

        self.material_strain = self.material_strain_last_loadstep + self.material_strain_incr

    def save_to_material_strain_last_loadstep(self, load_step_convergence):

        if load_step_convergence:
            self.material_strain_last_loadstep = self.material_strain
            self.material_stress_last_loadstep = self.material_stress

    '''
    material_stress :

         It is calculated from
         material_strain

         Function needed:
         calculate_material_stress()
    '''
    def get_material_stress(self):

        return self.material_stress

    def initialize_material_stress(self):

        self.material_stress = 0.

    def update_material_stress(self):

        self.material_stress = self.calculate_material_stress()


    '''
    material_tangent_modulus :

         Function needed:
         calculate_material_tangent_modulus()
    '''
    def get_material_tangent_modulus(self):

        return self.material_tangent_modulus

    def initialize_material_tangent_modulus(self):

        self.material_tangent_modulus = self.youngs_modulus

    def update_material_tangent_modulus(self):

        self.material_tangent_modulus = self.calculate_material_tangent_modulus()


    '''
    material_cosi :

         Function needed:
         calculate_material_cosi()
    '''
    def get_material_cosi(self):

        return self.material_cosi

    def initialize_material_cosi(self):

        self.material_cosi = 0.

    def update_material_cosi(self):

        E = self.youngs_modulus
        b = self.b

        eps_intersect = ((self.sgr-self.last_sgr) + E*b*self.last_epr - E*self.epr) /E /(b-1)
        self.cosi = abs(eps_intersect - self.last_epr)

    '''
    Update model parameters:

         ep0
         sg0
         epr
         sgr
         R
         cosi

    '''
    def update_model_parameters(self, nz):

        E = self.youngs_modulus
        b = self.b

        if self.last_loadstep_material_strain_incr == 0:
            if nz > 4:
                self.ep0 = -self.ep0
                self.sg0 = -self.sg0
        else:
            self.last_epr = self.epr
            self.last_sgr = self.sgr
            self.epr = self.material_strain
            self.sgr = self.material_stress
            self.direction *= -1
            self.ep0 = (E*self.epr - self.sgr + self.direction*self.yield_stress*(1-b)) /E /(1-b)
            self.sg0 = b*E*self.ep0 + self.direction*self.yield_stress*(1-b)
            self.update_material_cosi()
            self.R = self.R0 - self.a1*self.cosi /(self.a2 + self.cosi)
            # print_ = False
            # if print_:
            #     print('cosi= ', self.cosi)
            #     print('R= ', self.R)
            #     print('ep0', self.ep0)
            #     print('epr', self.epr)
            #     print('stress = ', self.material_stress)
            #     print('change_in_material_strain_incr', self.change_in_material_strain_incr)
            #     print('last_loadstep_material_strain_incr = ', self.last_loadstep_material_strain_incr)
            #     print('material_strain_incr = ', self.material_strain_incr)
            #     print('material_strain = ', self.material_strain)
            '''

            '''
            self.initialize_material_cosi()

    '''
    other:
    '''

    def get_material_PEEQ(self):

        return self.cosi

################################################################################
################################################################################
    '''
    The second block provides :
    functions which are needed in the first block.
    '''
################################################################################
    def check_reversal(self, do_reversal, j, nx, ny, nz):

        print_ = False

        if do_reversal and (j == 1):
            reversal = True
        else:
            reversal = False


        #if True:
            #print('warning: reversal')

        #reversal *= abs(self.material_strain_incr) > 1e-7
        #reversal *= abs(self.last_loadstep_material_strain_incr) > 1e-7

        return reversal




    def calculate_material_stress(self):

        b = self.b
        eps = self.material_strain
        eps_star = (eps - self.epr) /(self.ep0 - self.epr)
        '''dirty work'''
        if eps_star < 0:
            eps_star = 1e-10
        dg = (eps_star**self.R).real
        #print('dg = ', dg)
        #print('eps_star = ', eps_star)
        sg_star = b*eps_star + (1-b)*eps_star /(1+dg)**(1/self.R)

        return sg_star * (self.sg0 - self.sgr) + self.sgr

    def print_warning(self, nx,ny,nz):
        eps_star = (self.material_strain - self.epr) /(self.ep0 - self.epr)
        # if eps_star < 0 and abs(self.material_strain)>1e-8:
        #     print('section', nx)
        #     print('ny, nz = ', ny, ',', nz)
        #     print('change_in_material_strain_incr', self.change_in_material_strain_incr)
        #     print('warning ')
        #     print('material_strain', self.material_strain)
        #     print('ep0', self.ep0)
        #     print('epr', self.epr)
        #     print('eps_star', eps_star)



    def calculate_material_tangent_modulus(self):

        b = self.b
        eps = self.material_strain
        eps_star = (eps - self.epr) /(self.ep0 - self.epr)
        '''dirty work'''
        if eps_star < 0:
            eps_star = 1e-10
        dg = (eps_star**self.R).real
        Et = b + (1-b) * ((1+dg)**(-1/self.R) - dg * (1+dg)**(-1-1/self.R))
        Et *= (self.sg0 - self.sgr) /(self.ep0 - self.epr)
        #if isnan(Et):
            #print('warning ')
            #print('ep0 = ', self.ep0)
            #print('epr = ', self.epr)

        return Et
