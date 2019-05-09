import numpy as np

class IsotropicLinearHardenging(object):
    """
    Material IsotropicLinearHardenging provides the constitutive model for each point.

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

    def __init__(self, id, youngs_modulus, hardening_modulus, yield_stress):
        """
        Create a new material.

        Parameters
        ----------
        id : int

        youngs_modulus : float

        hardening_modulus : float

        yield_stress : float
        """

        self.id = id
        self.youngs_modulus = youngs_modulus
        self.hardening_modulus = hardening_modulus
        self.yield_stress = yield_stress

        self.change_in_material_strain_incr = 0.
        self.material_strain_incr = 0.
        self.material_strain_last_loadstep = 0.
        self.material_strain = 0.
        self.material_PEEQ_incr = 0.
        self.material_PEEQ_last_loadstep = 0.
        self.material_PEEQ = 0.
        self.material_stress_lasttime_j = 0.
        self.material_stress = 0.
        self.material_bstress = 0.
        self.material_tangent_modulus = youngs_modulus

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
    def update_change_in_material_strain_incr(self, change_in_fiber_strain_incr):

        self.change_in_material_strain_incr = change_in_fiber_strain_incr

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

        self.material_strain = self.material_strain_last_loadstep + self.material_strain_incr

    def save_to_material_strain_last_loadstep(self):

        self.material_strain_last_loadstep = self.material_strain

    '''
    material_stress :

         It is calculated from
         material_strain
         change_in_material_stress_incr
         change_in_material_PEEQ_incr
         (return mapping algorithm)

         Function needed:
         calculate_change_in_material_stress_incr()
    '''
    def get_material_stress(self):

        return self.material_stress

    def initialize_material_stress(self):

        self.material_stress = 0.

    def update_material_stress(self):

        self.material_stress += self.calculate_change_in_material_stress_incr()

    def save_to_material_stress_lasttime_j(self):

        self.material_stress_lasttime_j = self.material_stress

    '''
    material_tangent_modulus :

         It is calculated from

         Function needed:
         calculate_material_tangent_modulus()
    '''
    def get_material_tangent_modulus(self):

        return self.material_tangent_modulus

    def initialize_material_tangent_modulus(self):

        self.material_tangent_modulus = youngs_modulus

    def update_material_tangent_modulus(self):

        self.material_tangent_modulus = self.calculate_material_tangent_modulus()

    '''
    material_PEEQ_incr :
         IF i == 0 and j == 0: it initialized as zero
         Else: it is updated from change_in_material_PEEQ_incr

         Function needed:
         calculate_change_in_material_PEEQ_incr()
    '''
    def get_material_PEEQ_incr(self):

        return self.material_PEEQ_incr

    def initialize_material_PEEQ_incr(self):

        self.material_PEEQ_incr = 0.

    def update_material_PEEQ_incr(self):

        self.material_PEEQ_incr += self.calculate_change_in_material_PEEQ_incr()

    '''
    material_PEEQ :

         It is updated from material_PEEQ_incr
    '''
    def get_material_PEEQ(self):

        return self.material_PEEQ

    def initialize_material_PEEQ(self):

        self.material_PEEQ = 0.

    def update_material_PEEQ(self):

        self.material_PEEQ = self.material_PEEQ_last_loadstep + self.material_PEEQ_incr

    def save_to_material_PEEQ_last_loadstep(self):

        self.material_PEEQ_last_loadstep = self.material_PEEQ

    '''
    other:
    '''
    def get_material_id(self):

        return self.id

################################################################################
################################################################################
    '''
    The second block provides :
    functions which are needed in the first block.
    '''
################################################################################
    def calculate_change_in_material_PEEQ_incr(self):

        material_stress_lasttime_j = self.material_stress_lasttime_j
        yield_stress = self.yield_stress
        youngs_modulus = self.youngs_modulus
        hardening_modulus = self.hardening_modulus
        material_PEEQ = self.material_PEEQ
        change_in_material_strain_incr = self.change_in_material_strain_incr

        material_stress_predicted = material_stress_lasttime_j + youngs_modulus * change_in_material_strain_incr
        ep_state = abs(material_stress_predicted) > (yield_stress + hardening_modulus * material_PEEQ)
        delta_lamda = ep_state * (abs(material_stress_predicted) - yield_stress - hardening_modulus * material_PEEQ) / (youngs_modulus + hardening_modulus)

        return delta_lamda

    def calculate_change_in_material_stress_incr(self):

        material_stress_lasttime_j = self.material_stress_lasttime_j
        yield_stress = self.yield_stress
        youngs_modulus = self.youngs_modulus
        hardening_modulus = self.hardening_modulus
        material_PEEQ =self.material_PEEQ
        change_in_material_strain_incr = self.change_in_material_strain_incr

        material_stress_predicted = material_stress_lasttime_j + youngs_modulus * change_in_material_strain_incr
        ep_state = abs(material_stress_predicted) > (yield_stress + hardening_modulus * material_PEEQ)
        delta_lamda = ep_state * (abs(material_stress_predicted) - yield_stress - hardening_modulus * material_PEEQ) / (youngs_modulus + hardening_modulus)
        change_in_material_stress_incr = youngs_modulus * (change_in_material_strain_incr - ep_state * delta_lamda * ((material_stress_predicted > 0) * 2 - 1))
        return change_in_material_stress_incr

    def calculate_material_tangent_modulus(self):

        material_stress_lasttime_j = self.material_stress_lasttime_j
        yield_stress = self.yield_stress
        youngs_modulus = self.youngs_modulus
        hardening_modulus = self.hardening_modulus
        material_PEEQ =self.material_PEEQ
        change_in_material_strain_incr = self.change_in_material_strain_incr

        material_stress_predicted = material_stress_lasttime_j + youngs_modulus * change_in_material_strain_incr
        ep_state = abs(material_stress_predicted) > (yield_stress + hardening_modulus * material_PEEQ)
        material_tangent_modulus = (1 - ep_state) * youngs_modulus + ep_state * youngs_modulus * hardening_modulus / (youngs_modulus + hardening_modulus)

        return material_tangent_modulus
