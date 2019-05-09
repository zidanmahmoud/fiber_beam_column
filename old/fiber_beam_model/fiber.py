import numpy as np
from .menegotto_pinto_model import MenegottoPintoModel
class Fiber(object):
    """
    Fiber provides the stress, strain and tangent modulus for each fiber.

    Attributes
    ----------
    fiber_position : ndarray
        [y, z]
        y, z : float

    fiber_area : float

    material_type : a material type class

    fiber_strain_incr : float
    fiber_strain : float
    fiber_strain_last_loadstep : float

    fiber_stress: float

    fiber_tangent_stiffness: float

    fiber_PEEQ : float
    """

    def __init__(self, fiber_position = np.zeros(2), fiber_area = 0, material_type = MenegottoPintoModel(1, 20000, 580, 60, 20, 18.5, 0.0002)):
        """
        Parameters
        ----------
        fiber_position : ndarray
            [y, z]
            y : float
            z : float

        fiber_area : float

        material_type : Material
        """

        self.fiber_position = fiber_position
        self.fiber_area = fiber_area
        self.material_type = material_type

        self.change_in_fiber_strain_incr = 0.
        self.fiber_strain_incr = 0.
        self.fiber_strain = 0.
        self.fiber_strain_last_loadstep = 0.

        self.fiber_stress = 0.
        self.fiber_tangent_stiffness = material_type.get_material_tangent_modulus()
        self.fiber_PEEQ = 0.

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
    def get_material(self):
        return self.material_type
    '''
    Initialization:
    '''

    def initialize_at_beginning(self, nz):

        self.material_type.determin_direction(nz)

    def initialize_fiber_strain_incr(self):

        self.material_type.initialize_material_strain_incr()
        self.fiber_strain_incr = 0.


    def update_change_in_fiber_strain_incr(self, change_in_secdisp_incr):

        self.change_in_fiber_strain_incr = self.calculate_change_in_fiber_strain_incr(change_in_secdisp_incr)

    def get_change_in_fiber_strain_incr(self):

        return self.change_in_fiber_strain_incr
    '''
    material_type:
         connect material to the fiber.
    '''
    def get_material_type(self):

        return self.material_type

    def update_material_parameters_in_loadstep(self, load_step_convergence, j, print_, nx,ny,nz):


        self.material_type.update_change_in_material_strain_incr(self.change_in_fiber_strain_incr)
        #if (ny == 9 and nz==9) or (ny == 7 and nz== 7) or (ny == 0 and nz == 0):
            #if nx == -1.0:
                #print('change_in_fiber_strain_incr = ', self.change_in_fiber_strain_incr)
        self.material_type.update_material_strain_incr()
        self.material_type.update_material_strain(print_)
        #if print_:
            #print('1.material_strain', self.material_type.get_material_strain())
        if self.material_type.check_reversal(load_step_convergence, j, print_, nx,ny,nz):
            self.material_type.update_model_parameters(nz,print_)
            reversal = True
        else:
            reversal = False
        #if print_:
            #print('2.material_strain', self.material_type.get_material_strain())
        self.material_type.update_material_stress()
        self.material_type.update_material_tangent_modulus()
        self.material_type.print_warning(nx,ny,nz)
        # print_ = False
        # if print_:
        #     print('material_stress = ', self.material_type.get_material_stress())
        #     print('tangent_modulus = ', self.material_type.get_material_tangent_modulus())

        return reversal


    '''
    fiber_strain_incr :

         According to STEP(10.1)
         1. If i == 1 and j == 1: it is initialized to be zero.
         2. else: it is updated from change_in_fiber_strain_incr
                  which is interpolated from change_in_secdisp_incr

         Function needed:
    '''
    def get_fiber_strain_incr(self):

        return self.fiber_strain_incr

    def update_fiber_strain_incr(self):

        self.fiber_strain_incr += self.change_in_fiber_strain_incr

    '''
    fiber_strain :

         According to STEP(10.2)
         It is initialized when k == 1.
         It is updated from fiber_strain_incr.
    '''
    def get_fiber_strain(self):

        return self.fiber_strain

    def initialize_fiber_strain(self):

        self.fiber_strain = 0.

    def update_fiber_strain(self):

        self.fiber_strain = self.fiber_strain_last_loadstep + self.fiber_strain_incr

    def save_to_fiber_strain_last_loadstep(self, load_step_convergence):

        if load_step_convergence:
            self.material_type.save_to_last_loadstep_material_strain_incr(load_step_convergence)
            self.material_type.save_to_material_strain_last_loadstep(load_step_convergence)
            self.fiber_strain_last_loadstep = self.fiber_strain

    '''
    fiber_stress :

         According to STEP(11.1)
         It is updated from fiber_strain and fiber_PEEQ.

         Function needed:
         calculate_fiber_stress()
    '''
    def get_fiber_stress(self):

        return self.fiber_stress

    def initialize_fiber_stress(self):

        self.fiber_stress = 0.

    def update_fiber_stress(self):

        self.fiber_stress = self.calculate_fiber_stress()

    '''
    fiber_tangent_stiffness :

         According to STEP(11.2)
         It is updated from fiber_strain and fiber_PEEQ.

         Function needed:
         calculate_fiber_tangent_stiffness()
    '''
    def get_fiber_tangent_stiffness(self):

        return self.fiber_tangent_stiffness

    def initialize_fiber_tangent_stiffness(self):

        self.fiber_tangent_stiffness = self.youngs_modulus

    def update_fiber_tangent_stiffness(self):

        self.fiber_tangent_stiffness = self.calculate_fiber_tangent_stiffness()

    '''
    fiber_PEEQ :

         It is updated from fiber_strain and fiber_PEEQ.

         Function needed:
         calculate_fiber_PEEQ()
    '''
    def get_fiber_PEEQ(self):

        return self.fiber_PEEQ

    def initialize_fiber_PEEQ(self):

        self.fiber_tangent_stiffness = 0.

    def update_fiber_PEEQ(self):

        self.fiber_PEEQ = self.calculate_fiber_PEEQ()

    '''
    others:
    '''
    def get_fiber_position(self):

        return self.fiber_position

    def get_fiber_area(self):

        return self.fiber_area

################################################################################
################################################################################
    '''
    The second block provides :
    functions which are needed in the first block.
    '''
################################################################################
    def calculate_change_in_fiber_strain_incr(self, change_in_secdisp_incr):

        change_in_fiber_strain_incr = 0
        vector = np.zeros(3)
        fiber_position = self.fiber_position
        vector[0] = -fiber_position[0]
        vector[1] = fiber_position[1]
        vector[2] = 1.
        #if (vector[0] == 0.25 and vector[1] == 0.25) or (vector[0] == 0.45 and vector[1] == 0.45):
        '''
        print('vector = ', vector)
        print('change_in_secdisp_incr = ', change_in_secdisp_incr)
        '''
        change_in_fiber_strain_incr = vector @ change_in_secdisp_incr

        return change_in_fiber_strain_incr

    def calculate_fiber_stress(self):

        return self.material_type.get_material_stress()

    def calculate_fiber_tangent_stiffness(self):

        return self.material_type.get_material_tangent_modulus()

    def calculate_fiber_PEEQ(self):

        return self.material_type.get_material_PEEQ()
