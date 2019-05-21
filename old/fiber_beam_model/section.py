import numpy as np
from numpy import linalg as la


class Section(object):
    """
    Section provides the section force, displacement and constitutive matrix.

    Attributes
    ----------
    section_position : int

    section_shape : ndarray
        [b, h]

    fiber_list : ndarray
        Array of Fiber objects providing real stress and tangent modulus of each fiber.

    secforce_incr : ndarray
        3x1 array

    secforces : ndarray
        3x1 array

    secdisp_incr : ndarray
        3x1 array

    residual_secdisps : ndarray
        3x1 array

    sec_stiffness_matrix : ndarray
        3x3 array

    section_tolerance : float

    To compute
    ----------
    change_in_secdisp_incr: ndarray
        3x1 array

    change_in_secforce_incr: ndarray
        3x1 array

    sec_resisting_forces : ndarray
        3x1 array

    """

    def __init__(self, section_position, section_weight, section_shape, fiber_list, section_tolerance):
        """
        Create a new section.

        Parameters
        ----------
        section_position : int

        section_shape : ndarray
            [b, h]

        fiber_list : ndarray
            Array of Fiber objects providing real stress and tangent modulus of each fiber.

        section_tolerance : float

        """
        self.section_position = section_position
        self.section_weight = section_weight
        self.section_shape = section_shape
        self.fiber_list = fiber_list
        self.section_tolerance = section_tolerance

        self.change_in_secforce_incr = np.zeros(3)
        self.change_in_secdisp_incr = np.zeros(3)

        self.secforce_incr = np.zeros(3)
        self.secforces = np.zeros(3)
        self.secforces_last_loadstep = np.zeros(3)

        self.secdisp_incr = np.zeros(3)
        self.secdisp = np.zeros(3)
        self.secdisp_last_loadstep = np.zeros(3)
        self.residual_secdisps = np.zeros(3)
        self.sec_stiffness_matrix = np.zeros([3, 3])

################################################################################
################################################################################
    '''
    The first block is implementing functions for attributes
    which are referred in the steps :
    To get value
    To initialize value
    To update value
    '''
################################################################################
    def initialize_at_beginning(self):

        self.update_sec_stiffness_matrix()

        for fibers_row in self.fiber_list:
            for nz, fiber in enumerate(fibers_row):
                fiber.initialize_at_beginning(nz)

    def initialize_in_loadstep(self):

        self.initialize_secforce_incr()
        self.initialize_secdisp_incr()

        fiber_list = self.fiber_list
        for ny in range(len(fiber_list[:][0])):
            for nz in range(len(fiber_list[:][0])):
                single_fiber = fiber_list[ny][nz]
                single_fiber.initialize_fiber_strain_incr()
                fiber_list[ny][nz] = single_fiber
        self.fiber_list = fiber_list

    def update_change_in_secforce_incr(self, change_in_element_force_incr):

        self.change_in_secforce_incr = self.calculate_change_in_secforce_incr(change_in_element_force_incr)
        #print('step8. change_in_secforce_incr = ', self.change_in_secforce_incr)

    def update_change_in_secdisp_incr(self):


        self.change_in_secdisp_incr = self.calculate_change_in_secdisp_incr()
        #print('step9. change_in_secdisp_incr = ', self.change_in_secdisp_incr)



    def update_fiber_parameters_in_loadstep(self, load_step_convergence, j, i, k):

        fiber_list = self.fiber_list
        reverse = 0
        for ny in range(len(fiber_list[:][0])):
            for nz in range(len(fiber_list[:][0])):
                single_fiber = fiber_list[ny][nz]

                single_fiber.update_change_in_fiber_strain_incr(self.change_in_secdisp_incr)

                # if (ny == 0 and nz==0) or (ny == 14 and nz==14) or (ny == 1 and nz==1) or (ny == 13 and nz== 13):
                #     if self.section_position<0 and self.section_position>-1:
                #         #print('section_position', self.section_position)
                #         #print('ny, nz = ', ny, ',', nz)
                #         print_ = True
                #     elif self.section_position>0.5 and self.section_position<2:
                #         #print('section_position', self.section_position)
                #         #print('ny, nz = ', ny, ',', nz)
                #         print_ = True
                #     elif self.section_position>-2 and self.section_position<-0.9:
                #         #print('section_position', self.section_position)
                #         #print('ny, nz = ', ny, ',', nz)
                #         print_ = True
                #     else:
                #         print_ = False
                # else:
                #     print_ = False

                single_fiber.update_fiber_strain_incr()
                single_fiber.update_fiber_strain()
                fiber_reverse = single_fiber.update_material_parameters_in_loadstep(load_step_convergence, j, False, self.section_position, ny, nz)
                reverse += fiber_reverse
                single_fiber.update_fiber_stress()
                single_fiber.update_fiber_tangent_stiffness()
                single_fiber.update_fiber_PEEQ()
                fiber_list[ny][nz] = single_fiber

                if (ny == 0 and nz==5) or (ny == 9 and nz==9) or (ny == 0 and nz==0):
                    if self.section_position<0 and self.section_position>-0.9:
                        '''
                        print('step10. last_loadstep_material_strain_incr = ', (single_fiber.get_material()).get_last_loadstep_material_strain_incr())
                        print('step10. change_in_fiber_strain_incr = ', single_fiber.get_change_in_fiber_strain_incr())
                        print('step10. fiber_strain_incr = ', single_fiber.get_fiber_strain_incr())
                        print('step10. fiber_strain = ', single_fiber.get_fiber_strain())
                        print('step11. fiber_stress = ', single_fiber.get_fiber_stress())
                        print('step11. fiber_fiber_tangent_stiffness = ', single_fiber.get_fiber_tangent_stiffness())

                        print('step10. change_in_material_strain_incr = ', (single_fiber.get_material()).get_change_in_material_strain_incr())
                        print('step10. material_strain_incr = ', (single_fiber.get_material()).get_material_strain_incr())
                        print('step10. material_strain = ', (single_fiber.get_material()).get_material_strain())
                        '''
        self.fiber_list = fiber_list

        return reverse

    def save_to_last_loadstep(self, load_step_convergence):

        if load_step_convergence:
            fiber_list = self.fiber_list
            for ny in range(len(fiber_list[:][0])):
                for nz in range(len(fiber_list[:][0])):
                    single_fiber = fiber_list[ny][nz]
                    single_fiber.save_to_fiber_strain_last_loadstep(load_step_convergence)
                    fiber_list[ny][nz] = single_fiber
            self.fiber_list = fiber_list

            self.secforces_last_loadstep = self.secforces
            self.secdisp_last_loadstep = self.secdisp




    '''
    secforce_incr :

         According to STEP(8.1)
         1. If i == 1 and j == 1: it is initialized to be zero.
         2. else: it is updated from change_in_secforce_incr
                  which is interpolated from change_in_element_force_incr

         Function needed:
         calculate_change_in_secforce_incr(change_in_element_force_incr)
    '''
    def get_secforce_incr(self):

        return self.secforce_incr

    def initialize_secforce_incr(self):

        self.secforce_incr = np.zeros(3)

    def update_secforce_incr(self):

        self.secforce_incr += self.change_in_secforce_incr

    '''
    secforces :

         According to STEP(8.2)
         It is initialized when k == 1.
         It is updated from secforce_incr.
    '''
    def get_secforces(self):

        return self.secforces

    def initialize_secforces(self):

        self.secforces = np.zeros(3)

    def update_secforces(self):

        self.secforces = self.secforces_last_loadstep + self.secforce_incr


    '''
    secdisp_incr :

         According to STEP(9)
         If i == 1 and j == 1: it is initialized to be zeros.
         Else: it is updated from residual_secdisps and change_in_secforce_incr

         Function needed:
         get_residual_secdisps()
         calculate_change_in_secdisp_incr(change_in_element_force_incr)
    '''
    def get_secdisp_incr(self):

        return self.secdisp_incr

    def initialize_secdisp_incr(self):

        self.secdisp_incr = np.zeros(3)

    def update_secdisp_incr(self):

        self.secdisp_incr += self.change_in_secdisp_incr

    def get_secdisp(self):

        return self.secdisp

    def update_secdisp(self):

        self.secdisp = self.secdisp_last_loadstep + self.secdisp_incr

    '''
    sec_stiffness_matrix :

         According to STEP(12)
         It is updated as integration of fiber_tangent_stiffness over the section.

         Function needed:
         calculate_sec_stiffness_matrix()
    '''
    def get_sec_stiffness_matrix(self):

        return self.sec_stiffness_matrix

    def update_sec_stiffness_matrix(self):
        #print('update_sec_stiffness_matrix called')

        self.sec_stiffness_matrix = self.calculate_sec_stiffness_matrix()

    '''
    residual_secdisps :

         According to STEP(9) and (15)
         If j == 0: it is initialized as zeros.
         It is updated from sec_stiffness_matrix and secforces and sec_resisting_forces

         Function needed:
         get_sec_stiffness_matrix()
         get_secforces()
         calculate_sec_resisting_forces()
    '''
    def get_residual_secdisps(self):

        return self.residual_secdisps

    def initialize_residual_secdisps(self):

        self.residual_secdisps = np.zeros(3)

    def update_residual_secdisps(self):

        sec_flexibility_matrix = la.inv(self.sec_stiffness_matrix)
        secforces = self.secforces
        sec_resisting_forces = self.calculate_sec_resisting_forces()

        self.residual_secdisps = sec_flexibility_matrix @ (secforces - sec_resisting_forces)

    '''
    check_section_convergence :

         According to STEP(17.1)
         Norm of section_unbalanced_forces < section_tolerance
    '''
    def check_section_convergence(self):

        secforces = self.secforces
        sec_resisting_forces = self.calculate_sec_resisting_forces()
        section_unbalanced_forces = secforces - sec_resisting_forces

        if la.norm(section_unbalanced_forces) < self.section_tolerance:
            return True
        else:
            return False

################################################################################
################################################################################
    '''
    The second block provides :
    functions which are needed in the first block.
    '''
################################################################################
    def calculate_change_in_secforce_incr(self, change_in_element_force_incr):

        b_matrix = self.calculate_b_matrix(self.section_position)
        #print('b matrix = ', b_matrix)
        change_in_secforce_incr = b_matrix @ change_in_element_force_incr

        return change_in_secforce_incr

    def calculate_change_in_secdisp_incr(self):

        change_in_secforce_incr = self.change_in_secforce_incr
        sec_flexibility_matrix = la.inv(self.sec_stiffness_matrix)
        residual_secdisps = self.residual_secdisps
        #print('sec_stiffness_matrix = ', self.sec_stiffness_matrix)
        #print('sec_flexibility_matrix = ', sec_flexibility_matrix)
        #print('change_in_secforce_incr = ',change_in_secforce_incr)
        return residual_secdisps + sec_flexibility_matrix @ change_in_secforce_incr

    def calculate_sec_resisting_forces(self):

        sec_resisting_forces = np.zeros(3)
        vector = np.zeros(3)
        for ny, fiber_row in enumerate(self.fiber_list[:]):
            for nz, fiber in enumerate(fiber_row):
                fiber_position = fiber.get_fiber_position()
                fiber_stress = fiber.get_fiber_stress()
                fiber_area = fiber.get_fiber_area()
                vector[0] = -fiber_position[0]
                vector[1] = fiber_position[1]
                vector[2] = 1.
                sec_resisting_forces += fiber_stress * fiber_area * vector

        return sec_resisting_forces


    def calculate_sec_stiffness_matrix(self):
        #print('calculate_sec_stiffness_matrix called')
        sec_stiffness_matrix = np.zeros([3, 3])
        fiber_list = self.fiber_list
        vector = np.zeros(3)

        for ny in range(len(fiber_list[:][0])):
            for nz in range(len(fiber_list[:][0])):
                single_fiber = fiber_list[ny][nz]
                fiber_position = single_fiber.get_fiber_position()
                fiber_tangent_stiffness = single_fiber.get_fiber_tangent_stiffness()
                fiber_area = single_fiber.get_fiber_area()
                vector[0] = -fiber_position[0]
                vector[1] = fiber_position[1]
                vector[2] = 1.
                EA = fiber_tangent_stiffness * fiber_area
                sec_stiffness_matrix += EA * np.outer(vector, vector)
                fiber_list[ny][nz] = single_fiber
        self.fiber_list = fiber_list

        return sec_stiffness_matrix

    '''
    b_matrix:

         3x5 array
         Interpolation from nodal forces (5x1) to section forces(3x1).
    '''
    def calculate_b_matrix(self, cosi):

        b_matrix = np.zeros([3, 5])
        b_matrix[0, 0] = cosi /2 - 1/2
        b_matrix[0, 1] = cosi /2 + 1/2
        b_matrix[1, 2] = cosi /2 - 1/2
        b_matrix[1, 3] = cosi /2 + 1/2
        b_matrix[2, 4] = 1

        return b_matrix
