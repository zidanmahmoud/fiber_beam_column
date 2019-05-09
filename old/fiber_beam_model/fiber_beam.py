import numpy as np
import numpy.linalg as la
import sys
from .integral_scheme import GaussLobatto
import matplotlib.pyplot as plt

class FiberBeam(object):

    """
    Three dimensional beam element providing:
    element_global_stiffness_matrix
    element_global_internal_forces
    element_convergence

    Attention:
    This is the 3d beam element for infinitesimal non-torsional deformation.


    Attributes
    ----------
    id : int

    reference_x_a: ndarray
        First node location.

    reference_x_b: ndarray
        Second node location.

    reference_abt_vec : ndarray
        3x1 arbitary input vector
        Note: it must NOT be parallel with the undeformed beam element.

    element_disp_incr : ndarray
        5x1 Array

    element_local_stiffness_matrix : ndarray
        5x5 Array

    element_local_stiffness_matrix_lasttime_i : ndarray
        5x5 Array

    element_force_incr : ndarray
        5x1 Array

    element_resisting_forces : ndarray
        5x1 Array

    residual_element_disps : ndarray
        5x1 Array

    To calculate
    ------------
    reference_triad : ndarray
        3x3 Array of three unit vetors of initial triad.

    change_in_element_disp_incr: ndarray
        5x1 Array

    element_global_stiffness_matrix : ndarray
        12x12 Array

    element_global_internal_forces : ndarray
        12x1 Array

    reference_transform_matrix : ndarray
        12x5 Array

    element_convergence : bool

    """

    def __init__(self, id, Nd_a, Nd_b, vec, section_list):
        """
        Creat a new beam element.

        Attributes
        ----------
        id : int

        reference_x_a: ndarray
            First node location.

        reference_x_b: ndarray
            Second node location.

        reference_abt_vec : ndarray
            3x1 arbitary input vector
            Note: it must NOT be parallel with the undeformed beam element.

        """

        self.id = id
        self.Nd_a = Nd_a
        self.Nd_b = Nd_b
        self.reference_abt_vec = vec
        self.section_list = section_list

        self.change_in_element_disp_incr = np.zeros(5)
        self.change_in_element_force_incr = np.zeros(5)

        self.element_disp_incr = np.zeros(5)
        self.element_disp = np.zeros(5)
        self.element_disp_last_loadstep = np.zeros(5)

        self.element_local_stiffness_matrix = np.zeros([5, 5])
        self.element_local_stiffness_matrix_last_NR_i = np.zeros([5, 5])

        self.element_force_incr = np.zeros(5)
        self.element_resisting_forces = np.zeros(5)
        self.element_resisting_forces_last_loadstep = np.zeros(5)

        self.residual_element_disps = np.zeros(5)


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
    '''
    Initialization:
    '''
    def initialize_at_beginning(self):
        for section in self.section_list:
            section.initialize_at_beginning()

        self.update_element_local_stiffness_matrix()
        self.save_to_NRstep()

    def initialize_in_loadstep(self):

        self.initialize_element_disp_incr()
        self.initialize_element_force_incr()

        section_list = self.section_list
        for nx in range(len(section_list)):
            single_section = section_list[nx]
            single_section.initialize_in_loadstep()
            section_list[nx] = single_section
        self.section_list = section_list

    def update_change_in_element_disp_incr(self):
        #print('step6. change_in_element_disp_incr', self.calculate_change_in_element_disp_incr())
        self.change_in_element_disp_incr = self.calculate_change_in_element_disp_incr()

    def update_change_in_element_force_incr(self, j):
        #print('step6. change_in_element_force_incr', self.calculate_change_in_element_force_incr(j))
        self.change_in_element_force_incr = self.calculate_change_in_element_force_incr(j)

    def update_section_parameters_in_loadstep(self, load_step_convergence, j, i, k):

        reverse = 0


        # print("began section loop")
        for section in self.section_list:

            section.update_change_in_secforce_incr(self.change_in_element_force_incr)
            section.update_change_in_secdisp_incr()
            section.update_secforce_incr()
            section.update_secforces()
            section.update_secdisp_incr()
            section.update_secdisp() # Question: not in the thesis' steps!
            #print('step9. secdisp_incr = ', single_section.get_secdisp_incr())

            sec_reverse = section.update_fiber_parameters_in_loadstep(load_step_convergence, j, i, k)
            reverse += sec_reverse
            section.update_sec_stiffness_matrix()
            section.update_residual_secdisps()
            #print('step15. residual_secdisps = ', single_section.get_residual_secdisps())
            #print('\n\n')

        return reverse

    def save_to_NRstep(self):

        self.element_local_stiffness_matrix_last_NR_i = self.element_local_stiffness_matrix

    def update_element_disp(self):

        self.element_disp = self.element_disp_last_loadstep + self.element_disp_incr

    def get_element_disp(self):

        return self.element_disp


    def save_to_last_loadstep(self, load_step_convergence):

        if load_step_convergence:
            section_list = self.section_list
            for nx in range(len(section_list)):
                single_section = section_list[nx]
                single_section.save_to_last_loadstep(load_step_convergence)
                section_list[nx] = single_section
            self.section_list = section_list

            self.element_resisting_forces_last_loadstep = self.element_resisting_forces
            self.element_disp_last_loadstep = self.element_disp

    '''
    element_disp_incr :

         According to STEP(4)
         1. If i == 1: it is initialized to be zero.
         2. If i > 1: it is updated from global_disp_incr
         3. It does NOT change during loop j.
    '''
    def get_element_disp_incr(self):

        return self.element_disp_incr

    def initialize_element_disp_incr(self):

        self.element_disp_incr = np.zeros(5)

    def update_element_disp_incr(self):

        self.element_disp_incr += self.change_in_element_disp_incr

    '''
    element_local_stiffness_matrix :

         According to STEP(6.2)
         Must be saved as an attribute for the object

         Function needed:
         calculate_element_local_stiffness_matrix()
    '''
    def get_element_local_stiffness_matrix(self):

        return self.element_local_stiffness_matrix

    def update_element_local_stiffness_matrix(self):

        self.element_local_stiffness_matrix = self.calculate_element_local_stiffness_matrix()

    '''
    element_force_incr :

         According to STEP(6) and (7.1)
         1. If i == 1 and j == 1: it is initialized to be zero,
         2. If i > 1 and j == 1: it is updated from change_in_element_disp_incr.
         3. If j > 1: it is updated from residual_element_disps
    '''
    def get_element_force_incr(self):

        return self.element_force_incr

    def initialize_element_force_incr(self):

        self.element_force_incr = np.zeros(5)

    def update_element_force_incr(self):

        self.element_force_incr += self.change_in_element_force_incr

    '''
    element_resisting_forces :

         According to STEP(7.2)
         It is updated from element_force_incr
    '''
    def get_element_resisting_forces(self):

        return self.element_resisting_forces

    def initialize_element_resisting_forces(self):

        self.element_resisting_forces = np.zeros(5)

    def update_element_resisting_forces(self):

        self.element_resisting_forces = self.element_resisting_forces_last_loadstep + self.element_force_incr

    '''
    residual_element_disps :

         According to STEP(17.2)
         After checking element convergence
         It is updated as integration of residual_secdisps

         Function needed:
         calculate_residual_element_deformation()
    '''
    def get_residual_element_disps(self):

        return self.residual_element_disps

    def update_residual_element_disps(self):

        self.residual_element_disps = self.calculate_residual_element_deformation()

    '''
    element_convergence:

         Function needed:
         Section.check_section_convergence()
    '''
    def check_for_element_convergence(self):

        element_convergence = True

        section_list = self.section_list

        for nk in range(len(section_list)):
            single_section = section_list[nk]
            element_convergence *= single_section.check_section_convergence()

        return element_convergence
################################################################################
################################################################################
    '''
    The second block provides :
    functions which are needed in the first block.
    '''
################################################################################
    '''
    change_in_element_disp_incr:

         Get rid of rigid body motion with reference_transform_matrix.
         q = T.transpose * p
    '''
    def calculate_change_in_element_disp_incr(self):

        change_in_node_displacement_incr_a = self.Nd_a.get_change_in_node_displacement_incr()
        change_in_node_displacement_incr_b = self.Nd_b.get_change_in_node_displacement_incr()
        element_global_displacements = np.concatenate((change_in_node_displacement_incr_a, change_in_node_displacement_incr_b), axis = 0)
        reference_transform_matrix = self.get_reference_transform_matrix()

        return reference_transform_matrix.T @ element_global_displacements

    '''
    change_in_element_force_incr:

         calculate from change_in_element_disp_incr
    '''
    def calculate_change_in_element_force_incr(self, j):

        if j == 1:
            change_in_element_disp_incr = self.change_in_element_disp_incr
            print(change_in_element_disp_incr); input()
            element_local_stiffness_matrix_last_NR_i = self.element_local_stiffness_matrix_last_NR_i
            change_in_element_force_incr = element_local_stiffness_matrix_last_NR_i @ change_in_element_disp_incr
        elif j > 1:
            residual_element_disps = self.residual_element_disps
            print(residual_element_disps); input()
            element_local_stiffness_matrix = self.element_local_stiffness_matrix
            change_in_element_force_incr = -element_local_stiffness_matrix @ residual_element_disps

        return change_in_element_force_incr

    '''
    element_local_stiffness_matrix:

         5x5 array

         Function needed:
         Section.get_sec_stiffness_matrix()
    '''
    def calculate_element_local_stiffness_matrix(self):

        section_list = self.section_list

        element_local_flexibility_matrix = np.zeros([5, 5])
        reference_length = self.get_reference_length()

        for section in section_list:
            section_flexibility_matrix = la.inv(section.get_sec_stiffness_matrix())
            #print('step12. get_sec_stiffness_matrix = ', single_section.get_sec_stiffness_matrix())
            b_matrix = self.calculate_b_matrix(section.section_position)
            element_local_flexibility_matrix += reference_length / 2 * section.section_weight * (b_matrix.T @ section_flexibility_matrix @ b_matrix)

        return la.inv(element_local_flexibility_matrix)

    '''
    residual_element_disps:

         5x1 array

         Function needed:
         Section.get_residual_secdisps()
    '''
    def calculate_residual_element_deformation(self):

        section_list = self.section_list

        residual_element_disps = np.zeros(5)
        reference_length = self.get_reference_length()
        [gauss_points, weights] = GaussLobatto(len(section_list))

        for nk in range(len(section_list)):
            single_section = section_list[nk]
            residual_secdisps = single_section.get_residual_secdisps()
            b_matrix = self.calculate_b_matrix(gauss_points[nk])
            residual_element_disps += reference_length /2. * weights[nk] * (b_matrix.T @ residual_secdisps)

        #print('residual_element_disps = ', residual_element_disps)

        return residual_element_disps

################################################################################
################################################################################
    '''
    The third block provides :
    stiffness matrix and internal forces at global coordinates.
    '''
################################################################################
    '''
    element_global_stiffness_matrix:

         element stiffness matrix at global coordinates.
         12x12 array
    '''
    def calculate_element_global_stiffness_matrix(self):

        reference_transform_matrix = self.get_reference_transform_matrix()
        element_local_stiffness_matrix = self.element_local_stiffness_matrix
        return reference_transform_matrix @ element_local_stiffness_matrix @ reference_transform_matrix.T

    '''
    element_global_internal_forces:

         element internal forces at global coordinates.
         12x1 array
    '''
    def calculate_element_global_internal_forces(self):

        reference_transform_matrix = self.get_reference_transform_matrix()

        return reference_transform_matrix @ self.element_resisting_forces


################################################################################
################################################################################
    '''
    The forth block is providing functions
    for detailed calculation.
    '''
################################################################################
    '''
    reference_transform_matrix:

         12x5 array

         Force balance at reference configuration.
         reference_transform_matrix is defined as T * Q = P
         Thus for displacement:
         q = T.transpose * p
    '''
    def get_reference_transform_matrix(self):

        L = self.get_reference_length()
        reference_triad = self.get_reference_triad()
        e1 = reference_triad[:, 0]
        e2 = reference_triad[:, 1]
        e3 = reference_triad[:, 2]
        reference_transform_matrix = np.zeros([12, 5])

        # Forces of first node
        reference_transform_matrix[0:3, 0] = e2 /L
        reference_transform_matrix[0:3, 1] = e2 /L
        reference_transform_matrix[0:3, 2] = -e3 /L
        reference_transform_matrix[0:3, 3] = -e3 /L
        reference_transform_matrix[0:3, 4] = -e1
        # Moments of first node
        reference_transform_matrix[3:6, 0] = e3
        reference_transform_matrix[3:6, 2] = e2
        # Forces of second node
        reference_transform_matrix[6:9, :] = - reference_transform_matrix[0:3, :]
        # Moments of second node
        reference_transform_matrix[9:12, 1] = e3
        reference_transform_matrix[9:12, 3] = e2

        return reference_transform_matrix

    '''
    reference_triad:

         3 unit vectors composing the reference local coordinate.
    '''
    def get_reference_triad(self):

        reference_x_a = self.Nd_a.get_reference_location()
        reference_x_b = self.Nd_b.get_reference_location()

        e1 = (reference_x_b - reference_x_a) / (la.norm(reference_x_b - reference_x_a) + 1e-50)
        e2 = np.cross(self.reference_abt_vec, e1) / la.norm(np.cross(self.reference_abt_vec, e1))
        e3 = np.cross(e1, e2)

        return np.column_stack((e1, e2, e3))

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

    def get_reference_length(self):

        reference_x_a = self.Nd_a.get_reference_location()
        reference_x_b = self.Nd_b.get_reference_location()
        return la.norm(reference_x_b - reference_x_a)
