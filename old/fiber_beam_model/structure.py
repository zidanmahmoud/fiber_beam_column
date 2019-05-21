import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

class Structure(object):

    """
    Three dimensional structure providing:
    global_stiffness_matrix
    global_internal_forces
    convergence

    Attention:
    This is the 3d beam element for infinitesimal non-torsional deformation.


    Attributes
    ----------


    """

    def __init__(self, id, nodes_list, elements_list, element_node_map, length_incr, bd_conditions, solving_strategy, tolerance):
        """
        Creat a new structure.

        Parameters
        ----------
        id : int

        nodes_array: ndarray

        elements_array: ndarray

        node_element_map : ndarray

        solving_strategy : str
        """

        self.id = id
        self.nodes_list = nodes_list
        self.elements_list = elements_list
        self.element_node_map = element_node_map
        self.length_incr = length_incr #Question: what is length incr?
        self.bd_conditions = bd_conditions
        self.solving_strategy = solving_strategy
        self.tolerance = tolerance

        self.DoF = 6*len(nodes_list)

        self.change_in_structural_displacement_incr = np.zeros(self.DoF)
        self.structural_displacement_incr = np.zeros(self.DoF)
        self.structural_displacement = np.zeros(self.DoF)
        self.structural_displacement_last_loadstep = np.zeros(self.DoF)

        self.change_in_load_factor_incr = 0
        self.load_factor_incr = 0
        self.load_factor = 0
        self.load_factor_last_loadstep = 0

        self.structural_stiffness_matrix = np.zeros([self.DoF, self.DoF])
        self.structural_external_force = np.zeros(self.DoF)
        self.structural_resisting_force = np.zeros(self.DoF)
        self.structural_unbalanced_force = np.zeros(self.DoF)

        self.controled_displacement_id = self.bd_conditions.get_controled_displacement_id()
        self.load_vector = np.zeros(self.DoF)
        self.load_vector[self.controled_displacement_id] = 1

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
    def new_loading(self, new_length_incr, new_controled_displacement):

        self.bd_conditions.new_loading(new_controled_displacement)
        self.length_incr = new_length_incr

    '''
    Step 3.1:
        solve the global system of equations
        update the structure displacement increments
    '''
    def initialize_at_beginning(self):
        for element in self.elements_list:
            element.initialize_at_beginning()
        self.update_structural_stiffness_matrix()

    def initialize_in_loadstep(self):

        self.initialize_structural_displacement_incr()
        self.initialize_load_factor_incr()

        elements_list = self.elements_list
        for ne in range(len(elements_list)):
            single_element = elements_list[ne]
            single_element.initialize_in_loadstep()
            elements_list[ne] = single_element
        self.elements_list = elements_list

        nodes_list = self.nodes_list
        for nn in range(len(nodes_list)):
            single_node = nodes_list[nn]
            single_node.initialize_node_displacement_incr()
            nodes_list[nn] = single_node
        self.nodes_list = nodes_list

    '''
    Step 3.2:
        solve the global system of equations
        update the structure displacement increments
    '''
    def solver_for_displacement_control(self):

        structural_stiffness_matrix = self.bd_conditions.apply_dirichlet_boundary_condition_to_K(self.structural_stiffness_matrix)
        structural_unbalanced_force = self.structural_unbalanced_force
        load_vector = self.load_vector
        structural_displacement_incr = self.structural_displacement_incr

        modified_K = np.zeros([len(load_vector) + 1, len(load_vector) + 1])
        modified_K[0:len(load_vector), 0:len(load_vector)] = structural_stiffness_matrix
        modified_K[0:len(load_vector), len(load_vector)] = -load_vector
        modified_K[len(load_vector), 0:len(load_vector)] = -load_vector
        modified_F = np.zeros(len(load_vector) + 1)
        modified_F[0:len(load_vector)] = structural_unbalanced_force
        modified_F[len(load_vector)] = load_vector @ structural_displacement_incr - self.length_incr

        solved_change_in_incr = la.inv(modified_K) @ modified_F
        self.change_in_structural_displacement_incr = solved_change_in_incr[0:self.DoF]
        self.change_in_load_factor_incr = solved_change_in_incr[self.DoF]

    def update_node_parameters_in_loadstep(self):

        nodes_list = self.nodes_list
        change_in_structural_displacement_incr = self.change_in_structural_displacement_incr
        for index, node in enumerate(nodes_list):
            i = 6*index
            node.update_change_in_node_displacement_incr(change_in_structural_displacement_incr[i:i+6])
            node.update_node_displacement_incr()
            node.update_node_displacements()

    def update_element_parameters_in_loadstep(self, load_step_convergence, i, k):

        for element in self.elements_list:
            for j in range(1, 100+1):

                # print('j = ', j)
                # print('---------')
                element.update_change_in_element_disp_incr() # Question: shouldn't this be outside the j-loop?
                element.update_element_disp_incr()
                element.update_element_disp() # Question: this is not part of the thesis' steps..
                reverse = self.execute_element_state_determination(element, load_step_convergence , j, i, k)
                if element.check_for_element_convergence():
                    print('**** Element converged, j = ', j)
                    element.save_to_NRstep()
                    break
                else:
                    element.update_residual_element_disps()
                if j == 100:
                    print("WARNING: Element did not converge!!")

        return reverse

    def save_to_last_loadstep(self, load_step_convergence):

        if load_step_convergence:

            elements_list = self.elements_list
            for ne in range(len(elements_list)):
                single_element = elements_list[ne]
                single_element.save_to_last_loadstep(load_step_convergence)
                elements_list[ne] = single_element
            self.elements_list = elements_list

            nodes_list = self.nodes_list
            for nn in range(len(nodes_list)):
                single_node = nodes_list[nn]
                single_node.save_to_last_loadstep(load_step_convergence)
                nodes_list[nn] = single_node
            self.nodes_list = nodes_list

            self.structural_displacement_last_loadstep = self.structural_displacement
            self.load_factor_last_loadstep = self.load_factor




    '''
    structural_displacement_incr :

         According to STEP(3.3)
         1. If i == 1: it is initialized to be zero.
         2. If i > 1: it is updated
    '''
    def get_structural_displacement_incr(self):

        return self.structural_displacement_incr

    def initialize_structural_displacement_incr(self):

        self.structural_displacement_incr = np.zeros(self.DoF)

    def update_structural_displacement_incr(self):

        self.structural_displacement_incr += self.change_in_structural_displacement_incr

    def get_structural_displacement(self):

        return self.structural_displacement

    def update_structural_displacement(self):

        self.structural_displacement = self.structural_displacement_last_loadstep + self.structural_displacement_incr

    '''
    load_factor_incr :

         1. If i == 1: it is initialized to be zero.
         2. If i > 1: it is updated
    '''
    def get_load_factor_incr(self):

        return self.load_factor_incr

    def initialize_load_factor_incr(self):

        self.load_factor_incr = 0

    def update_load_factor_incr(self):

        self.load_factor_incr += self.change_in_load_factor_incr

    def get_load_factor(self):

        return self.load_factor

    def update_load_factor(self):

        self.load_factor = self.load_factor_last_loadstep + self.load_factor_incr

    '''
    structural_stiffness_matrix :

         According to STEP(3.2)

    '''
    def get_structural_stiffness_matrix(self):

        return self.structural_stiffness_matrix

    def update_structural_stiffness_matrix(self):

        self.structural_stiffness_matrix = self.calculate_structural_stiffness_matrix()

    '''
    structural_resisting_force :

         According to STEP(6.1)

    '''
    def get_structural_resisting_force(self):

        return self.structural_resisting_force

    def update_structural_resisting_force(self):

        self.structural_resisting_force = self.calculate_structural_resisting_force()

    '''
    structural_external_force :

         According to STEP(6.1)

    '''
    def get_structural_external_force(self):

        return self.structural_external_force

    def update_structural_external_force(self):

        self.structural_external_force = self.calculate_structural_external_force()

    '''
    structural_unbalanced_force :

         According to STEP(6.1)

    '''
    def get_structural_unbalanced_force(self):

        return self.structural_unbalanced_force

    def update_structural_unbalanced_force(self):

        self.structural_unbalanced_force = self.calculate_structural_unbalanced_force()

    def check_for_structural_convergence(self):
        if la.norm(self.structural_unbalanced_force) < self.tolerance:
            return True


################################################################################
################################################################################
    '''
    The second block provides :
    functions which are needed in the first block.
    '''
################################################################################
    '''
    structural_stiffness_matrix:

    '''
    def calculate_structural_stiffness_matrix(self):

        elements_list = self.elements_list
        element_node_map = self.element_node_map
        structural_stiffness_matrix = np.zeros([self.DoF, self.DoF])

        for ne in range(len(elements_list)):
            single_element = elements_list[ne]
            id1 = element_node_map[0, ne]
            id2 = element_node_map[1, ne]
            KE = single_element.calculate_element_global_stiffness_matrix()
            structural_stiffness_matrix[6*id1:6*(id1+1), 6*id1:6*(id1+1)] += KE[0:6, 0:6]
            structural_stiffness_matrix[6*id2:6*(id2+1), 6*id2:6*(id2+1)] += KE[6:12, 6:12]
            structural_stiffness_matrix[6*id1:6*(id1+1), 6*id2:6*(id2+1)] += KE[0:6, 6:12]
            structural_stiffness_matrix[6*id2:6*(id2+1), 6*id1:6*(id1+1)] += KE[6:12, 0:6]

        return KE

    def calculate_structural_resisting_force(self):

        element_node_map = self.element_node_map
        structural_resisting_force = np.zeros(self.DoF)

        for i, element in enumerate(self.elements_list):
            id1 = element_node_map[0, i]
            id2 = element_node_map[1, i]
            FE = element.calculate_element_global_internal_forces()
            structural_resisting_force[6*id1:6*(id1+1)] += FE[0:6]
            structural_resisting_force[6*id2:6*(id2+1)] += FE[6:12]

        return structural_resisting_force

    def calculate_structural_external_force(self):
        return self.load_factor * self.load_vector

    def calculate_structural_unbalanced_force(self):

        return self.bd_conditions.apply_dirichlet_boundary_condition_to_F(self.structural_external_force - self.structural_resisting_force)

    def execute_element_state_determination(self, single_element, load_step_convergence, j,  i, k):

        single_element.update_change_in_element_force_incr(j)
        single_element.update_element_force_incr()
        single_element.update_element_resisting_forces()
        print(single_element.element_resisting_forces)
        input()
        reverse = single_element.update_section_parameters_in_loadstep(load_step_convergence, j, i, k)
        single_element.update_element_local_stiffness_matrix()

        return reverse
