import numpy as np

class Node(object):
    """

    Three dimensional Node providing displacements.

    Attributes
    ----------
    id : int
        Unique ID.

    reference_x : ndarray
        Reference position
        3x1 array

    x : ndarray
        Actual position
        3x1 array

    node_displacements : ndarray
        6x1 array

    node_displacement_incr : ndarray
        6x1 array

    change_in_node_displacement_incr : ndarray
        6x1 array

    """

    def __init__(self, id, x):
        """

        Create a new node.

        Parameters
        ----------
        id : int or str
            Unique ID of the node.

        x : ndarray
            Initial position.

        """
        self.id = id
        self.x = x
        self.reference_x = x

        self.node_displacements = np.zeros(6)
        self.node_displacements_last_loadstep = np.zeros(6)
        self.node_displacement_incr = np.zeros(6)
        self.change_in_node_displacement_incr = np.zeros(6)

    # Change in displacement increments
    def get_change_in_node_displacement_incr(self):

        return self.change_in_node_displacement_incr

    def initialize_change_in_node_displacement_incr(self):

        self.change_in_node_displacement_incr = np.zeros(6)

    def update_change_in_node_displacement_incr(self, new_change_in_node_displacement_incr):

        self.change_in_node_displacement_incr = new_change_in_node_displacement_incr

    # Displacement increments
    def get_node_displacement_incr(self):

        return self.node_displacement_incr

    def initialize_node_displacement_incr(self):

        self.node_displacement_incr = np.zeros(6)

    def update_node_displacement_incr(self):

        self.node_displacement_incr += self.change_in_node_displacement_incr

    # Displacements
    def get_node_displacements(self):

        return node_displacements

    def initialize_node_displacements(self):

        self.node_displacements = np.zeros(6)

    def update_node_displacements(self):

        self.node_displacements = self.node_displacements_last_loadstep + self.node_displacement_incr

    def save_to_last_loadstep(self, load_step_convergence):
        if load_step_convergence:
            self.node_displacements_last_loadstep = self.node_displacements


    # Get position
    def get_reference_location(self):

        return self.reference_x

    def get_actual_location(self):

        return self.x
