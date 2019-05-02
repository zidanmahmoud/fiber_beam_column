"""This module only contains the Node class.

Author: Thomas Oberbichler
"""

import numpy as np

class Node(object):
    """Three dimensional Node providing Dofs for displacements.

    Attributes
    ----------
    reference_x : float
        Reference X coordinate.
    reference_y : float
        Reference Y coordinate.
    reference_z : float
        Reference Z coordinate.
    x : float
        Actual X coordinate.
    y : float
        Actual Y coordinate.
    z : float
        Actual Z coordinate.
    u : float
        Displacement in x direction.
    v : float
        Displacement in y direction.
    w : float
        Displacement in z direction.
    """

    def __init__(self, x, y, z):
        """Create a new node.

        Parameters
        ----------
        x : float
            Initial X coordinate of the node.
        y : float
            Initial Y coordinate of the node.
        z : float
            Initial Z coordinate of the node.
        """
        self._x = x
        self._y = y
        self._z = z
        self.reference_x = x
        self.reference_y = y
        self.reference_z = z

    @property
    def u(self):
        return self._x - self.reference_x

    @u.setter
    def u(self, value):
        self._x = self.reference_x + value

    @property
    def v(self):
        return self._y - self.reference_y

    @v.setter
    def v(self, value):
        self._y = self.reference_y + value

    @property
    def w(self):
        return self._z - self.reference_z

    @w.setter
    def w(self, value):
        self._z = self.reference_z + value

    def get_reference_location(self):
        """Location of the node in the reference configuration.

        Returns
        -------
        location : ndarray
            Numpy array containing the reference coordinates X, Y and Z.
        """
        _x = self.reference_x
        _y = self.reference_y
        _z = self.reference_z

        return np.array([_x, _y, _z], dtype=float)

    def get_actual_location(self):
        """Location of the node in the actual configuration.

        Returns
        -------
        location : ndarray
            Numpy array containing the actual coordinates X, Y and Z.
        """
        _x = self._x
        _y = self._y
        _z = self._z

        return np.array([_x, _y, _z], dtype=float)

    def get_displacement(self):
        """Displacement of the node in the actual configuration.

        Returns
        -------
        displacement : ndarray
            A numpy array containing the displacements u, v and w.
        """
        return self.get_reference_location() - self.get_actual_location()

    def get_dof_state(self, dof_type):
        """Get the current value of the given dof type.

        Parameters
        ----------
        dof_type : string
            Type of the dof.

        Returns
        -------
        value : float
            The current value of the dof type

        Raises
        ------
        AttributeError
            If `dof_type` does not exist.
        """
        if dof_type == 'u':
            return self.u
        if dof_type == 'v':
            return self.v
        if dof_type == 'w':
            return self.w

        raise AttributeError(f"Node has no dof of type \'{dof_type}\'")

    def set_dof_state(self, dof_type, value):
        """Update the node according to the value of the given dof type.

        Parameters
        ----------
        dof_type : string
            Type of the Dof.
        value : float
            The value of the given dof.

        Raises
        ------
        AttributeError
            If `dof_type` does not exist.
        """
        if dof_type == 'u':
            self.u = value
        elif dof_type == 'v':
            self.v = value
        elif dof_type == 'w':
            self.w = value
        else:
            raise AttributeError(f"Node has no dof of type \'{dof_type}\'")
