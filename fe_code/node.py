"""
This module only contains the Node class.
"""

import numpy as np


class Node:
    """Three dimensional Node providing Dofs for displacements.

    Parameters
    ----------
    x : float
        Initial X coordinate of the node.
    y : float
        Initial Y coordinate of the node.
    z : float
        Initial Z coordinate of the node.

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

    def __init__(self, node_id, x, y, z):
        self._id = node_id
        self._x = x
        self._y = y
        self._z = z
        self.reference_x = x
        self.reference_y = y
        self.reference_z = z

    @property
    def id(self):
        """ id """
        return self._id

    @property
    def u(self):
        """ Displacement in x direction. """
        return self._x - self.reference_x

    @u.setter
    def u(self, value):
        """ Displacement in y direction. """
        self._x = self.reference_x + value

    @property
    def v(self):
        """ Displacement in y direction. """
        return self._y - self.reference_y

    @v.setter
    def v(self, value):
        self._y = self.reference_y + value

    @property
    def w(self):
        """ Displacement in z direction. """
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
