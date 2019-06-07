"""
Module contains the class dof
"""


class DoF:
    """
    Degree of freedom

    Attributes
    ----------
    node_id : int

    type : char
        dof_type
    """

    def __init__(self, node_id, dof_type):
        self._node_id = node_id
        self._dof_type = dof_type

    @property
    def node_id(self):
        """ node id """
        return self._node_id

    @property
    def type(self):
        """ dof_type """
        return self._dof_type

    def __repr__(self):
        return f"({self._node_id}, {self._dof_type})"
