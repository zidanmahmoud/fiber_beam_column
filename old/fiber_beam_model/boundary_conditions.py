import numpy as np
from numpy import linalg as la

class BoundaryConditions(object):
    """
    BoundaryConditions provides von Neumann and Dirichlet boundary conditions.

    Attributes
    ----------
    id : int

    given_displacement : dict

    given_force : dict
    """

    def __init__(self, dirichlet_bcs, controled_displacement_id, given_force):

        self.dirichlet_bcs = dirichlet_bcs
        self.controled_displacement_id = controled_displacement_id
        self.given_force = given_force



    def find_constraint_dof_id(self):

        return [index for index,id in enumerate(self.dirichlet_bcs) if id is '0']


    def apply_dirichlet_boundary_condition_to_K(self, K):

        ids = self.find_constraint_dof_id()
        K[ids, :] = 0
        K[:, ids] = 0
        K[ids, ids] = 1

        return K

    def apply_dirichlet_boundary_condition_to_F(self, F):

        ids = self.find_constraint_dof_id()
        F[ids] = 0

        return F

    def get_controled_displacement_id(self):

        return self.controled_displacement_id

    def get_controled_displacement(self):

        return self.dirichlet_bcs[self.controled_displacement_id]

    def new_loading(self, new_controled_displacement):

        self.dirichlet_bcs[self.controled_displacement_id] = new_controled_displacement
