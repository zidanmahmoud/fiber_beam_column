import sys
sys.path.append('C:/Users/Fei/OngoingProjects/Force_based_fiber_beam_element/small_disp_fiberbeam/test_code_1/fiber_beam_model')
# import necessary modules
import numpy as np
from fiber_beam_model import *


dirichlet_bcs = ['0', '0', '0', '0', '0', '0', 'f', 'f', 0.1, 'f', 'f', 'f']
controled_displacement_id = 8
given_force = None
bd_conditions = BoundaryConditions(dirichlet_bcs, controled_displacement_id, given_force)

print(bd_conditions.get_controled_displacement())
print(bd_conditions.get_controled_displacement_id())
print(bd_conditions.apply_dirichlet_boundary_condition_to_K(np.ones([12, 12])))
print(bd_conditions.apply_dirichlet_boundary_condition_to_F(np.ones(12)))
