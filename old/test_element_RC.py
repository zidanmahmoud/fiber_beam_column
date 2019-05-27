import sys
sys.path.append('C:/Users/Fei/OngoingProjects/Force_based_fiber_beam_element/small_disp_fiberbeam/test_code_1/fiber_beam_model')
# import necessary modules
import numpy as np
from fiber_beam_model import *
from pprint import pprint
import integral_scheme as intg
import matplotlib.pyplot as plt

H = 1.
B = 1.
NSx = 4
NFz = 10
NFy = 10
z_st = H/2*(1/NFz - 1)
y_st = B/2*(1/NFy - 1)
dz = H/NFz
dy = B/NFy
[gauss_points, weights] = intg.GaussLobatto(NSx)
fiber_list = [[[Fiber() for nfz in range(NFz)] for nfy in range(NFy)] for nsx in range(NSx)]
section_list = []
for nsx in range(NSx):
    for nfz in range(NFz):
        for nfy in range(NFy):
            if (nfz == 2 or nfz == 7) and (nfy == 2 or nfy == 7):
                y = y_st + B/NFy*nfy
                z = z_st + H/NFz*nfz
                fiber_list[nsx][nfy][nfz] = Fiber([y, z], 1/NFz/NFy, MenegottoPintoModel(1, 29000, 580, 60, 20, 18.5, 0.0002))
            else:
                y = y_st + B/NFy*nfy
                z = z_st + H/NFz*nfz
                fiber_list[nsx][nfy][nfz] = Fiber([y, z], 1/NFz/NFy, KentParkModel(1, 41.37, 1, -0.0096, 128.205))
    section_list.append(Section(gauss_points[nsx], [1., 1.], fiber_list[nsx][:][:], 1e-5))


Nd_a = Node(1, np.array([0,0,0]))
Nd_b = Node(1, np.array([1,0,0]))
fbe = FiberBeam(1, Nd_a, Nd_b, np.array([0,0,1]), section_list)
Nd_a.update_change_in_node_displacement_incr(np.array([0,0,0,0,0.0005,0]))
Nd_b.update_change_in_node_displacement_incr(np.array([0,0,0,0,0.001,0]))
fbe.initialize_in_loadstep()
#fbe.update_section_parameters_in_loadstep(0)
print('step3. element_local_stiffness_matrix = ', fbe.get_element_local_stiffness_matrix())
fbe.update_change_in_element_disp_incr()
fbe.update_element_disp_incr()
for J in range(1, 100):
    print('#####################################################################')
    print('J = ', J)
    print('#####################################################################')
    fbe.update_change_in_element_force_incr(J)
    fbe.update_element_force_incr()
    print('step7. element_force_incr = ', fbe.get_element_force_incr())
    fbe.update_element_resisting_forces()
    print('step7. element_resisting_forces = ', fbe.get_element_resisting_forces())
    fbe.update_section_parameters_in_loadstep(J)
    fbe.update_element_local_stiffness_matrix()
    if fbe.check_for_element_convergence():
        print('J = ', J)
        print('converged')
        break
    else:
        fbe.update_residual_element_disps()
        print('step17. residual_element_disps = ', fbe.get_residual_element_disps())
    print('step16. element_local_stiffness_matrix = ', fbe.get_element_local_stiffness_matrix())
plt.grid(True)
plt.show()
