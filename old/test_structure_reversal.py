import sys
sys.path.append('C:/Users/Fei/OngoingProjects/Force_based_fiber_beam_element/small_disp_fiberbeam/test_code_1/fiber_beam_model')
# import necessary modules
import numpy as np
from fiber_beam_model import *
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
            y = y_st + B/NFy*nfy
            z = z_st + H/NFz*nfz
            fiber_list[nsx][nfy][nfz] = Fiber([y, z], 1/NFz/NFy, MenegottoPintoModel(1, 20000, 580, 60, 20, 18.5, 0.0002))
    section_list.append(Section(gauss_points[nsx], [1., 1.], fiber_list[nsx][:][:], 1e-7))
nodes_list = []
nodes_list.append(Node(0, np.array([0,0,0])))
nodes_list.append(Node(1, np.array([1,0,0])))
elements_list = []
fbe = FiberBeam(1, nodes_list[0], nodes_list[1], np.array([0,0,1]), section_list)
elements_list.append(fbe)
element_node_map = np.array([[0], [1]])
disp_incr = 0.0005
disp_max = 0.005
dirichlet_bcs = ['0', '0', '0', '0', '0', '0', 'f', 'f', disp_max, '0', 'f', 'f']
controled_displacement_id = 8
given_force = None
bd_conditions = BoundaryConditions(dirichlet_bcs, controled_displacement_id, given_force)
st = Structure(1, nodes_list, elements_list, element_node_map, disp_incr, bd_conditions, 'displacement_control', 1e-5)

st.initialize_at_beginning()
plt.plot(0, 0, 'yo')



k = 1
while (k <= 10):
    st.initialize_in_loadstep()
    print('#####################################################################')
    print('k = ', k)
    print('#####################################################################')
    for i in range(1,100):
        if i == 1:
            load_step_convergence = True
        else:
            load_step_convergence = False
        print('#####################################################################')
        print('i = ', i)
        print('#####################################################################')
        st.solver_for_displacement_control()
        st.update_change_in_structural_displacement_incr()
        st.update_change_in_load_factor_incr()
        st.update_node_parameters_in_loadstep()

        st.update_structural_displacement_incr()
        st.update_structural_displacement()
        st.update_load_factor_incr()
        st.update_load_factor()
        reverse = st.update_element_parameters_in_loadstep(load_step_convergence)
        if reverse:
            k-=1
            break

        st.update_structural_stiffness_matrix()
        st.update_structural_resisting_force()
        #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
        st.update_structural_external_force()
        #print('step18. structural_external_force = ', st.get_structural_external_force())
        st.update_structural_unbalanced_force()
        #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

        if st.check_for_structural_convergence():
            print('###########################################################################')
            print('NR converged at i = ', i)
            st.save_to_last_loadstep(True)
            break
    k += 1

st.new_loading(-0.002*disp_incr, 0)

while (k <= 12):
    st.initialize_in_loadstep()
    print('#####################################################################')
    print('k = ', k)
    print('#####################################################################')
    for i in range(1,100):
        if i == 1:
            load_step_convergence = True
        else:
            load_step_convergence = False
        print('#####################################################################')
        print('i = ', i)
        print('#####################################################################')
        st.solver_for_displacement_control()
        st.update_change_in_structural_displacement_incr()
        st.update_change_in_load_factor_incr()
        st.update_node_parameters_in_loadstep()

        st.update_structural_displacement_incr()
        st.update_structural_displacement()
        st.update_load_factor_incr()
        st.update_load_factor()
        reverse = st.update_element_parameters_in_loadstep(load_step_convergence)
        if reverse:
            k-=1
            break

        st.update_structural_stiffness_matrix()
        st.update_structural_resisting_force()
        #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
        st.update_structural_external_force()
        #print('step18. structural_external_force = ', st.get_structural_external_force())
        st.update_structural_unbalanced_force()
        #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

        if st.check_for_structural_convergence():
            print('###########################################################################')
            print('NR converged at i = ', i)
            st.save_to_last_loadstep(True)
            break
    k += 1
plt.grid(True)
plt.show()
