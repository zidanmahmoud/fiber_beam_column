import sys
sys.path.append("./fiber_beam_model")
# import necessary modules
import numpy as np
np.set_printoptions(linewidth=190)
import fiber_beam_model as fbm
import matplotlib.pyplot as plt
import time
start = time.time()


H = 8
B = 5.
NSx = 4
NFz = 15
NFy = 15
z_st = H/2*(1/NFz - 1)
y_st = B/2*(1/NFy - 1)
dz = H/NFz
dy = B/NFy
[gauss_points, weights] = fbm.GaussLobatto(NSx)
fiber_list = [[[fbm.Fiber() for nfz in range(NFz)] for nfy in range(NFy)] for nsx in range(NSx)]
section_list = []

for nsx in range(NSx):
    for nfz in range(NFz):
        for nfy in range(NFy):
            if (nfz == 1 or nfz == 13) and (nfy == 1 or nfy == 13):
                # Steel fibers
                y = y_st + B/NFy*nfy
                z = z_st + H/NFz*nfz
                fiber_list[nsx][nfy][nfz] = fbm.Fiber([y, z], H*B/NFz/NFy, fbm.MenegottoPintoModel(1, 29000, 0.0042*29000, 60, 20, 18.5, 0.0002))
            else:
                # Concrete fibers
                y = y_st + B/NFy*nfy
                z = z_st + H/NFz*nfz
                fiber_list[nsx][nfy][nfz] = fbm.Fiber([y, z], H*B/NFz/NFy, fbm.KentParkModel(1, 6.95, 1, -0.07, 770))
    section_list.append(fbm.Section(gauss_points[nsx], weights[nsx], [1., 1.], fiber_list[nsx][:][:], 0.05))

nodes_list = []
nodes_list.append(fbm.Node(0, np.array([0,0,0])))
nodes_list.append(fbm.Node(1, np.array([100,0,0])))
elements_list = []
fbe = fbm.FiberBeam(1, nodes_list[0], nodes_list[1], np.array([0,0,1]), section_list)
elements_list.append(fbe)
element_node_map = np.array([[0], [1]])
disp_incr = 0.4
disp_max = 0.005
step_n = 10
dirichlet_bcs = ['0', '0', '0', '0', '0', '0', 'f', 'f', disp_max, '0', 'f', 'f']
controled_displacement_id = 8
given_force = None
bd_conditions = fbm.BoundaryConditions(dirichlet_bcs, controled_displacement_id, given_force)
st = fbm.Structure(1, nodes_list, elements_list, element_node_map, disp_incr, bd_conditions, 'displacement_control', 0.05)

Fl = []
ul = []
Fl.append(0)
ul.append(0)
lamda_array = []


st.initialize_at_beginning()

k = 1
while (k <= step_n):
    print('k =', k)
    print("-----")
    for i in range(1,100+1):
        if i == 1 and k == 1:
            do_reversal = True
        else:
            do_reversal = False
        # print('i =', i)
        # print("-------")
        st.solver_for_displacement_control()
        st.update_node_parameters_in_loadstep()

        st.update_structural_displacement_incr()
        st.update_structural_displacement()
        st.update_load_factor_incr()
        st.update_load_factor()

        reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

        st.update_structural_stiffness_matrix()
        st.update_structural_resisting_force()
        #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
        st.update_structural_external_force()
        #print('step18. structural_external_force = ', st.get_structural_external_force())
        st.update_structural_unbalanced_force()
        #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

        if st.check_for_structural_convergence():
            print('**** NR converged at i =', i)
            st.save_to_last_loadstep(True)
            break
        if i == 99:
            print("WARNING: NewtonRaphson iteration did not converge!!")
    k += 1
    external_force = st.get_structural_external_force()
    structural_displacement = st.get_structural_displacement()
    Fl.append(external_force[8])
    ul.append(structural_displacement[8])




while (k <= step_n + 20 and k > step_n):
    if k == step_n + 1:
        st.new_loading(-0.01*disp_incr, 0)
    elif k > step_n + 1:
        st.new_loading(-disp_incr, 0)
    st.initialize_in_loadstep()
    print('k = ', k)
    print("-----")
    for i in range(1,100):
        if i == 1 and k == step_n + 1:
            do_reversal = True
        else:
            do_reversal = False
        # print('i = ', i)
        # print("-------")
        st.solver_for_displacement_control()
        st.update_node_parameters_in_loadstep()

        st.update_structural_displacement_incr()
        st.update_structural_displacement()
        st.update_load_factor_incr()
        st.update_load_factor()

        reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

        st.update_structural_stiffness_matrix()
        st.update_structural_resisting_force()
        #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
        st.update_structural_external_force()
        #print('step18. structural_external_force = ', st.get_structural_external_force())
        st.update_structural_unbalanced_force()
        #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

        if st.check_for_structural_convergence():
            print('NR converged at i = ', i)
            st.save_to_last_loadstep(True)
            break
        if i == 99:
            print("WARNING: NewtonRaphson iteration did not converge!!")
    k += 1
    external_force = st.get_structural_external_force()
    structural_displacement = st.get_structural_displacement()
    Fl.append(external_force[8])
    ul.append(structural_displacement[8])


while (k <= step_n + 45 and k > step_n + 20):
    if k == step_n + 21:
        st.new_loading(0.01*disp_incr, 0)
    elif k > step_n + 21:
        st.new_loading(disp_incr, 0)
    st.initialize_in_loadstep()
    print('k = ', k)
    print("-----")
    for i in range(1,100):
        if i == 1 and k == step_n + 21:
            do_reversal = True
        else:
            do_reversal = False
        # print('i = ', i)
        # print("-------")
        st.solver_for_displacement_control()
        st.update_node_parameters_in_loadstep()

        st.update_structural_displacement_incr()
        st.update_structural_displacement()
        st.update_load_factor_incr()
        st.update_load_factor()

        reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

        st.update_structural_stiffness_matrix()
        st.update_structural_resisting_force()
        #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
        st.update_structural_external_force()
        #print('step18. structural_external_force = ', st.get_structural_external_force())
        st.update_structural_unbalanced_force()
        #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

        if st.check_for_structural_convergence():
            print('NR converged at i = ', i)
            st.save_to_last_loadstep(True)
            break
        if i == 99:
            print("WARNING: NewtonRaphson iteration did not converge!!")
    k += 1
    external_force = st.get_structural_external_force()
    structural_displacement = st.get_structural_displacement()
    Fl.append(external_force[8])
    ul.append(structural_displacement[8])

# while (k <= step_n + 75 and k > step_n + 45):
#     if k == step_n + 46:
#         st.new_loading(-0.01*disp_incr, 0)
#     elif k > step_n + 46:
#         st.new_loading(-disp_incr, 0)
#     st.initialize_in_loadstep()
#     print('k = ', k)
#     print("-----")
#     for i in range(1,100):
#         if i == 1 and k == step_n + 46:
#             do_reversal = True
#         else:
#             do_reversal = False
#         # print('i = ', i)
#         # print("-------")
#         st.solver_for_displacement_control()
#         st.update_node_parameters_in_loadstep()

#         st.update_structural_displacement_incr()
#         st.update_structural_displacement()
#         st.update_load_factor_incr()
#         st.update_load_factor()

#         reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

#         st.update_structural_stiffness_matrix()
#         st.update_structural_resisting_force()
#         #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
#         st.update_structural_external_force()
#         #print('step18. structural_external_force = ', st.get_structural_external_force())
#         st.update_structural_unbalanced_force()
#         #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

#         if st.check_for_structural_convergence():
#             print('NR converged at i = ', i)
#             st.save_to_last_loadstep(True)
#             break
#         if i == 99:
#             print("WARNING: NewtonRaphson iteration did not converge!!")
#     k += 1
#     external_force = st.get_structural_external_force()
#     structural_displacement = st.get_structural_displacement()
#     Fl.append(external_force[8])
#     ul.append(structural_displacement[8])

# while (k <= step_n + 91 and k > step_n + 75):
#     if k == step_n + 76:
#         st.new_loading(0.01*disp_incr, 0)
#     elif k > step_n + 76:
#         st.new_loading(disp_incr, 0)
#     st.initialize_in_loadstep()
#     print('k = ', k)
#     print("-----")
#     for i in range(1,100):
#         if i == 1 and k == step_n + 76:
#             do_reversal = True
#         else:
#             do_reversal = False
#         print('i = ', i)
#         print("-------")
#         st.solver_for_displacement_control()
#         st.update_node_parameters_in_loadstep()

#         st.update_structural_displacement_incr()
#         st.update_structural_displacement()
#         st.update_load_factor_incr()
#         st.update_load_factor()

#         reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

#         st.update_structural_stiffness_matrix()
#         st.update_structural_resisting_force()
#         #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
#         st.update_structural_external_force()
#         #print('step18. structural_external_force = ', st.get_structural_external_force())
#         st.update_structural_unbalanced_force()
#         #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

#         if st.check_for_structural_convergence():
#             print('NR converged at i = ', i)
#             st.save_to_last_loadstep(True)
#             break
#         if i == 99:
#             print("WARNING: NewtonRaphson iteration did not converge!!")
#     k += 1
#     external_force = st.get_structural_external_force()
#     structural_displacement = st.get_structural_displacement()
#     Fl.append(external_force[8])
#     ul.append(structural_displacement[8])

# while (k <= step_n + 102 and k > step_n + 91):
#     if k == step_n + 92:
#         st.new_loading(-0.01*disp_incr, 0)
#     elif k > step_n + 92:
#         st.new_loading(-disp_incr, 0)
#     st.initialize_in_loadstep()
#     print('k = ', k)
#     print("-----")
#     for i in range(1,100):
#         if i == 1 and k == step_n + 92:
#             do_reversal = True
#         else:
#             do_reversal = False
#         # print('i = ', i)
#         # print("-------")
#         st.solver_for_displacement_control()
#         st.update_node_parameters_in_loadstep()

#         st.update_structural_displacement_incr()
#         st.update_structural_displacement()
#         st.update_load_factor_incr()
#         st.update_load_factor()

#         reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

#         st.update_structural_stiffness_matrix()
#         st.update_structural_resisting_force()
#         #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
#         st.update_structural_external_force()
#         #print('step18. structural_external_force = ', st.get_structural_external_force())
#         st.update_structural_unbalanced_force()
#         #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

#         if st.check_for_structural_convergence():
#             print('NR converged at i = ', i)
#             st.save_to_last_loadstep(True)
#             break
#         if i == 99:
#             print("WARNING: NewtonRaphson iteration did not converge!!")
#     k += 1
#     external_force = st.get_structural_external_force()
#     structural_displacement = st.get_structural_displacement()
#     Fl.append(external_force[8])
#     ul.append(structural_displacement[8])

# while (k <= step_n + 107 and k > step_n + 102):
#     if k == step_n + 103:
#         st.new_loading(0.01*disp_incr, 0)
#     elif k > step_n + 103:
#         st.new_loading(disp_incr, 0)
#     st.initialize_in_loadstep()
#     print('k = ', k)
#     print("-----")
#     for i in range(1,100):
#         if i == 1 and k == step_n + 103:
#             do_reversal = True
#         else:
#             do_reversal = False
#         # print('i = ', i)
#         # print("-------")
#         st.solver_for_displacement_control()
#         st.update_node_parameters_in_loadstep()

#         st.update_structural_displacement_incr()
#         st.update_structural_displacement()
#         st.update_load_factor_incr()
#         st.update_load_factor()

#         reverse = st.update_element_parameters_in_loadstep(do_reversal, i, k)

#         st.update_structural_stiffness_matrix()
#         st.update_structural_resisting_force()
#         #print('step18. structural_resisting_force = ', st.get_structural_resisting_force())
#         st.update_structural_external_force()
#         #print('step18. structural_external_force = ', st.get_structural_external_force())
#         st.update_structural_unbalanced_force()
#         #print('step18. structural_unbalanced_force = ', st.get_structural_unbalanced_force())

#         if st.check_for_structural_convergence():
#             print('NR converged at i = ', i)
#             st.save_to_last_loadstep(True)
#             break
#         if i == 99:
#             print("WARNING: NewtonRaphson iteration did not converge!!")
#     k += 1
#     external_force = st.get_structural_external_force()
#     structural_displacement = st.get_structural_displacement()
#     Fl.append(external_force[8])
#     ul.append(structural_displacement[8])


def mp(alist, num):
    for i in range(len(alist)):
        alist[i] *= num
    return alist


fig, ax = plt.subplots()
ax.plot(0, 0, 'yo')
ax.plot(mp(ul, 1/10000*3), mp(Fl, 100))
ax.set(
    xlabel='curvature (rad/in)',
    ylabel='moment (kip*in)'
)
ax.grid(True)
plt.savefig("blah.png")
plt.show()

end = time.time()
print(f"\n\nExecution time: {end-start} seconds.")
