"""
Example
"""

import plotting as p
from fe_code import io
from models.column import *
from disp_calc import *


def advance_in_load(structure, load_step):
    """ load stepping loop """

    if load_step < STEPS[0]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[1]:
        structure.controlled_dof_increment += -STEP
    elif load_step < STEPS[2]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[3]:
        structure.controlled_dof_increment += -STEP
    elif load_step < STEPS[4]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[5]:
        structure.controlled_dof_increment += -STEP
    elif load_step < STEPS[6]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[7]:
        structure.controlled_dof_increment += -STEP
    elif load_step < STEPS[8]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[9]:
        structure.controlled_dof_increment += -STEP
    elif load_step < STEPS[10]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[11]:
        structure.controlled_dof_increment += -STEP
    elif load_step < STEPS[12]:
        structure.controlled_dof_increment += STEP
    elif load_step < STEPS[13]:
        structure.controlled_dof_increment += -STEP
    else:
        structure.controlled_dof_increment += STEP


def solution_loop(structure, result_filename="results.dat"):
    max_nr_iterations = 10
    max_ele_iterations = 100

    structure.initialize()
    print(":: Initialized the solver ::")
    print("\n:: Starting solution loop ::")

    with open(result_filename, "w") as wfile:
        wfile.write("#DISPLACEMENT LOAD\n")
        wfile.write("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")

        for k in range(1, STEPS[-1]):
            print(f"\nLOAD STEP : {k}")
            advance_in_load(structure, k)

            for i in range(1, max_nr_iterations + 1):
                convergence, residual = structure.solve_NR_iteration(max_ele_iterations)
                if convergence:
                    print(f"NR converged with {i} iteration(s). Residual = {residual}")
                    break
                if i == max_nr_iterations:
                    io.warning(f"Newton-Raphson did not converge {max_nr_iterations} iterations")
                    io.warning("FATAL ERROR: The solution is unstable")
                    break
            if i == max_nr_iterations:
                break

            structure.finalize_load_step()

            displacements = structure.get_displacements()
            loads = structure.get_forces()
            for disp in displacements:
                wfile.write(str(disp)+" ")
            for load in loads:
                wfile.write(str(load)+" ")
            wfile.write("\n")

    print("\n:: Finished solution loop ::")


if __name__ == "__main__":
    STEP = 0.1
    STEPS = calculate_loadsteps(STEP)
    STRUCTURE = model1_4()
    # p.plot_disctrized_2d(STRUCTURE.get_element(1).get_section(1))

    solution_loop(STRUCTURE)

    # import matplotlib.pyplot as plt
    # data = np.loadtxt("results.dat", usecols=(8, 14))
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(data[:,0], -data[:,1])
    # ax.grid()
    # plt.show()
