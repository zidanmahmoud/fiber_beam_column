"""
Example
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from fe_code import io, Structure, MenegottoPinto, KentPark

PLOT_FLAG = True
SAVE_PLOT = True

# == MODELING PARAMETERS
LENGTH = 100.0
WIDTH = 5.0
HEIGHT = 8.0
NO_SECTIONS = 4
NO_FIBERS_Y = 15
NO_FIBERS_Z = 15


def model_structure():
    """ initiate the structural model """
    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, LENGTH, 0.0, 0.0)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(NO_SECTIONS):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    fiber_area = (WIDTH / NO_FIBERS_Y) * (HEIGHT / NO_FIBERS_Z)
    y_st = WIDTH / 2 * (1 / NO_FIBERS_Y - 1)
    z_st = HEIGHT / 2 * (1 / NO_FIBERS_Z - 1)
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(NO_FIBERS_Y):
            # y = width/no_fibers_y * (i + 0.5)
            y = y_st + WIDTH / NO_FIBERS_Y * i
            for j in range(NO_FIBERS_Z):
                z = z_st + HEIGHT / NO_FIBERS_Z * j
                # z = height/no_fibers_z * (j + 0.5)
                if i in (1, 13) and j in (1, 13):
                    section.add_fiber(
                        counter, y, z, i, j, fiber_area,
                        MenegottoPinto(29000, 0.0042 * 29000, 60, 20, 18.5, 0.0002),
                    )
                else:
                    section.add_fiber(
                        counter, y, z, i, j, fiber_area, KentPark(6.95, 1, -0.07, 770)
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    return stru


def add_solution_parameters(structure):
    """ add tolerance values and boundary conditions """
    # CONVERGENCE TOLERANCE VALUES
    structure.tolerance = 1e-10
    structure.set_section_tolerance(0.05)

    # BOUNDARY CONDITIONS
    structure.set_controled_dof(2, "w")
    structure.add_neumann_condition(2, "w", 1.0)
    structure.add_dirichlet_condition(1, "uvwxyz", 0)
    structure.add_dirichlet_condition(2, "x", 0)
    structure.add_dirichlet_condition(2, "z", 0.005)  # useless
    print("Added the boundary conditions.")


def advance_in_load(structure, load_step):
    """ load stepping loop """

    if load_step < 10 + 1:
        structure.controled_dof_increment = 0.4

    elif load_step == 10 + 1:
        structure.reverse_all_fibers()
        structure.controled_dof_increment = -0.01 * 0.4

    elif load_step < 30 + 1:
        structure.controled_dof_increment = -0.4

    elif load_step == 30 + 1:
        structure.reverse_all_fibers()
        structure.controled_dof_increment = 0.01 * 0.4

    elif load_step < 55 + 1:
        structure.controled_dof_increment = 0.4

    elif load_step == 55 + 1:
        structure.reverse_all_fibers()
        structure.controled_dof_increment = -0.01 * 0.4

    elif load_step < 85 + 1:
        structure.controled_dof_increment = -0.4

    elif load_step == 85 + 1:
        structure.reverse_all_fibers()
        structure.controled_dof_increment = 0.01 * 0.4

    elif load_step < 101 + 1:
        structure.controled_dof_increment = 0.4

    elif load_step == 101 + 1:
        structure.reverse_all_fibers()
        structure.controled_dof_increment = -0.01 * 0.4

    elif load_step < 112 + 1:
        structure.controled_dof_increment = -0.4

    elif load_step == 112 + 1:
        structure.reverse_all_fibers()
        structure.controled_dof_increment = 0.01 * 0.4

    else:
        structure.controled_dof_increment = 0.4


def initiate_plot():
    """ setup the plot """
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.set(title="Moment-Curvature Plot", xlabel="Curvature [rad/in]", ylabel="Moment [kip*in]")
    axes.plot(0, 0, "yo")
    line, = axes.plot(0, 0, "bo-", mfc="none")
    axes.set_autoscalex_on(True)
    axes.set_autoscaley_on(True)
    axes.grid(True)
    return fig, axes, line


def initiate_plot_3d(structure):
    """ setup the plot """
    fig = plt.figure()
    axes = fig.add_subplot(111, projection=Axes3D.name)
    x = list()
    y = list()
    z = list()
    for node in structure.nodes:
        loc = node.get_reference_location()
        x.append(loc[0])
        y.append(loc[1])
        z.append(loc[2])
    line, = axes.plot(x, y, z, "bo-")
    axes.set(xlabel="< x >", ylabel="< y >", zlabel="< z >")
    axes.view_init(elev=0, azim=-90)
    axes.set_xlim3d(min(x) - 1, max(x) + 1)
    axes.set_ylim3d(min(y) - 1, max(y) + 1)
    axes.set_zlim3d(-4, 4)
    return fig, axes, line


def update_plot(axes, line, x, y):
    """ update the plot dynamically """
    line.set_xdata(x)
    line.set_ydata(y)
    axes.relim()
    axes.autoscale_view()
    axes.get_figure().canvas.draw()
    plt.pause(1e-20)


def update_plot_3d(axes, line, structure):
    """ update the plot dynamically """
    x = list()
    y = list()
    z = list()
    for node in structure.nodes:
        loc = node.get_actual_location()
        x.append(loc[0])
        y.append(loc[1])
        z.append(loc[2])
    line.set_xdata(x)
    line.set_ydata(y)
    line.set_3d_properties(z)
    axes.set_xlim3d(min(x) - 1, max(x) + 1)
    axes.set_ylim3d(min(y) - 1, max(y) + 1)
    axes.set_zlim3d(-4, 4)
    axes.set_top_view()
    axes.get_figure().canvas.draw()
    plt.pause(1e-20)


def main():
    """main function"""

    stru = model_structure()
    add_solution_parameters(stru)

    max_nr_iterations = 100
    max_ele_iterations = 100

    stru.initialize()
    print(":: Initialized the solver ::")

    load = [0]
    disp = [0]
    print("\n:: Starting solution loop ::")

    if PLOT_FLAG:
        fig, axes, line = initiate_plot()
        # fig3d, axes3d, line3d = initiate_plot_3d(stru)

    for k in range(1, 117 + 1):
        print(f"\nLOAD STEP : {k}")
        advance_in_load(stru, k)

        for i in range(1, max_nr_iterations + 1):
            stru.solve(max_ele_iterations)
            convergence, residual = stru.check_nr_convergence()
            if convergence:
                print(f"NR converged with {i} iteration(s). Residual = {residual}")
                break
            if i == max_nr_iterations:
                io.warning(f"Newton-Raphson did not converge {max_nr_iterations} iterations")

        stru.finalize_load_step()

        if PLOT_FLAG:
            load.append(100 * stru.converged_load_factor)
            disp.append(1 / 10000 * 3 * stru.converged_controled_dof)
            update_plot(axes, line, disp, load)
            # update_plot_3d(axes3d, line3d, stru)

    print("\n:: Finished solution loop ::")

    if PLOT_FLAG and not SAVE_PLOT:
        fig.show()
    if PLOT_FLAG and SAVE_PLOT:
        fig.savefig("moment_curvature_reference.png")


if __name__ == "__main__":
    main()
