"""
Example
"""

from fe_code import io, Structure, MenegottoPinto, KentPark
from plotting import plot_disctrized_2d, initiate_plot, update_plot

from disp_calc import calculate_loadsteps

STEP = 0.4
STEPS = calculate_loadsteps(STEP)

PLOT_FLAG = True
SAVE_PLOT = True

# == MODELING PARAMETERS
LENGTH = 100.0
WIDTH = 4.83
HEIGHT = 8.0
NO_SECTIONS = 3


def model_structure():
    """ initiate the structural model """

    NO_FIBERS_Y = 15
    NO_FIBERS_Z = 15

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
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPinto(29000, 0.0042, 60, 20, 18.5, 0.0002),
                    )
                else:
                    section.add_fiber(counter, y, z, fiber_area, KentPark(6.95, 1, 770, 0.0027))
                counter += 1
    print(f"Added {counter - 1} fibers.")

    return stru


def model_structure_2():
    """ initiate the structural model """

    NO_SECTIONS = 3
    NO_FIBERS_Y = 10
    NO_FIBERS_Z = 16

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
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(NO_FIBERS_Y):
            y = WIDTH * (-0.5 + (i + 0.5) / NO_FIBERS_Y)
            for j in range(NO_FIBERS_Z):
                z = HEIGHT * (-0.5 + (j + 0.5) / NO_FIBERS_Z)
                if i in (2, 7) and j in (2, 13):
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPinto(29000, 0.0042, 48.4, 20, 18.5, 0.0002),
                    )
                # elif i < 2 or i > 7 or j < 2 or j > 13:
                #     section.add_fiber(counter, y, z, fiber_area, KentPark.eu(6.95, 1, 0.00292, 0.0027))
                else:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentPark(6.95, 1.0763, 40.77, 0.002)
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    return stru


def model_structure_3():
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
    counter = 1

    # == steel
    y = (WIDTH / 2) - 0.79
    z = (HEIGHT / 2) - 1.0
    y_steel = [-y, y, -y, y]
    z_steel = [-z, -z, z, z]
    from math import pi, sqrt
    import numpy as np

    area = pi * (0.5 / 2) ** 2
    w = sqrt(area)
    h = sqrt(area)
    for section in stru.get_element(1).sections:
        for y, z in zip(y_steel, z_steel):
            section.add_fiber(
                counter, y, z, area, MenegottoPinto(29000, 0.0042, 48.4, 20, 18.5, 0.0002), w, h
            )
            counter += 1

    # == confined concrete
    y = y / 2
    z = z - z / 12
    w = (WIDTH / 2) - 0.79
    h = ((HEIGHT / 2) - 1.0) / 6
    area = w * h
    y_confined = np.tile(np.linspace(-y, y, 2), 12)
    z_confined = np.repeat(np.linspace(-z, z, 12), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(y_confined, z_confined):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # == unconfined sides
    y = (WIDTH / 2) - (0.79 / 2)
    z = (HEIGHT / 2) - 1.0
    z = z - z / 10
    w = 0.79
    h = ((HEIGHT / 2) - 1.0) / 5
    area = w * h
    y_unconfined_sides = np.tile(np.linspace(-y, y, 2), 10)
    z_unconfined_sides = np.repeat(np.linspace(-z, z, 10), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(y_unconfined_sides, z_unconfined_sides):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    # == unconfined bottom
    y = (WIDTH / 2) / 2
    zmin = (HEIGHT / 2) - ((1.0 / 4) / 2)
    zmax = (HEIGHT / 2) - 1.0 + ((1.0 / 4) / 2)
    w = WIDTH / 2
    h = 1.0 / 4
    area = w * h
    y_unconfined_bottom = np.tile(np.linspace(-y, y, 2), 4)
    z_unconfined_bottom = np.repeat(np.linspace(-zmin, -zmax, 4), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(y_unconfined_bottom, z_unconfined_bottom):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    # == unconfined top
    y_unconfined_top = np.tile(np.linspace(-y, y, 2), 4)
    z_unconfined_top = np.repeat(np.linspace(zmin, zmax, 4), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(y_unconfined_top, z_unconfined_top):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    print(f"Added {counter - 1} fibers.")

    return stru


def add_solution_parameters(structure):
    """ add tolerance values and boundary conditions """
    # CONVERGENCE TOLERANCE VALUES
    structure.tolerance = 0.05
    structure.set_section_tolerance(0.05)

    # BOUNDARY CONDITIONS
    structure.set_controled_dof(2, "w")
    structure.add_dirichlet_condition(1, "uvwxyz", 0)
    structure.add_dirichlet_condition(2, "x", 0)
    structure.add_dirichlet_condition(2, "z", 0.005)  # useless
    print("Added the boundary conditions.")


def advance_in_load(structure, load_step):
    """ load stepping loop """

    if load_step < STEPS[0]:
        structure.controled_dof_increment = STEP

    elif load_step < STEPS[1]:
        structure.controled_dof_increment = -STEP

    elif load_step < STEPS[2]:
        structure.controled_dof_increment = STEP

    elif load_step < STEPS[3]:
        structure.controled_dof_increment = -STEP

    elif load_step < STEPS[4]:
        structure.controled_dof_increment = STEP

    elif load_step < STEPS[5]:
        structure.controled_dof_increment = -STEP

    else:
        structure.controled_dof_increment = STEP


def solution_loop(structure):
    max_nr_iterations = 100
    max_ele_iterations = 100

    structure.initialize()
    print(":: Initialized the solver ::")

    load = [0]
    disp = [0]
    print("\n:: Starting solution loop ::")

    if PLOT_FLAG:
        fig, axes, line = initiate_plot()

    for k in range(1, STEPS[-1]):
        print(f"\nLOAD STEP : {k}")
        advance_in_load(structure, k)

        for i in range(1, max_nr_iterations + 1):
            structure.solve(max_ele_iterations)
            convergence, residual = structure.check_nr_convergence()
            if convergence:
                print(f"NR converged with {i} iteration(s). Residual = {residual}")
                break
            if i == max_nr_iterations:
                io.warning(f"Newton-Raphson did not converge {max_nr_iterations} iterations")

        structure.finalize_load_step()

        if PLOT_FLAG:
            load.append(100 * structure.converged_load_factor)
            disp.append(1 / 10000 * 3 * structure.converged_controled_dof)
            update_plot(axes, line, disp, load)

    print("\n:: Finished solution loop ::")
    if PLOT_FLAG and not SAVE_PLOT:
        fig.show()
    if PLOT_FLAG and SAVE_PLOT:
        fig.savefig("moment_curvature.png")


if __name__ == "__main__":
    stru = model_structure_3()
    plot_disctrized_2d(stru.get_element(1))
    add_solution_parameters(stru)
    solution_loop(stru)
