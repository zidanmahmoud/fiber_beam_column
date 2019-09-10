"""
Example
"""

import plotting as p
from fe_code import io
from models.column import model1, model2, model3
from disp_calc import calculate_loadsteps


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


def solution_loop(structure, *plot_args, **plot_kwargs):
    max_nr_iterations = 10
    max_ele_iterations = 100

    structure.initialize()
    print(":: Initialized the solver ::")

    load = [0]
    disp = [0]
    print("\n:: Starting solution loop ::")

    if PLOT_FLAG:
        fig, axes, line = p.initiate_plot(*plot_args, **plot_kwargs)

    stresses = []
    strains = []

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
                io.warning("FATAL ERROR: The solution is unstable")
                break

        structure.finalize_load_step()
        # load.append(100 * structure.converged_load_factor)
        # disp.append(1 / 10000 * 3 * structure.converged_controled_dof)
        load.append(structure.get_force((1, "y")))
        disp.append(-1/40 * structure.get_dof_value((2, "y")))

        # fiber = structure.get_element(1).get_section(1).get_fiber(279)
        # stresses.append(fiber.stress)
        # strains.append(fiber.strain)

        if PLOT_FLAG:
            p.update_plot(axes, line, disp, load)

    print("\n:: Finished solution loop ::")
    if PLOT_FLAG and not SAVE_PLOT:
        p.keep_plot()
    if PLOT_FLAG and SAVE_PLOT:
        fig.savefig("moment_curvature.png")

    return disp, load, stresses, strains


if __name__ == "__main__":
    PLOT_FLAG = True
    SAVE_PLOT = False
    STEP = 0.4
    STEPS = calculate_loadsteps(STEP)

    stru = model3()
    p.plot_disctrized_2d(stru.get_element(1))
    d, l, sig, eps = solution_loop(stru, "-o", color="blue", mfc="none")
    # p.custom_2d_plot(eps, sig, "blue")
