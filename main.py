"""main"""
import fe_code
from fe_code import debug, warning

# == MODELING PARAMETERS
LENGTH = 100
WIDTH = 5
HEIGHT = 8
NO_SECTIONS = 4
NO_FIBERS_Y = 15
NO_FIBERS_Z = 15


def main():
    """main function"""

    # STRUCTURE INITIALIZATION
    stru = fe_code.Structure()
    debug("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0, 0, 0)
    stru.add_node(2, LENGTH, 0, 0)
    debug(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    debug(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(NO_SECTIONS):
        stru.get_element(1).add_section(i + 1)
    debug(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

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
                        fe_code.MenegottoPintoModel(29000, 0.0042 * 29000, 60, 20, 18.5, 0.0002),
                    )
                else:
                    section.add_fiber(
                        counter, y, z, i, j, fiber_area, fe_code.KentParkModel(6.95, 1, -0.07, 770)
                    )
                counter += 1
    debug(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.tolerance = 0.05
    for section in stru.get_element(1).sections:
        section.tolerance = 0.05

    # BOUNDARY CONDITIONS
    controled_dof = (2, "w")
    stru.add_neumann_condition(controled_dof[0], controled_dof[1], 1)
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "x", 0)
    # stru.add_dirichlet_condition(2, "w", 0.005)
    debug("Added the boundary conditions.")

    max_nr_iterations = 100
    max_ele_iterations = 100

    stru.initialize()
    debug(":: Initialized the solver ::")

    debug("\n:: Starting solution loop ::")
    for k in range(1, 10 + 1):
        debug(f"\nLOAD STEP : {k}")

        for i in range(1, max_nr_iterations + 1):

            stru.solve(max_ele_iterations)
            if stru.check_nr_convergence():
                debug(f"NR converged with {i} iteration(s).")
                break
            if i == max_nr_iterations:
                warning(f"Newton-Raphson did not converge {max_nr_iterations} iterations")

        stru.finalize_load_step()


if __name__ == "__main__":
    main()
