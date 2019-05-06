"""main"""
import fe_code

def debug(message, *args, **kwargs):
    print(message, *args, **kwargs)

def warning(message, *args, **kwargs):
    print('\33[93m'+"WARNING: "+message.upper()+'\33[0m', *args, **kwargs)

def main():
    """main function"""

    # == MODELING PARAMETERS
    length = 100
    width = 5
    height = 8
    no_sections = 4
    no_fibers_y = 15
    no_fibers_z = 15

    # STRUCTURE INITIALIZATION
    stru = fe_code.Structure()

    # NODES
    stru.add_node(1, 0, 0, 0)
    stru.add_node(2, length, 0, 0)

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i+1)

    # FIBERS
    concrete_material = fe_code.KentParkModel(6.95, 1, -0.07, 770)
    steel_material = fe_code.MenegottoPintoModel(29000, 0.0042*29000, 60, 20, 18.5, 0.0002)
    fiber_area = (width/no_fibers_y) * (height/no_fibers_z)
    y_st = width/2*(1/no_fibers_y - 1)
    z_st = height/2*(1/no_fibers_z - 1)
    for section in stru.get_element(1).sections:
        counter = 1
        for i in range(no_fibers_y):
            # y = width/no_fibers_y * (i + 0.5)
            y = y_st + width/no_fibers_y*i
            for j in range(no_fibers_z):
                z = z_st + height/no_fibers_z*j
                # z = height/no_fibers_z * (j + 0.5)
                if i in (1, 13) and j in (1, 13):
                    section.add_fiber(counter, y, z, fiber_area, steel_material)
                else:
                    section.add_fiber(counter, y, z, fiber_area, concrete_material)
                counter += 1

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

    max_nr_iterations = 100
    max_ele_iterations = 100

    for k in range(1, 10+1):
        debug(f"\nLOAD STEP : {k}")

        for i in range(1, max_nr_iterations+1):

            stru.solve(max_ele_iterations)
            if stru.check_nr_convergence():
                debug(f"NR converged with {i} iteration(s).")
                break
            if i == max_nr_iterations:
                warning(f"Newton-Raphson did not converge {max_nr_iterations} iterations")

        stru.finalize_load_step()

        stru.load_factor += 0.5


if __name__ == "__main__":
    main()
