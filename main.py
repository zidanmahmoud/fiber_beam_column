"""main"""
import fe_code

def main():
    """main function"""
    length = 100
    width = 5
    height = 8
    no_sections = 4
    no_fibers_y = 15
    no_fibers_z = 15

    stru = fe_code.Structure()
    stru.add_node(1, 0, 0, 0)
    stru.add_node(2, length, 0, 0)

    stru.add_fiber_beam_element(1, 1, 2)

    for i in range(no_sections):
        stru.get_element(1).add_section(i+1)

    concrete_material = fe_code.KentParkModel(6.95, 1, -0.07, 770)
    steel_material = fe_code.MenegottoPintoModel(29000, 0.0042*29000, 60, 20, 18.5, 0.0002)

    fiber_area = (width/no_fibers_y) * (height/no_fibers_z)
    for section in stru.get_element(1).sections:
        counter = 1
        for i in range(no_fibers_y):
            y = width/no_fibers_y * (i + 0.5)
            for j in range(no_fibers_z):
                z = height/no_fibers_z * (j + 0.5)
                if i in (1, 13) and j in (1, 13):
                    section.add_fiber(counter, y, z, fiber_area, steel_material)
                else:
                    section.add_fiber(counter, y, z, fiber_area, concrete_material)
                counter += 1

    stru.tolerance = 0.05
    for section in stru.get_element(1).sections:
        section.tolerance = 0.05

    stru.add_dirichlet_condition(1, "uvwxyz", 0.05)
    stru.add_dirichlet_condition(2, "w", 0.04)
    stru.add_dirichlet_condition(2, "x", 0)

    # stru.initialize()
    # stru.calculate_stiffness_matrix()


if __name__ == "__main__":
    main()
