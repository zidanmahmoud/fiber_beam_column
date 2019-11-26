import numpy as np
from math import sqrt, pi

from fe_code import io, Structure, MenegottoPinto, KentPark


def model1():
    """ initiate the structural model """

    length = 100
    width = 5
    height = 8

    no_fibers_y = 15
    no_fibers_z = 15
    no_sections = 4

    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, length, 0.0, 0.0)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    w = width / no_fibers_y
    h = height / no_fibers_z
    fiber_area = w * h
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(no_fibers_y):
            # y = width / no_fibers_y * (i + 0.5)
            y = 0.5 * (w - width) + i * w
            for j in range(no_fibers_z):
                # z = height / no_fibers_z * (j + 0.5)
                z = 0.5 * (h - height) + j * h
                # section.add_fiber(
                #         counter, y, z, fiber_area, KentPark(6.95, 1, 770, 0.0027), w, h
                #     )
                if i in (1, no_fibers_y-2) and j in (1, no_fibers_z-2):
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPinto(29000, 0.0042, 60, 20, 18.5, 0.0002),
                        w,
                        h,
                    )
                else:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentPark(6.95, 1, 770, 0.0027), w, h
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.tolerance = 1e-6
    stru.set_section_tolerance(1e-6)

    # BOUNDARY CONDITIONS
    stru.set_controlled_dof_dof(2, "w")
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "vxz", 0)
    stru.add_neumann_condition(2, "w", 1.0)
    print("Added the boundary conditions.")

    return stru


def model1c():
    """ initiate the structural model """

    length = 100
    width = 5
    height = 8

    no_fibers_y = 15
    no_fibers_x = 20
    no_sections = 4

    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, 0.0, 0.0, length)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    w = width / no_fibers_y
    h = height / no_fibers_x
    fiber_area = w * h
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(no_fibers_y):
            y = width / no_fibers_y * (i + 0.5)
            for j in range(no_fibers_x):
                z = height / no_fibers_x * (j + 0.5)
                if i in (1, 13) and j in (1, 18):
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPinto(29000, 0.0042, 60, 20, 18.5, 0.0002),
                        w,
                        h,
                    )
                else:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentPark(6.95, 1, 770, 0.0027), w, h
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.tolerance = 0.05
    stru.set_section_tolerance(0.05)

    # BOUNDARY CONDITIONS
    stru.set_controlled_dof_dof(2, "u")
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "vxz", 0)
    print("Added the boundary conditions.")

    return stru


def model2():
    """ initiate the structural model """

    length = 100
    width = 5
    height = 8

    no_fibers_y = 15
    no_fibers_z = 20
    no_sections = 4

    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, length, 0.0, 0.0)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    w = width / no_fibers_y
    h = height / no_fibers_z
    fiber_area = w * h
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(no_fibers_y):
            y = width * (-0.5 + (i + 0.5) / no_fibers_y)
            for j in range(no_fibers_z):
                z = height * (-0.5 + (j + 0.5) / no_fibers_z)
                if i in (1, 13) and j in (1, 18):
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPinto(29000, 0.0042, 48.4, 20, 18.5, 0.0002),
                        w,
                        h,
                    )
                elif i < 1 or i > 13 or j < 1 or j > 18:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h
                    )
                else:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentPark.eu(6.95, 1.0763, 0.03810, 0.0027), w, h
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.tolerance = 0.05
    stru.set_section_tolerance(0.05)

    # BOUNDARY CONDITIONS
    stru.set_controlled_dof_dof(2, "w")
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "x", 0)
    print("Added the boundary conditions.")

    return stru


def model3():
    """ initiate the structural model """

    length = 100
    width = 4 + 15/16
    height = 8.0
    cover_y = 1.0
    cover_z = 1.0

    no_sections = 3

    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, length, 0.0, 0.0)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    counter = 1

    # == confined concrete
    no_y = 2
    no_z = 12
    y = width / 2 - cover_y
    y = y - y / no_y
    z = height / 2 - cover_z
    z = z - z / no_z
    w = ((width / 2) - cover_y) / (no_y / 2)
    h = ((height / 2) - cover_z) / (no_z / 2)
    area = w * h
    ys = np.tile(np.linspace(-y, y, no_y), no_z)
    zs = np.repeat(np.linspace(-z, z, no_z), no_y)
    for section in stru.get_element(1).sections:
        for y, z in zip(ys, zs):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # == unconfined sides
    y = (width / 2) - (cover_y / 2)
    z = (height / 2) - cover_z
    z = z - z / 10
    w = cover_y
    h = ((height / 2) - cover_z) / 5
    area = w * h
    ys = np.tile(np.linspace(-y, y, 2), 10)
    zs = np.repeat(np.linspace(-z, z, 10), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(ys, zs):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    # == unconfined bottom
    y = (width / 2) / 2
    zmin = (height / 2) - ((cover_z / 4) / 2)
    zmax = (height / 2) - cover_z + ((cover_z / 4) / 2)
    w = width / 2
    h = cover_z / 4
    area = w * h
    ys = np.tile(np.linspace(-y, y, 2), 4)
    zs = np.repeat(np.linspace(-zmin, -zmax, 4), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(ys, zs):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    # == unconfined top
    ys = np.tile(np.linspace(-y, y, 2), 4)
    zs = np.repeat(np.linspace(zmin, zmax, 4), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(ys, zs):
            section.add_fiber(counter, y, z, area, KentPark.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    # == steel
    y = (width / 2) - cover_y
    z = (height / 2) - cover_z
    ys = [-y, y, -y, y]
    zs = [-z, -z, z, z]
    area = pi * (0.5 / 2) ** 2
    w = sqrt(area)
    h = sqrt(area)
    for section in stru.get_element(1).sections:
        for y, z in zip(ys, zs):
            section.add_fiber(
                counter, y, z, area, MenegottoPinto(29000, 0.0042, 48.4, 20, 18.5, 0.0002), w, h
            )
            counter += 1

    print(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.tolerance = 0.05
    stru.set_section_tolerance(0.05)

    # BOUNDARY CONDITIONS
    stru.set_controlled_dof_dof(2, "w")
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "vxz", 0)
    print("Added the boundary conditions.")

    return stru


def model4():
    """ initiate the structural model """

    length = 100
    width = 4 + 15/16
    height = 8.0
    topCover = 1
    bottomCover = 1
    sidesCover = 1
    topNumberOfSteelRebars = 2
    bottomNumberOfSteelRebars = 2
    confinedConcrete = 1
    sideConcrete = 1
    topConcrete = 1
    bottomConcrete = 1
    topBarsDia = sqrt(pi * (0.5 / 2) ** 2)
    bottomBarsDia = sqrt(pi * (0.5 / 2) ** 2)

    no_sections = 4

    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, length, 0.0, 0.0)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    counter = 1

    # bottom steel
    val = (width-2*sidesCover-bottomBarsDia) / 2
    ys = np.linspace(-val, val, bottomNumberOfSteelRebars)
    z = height/2 - bottomCover - bottomBarsDia/2
    w = bottomBarsDia
    h = bottomBarsDia
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(
                counter, y, -z, w*h, MenegottoPinto(29000, 0.0042, 48.4, 20, 18.5, 0.0002), w, h
            )
            counter += 1
    spacing = (width-2*sidesCover-bottomNumberOfSteelRebars*bottomBarsDia)/(bottomNumberOfSteelRebars-1)
    val = width/2 - sidesCover - bottomBarsDia - spacing/2
    ys = np.linspace(-val, val, bottomNumberOfSteelRebars-1)
    h = bottomBarsDia
    w = spacing
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(counter, y, -z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # top steel
    val = (width-2*sidesCover-topBarsDia) / 2
    ys = np.linspace(-val, val, topNumberOfSteelRebars)
    z = height/2 - topCover - topBarsDia/2
    w = topBarsDia
    h = topBarsDia
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(
                counter, y, z, w*h, MenegottoPinto(29000, 0.0042, 48.4, 20, 18.5, 0.0002), w, h
            )
            counter += 1
    spacing = (width-2*sidesCover-topNumberOfSteelRebars*topBarsDia)/(topNumberOfSteelRebars-1)
    val = width/2 - sidesCover - topBarsDia - spacing/2
    ys = np.linspace(-val, val, topNumberOfSteelRebars-1)
    w = spacing
    h = topBarsDia
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(counter, y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # confined concrete
    total = height-bottomCover-topCover-topBarsDia-bottomBarsDia
    area = (width/2 - sidesCover) * (total/confinedConcrete)
    val1 = height/2 - bottomCover - bottomBarsDia - total/confinedConcrete/2
    val2 = height/2 - topCover - topBarsDia - total/confinedConcrete/2
    zs = np.linspace(-val1, val2, confinedConcrete)
    y = (width/2 - sidesCover)/2
    w = width/2 - sidesCover
    h = total/confinedConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            section.add_fiber(counter, -y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1
            section.add_fiber(counter, y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # sides concrete
    total = height-bottomCover-topCover
    area = (sidesCover) * (total/sideConcrete)
    val1 = height/2 - bottomCover - total/sideConcrete/2
    val2 = height/2 - topCover - total/sideConcrete/2
    zs = np.linspace(-val1, val2, sideConcrete)
    y = width/2 - sidesCover/2
    w = sidesCover
    h = total/sideConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            section.add_fiber(counter, -y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1
            section.add_fiber(counter, y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # top concrete
    val1 = height/2 - topCover + topCover/topConcrete/2
    val2 = height/2 - topCover/topConcrete/2
    zs = np.linspace(val1, val2, topConcrete)
    y = width/4
    w = width/2
    h = topCover/topConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            section.add_fiber(counter, -y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1
            section.add_fiber(counter, y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    # bottom concrete
    val1 = height/2 - bottomCover + bottomCover/bottomConcrete/2
    val2 = height/2 - bottomCover/bottomConcrete/2
    zs = np.linspace(-val1, -val2, bottomConcrete)
    y = width/4
    w = width/2
    h = bottomCover/bottomConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            section.add_fiber(counter, -y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1
            section.add_fiber(counter, y, z, w*h, KentPark.eu(6.95, 1, 0.03810, 0.0027), w, h)
            counter += 1

    print(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.set_tolerance(1e-9)
    stru.set_section_tolerance(1e-9)

    # BOUNDARY CONDITIONS
    stru.set_controlled_dof(2, "w")
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "vxz", 0)
    stru.add_neumann_condition(2, "w", 1.0)
    print("Added the boundary conditions.")

    return stru


def model5():
    """ initiate the structural model """

    length = 71
    width = 9
    height = 16
    topCover = 2
    bottomCover = 1
    sidesCover = 0.75
    topNumberOfSteelRebars = 4
    bottomNumberOfSteelRebars = 3
    confinedConcrete = 5
    sideConcrete = 1
    topConcrete = 1
    bottomConcrete = 1
    topBarsDia = sqrt(pi * (0.75 / 2) ** 2)
    bottomBarsDia = sqrt(pi * (0.625 / 2) ** 2)
    # topBarsDia = 0.75
    # bottomBarsDia = 0.625

    no_sections = 2

    # STRUCTURE INITIALIZATION
    stru = Structure()
    print("Constructed an empty stucture.")

    # NODES
    stru.add_node(1, 0.0, 0.0, 0.0)
    stru.add_node(2, length, 0.0, 0.0)
    print(f"Added {len(stru.nodes)} nodes.")

    # ELEMENTS
    stru.add_fiber_beam_element(1, 1, 2)
    print(f"Added {len(stru.elements)} elements.")

    # SECTIONS
    for i in range(no_sections):
        stru.get_element(1).add_section(i + 1)
    print(f"Added {sum([len(element.sections) for element in stru.elements])} sections.")

    # FIBERS
    counter = 1

    # bottom steel
    val = (width-2*sidesCover-bottomBarsDia) / 2
    ys = np.linspace(-val, val, bottomNumberOfSteelRebars)
    z = height/2 - bottomCover - bottomBarsDia/2
    w = bottomBarsDia
    h = bottomBarsDia
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(
                counter, y, z, w*h, MenegottoPinto(29000, 0.0085, 66.5, 20, 18.5, 0.0002), w, h
            )
            counter += 1
    spacing = (width-2*sidesCover-bottomNumberOfSteelRebars*bottomBarsDia)/(bottomNumberOfSteelRebars-1)
    val = width/2 - sidesCover - bottomBarsDia - spacing/2
    ys = np.linspace(-val, val, bottomNumberOfSteelRebars-1)
    h = bottomBarsDia
    w = spacing
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(counter, y, z, w*h, KentPark.eu(5.43, 1, 0.069, 0.00214), w, h)
            counter += 1

    # top steel
    val = (width-2*sidesCover-topBarsDia) / 2
    ys = np.linspace(-val, val, topNumberOfSteelRebars)
    z = height/2 - topCover - topBarsDia/2
    w = topBarsDia
    h = topBarsDia
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(
                counter, y, -z, w*h, MenegottoPinto(29000, 0.0085, 66.5, 20, 18.5, 0.0002), w, h
            )
            counter += 1
    spacing = (width-2*sidesCover-topNumberOfSteelRebars*topBarsDia)/(topNumberOfSteelRebars-1)
    val = width/2 - sidesCover - topBarsDia - spacing/2
    ys = np.linspace(-val, val, topNumberOfSteelRebars-1)
    w = spacing
    h = topBarsDia
    for section in stru.get_element(1).sections:
        for y in ys:
            section.add_fiber(counter, y, -z, w*h, KentPark.eu(5.43, 1, 0.069, 0.00214), w, h)
            counter += 1

    # confined concrete
    total = height-bottomCover-topCover-topBarsDia-bottomBarsDia
    area = (width/2 - sidesCover) * (total/confinedConcrete)
    val1 = height/2 - bottomCover - bottomBarsDia - total/confinedConcrete/2
    val2 = height/2 - topCover - topBarsDia - total/confinedConcrete/2
    zs = np.linspace(-val1, val2, confinedConcrete)
    y = (width/2 - sidesCover)/2
    w = width/2 - sidesCover
    h = total/confinedConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            section.add_fiber(counter, -y, -z, w*h, KentPark.eu(5.43, 1, 0.069, 0.00265), w, h)
            counter += 1
            section.add_fiber(counter, y, -z, w*h, KentPark.eu(5.43, 1, 0.069, 0.00265), w, h)
            counter += 1

    # sides concrete
    total = height-bottomCover-topCover
    area = (sidesCover) * (total/sideConcrete)
    val1 = height/2 - bottomCover - total/sideConcrete/2
    val2 = height/2 - topCover - total/sideConcrete/2
    zs = np.linspace(-val1, val2, sideConcrete)
    y = width/2 - sidesCover/2
    w = sidesCover
    h = total/sideConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            section.add_fiber(counter, -y, -z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            counter += 1
            section.add_fiber(counter, y, -z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            counter += 1

    # top concrete
    val1 = height/2 - topCover + topCover/topConcrete/2
    val2 = height/2 - topCover/topConcrete/2
    zs = np.linspace(val1, val2, topConcrete)
    y = width/4
    # w = width/2
    w = width
    h = topCover/topConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            # section.add_fiber(counter, -y, z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            # counter += 1
            # section.add_fiber(counter, y, z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            # counter += 1
            section.add_fiber(counter, 0, -z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            counter += 1

    # bottom concrete
    val1 = height/2 - bottomCover + bottomCover/bottomConcrete/2
    val2 = height/2 - bottomCover/bottomConcrete/2
    zs = np.linspace(-val1, -val2, bottomConcrete)
    y = width/4
    # w = width/2
    w = width
    h = bottomCover/bottomConcrete
    for section in stru.get_element(1).sections:
        for z in zs:
            # section.add_fiber(counter, -y, z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            # counter += 1
            # section.add_fiber(counter, y, z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            # counter += 1
            section.add_fiber(counter, 0, -z, w*h, KentPark.eu(5.07, 1, 0.003, 0.002), w, h)
            counter += 1

    print(f"Added {counter - 1} fibers.")

    # CONVERGENCE TOLERANCE VALUES
    stru.tolerance = 1e-9
    stru.set_section_tolerance(1e-9)

    # BOUNDARY CONDITIONS
    stru.set_controlled_dof_dof(2, "w")
    stru.add_dirichlet_condition(1, "uvwxyz", 0)
    stru.add_dirichlet_condition(2, "vxz", 0)
    print("Added the boundary conditions.")

    return stru
