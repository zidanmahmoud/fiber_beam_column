import numpy as np
from math import sqrt, pi

from fe_code import io, Structure, MenegottoPinto, KentPark, KentParkMod, MenegottoPintoMod


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
    fiber_area = (width / no_fibers_y) * (height / no_fibers_z)
    w = width / no_fibers_y
    h = height / no_fibers_z
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(no_fibers_y):
            y = width / no_fibers_y * (i + 0.5)
            for j in range(no_fibers_z):
                z = height / no_fibers_z * (j + 0.5)
                if i in (1, 13) and j in (1, 13):
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPintoMod(29000, 0.0042, 60, 20, 18.5, 0.0002),
                        w,
                        h,
                    )
                else:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentParkMod(6.95, 1, 770, 0.0027), w, h
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    return stru


def model2():
    """ initiate the structural model """

    length = 100
    width = 5
    height = 8

    no_fibers_y = 10
    no_fibers_z = 16
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
    fiber_area = (width / no_fibers_y) * (height / no_fibers_z)
    w = sqrt(fiber_area)
    h = sqrt(fiber_area)
    counter = 1
    for section in stru.get_element(1).sections:
        for i in range(no_fibers_y):
            y = width * (-0.5 + (i + 0.5) / no_fibers_y)
            for j in range(no_fibers_z):
                z = height * (-0.5 + (j + 0.5) / no_fibers_z)
                if i in (2, 7) and j in (2, 13):
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        MenegottoPintoMod(29000, 0.0042, 48.4, 20, 18.5, 0.0002),
                        w,
                        h,
                    )
                elif i < 2 or i > 7 or j < 2 or j > 13:
                    section.add_fiber(
                        counter, y, z, fiber_area, KentParkMod.eu(6.95, 1, 0.00292, 0.0027), w, h
                    )
                else:
                    section.add_fiber(
                        counter,
                        y,
                        z,
                        fiber_area,
                        KentParkMod.eu(6.95, 1.0763, 0.03810, 0.0027),
                        w,
                        h,
                    )
                counter += 1
    print(f"Added {counter - 1} fibers.")

    return stru


def model3():
    """ initiate the structural model """

    length = 100
    width = 4.83
    height = 8.0
    cover_y = 0.79
    cover_z = 1.0

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
            section.add_fiber(counter, y, z, area, KentParkMod.eu(6.95, 1, 0.03810, 0.0027), w, h)
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
            section.add_fiber(counter, y, z, area, KentParkMod.eu(6.95, 1, 0.00292, 0.0027), w, h)
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
            section.add_fiber(counter, y, z, area, KentParkMod.eu(6.95, 1, 0.00292, 0.0027), w, h)
            counter += 1

    # == unconfined top
    ys = np.tile(np.linspace(-y, y, 2), 4)
    zs = np.repeat(np.linspace(zmin, zmax, 4), 2)
    for section in stru.get_element(1).sections:
        for y, z in zip(ys, zs):
            section.add_fiber(counter, y, z, area, KentParkMod.eu(6.95, 1, 0.00292, 0.0027), w, h)
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
                counter, y, z, area, MenegottoPintoMod(29000, 0.0042, 48.4, 20, 18.5, 0.0002), w, h
            )
            counter += 1

    print(f"Added {counter - 1} fibers.")

    return stru
