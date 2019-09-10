import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import sqrt

from fe_code import KentPark, MenegottoPinto


def plot_disctrized_2d(beam):

    fig = plt.figure()
    axes = fig.add_subplot(111)
    for fiber in beam.get_section(1).fibers:
        y = -fiber.direction[0]
        z = fiber.direction[1]
        side = sqrt(fiber.area)
        if isinstance(fiber._material, KentPark):
            fcolor = "grey"
        elif isinstance(fiber._material, MenegottoPinto):
            fcolor = "black"
        else:
            raise TypeError(f"Unknown fiber material law: {fiber._material}")
        ecolor = "black"
        rec = patches.Rectangle(
            xy=(y - fiber.w / 2, z - fiber.h / 2),
            width=fiber.w,
            height=fiber.h,
            facecolor=fcolor,
            edgecolor=ecolor,
        )
        axes.add_patch(rec)
    axes.autoscale()
    axes.set(aspect="equal")
    plt.show()


def custom_2d_plot(x, y, *args, **kwargs):
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.plot(x, y, *args, **kwargs)
    axes.grid()
    plt.show()
