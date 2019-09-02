import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def initiate_plot():
    """ setup the plot """
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.set(title="Moment-Curvature Plot", xlabel="Curvature [rad/in]", ylabel="Moment [kip*in]")
    axes.plot(0, 0, "yo")
    line, = axes.plot(0, 0, "blue")
    # axes.set_autoscalex_on(True)
    # axes.set_autoscaley_on(True)
    axes.set(
        xlim=[-0.002, 0.002],
        ylim=[-160, 160],
        xticks=[-0.002, -0.001, 0, 0.001, 0.002],
        yticks=[-160, -120, -80, -40, 0, 40, 80, 120, 160],
    )
    axes.grid(True)
    return fig, axes, line


def update_plot(axes, line, x, y):
    """ update the plot dynamically """
    line.set_xdata(x)
    line.set_ydata(y)
    # axes.relim()
    # axes.autoscale_view()
    axes.get_figure().canvas.draw()
    plt.pause(1e-20)


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
    axes.set_zlim3d(-6, 6)
    return fig, axes, line


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
    axes.set_zlim3d(-6, 6)
    axes.get_figure().canvas.draw()
    plt.pause(1e-20)
