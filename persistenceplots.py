import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import gudhi

gudhi.persistence_graphical_tools._gudhi_matplotlib_use_tex=False


def draw_point_cloud(data, e_radius=None, save_file=None, holes=None):
    """
    Draws point cloud in either 2 or 3 dimensional space

    :param data: Point cloud as numpy array
    :param e_radius: Radius of epsilon balls. Use None to not include them
    :param save_file: Save file location (without .png) to save plot, Also acts as title
    :return: figure, wont open window of figure, so you must do plt.show()
    """
    fig = plt.figure()

    if data.shape[1] == 2:
        x = data.T[0]
        y = data.T[1]

        ax = fig.add_subplot()

        ax.scatter(x, y)

        if e_radius is not None:
            theta = np.linspace(0, 2 * np.pi, 150)
            for xp, yp in zip(x, y):
                ax.plot(e_radius * np.cos(theta) + xp, e_radius * np.sin(theta) + yp, color='black')

    elif data.shape[1] == 3:
        x, y, z = data.transpose()

        ax = fig.add_subplot(projection='3d')

        ax.scatter(x, y, z)

        if e_radius is not None:
            phi, theta = np.mgrid[0:2 * np.pi:10j, 0:np.pi:10j]
            for xp, yp, zp in zip(x, y, z):
                ax.plot_surface(e_radius * np.cos(phi) * np.sin(theta) + xp, e_radius * np.sin(phi) * np.sin(theta) + yp,
                                e_radius * np.cos(theta) + zp, alpha=0.5)

    else:
        raise ValueError("Data must have dimension 2 or 3 to create plot.")

    if holes is not None:
        for hole in holes:
            for face in hole:
                if data.shape[1] == 2:
                    ax.plot([face[0][0], face[1][0]], [face[0][1], face[1][1]], color='black')
                elif data.shape[1] == 3:
                    ax.add_collection3d(Poly3DCollection(np.array(face), color='black', alpha=.5))

    if save_file is not None:
        ax.set_title(save_file)
        fig.savefig(save_file)
    return fig


def plot_barcodes(barcodes):
    if len(barcodes[0]) > 2:
        barcodes = [bc[:2] for bc in barcodes]

    ax = gudhi.plot_persistence_barcode(barcodes, legend=True)
    return ax


def plot_persistence_diagram(barcodes):
    if len(barcodes[0]) > 2:
        barcodes = [bc[:2] for bc in barcodes]

    ax = gudhi.plot_persistence_diagram(barcodes, legend=True)
    return ax