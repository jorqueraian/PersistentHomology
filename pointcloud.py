import numpy as np
import random
import gudhi
from simplex import Simplex
from filtered_simplical_complex import FilteredSimplicialComplex
import os


gudhi.persistence_graphical_tools._gudhi_matplotlib_use_tex=False


def generate_point_cloud_from_sphere(num_points=100, radius=1, dim=2, error=.2):
    """
    Generates a point cloud from any dimensional sphere

    :param num_points: Number of points in generated point cloud
    :param radius: radius of sphere that points are picked from
    :param dim: dimension of the space in which the sphere will be created. dim=2 a circle, dim=3 a sphere and so on
    :param error: this specifies uncertainly in the points. With non-zero error points will be selected off of sphere with a margin of error
    :return: numpy array of the points in the point cloud
    """
    def apply_error(val):
        return val + random.uniform(-1 * error, error)
    # First generate random array of size (dim, num_points) and normalize them so they fall on unit ball
    random_points = np.random.normal(size=(dim, num_points))
    random_points /= np.linalg.norm(random_points, axis=0)

    # Return the list of random points with an error applied
    return radius * np.array([[apply_error(p) for p in row] for row in random_points]).T


def get_filtered_simplicial_complex_from_point_cloud(point_cloud, max_eps, max_dim, field):
    rips_complex = gudhi.RipsComplex(points=point_cloud, max_edge_length=2*max_eps)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dim)
    simplices = []

    for s in simplex_tree.get_filtration():
        simplices.append(Simplex(s[0], s[1]))

    return FilteredSimplicialComplex(simplices, field)


def save_point_cloud_to_csv(point_cloud, csv_loc):
    if os.path.exists(csv_loc):
        os.remove(csv_loc)
    num_points = point_cloud.shape[0]
    dim = point_cloud.shape[1]

    with open(csv_loc, "a") as f:
        for i in range(num_points):
            for j in range(dim-1):
                f.write(f"{point_cloud[i, j]},")
            f.write(f"{point_cloud[i, j+1]}")
            f.write(f"\n")


def load_point_cloud_from_csv(csv_loc):
    return np.genfromtxt(csv_loc, delimiter=',')


def gudhi_barcodes(point_cloud, max_eps, max_dim):
    rips_complex = gudhi.RipsComplex(points=point_cloud, max_edge_length=2 * max_eps)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dim)

    return simplex_tree.persistence()
