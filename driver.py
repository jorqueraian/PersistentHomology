import fieldmath as fm
import pointcloud
import persistenceplots
import ph_algo as ph

if __name__ == '__main__':
    # Generate point cloud
    dimension = 2
    points = pointcloud.generate_point_cloud_from_sphere(num_points=20, radius=1, dim=dimension, error=0.2)

    # Plot point cloud
    persistenceplots.draw_point_cloud(points)
    persistenceplots.plt.show()

    # Create Filtered Simplicial complex
    fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud=points, max_eps=1, max_dim=dimension, field=fm.Zp(3))

    # print Gudhi Barcodes
    print(pointcloud.gudhi_barcodes(point_cloud=points, max_eps=1, max_dim=dimension))

    # Print barcodes from PH algo
    barcodes = ph.involuted_ph_algo_clearing(fsc, dimension, True)
    print([bc[:2] for bc in barcodes])

    # Plot point cloud with holes and voids
    holes = [[[points[point] for point in face] for face in barcode[2]] for barcode in barcodes if barcode[0] > 0]
    persistenceplots.draw_point_cloud(points, holes=holes)
    persistenceplots.plt.show()

    # Plot barcodes
    persistenceplots.plot_barcodes(barcodes)
    persistenceplots.plt.show()
