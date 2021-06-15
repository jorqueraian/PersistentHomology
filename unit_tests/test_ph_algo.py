from unittest import TestCase

import fieldmath as fm
import pointcloud
import ph_algo

TEST_POINT_CLOUD = "point_cloud_test.csv"


class Test(TestCase):
    def test_standard_algo(self):
        point_cloud = pointcloud.load_point_cloud_from_csv(TEST_POINT_CLOUD)
        fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud, 2, 2, fm.Zp(3))
        bmap = fsc.create_boundary_map()
        barcodes = ph_algo.get_bar_codes(bmap, 2, cohomology=False, get_representatives=True)

        # A few comments: Ripserer computes a different rep for the 1-hole and doesnt compute reps for the 0-holes
        assert barcodes == [[1, (1.1129313486332384, 1.8578077893119918), [(2, 6), (6, 7), (1, 5), (0, 1), (2, 3), (0, 4), (5, 7), (3, 4)]], [0, (0.0, float("inf")), [(0,)]], [0, (0.0, 0.9081243849936497), [(2,), (3,)]], [0, (0.0, 0.7798078937472834), [(0,), (1,)]], [0, (0.0, 0.7640271346433212), [(0,), (4,)]], [0, (0.0, 0.6303941484701829), [(1,), (2,)]], [0, (0.0, 0.609684413223288), [(1,), (5,)]], [0, (0.0, 0.509522144845916), [(5,), (9,)]], [0, (0.0, 0.3700465613624214), [(6,), (7,)]], [0, (0.0, 0.35456706909787444), [(2,), (6,)]], [0, (0.0, 0.193046610192155), [(0,), (8,)]]]

    def test_dual_algorithm(self):
        point_cloud = pointcloud.load_point_cloud_from_csv(TEST_POINT_CLOUD)
        fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud, 2, 2, fm.Zp(3))
        cbmap = fsc.create_coboundary_map()
        barcodes = ph_algo.get_bar_codes(cbmap, 2, cohomology=True, get_representatives=True)

        # this seems wrong. Ripserer finds the rep to be [(3, 4), (3, 8), (2, 4), (0, 3), (4, 6)]  (this is adjusting for 0 indexing)
        self.fail("this seems wrong. Ripserer finds the rep to be [(3, 4), (3, 8), (2, 4), (0, 3), (4, 6)]  (this is adjusting for 0 indexing)")

        assert barcodes == []  # my algorith finds this: [[1, (1.1129313486332384, 1.8578077893119918), [(4, 6), (2, 4), (3, 4)]], [0, (0.0, float("inf")), []], [0, (0.0, 0.9081243849936497), [(3,), (2,), (1,)]], [0, (0.0, 0.7798078937472834), [(1,)]], [0, (0.0, 0.7640271346433212), [(8,), (4,)]], [0, (0.0, 0.6303941484701829), [(9,), (7,), (5,), (2,)]], [0, (0.0, 0.609684413223288), [(5,)]], [0, (0.0, 0.509522144845916), [(9,)]], [0, (0.0, 0.3700465613624214), [(7,), (6,)]], [0, (0.0, 0.35456706909787444), [(6,)]], [0, (0.0, 0.193046610192155), [(8,)]]]

    def test_involuted_ph_algo(self):
        # I unfortunitly can not say with 100% certainty that this is in fact working but it seems to be
        point_cloud = pointcloud.load_point_cloud_from_csv(TEST_POINT_CLOUD)
        fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud, 2, 2, fm.Zp(3))
        barcodes = ph_algo.involuted_ph_algo(fsc, 2, True)

        # again this isnt the output from ripserer but its still a rep for the hole i guess
        assert barcodes == [[1, (1.1129313486332384, 1.8578077893119918), [(2, 6), (6, 7), (1, 5), (0, 1), (2, 3), (0, 4), (5, 7), (3, 4)]], [0, (0.0, float("inf")), [(0,)]], [0, (0.0, 0.9081243849936497), [(2,), (3,)]], [0, (0.0, 0.7798078937472834), [(0,), (1,)]], [0, (0.0, 0.7640271346433212), [(0,), (4,)]], [0, (0.0, 0.6303941484701829), [(1,), (2,)]], [0, (0.0, 0.609684413223288), [(1,), (5,)]], [0, (0.0, 0.509522144845916), [(5,), (9,)]], [0, (0.0, 0.3700465613624214), [(6,), (7,)]], [0, (0.0, 0.35456706909787444), [(2,), (6,)]], [0, (0.0, 0.193046610192155), [(0,), (8,)]]]

    def test_involuted_on_small_example(self):
        from simplex import Simplex
        from filtered_simplical_complex import FilteredSimplicialComplex
        filtered_simplicial_complex = FilteredSimplicialComplex(
            [Simplex([0], 0), Simplex([1], 0),
             Simplex([2], 1), Simplex([3], 1), Simplex([0, 1], 1), Simplex([1, 2], 1),
             Simplex([2, 3], 2), Simplex([0, 3], 2),
             Simplex([0, 2], 3),
             Simplex([0, 1, 2], 4),
             Simplex([0, 2, 3], 5)], fm.Zp(3))

        print(ph_algo.involuted_ph_algo(filtered_simplicial_complex, 2, True))

    def test_involuted_on_small_example_2(self):
        from simplex import Simplex
        from filtered_simplical_complex import FilteredSimplicialComplex
        filtered_simplicial_complex = FilteredSimplicialComplex(
            [Simplex([0], 0), Simplex([1], 0), Simplex([2], 0), Simplex([3], 0), Simplex([4], 0),
             Simplex([0, 1], 1), Simplex([1, 2], 1), Simplex([0, 2], 1), Simplex([0, 3], 1), Simplex([0, 4], 1), Simplex([3, 4], 1),
             Simplex([0, 1, 2], 2), Simplex([0, 3, 4], 2)], fm.Zp(3))

        cmap = filtered_simplicial_complex.create_coboundary_map()
        print(ph_algo.get_bar_codes(cmap, 3, cohomology=True, get_representatives=True))

    def test_reduction_alg_bmap(self):
        # Results:
        # 2-dim 20 points
        # Number Of Simplices: 1350
        # Average Time for 'normal' alg: 9.86622365474701
        # Average Time for 'clearing' alg: 9.781163196563721
        # Average Time for 'clearing2' alg: 10.46381893634796

        # 3-dim 12 points
        # Number Of Simplices: 793
        # Average Time for 'normal' alg: 3.5520474529266357
        # Average Time for 'clearing' alg: 2.4165234184265136
        # Average Time for 'clearing2' alg: 2.8779487562179566

        import time
        dim = 3
        point_cloud = pointcloud.generate_point_cloud_from_sphere(num_points=12, radius=1, dim=dim, error=.2)
        # point_cloud = pointcloud.load_point_cloud_from_csv(TEST_POINT_CLOUD)
        fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud=point_cloud, max_eps=2,
                                                                          max_dim=dim,
                                                                          field=fm.Zp(11))
        print(f"Number Of Simplices: {fsc.num_simplices}")
        bmap = fsc.create_boundary_map()
        bmap_2 = bmap.copy()
        bmap_3 = bmap.copy()

        avg_time_1 = 0
        avg_time_2 = 0
        avg_time_3 = 0

        its = 50

        for _ in range(its):
            start_time = time.time()
            one, _ = ph_algo.reduce_boundary_map(bmap)
            end_time = time.time()
            avg_time_1 += end_time - start_time
        print(f"Average Time for 'normal' alg: {avg_time_1 / its}")

        for _ in range(its):
            start_time = time.time()
            two = ph_algo.reduce_full_boundary_map_with_clearing(bmap_2, dim)
            end_time = time.time()
            avg_time_2 += end_time - start_time
            assert one == two
        print(f"Average Time for 'clearing' alg: {avg_time_2 / its}")

        for _ in range(its):
            start_time = time.time()
            three = ph_algo.reduce_full_boundary_map_with_clearing_hom(bmap_3, dim)
            end_time = time.time()
            avg_time_3 += end_time-start_time
            assert two == three
        print(f"Average Time for 'clearing2' alg: {avg_time_3/its}")

    def test_reduction_alg_cbmap(self):
        # Results over 50 trials:

        # 2-dim 20 points
        # Number Of Simplices: 1350
        # Average Time for 'normal' alg: 4.155983452796936
        # Average Time for 'clearing' alg: 4.051275854110718
        # Average Time for 'clearing2' alg: 3.672375521659851

        # 3-dim 12 points
        # Number Of Simplices: 793
        # Average Time for 'normal' alg: 1.7547970008850098
        # Average Time for 'clearing' alg: 1.7041664409637451
        # Average Time for 'clearing2' alg: 1.7190111446380616

        # 4-dim 10 points
        # Number Of Simplices: 637
        # Average Time for 'normal' alg: 1.8635532236099244
        # Average Time for 'clearing' alg: 1.9906368350982666
        # Average Time for 'clearing2' alg: 1.5680986881256103

        import time
        dim = 4
        point_cloud = pointcloud.generate_point_cloud_from_sphere(num_points=10, radius=1, dim=dim, error=.2)
        # point_cloud = pointcloud.load_point_cloud_from_csv(TEST_POINT_CLOUD)
        fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud=point_cloud, max_eps=2, max_dim=dim,
                                                                          field=fm.Zp(11))
        print(f"Number Of Simplices: {fsc.num_simplices}")
        bmap = fsc.create_coboundary_map()
        bmap_2 = bmap.copy()
        bmap_3 = bmap.copy()

        avg_time_1 = 0
        avg_time_2 = 0
        avg_time_3 = 0

        its = 50

        for _ in range(its):
            start_time = time.time()
            one, _ = ph_algo.reduce_boundary_map(bmap)
            end_time = time.time()
            avg_time_1 += end_time - start_time
        print(f"Average Time for 'normal' alg: {avg_time_1 / its}")

        for _ in range(its):
            start_time = time.time()
            two = ph_algo.reduce_full_boundary_map_with_clearing(bmap_2, dim)
            end_time = time.time()
            avg_time_2 += end_time - start_time
            assert one == two
        print(f"Average Time for 'clearing' alg: {avg_time_2 / its}")

        for _ in range(its):
            start_time = time.time()
            three = ph_algo.reduce_full_boundary_map_with_clearing_cohom(bmap_3)
            end_time = time.time()
            avg_time_3 += end_time-start_time
            assert two == three
        print(f"Average Time for 'clearing2' alg: {avg_time_3/its}")

    def test_involuted_alg_timing(self):
        # Results over 50 trials

        # 2-dim 20 points
        # Number Of Simplices: 1350
        # Average Time for 'normal' alg: 2.0548864507675173
        # Average Time for 'clearing' alg: 0.680950379371643

        # 3-dim 12 points
        # Number Of Simplices: 793
        # Average Time for 'normal' alg: 1.4154572200775146
        # Average Time for 'clearing' alg: 0.26904468059539793

        # 4-dim 10 points
        # Number Of Simplices: 637
        # Average Time for 'normal' alg: 0.8057827091217041
        # Average Time for 'clearing' alg: 0.34567726135253907

        import time
        dim = 2
        point_cloud = pointcloud.generate_point_cloud_from_sphere(num_points=20, radius=1, dim=dim, error=.2)
        # point_cloud = pointcloud.load_point_cloud_from_csv(TEST_POINT_CLOUD)
        fsc = pointcloud.get_filtered_simplicial_complex_from_point_cloud(point_cloud=point_cloud, max_eps=2,
                                                                          max_dim=dim,
                                                                          field=fm.Zp(11))
        print(f"Number Of Simplices: {fsc.num_simplices}")

        avg_time_1 = 0
        avg_time_2 = 0
        avg_time_3 = 0

        its = 50

        for _ in range(its):
            start_time = time.time()
            ph_algo.involuted_ph_algo(fsc=fsc, max_dim=dim, reps=True)
            end_time = time.time()
            avg_time_1 += end_time - start_time
        print(f"Average Time for 'normal' alg: {avg_time_1 / its}")

        for _ in range(its):
            start_time = time.time()
            ph_algo.involuted_ph_algo_clearing(fsc=fsc, max_dim=dim, reps=True)
            end_time = time.time()
            avg_time_2 += end_time - start_time
        print(f"Average Time for 'clearing' alg: {avg_time_2 / its}")
