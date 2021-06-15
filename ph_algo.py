import matrixmath as mm


def reduce_boundary_map(b_map):
    lin_com = [[c] for c in range(b_map.columns)]
    rb_map = mm.Matrix(b_map[:, :], b_map.f, b_map.row_bases, b_map.column_bases)
    l_arr = [-1] * rb_map.rows
    for i in range(rb_map.columns):
        while rb_map.low(i) != -1 and l_arr[rb_map.low(i)] != -1:
            j = l_arr[rb_map.low(i)]
            lam = rb_map.f.multiply(rb_map[rb_map.low(i), i], rb_map.f.reciprocal(rb_map[rb_map.low(j), j]))
            rb_map.add_columns(j, i, rb_map.f.negate(lam))
            lin_com[i].append(j)
        if rb_map.low(i) != -1:
            l_arr[rb_map.low(i)] = i
    return rb_map, lin_com


def reduce_full_boundary_map_with_clearing(b_map, dim):
    rb_map = mm.Matrix(b_map[:, :], b_map.f, b_map.row_bases, b_map.column_bases)
    l_arr = [-1] * rb_map.rows
    for d in range(dim, -1, -1):  # I guess you can only go to dim 1, as dim 0 are already trivial columns
        for i in range(rb_map.columns):
            if rb_map.column_bases[i].dimension == d:
                while rb_map.low(i) != -1 and l_arr[rb_map.low(i)] != -1:
                    j = l_arr[rb_map.low(i)]
                    lam = rb_map.f.multiply(rb_map[rb_map.low(i), i], rb_map.f.reciprocal(rb_map[rb_map.low(j), j]))
                    rb_map.add_columns(j, i, rb_map.f.negate(lam))
                if rb_map.low(i) != -1:
                    l_arr[rb_map.low(i)] = i
                    if rb_map.low(i) != -1:
                        rb_map.clear_column(rb_map.low(i))
                        # rb_map[:, rb_map.low(i)] = [rb_map.f.zero]*rb_map.rows
    return rb_map


def reduce_full_boundary_map_with_clearing_hom(b_map, dim):
    rb_map = mm.Matrix(b_map[:, :], b_map.f, b_map.row_bases, b_map.column_bases)
    l_arr = [-1] * rb_map.rows
    for d in range(dim, 0, -1):  # I guess you can only go to dim 1, as dim 0 are already trivial columns
        for i in range(rb_map.columns):
            if rb_map.column_bases[i].dimension == d:
                while rb_map.low(i) != -1 and l_arr[rb_map.low(i)] != -1:
                    j = l_arr[rb_map.low(i)]
                    lam = rb_map.f.multiply(rb_map[rb_map.low(i), i], rb_map.f.reciprocal(rb_map[rb_map.low(j), j]))
                    rb_map.add_columns(j, i, rb_map.f.negate(lam))
                if rb_map.low(i) != -1:
                    l_arr[rb_map.low(i)] = i
                    if rb_map.low(i) != -1:
                        rb_map.clear_column(rb_map.low(i))
                        # rb_map[:, rb_map.low(i)] = [rb_map.f.zero]*rb_map.rows
    return rb_map


def reduce_full_boundary_map_with_clearing_cohom(b_map):
    rb_map = mm.Matrix(b_map[:, :], b_map.f, b_map.row_bases, b_map.column_bases)
    l_arr = [-1] * rb_map.rows
    for i in range(rb_map.columns):
        while rb_map.low(i) != -1 and l_arr[rb_map.low(i)] != -1:
            j = l_arr[rb_map.low(i)]
            lam = rb_map.f.multiply(rb_map[rb_map.low(i), i], rb_map.f.reciprocal(rb_map[rb_map.low(j), j]))
            rb_map.add_columns(j, i, rb_map.f.negate(lam))
        if rb_map.low(i) != -1:
            l_arr[rb_map.low(i)] = i
            if rb_map.low(i) != -1:
                rb_map.clear_column(rb_map.low(i))
                # rb_map[:, rb_map.low(i)] = [rb_map.f.zero]*rb_map.rows
    return rb_map


def reduce_p_boundary_and_get_clearing(b_map):
    lin_com = [[c] for c in range(b_map.columns)]
    clearing = []
    rb_map = mm.Matrix(b_map[:, :], b_map.f, b_map.row_bases, b_map.column_bases)
    l_arr = [-1] * rb_map.rows
    for i in range(rb_map.columns):
        while rb_map.low(i) != -1 and l_arr[rb_map.low(i)] != -1:
            j = l_arr[rb_map.low(i)]
            lam = rb_map.f.multiply(rb_map[rb_map.low(i), i], rb_map.f.reciprocal(rb_map[rb_map.low(j), j]))
            rb_map.add_columns(j, i, rb_map.f.negate(lam))
            lin_com[i].append(j)
        if rb_map.low(i) != -1:
            l_arr[rb_map.low(i)] = i
            clearing.append(rb_map.low(i))
    return rb_map, clearing, lin_com


def get_bar_codes(b_map, max_dim, cohomology=False, get_representatives=True):
    """ This is only for standard and dual algorithms. And is very slow """
    # Reduce boundary/coboundary map
    reduced_b_map, lin_combs = reduce_boundary_map(b_map)
    field = reduced_b_map.f
    persistence_pairs = []
    barcodes = []
    paired = [False] * reduced_b_map.columns
    for i in range(reduced_b_map.columns):
        if reduced_b_map.low(i) != -1:
            paired[reduced_b_map.low(i)] = True
            paired[i] = True
            persistence_pairs.append([reduced_b_map.low(i), i])

    # Add essential simplicies
    for i, p in enumerate(paired):
        if not p and reduced_b_map.row_bases[i].dimension < max_dim:
            persistence_pairs.append([i, None])

    # compute reps
    if get_representatives:
        representatives = []
        for pp in persistence_pairs:
            if pp[1] is not None:
                if cohomology:
                    # Cohomology persistent pairs
                    # Honesty have no clue if this is correct. Its not. or its probably not
                    representatives.append(
                        [reduced_b_map.row_bases[i].simplex for i, v in enumerate(reduced_b_map[pp[0], :])  # [pp[0], :]
                         if not field.equals(v, field.zero)])
                else:
                    # Homology persistent pairs
                    representatives.append(
                        [reduced_b_map.row_bases[i].simplex for i, v in enumerate(reduced_b_map[:, pp[1]])
                         if not field.equals(v, field.zero)])
            else:
                if cohomology:
                    representatives.append([])
                else:
                    # Rep cycles for essential Simplices
                    representatives.append([])
                    rep = lin_combs[pp[0]]
                    for i in rep:
                        if i != pp[0]:
                            for j in lin_combs[i]:
                                if i != j:
                                    rep += [j]
                        representatives[-1].append(reduced_b_map.column_bases[i].simplex)

    # Now put everything together
    for i, x in enumerate(persistence_pairs):
        bc_dim = reduced_b_map.row_bases[x[1]].dimension if cohomology and x[1] else reduced_b_map.row_bases[x[0]].dimension
        if x[1] is not None:
            if cohomology:
                bc = tuple([reduced_b_map.row_bases[x[1]].filtration, reduced_b_map.row_bases[x[0]].filtration])
            else:
                bc = tuple([reduced_b_map.row_bases[x[0]].filtration, reduced_b_map.row_bases[x[1]].filtration])
        else:
            bc = tuple([reduced_b_map.row_bases[x[0]].filtration, float("inf")])
        barcode = [bc_dim, bc]

        if get_representatives:
            barcode.append(representatives[i])
        barcodes.append(barcode)

    return sort_and_clean_barcodes(barcodes)


def get_cohomology_pairs_from_reduced_boundary_map(reduced_b_map):
    # helper function for involuted_ph_algo
    cohomology_simplex_pairs = []
    for i in range(reduced_b_map.columns):
        if reduced_b_map.low(i) != -1:
            cohomology_simplex_pairs.append([reduced_b_map.low(i), reduced_b_map.row_bases[reduced_b_map.low(i)], i,
                                             reduced_b_map.column_bases[i]])
    return cohomology_simplex_pairs


def involuted_ph_algo(fsc, max_dim, reps=True):
    barcodes = []
    paired = [False] * fsc.num_simplicies

    # we store these prev values because the essential p-simplices are computed on the p+1 iteration so we end up
    # computing $D_p$ on the p+1 iteration
    prev_cohomology_pairs = []
    prev_hom_death_simplicies = []

    for p in range(max_dim + 2):
        hom_death_simplicies = []
        essential_simplicies = []
        if p < max_dim + 1:
            # compute reduced coboundary matrix and get cohomology pairs
            # check to see if the lin_com are the cohom cocycles??
            cb_map, _ = reduce_boundary_map(fsc.create_p_coboundary_map(p))
            cohomology_pairs = get_cohomology_pairs_from_reduced_boundary_map(cb_map)

            # compute the homology death p-simplicies
            for i, birth_simplex, j, death_simplex in cohomology_pairs:
                paired[fsc.simplex_to_index_map_inverted[cb_map.row_bases[i]]] = True
                paired[fsc.simplex_to_index_map_inverted[cb_map.column_bases[j]]] = True
                hom_death_simplicies.append(birth_simplex)
                # I dont compute cocycle reps. But you theoretically can here

            # Compute essential (p-1)-simplices
            for i, val in enumerate(paired):
                if not val and fsc.index_to_simplex_map_inverted[i].dimension < p:   # Big unsure if this will work
                    essential_simplicies.append(fsc.index_to_simplex_map_inverted[i])
                    paired[i] = True

        # at this point we want to compute D_{p-1}
        if len(prev_hom_death_simplicies) != 0 or len(essential_simplicies) != 0:
            b_minor_map = fsc.construct_boundary_map_from_simplicies(prev_hom_death_simplicies+essential_simplicies)
            reduced_minor_b_matrix, hom_lin_combs = reduce_boundary_map(b_minor_map)

            cycles = []
            ess_cycles = []
            # Compute homology rep cycles
            for col in range(reduced_minor_b_matrix.columns):
                if reduced_minor_b_matrix.low(col) != -1:
                    cycles.append([])
                    for i, v in enumerate(reduced_minor_b_matrix[:, col]):
                        if not fsc.field.equals(v, fsc.field.zero):
                            cycles[-1].append(reduced_minor_b_matrix.row_bases[i])
                else:
                    ess_cycles.append([])
                    rep = hom_lin_combs[col]
                    for i in rep:
                        if i != col:
                            for j in hom_lin_combs[i]:
                                if i != j:
                                    rep += [j]
                        ess_cycles[-1].append(reduced_minor_b_matrix.column_bases[i])

            # Put everything together
            assert len(prev_cohomology_pairs) == len(cycles), "something went wrong with cohomology pairs and cycle reps!"
            for pair, cycle in zip(prev_cohomology_pairs, cycles):
                if reps:
                    barcode = [p - 2, tuple([pair[3].filtration, pair[1].filtration]), [s.simplex for s in cycle]]
                else:
                    barcode = [p - 2, tuple([pair[3].filtration, pair[1].filtration])]
                barcodes.append(barcode)

            assert len(essential_simplicies) == len(ess_cycles), "something went wrong with essential simplices and cycles!"
            for simplex, cycle in zip(essential_simplicies, ess_cycles):
                if reps:
                    barcode = [p - 1, tuple([simplex.filtration, float("inf")]), [s.simplex for s in cycle]]
                else:
                    barcode = [p - 1, tuple([simplex.filtration, float("inf")])]
                barcodes.append(barcode)
        prev_cohomology_pairs = cohomology_pairs
        prev_hom_death_simplicies = hom_death_simplicies

    return sort_and_clean_barcodes(barcodes)


def involuted_ph_algo_clearing(fsc, max_dim, reps=True):
    barcodes = []
    paired = [False] * fsc.num_simplicies

    # we store these prev values because the essential p-simplices are computed on the p+1 iteration so we end up
    # computing $D_p$ on the p+1 iteration
    prev_cohomology_pairs = []
    prev_hom_death_simplicies = []
    clearing = []

    for p in range(max_dim + 2):
        hom_death_simplicies = []
        essential_simplicies = []
        if p < max_dim + 1:
            # compute reduced coboundary matrix and get cohomology pairs
            # check to see if the lin_com are the cohom cocycles??
            cb_map = fsc.create_p_coboundary_map(p)
            for c in clearing:
                cb_map[:, c] = [cb_map.f.zero] * cb_map.rows
            cb_map, clearing, _ = reduce_p_boundary_and_get_clearing(cb_map)
            cohomology_pairs = get_cohomology_pairs_from_reduced_boundary_map(cb_map)

            # compute the homology death p-simplicies
            for i, birth_simplex, j, death_simplex in cohomology_pairs:
                paired[fsc.simplex_to_index_map_inverted[cb_map.row_bases[i]]] = True
                paired[fsc.simplex_to_index_map_inverted[cb_map.column_bases[j]]] = True
                hom_death_simplicies.append(birth_simplex)
                # I dont compute cocycle reps. But you theoretically can here

            # Compute essential (p-1)-simplices
            for i, val in enumerate(paired):
                if not val and fsc.index_to_simplex_map_inverted[i].dimension < p:   # Big unsure if this will work
                    essential_simplicies.append(fsc.index_to_simplex_map_inverted[i])
                    paired[i] = True

        # at this point we want to compute D_{p-1}
        if len(prev_hom_death_simplicies) != 0 or len(essential_simplicies) != 0:
            b_minor_map = fsc.construct_boundary_map_from_simplicies(prev_hom_death_simplicies+essential_simplicies)
            reduced_minor_b_matrix, hom_lin_combs = reduce_boundary_map(b_minor_map)

            cycles = []
            ess_cycles = []
            # Compute homology rep cycles
            for col in range(reduced_minor_b_matrix.columns):
                if reduced_minor_b_matrix.low(col) != -1:
                    cycles.append([])
                    for i, v in enumerate(reduced_minor_b_matrix[:, col]):
                        if not fsc.field.equals(v, fsc.field.zero):
                            cycles[-1].append(reduced_minor_b_matrix.row_bases[i])
                else:
                    ess_cycles.append([])
                    rep = hom_lin_combs[col]
                    for i in rep:
                        if i != col:
                            for j in hom_lin_combs[i]:
                                if i != j:
                                    rep += [j]
                        ess_cycles[-1].append(reduced_minor_b_matrix.column_bases[i])

            # Put everything together
            assert len(prev_cohomology_pairs) == len(cycles), "something went wrong with cohomology pairs and cycle reps!"
            for pair, cycle in zip(prev_cohomology_pairs, cycles):
                if reps:
                    barcode = [p - 2, tuple([pair[3].filtration, pair[1].filtration]), [s.simplex for s in cycle]]
                else:
                    barcode = [p - 2, tuple([pair[3].filtration, pair[1].filtration])]
                barcodes.append(barcode)

            assert len(essential_simplicies) == len(ess_cycles), "something went wrong with essential simplices and cycles!"
            for simplex, cycle in zip(essential_simplicies, ess_cycles):
                if reps:
                    barcode = [p - 1, tuple([simplex.filtration, float("inf")]), [s.simplex for s in cycle]]
                else:
                    barcode = [p - 1, tuple([simplex.filtration, float("inf")])]
                barcodes.append(barcode)
        prev_cohomology_pairs = cohomology_pairs
        prev_hom_death_simplicies = hom_death_simplicies

    return sort_and_clean_barcodes(barcodes)


def sort_and_clean_barcodes(barcodes):
    barcodes = [[bc[0], tuple([bc[1][0], bc[1][1]])] + bc[2:] for bc in barcodes if bc[1][1] is None or bc[1][0] != bc[1][1]]
    barcodes.sort(key=lambda x: x[0], reverse=True)

    prev = 0
    for i in range(1, len(barcodes) + 1):
        if i == len(barcodes) or barcodes[i][0] != barcodes[i - 1][0]:
            barcodes[prev: i] = sorted(barcodes[prev: i], key=lambda x: x[1][1] - x[1][0], reverse=True)
            prev = i

    return barcodes


def betti_numbers(barcodes, filtration_val):
    # TODO i guess
    betti_nums = [filtration_val] + [0] * barcodes[0][0]

    for bc in barcodes:
        pass
