import matrixmath as mm


class FilteredSimplicialComplex(object):
    __complex = []
    __simplex_to_index_map = None
    __index_to_simplex_map = None
    __simplex_to_index_map_inverted = None
    __index_to_simplex_map_inverted = None
    field = None

    def __init__(self, simplices, field):
        simplices.sort(key=lambda x: x.filtration)
        prev = 0
        for i in range(1, len(simplices) + 1):
            if i == len(simplices) or simplices[i][0] != simplices[i - 1][0]:
                simplices[prev: i] = sorted(simplices[prev: i], key=lambda x: x.dimension)
                prev = i
        self.__complex = simplices
        self.field = field

    def create_boundary_map(self):
        n = self.num_simplices
        b_map = mm.create_matrix(n, n, self.field, None, self.__complex, self.__complex)

        for i, s in enumerate(self):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = self.__complex.index(b_s)
                b_map[b_s_j, i] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def create_coboundary_map(self):
        n = self.num_simplices
        b_map = mm.create_matrix(n, n, self.field, None, self.__complex[::-1], self.__complex[::-1])

        for i, s in enumerate(self):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = self.num_simplices - self.__complex.index(b_s) - 1
                b_map[n - i - 1, b_s_j] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def create_p_boundary_map(self, p):
        p_simplices = self.get_p_simplices(p)
        return self.construct_boundary_map_from_simplices(p_simplices)

    def create_p_coboundary_map(self, p):
        pmo_simplices = self.get_p_simplices(p - 1)
        return self.construct_coboundary_map_from_simplices(pmo_simplices, p)

    def construct_boundary_map_from_simplices(self, simplices):
        # All simplices in simplices need to be the same dim
        pmo_simplices = self.get_p_simplices(simplices[0].dimension - 1) if simplices[0].dimension > 0 else [None]

        b_map = mm.create_matrix(max(1, len(pmo_simplices)), len(simplices), self.field, None, pmo_simplices, simplices)

        for i, s in enumerate(simplices):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = pmo_simplices.index(b_s)
                b_map[b_s_j, i] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def construct_coboundary_map_from_simplices(self, simplices, p, revered=False):
        # All simplices in simplices need to be the same dim
        p_simplices = self.get_p_simplices(p)
        if revered:
            simplices = simplices[::-1]

        b_map = mm.create_matrix(len(p_simplices), len(simplices), self.field, None, p_simplices[::-1], simplices[::-1])

        for i, s in enumerate(p_simplices):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = len(simplices) - (simplices + [b_s]).index(b_s) - 1
                if b_s_j != -1:
                    b_map[len(p_simplices) - i - 1, b_s_j] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def get_p_simplices(self, p):
        return [s for s in self if s.dimension == p]

    @property
    def simplex_to_index_map(self):
        if self.__simplex_to_index_map is None:
            self.__maps()
        return self.__simplex_to_index_map

    @property
    def index_to_simplex_map(self):
        if self.__index_to_simplex_map is None:
            self.__maps()
        return self.__index_to_simplex_map

    @property
    def simplex_to_index_map_inverted(self):
        if self.__simplex_to_index_map_inverted is None:
            self.__inverted_maps()
        return self.__simplex_to_index_map_inverted

    @property
    def index_to_simplex_map_inverted(self):
        if self.__index_to_simplex_map_inverted is None:
            self.__inverted_maps()
        return self.__index_to_simplex_map_inverted

    def __maps(self):
        self.__simplex_to_index_map = {}
        self.__index_to_simplex_map = [None]*self.num_simplices
        for i, s in enumerate(self):
            self.__simplex_to_index_map[s] = i
            self.__index_to_simplex_map[i] = s

    def __inverted_maps(self):
        self.__simplex_to_index_map_inverted = {}
        self.__index_to_simplex_map_inverted = [None]*self.num_simplices
        n = self.num_simplices - 1
        for i, s in enumerate(self):
            self.__simplex_to_index_map_inverted[s] = n-i
            self.__index_to_simplex_map_inverted[n-i] = s

    @property
    def num_simplices(self):
        return len(self.__complex)

    def get_filtration(self):
        return self.__complex

    def __getitem__(self, key):
        return self.__complex[key]
