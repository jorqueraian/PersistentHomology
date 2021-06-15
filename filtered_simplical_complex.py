import matrixmath as mm


class FilteredSimplicialComplex(object):
    __complex = []
    __simplex_to_index_map = None
    __index_to_simplex_map = None
    __simplex_to_index_map_inverted = None
    __index_to_simplex_map_inverted = None
    field = None

    def __init__(self, simplicies, field):
        simplicies.sort(key=lambda x: x.filtration)
        prev = 0
        for i in range(1, len(simplicies) + 1):
            if i == len(simplicies) or simplicies[i][0] != simplicies[i - 1][0]:
                simplicies[prev: i] = sorted(simplicies[prev: i], key=lambda x: x.dimension)
                prev = i
        self.__complex = simplicies
        self.field = field

    def create_boundary_map(self):
        n = self.num_simplicies
        b_map = mm.create_matrix(n, n, self.field, None, self.__complex, self.__complex)

        for i, s in enumerate(self):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = self.__complex.index(b_s)
                b_map[b_s_j, i] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def create_coboundary_map(self):
        n = self.num_simplicies
        b_map = mm.create_matrix(n, n, self.field, None, self.__complex[::-1], self.__complex[::-1])

        for i, s in enumerate(self):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = self.num_simplicies - self.__complex.index(b_s) - 1
                b_map[n - i - 1, b_s_j] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def create_p_boundary_map(self, p):
        p_simplicies = self.get_p_simplicies(p)
        return self.construct_boundary_map_from_simplicies(p_simplicies)

    def create_p_coboundary_map(self, p):
        pmo_simplicies = self.get_p_simplicies(p-1)
        return self.construct_coboundary_map_from_simplicies(pmo_simplicies, p)

    def construct_boundary_map_from_simplicies(self, simplicies):
        # All simplicies in simplicies need to be the same dim
        pmo_simplicies = self.get_p_simplicies(simplicies[0].dimension - 1) if simplicies[0].dimension > 0 else [None]

        b_map = mm.create_matrix(max(1, len(pmo_simplicies)), len(simplicies), self.field, None, pmo_simplicies, simplicies)

        for i, s in enumerate(simplicies):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = pmo_simplicies.index(b_s)
                b_map[b_s_j, i] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def construct_coboundary_map_from_simplicies(self, simplicies, p, revered=False):
        # All simplicies in simplicies need to be the same dim
        p_simplicies = self.get_p_simplicies(p)
        if revered:
            simplicies = simplicies[::-1]

        b_map = mm.create_matrix(len(p_simplicies), len(simplicies), self.field, None, p_simplicies[::-1], simplicies[::-1])

        for i, s in enumerate(p_simplicies):
            b_map_s = s.boundary()
            for j, b_s in enumerate(b_map_s):
                b_s_j = len(simplicies) - (simplicies+[b_s]).index(b_s) - 1
                if b_s_j != -1:
                    b_map[len(p_simplicies) - i - 1, b_s_j] = (self.field.negate(self.field.one) if j % 2 != 0 else self.field.one)
        return b_map

    def get_p_simplicies(self, p):
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
        self.__index_to_simplex_map = [None]*self.num_simplicies
        for i, s in enumerate(self):
            self.__simplex_to_index_map[s] = i
            self.__index_to_simplex_map[i] = s

    def __inverted_maps(self):
        self.__simplex_to_index_map_inverted = {}
        self.__index_to_simplex_map_inverted = [None]*self.num_simplicies
        n = self.num_simplicies - 1
        for i, s in enumerate(self):
            self.__simplex_to_index_map_inverted[s] = n-i
            self.__index_to_simplex_map_inverted[n-i] = s

    @property
    def num_simplicies(self):
        return len(self.__complex)

    def get_filtration(self):
        return self.__complex

    def __getitem__(self, key):
        return self.__complex[key]
