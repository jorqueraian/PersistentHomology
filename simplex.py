class Simplex(object):
    __data_tuple = None
    __filtration = 0

    def __init__(self, vertices, filtration=0):
        assert isinstance(vertices, list)
        vertices.sort()

        self.__data_tuple = tuple(vertices)
        self.__filtration = filtration

    def boundary(self):
        b = []
        if self.dimension < 1:
            return []
        for v in range(self.dimension + 1):
            b.append(self[:v] + self[v + 1:])
        return b

    def __getitem__(self, key):
        return self.__data_tuple[key]

    def __hash__(self):
        return hash(self.__data_tuple)

    def __eq__(self, other):
        return self[:] == other[:]

    def __str__(self):
        return str(self.simplex)

    @property
    def dimension(self):
        return len(self.__data_tuple) - 1

    @property
    def filtration(self):
        return self.__filtration

    @property
    def simplex(self):
        return self.__data_tuple
