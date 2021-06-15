import fieldmath
import numpy as np


class Matrix:
    # This class adds a few features to the np.array
    __data = None

    def __init__(self, matrix_values, field, row_bases=None, column_bases=None):
        """ initializes Matrix from list like param matrix_values over specified field """
        if not isinstance(field, fieldmath.Field):
            raise TypeError("Unknown Field")

        self.f = field
        self.data_type = self.f.data_type
        self.__data = np.array(matrix_values, dtype=self.data_type)
        self.rows, self.columns = self.__data.shape

        self.row_bases = row_bases
        self.column_bases = column_bases

        if self.row_bases is not None or self.column_bases is not None:
            assert self.row_bases is not None and self.column_bases is not None
            assert len(self.row_bases) == self.rows and len(self.column_bases) == self.columns

        self.__low_vals = [None]*self.columns

    def __getitem__(self, key):
        return self.__data[key]

    def __setitem__(self, key, value):
        if isinstance(key[1], slice):
            self.__low_vals[key[1]] = [None]*((key[1].stop-key[1].start)//key[1].step)
        elif isinstance(key[1], int):
            self.__low_vals[key[1]] = None
        self.__data[key] = value

    def __mul__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError()
        if self.f != other.f:
            raise Exception("Fields do not align.")

        if self.columns != other.rows:
            raise Exception("Can not multiple matrices, inner dimensions do no align.")

        # Is this right? idk
        result = create_matrix(self.rows, other.columns, self.f, None, row_bases=other.row_bases, column_bases=self.column_bases)
        for r in range(result.rows):
            for c in range(result.columns):
                val = self.f.zero
                for i in range(self.columns):
                    val = self.f.add(val, self.f.multiply(self[r, i], other[i, c]))
                result[r, c] = val
        return result

    def __str__(self):
        delim = 5
        matrix_str = "    "
        if self.column_bases is not None:
            col_bases_str = [str(s) for s in self.column_bases]
            row_bases_str = [str(s) for s in self.row_bases]
            delim = max([len(cbs) for cbs in col_bases_str])
            row_bases_delim = max([len(rbs) for rbs in row_bases_str]) + 1

            matrix_str += " "*row_bases_delim + " ".join(str(val) + " " * (delim - len(val)) for val in col_bases_str) + " \n    "
        for (i, row) in enumerate(self):
            if i > 0:
                matrix_str += " \n    "
            if self.column_bases is not None:
                matrix_str += row_bases_str[i] + ": " + " "*(row_bases_delim-len(row_bases_str[i]))
            matrix_str += " ".join(str(val) + " "*(delim-len(str(val))) for val in row)
        return matrix_str + " \n"

    def __eq__(self, other):
        if self.f != other.f:
            return False
        if self.__data.shape != other.__data.shape:
            return False

        for ai, bi in zip(self.__data.flat, other.__data.flat):
            if ai != bi:
                return False
        if self.row_bases != other.row_bases:
            return False
        if self.column_bases != other.column_bases:
            return False
        return True

    def copy(self):
        result = self.__class__(np.array(self.__data), self.f, row_bases=self.row_bases, column_bases=self.column_bases)
        return result

    def low(self, c_ind):
        if self.__low_vals[c_ind] is not None:
            return self.__low_vals[c_ind]
        for j in range(self.rows-1, -1, -1):
            if self[j, c_ind] != self.f.zero:
                self.__low_vals[c_ind] = j
                return j
        self.__low_vals[c_ind] = -1
        return -1

    def anti_transpose(self):
        if self.row_bases is not None and self.column_bases is not None:
            col_bases = self.row_bases[::-1]
            row_bases = self.column_bases[::-1]
            result = self.__class__(self.__data[::-1, ::-1].T, self.f, row_bases=row_bases, column_bases=col_bases)
        else:
            result = self.__class__(self.__data[::-1, ::-1].T, self.f)
        return result

    def remove_column(self, c):
        # if c < 0 or c >= self.columns:
        #     raise IndexError("Index for c out of bounds")
        self.__low_vals.pop(c)
        if self.column_bases is not None:
            self.column_bases.pop(c)
        self.columns -= 1
        np.delete(self.__data, c, 1)

    def add_columns(self, c1, c2, multiplier):
        """ c2 = c2 + c1 * multiplier """
        # if c1 < 0 or c2 < 0 or c1 >= self.columns or c2 >= self.columns:
        #     raise IndexError("Index for c1 or c2 out of bounds")
        max_row = max(self.low(c1), self.low(c2))
        new_low = None
        for r in range(max_row, -1, -1):
            self[r, c2] = self.f.add(self[r, c2], self.f.multiply(multiplier, self[r, c1]))
            if new_low is None and not self.f.equals(self[r, c2], self.f.zero):
                new_low = r
        self.__low_vals[c2] = new_low

    def clear_column(self, c):
        # if c < 0 or c >= self.columns:
        #     raise IndexError("Index for c out of bounds")
        self[:, c] = [self.f.zero]*self.rows
        self.__low_vals[c] = -1


def create_matrix(num_rows, num_columns, field, init_val=None, row_bases=None, column_bases=None):
    if init_val is None:
        init_val = field.zero
    return Matrix([[init_val]*num_columns]*num_rows, field, row_bases=row_bases, column_bases=column_bases)
