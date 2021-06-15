import numpy as np  # This is really only needed for np.int32


class Field:
    """
    Abstraction for fields.
    A field must satisfy a list of axioms. A list can be found here https://mathworld.wolfram.com/FieldAxioms.html
    In short this means for every field we must have a definition for addition, subtraction, multiplication, and
    division. To satisfy both subtraction and division we will define additive and multiplicative inverses. This also
    means we have to have a zero and an one element in our field, and these cannot be equal. The following class will
    act as a base class for field implementation and will be an object of operations such that for a field f. we can
    use the operations as follows
    ex. f.add(x, y)
    """

    @property
    def zero(self):
        """ additive zero element of field. x + 0 = 0 + x = x"""
        raise AssertionError("Not implemented")

    @property
    def one(self):
        """ multiplicative one element of field. x * 1 = 1 * x = x"""
        raise AssertionError("Not implemented")

    def equals(self, x, y):
        """ Checks is x and y are the same element of field"""
        raise AssertionError("Not implemented")

    def add(self, x, y):
        """ adds x and y"""
        raise AssertionError("Not implemented")

    def negate(self, x):
        """ additive inverse of the element x"""
        raise AssertionError("Not implemented")

    def subtract(self, x, y):
        """ acts as helper function to computer x + -y"""
        raise AssertionError("Not implemented")

    def multiply(self, x, y):
        """ multiply x and y"""
        raise AssertionError("Not implemented")

    def reciprocal(self, x):
        """ Multiplicative inverse element of field of the element x"""
        raise AssertionError("Not implemented")

    def divide(self, x, y):
        """ acts as helper function to computer x * y^{-1}"""
        raise AssertionError("Not implemented")

    def __eq__(self, other):
        """ Checks if two instances of field are the same field"""
        raise AssertionError("Not implemented")

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def data_type(self):
        """ underlying data type. This is only needed for integration into Matrix class. """
        raise AssertionError("Not implemented")


class Zp(Field):
    def __init__(self, p):
        """ Creates Z_p or Z/pZ where p is prime"""
        self.p = self.data_type(p)

    @property
    def data_type(self):
        return np.int32

    def _is_valid(self, x):
        """ x must be of type np.int32 and between 0 and p. I removed the validity check here to make things faster."""
        return x

    @property
    def zero(self):
        """ additive identity element of field"""
        return self.data_type(0)

    @property
    def one(self):
        """ multiplicative identity element of field"""
        return self.data_type(1)

    def equals(self, x, y):
        """ Checks is x and y are the same element of field"""
        return self._is_valid(x) == self._is_valid(y)

    def add(self, x, y):
        return (self._is_valid(x) + self._is_valid(y)) % self.p

    def negate(self, x):
        """ additive inverse element of field of the element x """
        return (-1 * self._is_valid(x)) % self.p

    def subtract(self, x, y):
        return (self._is_valid(x) - self._is_valid(y)) % self.p

    def multiply(self, x, y):
        return (self._is_valid(x) * self._is_valid(y)) % self.p

    def reciprocal(self, x):
        """ Multiplicative inverse element of field of the element x """
        def gcd_extended(a, b):
            if a == 0:
                return b, 0, 1

            gcd, x1, y1 = gcd_extended(b % a, a)

            xp = y1 - (b // a) * x1
            yp = x1

            return gcd, xp, yp

        one, x_inv, p = gcd_extended(self._is_valid(x), self.p)
        assert one == 1, "Extended Euclidean algorithm failed"
        return self.data_type(x_inv % self.p)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.p == other.p:
                return True
        return False
