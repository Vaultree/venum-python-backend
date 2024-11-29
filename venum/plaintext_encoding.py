from sympy import Poly
from sympy.abc import x
from typing import Iterable

from abc import ABC


class Encoder(ABC):
    """
    Interface for encoding and decoding plaintexts.
    """

    def encode(self, message: Iterable[int]) -> Poly:
        pass

    def decode(self, poly: Poly) -> Iterable[int]:
        pass


class PolynomialEncoder(Encoder):
    """
    Encodes messages as a polynomials by mapping each element to a
    coefficient of the polynomial in increasing order of degree.
    """

    def __init__(self, dist):
        """
        Initializes the PolynomialEncoder with the given GLWE distribution.
        """

        self.dist = dist

    def encode(self, message: Iterable[int]) -> Poly:
        """
        Encodes the given message as a polynomial by mapping each element to a
        coefficient of the polynomial in increasing order of degree.

        Args:
        - message: The message to encode.

        Returns:
        - A polynomial representing the message.
        """

        return Poly(reversed(message), x, domain=self.dist.plaintext_ring)

    def decode(self, poly: Poly) -> Iterable[int]:
        """
        Decodes the given polynomial by extracting the coefficients and
        returning them as an iterable.

        Args:
        - poly: The polynomial to decode.

        Returns:
        - The message represented by the polynomial.
        """

        q = self.dist.params.ciphertext_modulus
        p0 = self.dist.params.plaintext_modulus
        p1 = self.dist.params.noise_modulus
        p0p1 = p0 * p1
        k = (q // (2 * p0p1)) * p0p1

        coeffs = reversed(poly.all_coeffs())

        coeffs = list((((coef + k) % q) % p0p1) % p0
                      for coef in coeffs)
        if len(coeffs) < self.dist.params.dimension:
            coeffs += [0] * (self.dist.params.dimension - len(coeffs))
        return coeffs


class BatchEncoder(Encoder):
    def __init__(self, dist):
        self.dist = dist

    def encode(self, message: Iterable[int]) -> Poly:
        raise NotImplementedError

    def decode(self, poly: Poly) -> Iterable[int]:
        raise NotImplementedError
