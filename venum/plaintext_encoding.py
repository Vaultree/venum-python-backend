from sympy import Poly
from sympy.abc import x
from typing import Iterable
from abc import ABC


class Encoder(ABC):
    def encode(self, message: Iterable[int]) -> Poly:
        pass

    def decode(self, poly: Poly) -> Iterable[int]:
        pass


class PolynomialEncoder(Encoder):
    def __init__(self, dist):
        self.dist = dist

    def encode(self, message: Iterable[int]) -> Poly:
        return Poly(message, x, domain=self.dist.cipher_ring)

    def decode(self, poly: Poly) -> Iterable[int]:
        return poly.coeffs()


class BatchEncoder(Encoder):
    def __init__(self, dist):
        self.dist = dist

    def encode(self, message: Iterable[int]) -> Poly:
        pass

    def decode(self, poly: Poly) -> Iterable[int]:
        pass
