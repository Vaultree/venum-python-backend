from sympy import Poly
from sympy.abc import x

from typing import Iterable


def nth_digit(number, radix, n):
    """
    Return the n-th digit of a number in a given radix.

    Args:
    - number: int, the number to extract the digit from.
    - radix: int, the base of the number system.
    - n: int, the index of the digit to extract.j

    Returns:
    - int, the n-th digit of the number in the given radix.

    Raises:
    - ValueError: if number is negative, radix is less than 2, or n is negative
    """

    if number < 0:
        raise ValueError("Number must be non-negative.")
    if radix < 2:
        raise ValueError("Radix must be at least 2.")
    if n < 0:
        raise ValueError("Index n must be non-negative.")

    for _ in range(n):
        number //= radix

    return number % radix


def radix_decompose_poly(poly: Poly, radix: int,
                         num_components: int, domain) -> Iterable[Poly]:
    """
    Decompose a polynomial into its components in a given radix.
    The components are obtained by extracting the digits of the coefficients
    and constructing a new polynomial from them.

    Args:
    - poly: Poly, the polynomial to decompose.
    - radix: int, the base of the number system.
    - num_components: int, the number of components to extract.
    - domain: Domain, the domain of the polynomial.

    Returns:
    - Iterable[Poly], the components of the polynomial.
    """

    for n in range(num_components):
        coeffs = (nth_digit(coef, radix, n)
                  for coef in reversed(poly.all_coeffs()))
        decomposed = Poly(coeffs, x, domain=domain)
        yield decomposed
