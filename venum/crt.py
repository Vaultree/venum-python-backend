from .rns import Rns
from .logging import logger

from sympy import Poly
from sympy.abc import x


class CrtEncoder:
    def __init__(self, basis, plaintext_ring):
        """
        CRT encoder for encoding and decoding polynomials
        with respect to a pair of moduli.
        CRT encoding is described on https://eprint.iacr.org/2024/1105.pdf.

        Args:
        - basis: a pair of moduli for the CRT encoding
        - plaintext_ring: the ring of the plaintext polynomials
        """

        if len(basis) != 2:
            raise ValueError("CRT encoding requires two moduli")
        self.basis = basis
        self.plaintext_ring = plaintext_ring

    def _encode_coef(self, message, noise):
        return Rns(self.basis, [message, noise]).to_int()

    def _decode_coef(self, value):
        return self.basis.to_rns(value)

    def _normalized_coeffs(self, message: Poly, noise: Poly):
        message_coeffs = (message.all_coeffs()
                          if not message.is_zero
                          else [0] * len(noise.coeffs()))
        noise_coeffs = (noise.all_coeffs()
                        if not noise.is_zero
                        else [0] * len(message.coeffs()))
        diff_len = len(message_coeffs) - len(noise_coeffs)
        if diff_len > 0:
            noise_coeffs += [0] * diff_len
        elif diff_len < 0:
            message_coeffs += [0] * -diff_len
        logger.debug(f'message_coeffs: {message_coeffs}')
        logger.debug(f'noise_coeffs: {noise_coeffs}')
        return zip(message_coeffs, noise_coeffs)

    def encode(self, message: Poly, noise: Poly):
        """
        Encode a message and noise polynomial into a single polynomial
        using the CRT encoding.

        Args:
        - message: the message polynomial
        - noise: the noise polynomial

        Returns:
        - a CRT-encoded polynomial
        """

        logger.debug(f'CRT encoding message: {message} with noise: {noise}')
        msg_noise_pairs = self._normalized_coeffs(message, noise)
        coefs = (self._encode_coef(msg_coef, noise_coef)
                 for (msg_coef, noise_coef) in msg_noise_pairs)
        return Poly(coefs, x, domain=self.plaintext_ring)

    def _encode_with_zero(self, poly: Poly, component: int):
        zero = Poly(0, x, domain=self.plaintext_ring)
        if component == 0:
            return self.encode(poly, zero)
        elif component == 1:
            return self.encode(zero, poly)
        else:
            raise ValueError("component must be 0 or 1")

    def encode_pure_message(self, message: Poly):
        """
        Encode a message polynomial with zero noise.

        Args:
        - message: the message polynomial

        Returns:
        - a CRT-encoded polynomial
        """

        return self._encode_with_zero(message, 0)

    def encode_pure_noise(self, noise: Poly):
        """
        Encode a noise polynomial with zero message.

        Args:
        - noise: the noise polynomial

        Returns:
        - a CRT-encoded polynomial
        """

        return self._encode_with_zero(noise, 1)

    def decode(self, poly: Poly):
        """
        Decode a CRT-encoded polynomial into its message and noise components.
        """

        logger.debug(f'CRT decoding polynomial: {poly}')
        return [self._decode_coef(coeff)
                for coeff in poly.coeffs()]
