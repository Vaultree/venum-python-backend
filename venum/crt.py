from .rns import Rns
from .logging import logger

from sympy import Poly
from sympy.abc import x


class CrtEncoder:
    def __init__(self, basis, cipher_ring):
        if len(basis) != 2:
            raise ValueError("CRT encoding requires two moduli")
        self.basis = basis
        self.cipher_ring = cipher_ring

    def _encode_coef(self, message, noise):
        return Rns.from_residues(self.basis, [message, noise]).to_int()

    def _decode_coef(self, value):
        return self.basis.to_rns(value)

    def encode(self, message: Poly, noise: Poly):
        logger.debug(f'CRT encoding message: {message} with noise: {noise}')
        message_coeffs = (message.coeffs()
                          if not message.is_zero
                          else [0] * len(noise.coeffs()))
        noise_coeffs = (noise.coeffs()
                        if not noise.is_zero
                        else [0] * len(message.coeffs()))
        msg_noise_pairs = zip(message_coeffs, noise_coeffs)
        coefs = (self._encode_coef(msg_coef, noise_coef)
                 for (msg_coef, noise_coef) in msg_noise_pairs)
        return Poly(coefs, x, domain=self.cipher_ring)

    def decode(self, poly: Poly):
        logger.debug(f'CRT decoding polynomial: {poly}')
        return [self._decode_coef(coeff)
                for coeff in poly.coeffs()]
