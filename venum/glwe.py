from .rns import RnsBasis
from .crt import CrtEncoder
from .logging import logger

from sympy import Poly, GF
from sympy.abc import x
from sympy.polys.specialpolys import random_poly

import random
from dataclasses import dataclass


@dataclass
class EncryptionParameters:
    """
    Parameters for the encryption scheme.
    """

    dimension: int
    ciphertext_modulus: int
    plaintext_modulus: int
    noise_modulus: int
    seed: int = None

    def __post_init__(self):
        if (self.plaintext_modulus * self.noise_modulus
                >= self.ciphertext_modulus):
            raise ValueError(
                'Invalid parameters: plaintext_modulus * noise_modulus '
                '>= ciphertext_modulus')


class GlweSample:
    def __init__(self, mask, body):
        self.mask = mask
        self.body = body

    def __repr__(self):
        return f'GlweSample(mask={self.mask}, body={self.body})'


class GlweDistribution:
    def __init__(self, params: EncryptionParameters):
        if params.seed is not None:
            random.seed(params.seed)
            logger.warning(f"Setting random seed to {params.seed}")
        self.params = params
        self.plaintext_ring = GF(params.plaintext_modulus, symmetric=False)
        self.cipher_ring = GF(params.ciphertext_modulus, symmetric=False)
        crt_basis = RnsBasis(
            [self.params.plaintext_modulus, self.params.noise_modulus])
        self.crt_encoder = CrtEncoder(crt_basis, self.cipher_ring)

    def sample_polynomial(self, modulus=None):
        modulus = modulus or self.params.ciphertext_modulus
        poly = random_poly(x, self.params.dimension, 0,
                           modulus - 1).as_poly(domain=self.cipher_ring)
        return poly

    def sample_mask(self):
        return self.sample_polynomial()

    def sample_noise(self):
        # WARN: noise should be sampled from a subgaussian distribution
        return self.sample_polynomial(self.params.noise_modulus)

    def sample_zero_secret(self, secret: Poly):
        mask = self.sample_mask()
        noise = self.sample_noise()
        zero = Poly(0, x, domain=self.cipher_ring)
        crt_noise = self.crt_encoder.encode(zero, noise)

        # These operations happen in Z_q where q == ciphertext_modulus
        # Sympy's Poly class handles the modulo operation
        body = mask * secret + crt_noise
        return GlweSample(-mask, body)
