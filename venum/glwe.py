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
    """
    A sample from the GLWE distribution.

    The sample is a pair of polynomials (mask, body).
    """

    def __init__(self, mask, body):
        self.mask = mask
        self.body = body

    @staticmethod
    def _compute_mask(mask, noise, u, poly_modulus):
        new_mask = mask * u + noise
        new_mask = new_mask % poly_modulus
        return new_mask

    @staticmethod
    def _compute_body(body, noise, message, u, poly_modulus):
        new_body = body * u + message + noise
        new_body = new_body % poly_modulus
        return new_body

    @classmethod
    def compute_sample(cls, mask, mask_noise, body, body_noise,
                       message, u, poly_modulus):
        """
        Compute a new sample from the given parameters. The new sample
        corresponds to an encryption of the given message.
        """

        new_mask = cls._compute_mask(mask, mask_noise, u,
                                     poly_modulus)
        new_body = cls._compute_body(body, body_noise, message, u,
                                     poly_modulus)
        return cls(mask=new_mask, body=new_body)

    @classmethod
    def _compute_zero_sample(
            cls, mask: Poly, secret: Poly,
            crt_noise: Poly, poly_modulus: Poly):
        body = (mask * secret +
                crt_noise.set_domain(poly_modulus.domain)) % poly_modulus
        return cls(mask=-mask, body=body)

    def __repr__(self):
        return f'GlweSample(mask={self.mask}, body={self.body})'


class GlweDistribution:
    def __init__(self, params: EncryptionParameters):
        """
        Initialize the GLWE distribution with the given parameters.

        Args:
        - params: the encryption parameters.
        """

        if params.seed is not None:
            random.seed(params.seed)
            logger.warning(f"Setting random seed to {params.seed}")
        self.params = params
        self.plaintext_ring = GF(params.plaintext_modulus, symmetric=False)
        self.cipher_ring = GF(params.ciphertext_modulus, symmetric=False)
        self.poly_modulus = Poly(
            x ** params.dimension + 1, x, domain=self.cipher_ring)
        crt_basis = RnsBasis(
            [self.params.plaintext_modulus, self.params.noise_modulus])
        self.crt_encoder = CrtEncoder(crt_basis, self.plaintext_ring)

    def sample_polynomial(self, modulus=None):
        """
        Sample a polynomial with coefficients in the given modulus.

        Args:
        - modulus: the modulus for the coefficients. If None, the
            ciphertext modulus is used.

        Returns:
        - a polynomial with coefficients in the given modulus.
        """

        modulus = modulus or self.params.ciphertext_modulus
        degree = self.params.dimension - 1
        return random_poly(
            x,
            n=degree,
            inf=0,
            sup=modulus - 1
        ).as_poly(domain=self.cipher_ring)

    def sample_mask(self):
        """
        Sample a mask polynomial.

        Returns:

        - a polynomial with coefficients in the ciphertext modulus.
        """

        return self.sample_polynomial()

    def sample_noise(self):
        """
        Sample a noise polynomial.

        Returns:

        - a polynomial with coefficients in the noise modulus.
        """

        # FIX: for security reasons, noise should not be sampled from
        # a uniform distribution.
        return self.sample_polynomial(self.params.noise_modulus)

    def sample_crt_noise(self):
        """
        Sample a CRT noise polynomial.

        Returns:

        - a CRT-encoded noise polynomial.
        """

        noise = self.sample_noise()
        logger.debug(f"Sampled CRT noise: {noise}")
        crt_noise = self.crt_encoder.encode_pure_noise(noise)
        logger.debug(f"CRT noise: {crt_noise}")
        return crt_noise

    def sample_zero_secret(self, secret: Poly):
        """
        Produces a random GLWE sample corresponding to an encryption of
        zero message.

        Args:
        - secret: the secret polynomial to use for the encryption.

        Returns:
        - a GLWE sample corresponding to an encryption of zero.
        """

        mask = self.sample_mask()
        crt_noise = self.sample_crt_noise()
        return GlweSample._compute_zero_sample(
            mask, secret, crt_noise, self.poly_modulus)
