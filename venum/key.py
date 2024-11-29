from .glwe import GlweDistribution, GlweSample
from .logging import logger

from sympy import Poly

import math
from typing import Iterable, Tuple


class SecretKey:
    """
    A secret key for the scheme.

    Attributes:
    - dist: The GLWE distribution used to generate the secret key.
    - secret_poly: a secret polynomial.
    """

    def __init__(self, dist: GlweDistribution, secret_poly: Poly):
        self._dist = dist
        self.secret_poly = secret_poly

    def __repr__(self):
        return f'SecretKey({self.secret_poly})'

    @classmethod
    def rand(cls, dist: GlweDistribution, modulus=None) -> 'SecretKey':
        """
        Generates a random secret key.

        Args:
        - dist: A GLWE distribution to use.
        - modulus: The modulus to use for the secret key. If None, the modulus
          of the distribution is used.

        Returns:
        - A secret key.
        """

        secret = dist.sample_polynomial(modulus)
        logger.debug(f"Generating random secret key from: {secret}")
        return cls(dist, secret)

    @property
    def dist(self):
        return self._dist


class PublicKey:
    """
    A public key for the scheme.

    Attributes:
    - glwe_sample: A GLWE sample representing the public key.
    """

    def __init__(self, glwe_sample: GlweSample):
        self.glwe_sample = glwe_sample

    def __repr__(self):
        return f'PublicKey(mask={self.glwe_sample.mask}, '
        f'body={self.glwe_sample.body})'

    @classmethod
    def from_secret_key(cls, secret_key: SecretKey) -> 'PublicKey':
        """
        Generates a public key from a secret key.

        Args:
        - secret_key: The secret key to derive the public key from.

        Returns:
        - A public key.
        """

        sample = secret_key.dist.sample_zero_secret(secret_key.secret_poly)
        logger.debug(f"Generating public key with sample: {sample}")
        return cls(sample)


def gen_key_pair(dist: GlweDistribution,
                 modulus=None) -> Tuple[SecretKey, PublicKey]:
    """
    Generates a key pair.

    Args:
    - dist: The GLWE distribution to use.
    - modulus: The modulus to use for the secret key. If None, the modulus of
      the distribution is used.

    Returns:
    - A tuple (sk, pk) where sk is the secret key and pk is the public key.
    """

    sk = SecretKey.rand(dist, modulus)
    pk = PublicKey.from_secret_key(sk)
    return sk, pk


class RelinKey:
    """
    A relinearization key for the scheme. Used to normalize ciphertexts.

    Attributes:
    - aux_keys: A list of auxiliary keys as described in the reference paper.
    - base: The base/radix used to generate the auxiliary keys.
    """

    def __init__(self, aux_keys: Iterable[GlweSample], base: int):
        self.aux_keys = aux_keys
        self.base = base

    @staticmethod
    def _compute_aux_keys(sk: SecretKey, base: int) -> Iterable[GlweSample]:
        digit_count = math.log(sk.dist.params.ciphertext_modulus, base)
        digit_count = math.ceil(digit_count)
        aux_keys = []
        sk2 = sk.secret_poly ** 2
        for i in range(digit_count):
            mask = sk.dist.sample_mask()
            crt_noise = (sk.dist.sample_crt_noise()
                         .set_domain(sk.dist.cipher_ring))
            masked_secret = mask * sk.secret_poly
            masked_secret = masked_secret % sk.dist.poly_modulus
            noisy_secret = masked_secret + crt_noise
            noisy_secret = noisy_secret % sk.dist.poly_modulus
            message = base ** i * sk2
            message = message % sk.dist.poly_modulus
            body = (noisy_secret + message) % sk.dist.poly_modulus
            aux_keys.append(GlweSample(mask=-mask, body=body))
        return aux_keys

    @classmethod
    def from_secret_key(cls, secret_key: SecretKey,
                        base: int = 2) -> 'RelinKey':
        """
        Generates a relinearization key from a secret key.

        Args:
        - secret_key: The secret key to derive the relinearization key from.
        - base: The base to use for the relinearization key. Defaults to 2.

        Returns:
        - A relinearization key.
        """

        aux_keys = cls._compute_aux_keys(secret_key, base)
        return cls(aux_keys, base)

    def digit_count(self):
        """
        The number of parts in the relinearization key decomposition.
        """

        return len(self.aux_keys)
