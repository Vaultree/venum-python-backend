from .glwe import GlweDistribution, GlweSample
from .logging import logger

from sympy import Poly


class Key:
    pass


class SecretKey(Key):
    def __init__(self, dist: GlweDistribution, secret_poly: Poly):
        self._dist = dist
        self.secret_poly = secret_poly

    def __repr__(self):
        return f'SecretKey({self.secret_poly})'

    @classmethod
    def rand(cls, dist: GlweDistribution, modulus=None):
        secret = dist.sample_polynomial(modulus)
        logger.debug(f"Generating random secret key from: {secret}")
        return cls(dist, secret)

    @property
    def dist(self):
        return self._dist


class PublicKey(Key):
    def __init__(self, glwe_sample: GlweSample):
        self.glwe_sample = glwe_sample

    def __repr__(self):
        return f'PublicKey(mask={self.glwe_sample.mask}, '
        f'body={self.glwe_sample.body})'

    @classmethod
    def from_secret_key(cls, secret_key: SecretKey):
        sample = secret_key.dist.sample_zero_secret(secret_key.secret_poly)
        logger.debug(f"Generating public key with sample: {sample}")
        return cls(sample)


def gen_key_pair(dist: GlweDistribution, modulus=None):
    sk = SecretKey.rand(dist, modulus)
    pk = PublicKey.from_secret_key(sk)
    return sk, pk


class RelinKey(Key):
    def __init__(self):
        raise NotImplementedError

    @classmethod
    def from_secret_key(cls, secret_key: SecretKey):
        raise NotImplementedError
