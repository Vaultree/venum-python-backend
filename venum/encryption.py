from .logging import logger
from .glwe import GlweSample
from .key import SecretKey, PublicKey

from typing import Iterable

from sympy import Poly
from sympy.abc import x


class Cipher:
    def __init__(self, glwe_sample):
        self.glwe_sample = glwe_sample

    def __repr__(self):
        return f'Cipher(mask={self.glwe_sample.mask}, '
        f'body={self.glwe_sample.body})'


class Encryptor:
    def __init__(self, dist, encoder_cls):
        self._dist = dist
        self.plaintext_encoder = encoder_cls(dist)

    @property
    def dist(self):
        return self._dist

    def encrypt(self, pk: PublicKey, message: Iterable[int],
                plaintext_encoder=None) -> Cipher:
        logger.debug(f'Encrypting message: {message}')
        plaintext_encoder = plaintext_encoder or self.plaintext_encoder
        message = plaintext_encoder.encode(message)
        logger.debug(f'encoded message: {message}')
        zero = Poly(0, x, domain=self.dist.cipher_ring)

        crt_message = self.dist.crt_encoder.encode(message, zero)
        logger.debug(f'crt_message: {crt_message}')

        noise1 = self.dist.sample_noise()
        noise2 = self.dist.sample_noise()
        crt_noise1 = self.dist.crt_encoder.encode(zero, noise1)
        logger.debug(f'crt_noise1: {crt_noise1}')
        crt_noise2 = self.dist.crt_encoder.encode(zero, noise2)
        logger.debug(f'crt_noise2: {crt_noise2}')

        u = self.dist.sample_polynomial(modulus=2)
        logger.debug(f'u: {u}')

        # These operations happen in Z_q where q == ciphertext_modulus
        # Sympy's Poly class handles the modulo operation
        logger.debug(f"{pk.glwe_sample}")
        mask = pk.glwe_sample.mask * u + crt_noise2
        logger.debug(f"mask: {mask}")
        body = pk.glwe_sample.body * u + crt_noise1 + crt_message
        logger.debug(f"body: {body}")
        return Cipher(GlweSample(mask, body))

    def decrypt(self, sk: SecretKey, cipher: Cipher) -> Iterable[int]:
        # This operation happen in Z_q where q == ciphertext_modulus
        # Sympy's Poly class handles the modulo operation
        logger.debug(f"{cipher}")
        crt_message = (cipher.glwe_sample.body +
                       cipher.glwe_sample.mask * sk.secret_poly)
        logger.debug(f"{crt_message}")
        noisy_message = self.dist.crt_encoder.decode(crt_message)
        logger.debug(f"{noisy_message}")

        message_poly = Poly.from_list([rns[0] for rns in noisy_message],
                                      x, domain=self.dist.plaintext_ring)

        logger.debug(f"{message_poly}")
        return self.plaintext_encoder.decode(message_poly)
