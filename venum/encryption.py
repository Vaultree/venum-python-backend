from .logging import logger
from .glwe import GlweSample, GlweDistribution
from .key import SecretKey, PublicKey, RelinKey
from .numeric import radix_decompose_poly

from sympy import Poly
from sympy.abc import x

from typing import Iterable


class Cipher:
    """
    A class representing a GLWE ciphertext.

    Attributes:
    - glwe_sample: A GlweSample object representing the ciphertext.
    """

    def __init__(self, glwe_sample):
        self.glwe_sample = glwe_sample

    def __repr__(self):
        return f'Cipher(mask={self.glwe_sample.mask}, '
        f'body={self.glwe_sample.body})'


class Encryptor:
    """
    A class for handling encryption and decryption of messages.
    """

    def __init__(self, dist: GlweDistribution, plaintext_encoder):
        """
        Initializes an Encryptor object.

        Args:
        - dist: A GlweDistribution object representing the distribution
          used for encryption.
        - plaintext_encoder: An object that encodes and decodes messages
          according to the `venum.plaintext_encoding.Encoder` interface.
        """

        self._dist = dist
        self.plaintext_encoder = plaintext_encoder

    @property
    def dist(self):
        return self._dist

    def encrypt(self, pk: PublicKey, message: Iterable[int],
                plaintext_encoder=None) -> Cipher:
        """
        Encrypts a message.

        Args:
        - pk: A PublicKey object representing the public key.
        - message: An iterable of integers representing the message.
        - plaintext_encoder: An object that encodes and decodes messages
          according to the `venum.plaintext_encoding.Encoder` interface.
          If None, the default encoder is used.

        Returns:
        - A Cipher object representing the encrypted message.
        """

        logger.debug(f'Encrypting message: {message}')

        plaintext_encoder = plaintext_encoder or self.plaintext_encoder

        message = plaintext_encoder.encode(message)
        logger.debug(f'encoded message: {message}')

        crt_message = self.dist.crt_encoder.encode_pure_message(
            message).set_domain(self.dist.cipher_ring)
        logger.debug(f'crt_message: {crt_message}')

        crt_noise1 = (self.dist.sample_crt_noise()
                      .set_domain(self.dist.cipher_ring))
        crt_noise2 = (self.dist.sample_crt_noise()
                      .set_domain(self.dist.cipher_ring))

        u = self.dist.sample_polynomial(modulus=2)
        logger.debug(f'sampled u: {u}')

        logger.debug(f"using public key: {pk.glwe_sample}")
        sample = GlweSample.compute_sample(
            mask=pk.glwe_sample.mask,
            mask_noise=crt_noise2,
            body=pk.glwe_sample.body,
            body_noise=crt_noise1,
            message=crt_message,
            u=u,
            poly_modulus=self.dist.poly_modulus,
        )
        return Cipher(sample)

    def decrypt(self, sk: SecretKey, cipher: Cipher) -> Iterable[int]:
        """
        Decrypts a ciphertext.

        Args:
        - sk: A SecretKey object representing the secret key.
        - cipher: A Cipher object representing the ciphertext.

        Returns:
        - An iterable of integers representing the decrypted message.
        """

        logger.debug(f"{cipher}")
        cipher_mask = cipher.glwe_sample.mask
        cipher_body = cipher.glwe_sample.body

        crt_message = (cipher_body + cipher_mask * sk.secret_poly)
        crt_message = crt_message % self.dist.poly_modulus
        logger.debug(f"{crt_message}")
        noisy_message = self.dist.crt_encoder.decode(crt_message)
        logger.debug(f"{noisy_message}")

        message_poly = Poly.from_list([rns[0] for rns in noisy_message],
                                      x, domain=self.dist.plaintext_ring)

        logger.debug(f"{message_poly}")
        return self.plaintext_encoder.decode(message_poly)


class Rank2Cipher:
    """
    A class representing a rank-2 ciphertext, representing a non-normalized
    Cipher containing squared terms. Usually produced as an intermediate step
    during homomorphic multiplication.

    Attributes:
    - constant: A Poly object representing the constant term over the secret.
    - linear: A Poly object representing the linear term over the secret.
    - quadratic: A Poly object representing the quadratic term over the secret.
    """

    def __init__(self, constant: Poly, linear: Poly, quadratic: Poly):
        self.constant = constant
        self.linear = linear
        self.quadratic = quadratic

    def relinearize(self, relin_key: RelinKey, poly_modulus) -> Cipher:
        """
        Relinearizes the rank-2 ciphertext into a normalized Cipher.

        Args:
        - relin_key: A RelinKey object representing the relinearization key.
        - poly_modulus: A Poly object representing the polynomial modulus.

        Returns:
        - A Cipher object representing the relinearized ciphertext.
        """

        cipher_ring = poly_modulus.domain
        quad_decomposed = radix_decompose_poly(
            poly=self.quadratic,
            radix=relin_key.base,
            num_components=relin_key.digit_count(),
            domain=cipher_ring
        )
        mask = Poly([0], x, domain=cipher_ring)
        body = Poly([0], x, domain=cipher_ring)
        for aux_key, component in zip(relin_key.aux_keys,
                                      quad_decomposed):
            mask += aux_key.mask * component
            mask = mask % poly_modulus
            body += aux_key.body * component
            body = body % poly_modulus
        mask += self.linear
        body += self.constant
        mask = mask % poly_modulus
        body = body % poly_modulus
        return Cipher(GlweSample(mask=mask, body=body))
