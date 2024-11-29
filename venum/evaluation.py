from .logging import logger
from .glwe import GlweDistribution, GlweSample
from .key import RelinKey
from .encryption import Cipher, Rank2Cipher


class Evaluator:
    """
    Evaluator class for performing arithmetic operations on encrypted data.

    Attributes:
    - dist: GlweDistribution, the distribution used for encryption
    - relin_key: RelinKey, the relinearization key used for homomorphic
      multiplication. If not provided, multiplication will raise an error.
    """

    def __init__(self, dist: GlweDistribution, relin_key: RelinKey = None):
        logger.debug(f"Initializing Evaluator with {dist} and {relin_key}")
        self._dist = dist
        self.relin_key = relin_key

    @property
    def dist(self):
        return self._dist

    def add(self, lhs: Cipher, rhs: Cipher):
        """
        Add two ciphertexts together.

        Args:
        - lhs: Cipher, the left-hand side of the addition
        - rhs: Cipher, the right-hand side of the addition

        Returns:
        - Cipher, the sum of the two ciphertexts
        """

        logger.debug(f"Adding {lhs} and {rhs}")
        mask = lhs.glwe_sample.mask + rhs.glwe_sample.mask
        body = lhs.glwe_sample.body + rhs.glwe_sample.body
        return Cipher(GlweSample(mask=mask, body=body))

    def sub(self, lhs: Cipher, rhs: Cipher):
        """
        Subtract one ciphertext from another.

        Args:
        - lhs: Cipher, the left-hand side of the subtraction
        - rhs: Cipher, the right-hand side of the subtraction

        Returns:
        - Cipher, the difference of the two ciphertexts
        """

        logger.debug(f"Subtracting {lhs} and {rhs}")
        mask = lhs.glwe_sample.mask - rhs.glwe_sample.mask
        body = lhs.glwe_sample.body - rhs.glwe_sample.body
        return Cipher(GlweSample(mask=mask, body=body))

    def _compute_rank2_product(self, lhs: GlweSample,
                               rhs: GlweSample) -> Rank2Cipher:
        logger.debug(f"Computing rank 2 product of {lhs} and {rhs}")
        constant = lhs.body * rhs.body
        constant = constant % self.dist.poly_modulus

        linear = lhs.body * rhs.mask + lhs.mask * rhs.body
        linear = linear % self.dist.poly_modulus

        quadratic = lhs.mask * rhs.mask
        quadratic = quadratic % self.dist.poly_modulus

        return Rank2Cipher(constant, linear, quadratic)

    def mul(self, lhs: Cipher, rhs: Cipher):
        """
        Multiply two ciphertexts together.

        Args:
        - lhs: Cipher, the left-hand side of the multiplication
        - rhs: Cipher, the right-hand side of the multiplication

        Returns:
        - Cipher, the product of the two ciphertexts
        """

        if self.relin_key is None:
            raise ValueError("No relinearization key provided")

        # FIXME: either multiplication or relinearization needs
        # debugging.
        raise NotImplementedError(
            "Multiplication support is not yet implemented")

        logger.debug(f"Multiplying {lhs} and {rhs}")
        rank2 = self._compute_rank2_product(lhs.glwe_sample, rhs.glwe_sample)
        return rank2.relinearize(self.relin_key, self.dist.poly_modulus)
