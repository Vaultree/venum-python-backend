from venum.glwe import EncryptionParameters, GlweDistribution
from venum.encryption import Encryptor
from venum.plaintext_encoding import PolynomialEncoder
from venum.key import gen_key_pair, RelinKey
from venum.evaluation import Evaluator

from sympy import Poly
from sympy.abc import x
import pytest


@pytest.mark.skip(reason="multiplication needs fixing")
@pytest.mark.parametrize(
    "input",
    [
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=1400472361734830353,
                plaintext_modulus=12289,
                noise_modulus=3,
            ),
            "lhs": [0, 0, 0, 0],
            "rhs": [0, 0, 0, 0],
        },
    ])
def test_multiplication(input):
    params, lhs, rhs = input["params"], input["lhs"], input["rhs"]

    dist = GlweDistribution(params)

    expected = (Poly(reversed(lhs), x, domain=dist.plaintext_ring) *
                Poly(reversed(rhs), x, domain=dist.plaintext_ring) %
                dist.poly_modulus.set_domain(dist.plaintext_ring))
    expected = list(reversed(expected.all_coeffs()))

    sk, pk = gen_key_pair(dist)
    encryptor = Encryptor(dist, PolynomialEncoder(dist))
    lhs_cipher = encryptor.encrypt(pk, lhs)
    rhs_cipher = encryptor.encrypt(pk, rhs)
    relin_key = RelinKey.from_secret_key(sk)
    eval = Evaluator(dist, relin_key)
    cipher_result = eval.mul(lhs_cipher, rhs_cipher)
    decrypted = encryptor.decrypt(sk, cipher_result)
    assert decrypted == expected
