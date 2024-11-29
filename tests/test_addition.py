from venum.glwe import EncryptionParameters, GlweDistribution
from venum.encryption import Encryptor
from venum.plaintext_encoding import PolynomialEncoder
from venum.key import gen_key_pair
from venum.evaluation import Evaluator

import pytest


@pytest.mark.parametrize(
    "input",
    [
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=383,
                plaintext_modulus=127,
                noise_modulus=3,
                seed=0
            ),
            "lhs": [1, 2, 3, 4],
            "rhs": [5, 6, 7, 8],
        },
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=12289,
                plaintext_modulus=127,
                noise_modulus=3,
                seed=1
            ),
            "lhs": [1, 2, 3, 4],
            "rhs": [5, 6, 7, 8],
        },
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=383,
                plaintext_modulus=127,
                noise_modulus=3,
            ),
            "lhs": [1, 2, 3, 4],
            "rhs": [5, 6, 7, 8],
        },
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=12289,
                plaintext_modulus=127,
                noise_modulus=3,
            ),
            "lhs": [1, 2, 3, 4],
            "rhs": [5, 6, 7, 8],
        },
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
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=1400472361734830353,
                plaintext_modulus=12289,
                noise_modulus=3,
            ),
            "lhs": [10001, 10002, 10003, 10004],
            "rhs": [4, 3, 2, 1],
        },
    ])
def test_addition(input):
    params, lhs, rhs = input["params"], input["lhs"], input["rhs"]
    expected = [x + y for x, y in zip(lhs, rhs)]

    dist = GlweDistribution(params)
    sk, pk = gen_key_pair(dist)
    encryptor = Encryptor(dist, PolynomialEncoder(dist))
    lhs_cipher = encryptor.encrypt(pk, lhs)
    rhs_cipher = encryptor.encrypt(pk, rhs)
    eval = Evaluator(dist)
    cipher_result = eval.add(lhs_cipher, rhs_cipher)
    decrypted = encryptor.decrypt(sk, cipher_result)
    assert decrypted == expected
