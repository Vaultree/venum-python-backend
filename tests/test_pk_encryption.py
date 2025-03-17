from venum.glwe import EncryptionParameters, GlweDistribution
from venum.encryption import Encryptor
from venum.plaintext_encoding import PolynomialEncoder
from venum.key import gen_key_pair

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
            "message": [1, 2, 3, 4]
        },
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=12289,
                plaintext_modulus=127,
                noise_modulus=3,
                seed=1
            ),
            "message": [5, 6, 7, 8]
        },
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=383,
                plaintext_modulus=127,
                noise_modulus=3,
            ),
            "message": [5, 6, 7, 8]
        },
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=12289,
                plaintext_modulus=127,
                noise_modulus=3,
            ),
            "message": [1, 2, 3, 4]
        },
    ])
def test_encrypt_decrypt(input):
    params, message = input["params"], input["message"]
    dist = GlweDistribution(params)
    sk, pk = gen_key_pair(dist)
    encryptor = Encryptor(dist, PolynomialEncoder(dist))
    cipher = encryptor.encrypt(pk, message)
    decrypted = encryptor.decrypt(sk, cipher)
    assert decrypted == message
