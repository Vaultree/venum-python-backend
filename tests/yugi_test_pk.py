import sys
from tests.yugi_test_enc_sk import vetor_aleatorio
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
        {
            "params": EncryptionParameters(
                dimension=4,
                ciphertext_modulus=12289,
                plaintext_modulus=127,
                noise_modulus=3,
            ),
            "message": [0, 0, 0, 4]
        },
    ])
def test_encrypt_decrypt2(input):
    params, message = input["params"], input["message"]

    params.dimension = 8
    params.ciphertext_modulus=1400472361734830353
    # params.ciphertext_modulus = gerar_primo(2**60, 2**61, params.dimension)
    params.plaintext_modulus = 65537

    total = 1000
    quant_erros = 0
    for ct in range(total):
        message = vetor_aleatorio(params.dimension, 1, params.plaintext_modulus-1)
        print("MESSAGE: ", message)
    
        dist = GlweDistribution(params)
        sk, pk = gen_key_pair(dist)
        encryptor = Encryptor(dist, PolynomialEncoder(dist))
        cipher = encryptor.encrypt(pk, message)
        decrypted = encryptor.decrypt(sk, cipher)
        
        print("Decifragem: ", decrypted)
        print("-" * 50)
        
        assert decrypted == message
        if decrypted != message:
            quant_erros += 1        
    
    print("=" * 80)
    print("Quantidade de erros: ", quant_erros)
    sys.exit(1)
