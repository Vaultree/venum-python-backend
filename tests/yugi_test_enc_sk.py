from venum.glwe import EncryptionParameters, GlweDistribution
from venum.encryption import Encryptor
from venum.plaintext_encoding import PolynomialEncoder
from venum.key import gen_key_pair

import pytest
import random

def vetor_aleatorio(n, min_val, max_val):
    """
    Gera um vetor aleatório com n elementos.
    
    Parâmetros:
    - n (int): número de elementos do vetor.
    - min_val (float): valor mínimo do intervalo.
    - max_val (float): valor máximo do intervalo.
    
    Retorna:
    - list: vetor com n números aleatórios entre min_val e max_val.
    """
    return [random.randint(min_val, max_val) for _ in range(n)]

@pytest.mark.parametrize(
    "input",
    [
        {
            "params": EncryptionParameters(
                 dimension=4,
                ciphertext_modulus=12289,
                plaintext_modulus=127,
                noise_modulus=3,
                seed=0
            ),
            "message": [1, 2, 3, 4]
        },
    ])
def test_encrypt_decrypt(input):
    
    params, message = input["params"], input["message"]
    
    for i in range(1000):
        message = vetor_aleatorio(4, 0, 126)
        print(message)
        
        dist = GlweDistribution(params)
        sk, _ = gen_key_pair(dist)
        encryptor = Encryptor(dist, PolynomialEncoder(dist))
        cipher = encryptor.encrypt(sk, message)
        decrypted = encryptor.decrypt(sk, cipher)
        assert decrypted == message
