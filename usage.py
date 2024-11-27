from venum.glwe import EncryptionParameters, GlweDistribution
from venum.encryption import Encryptor
from venum.plaintext_encoding import PolynomialEncoder
from venum.key import gen_key_pair, SecretKey


def test_sympy():
    # params = EncryptionParameters(
    #     dimension=4,
    #     ciphertext_modulus=12289,
    #     plaintext_modulus=127,
    #     noise_modulus=3,
    #     seed=1
    # )
    params = EncryptionParameters(
        dimension=2,
        ciphertext_modulus=23,
        plaintext_modulus=7,
        noise_modulus=3,
        seed=1
    )
    dist = GlweDistribution(params)
    sk, pk = gen_key_pair(dist)
    encryptor = Encryptor(dist, PolynomialEncoder)
    cipher = encryptor.encrypt(pk, [1, 1])
    message = encryptor.decrypt(sk, cipher)
    print(message)


test_sympy()
