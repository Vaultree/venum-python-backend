import sys
import time
from venum.ntt import *
import pytest
import random
from venum.ntt_func import *

@pytest.mark.parametrize(
    "input",
    [
        {
        }
    ])
def test_scheme(input):
    # print("\n", "-" * 80)
    # n = 8  # Dimension of the polynomial
    # q = generate_modulus(2**60, 2**63, n)
    
    # # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    # psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    
    # # Example of n-dimensional polynomials
    # # Polynomial A(x) = 1 + 2x + 3x^2 + 4x^3 ...
    # # Polynomial B(x) = 8 + 7x + 6x^2 + 5x^3 ...
    # a = generate_random_vector(n, 0, q - 1)
    # b = generate_random_vector(n, 0, q - 1) 

    # print("Polynomial A:", a)
    # print("Polynomial B:", b)
    # # print("q...........:", q)
    
    # # Multiplication of polynomials using NTT/INTT
    # result = polymul_ntt(a, b, q, psi_rev, psi_inv_rev, n_inv, bar)
    # print("Result of multiplication: ", result)
    
    # # check if the result is correct
    # assert multiply_poly_mod(a, b, q) == result
    
    # --------------------------------------------------------------
    
    # Main parameters of scheme
    print("\n", "-" * 80)
    n = 8  # Dimension of the polynomial
    q = generate_modulus(2**62, 2**63, n)
    p1 = 65537  # space cleartext
    p2 = 3      # space noise
    
    print("q...........:", q)
    print("p1..........:", p1)
    print("p2..........:", p2)
    print("n...........:", n)
    
    # Generate a secret key
    sk = create_sk(n, 0, 1, q)
    print("Secret key:", sk)
    
    # Generate the required parameters for NTT/INTT (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    params = generate_parameters(n, q)

    # Generate the public key
    pk = generate_pk(sk, q, p1, p2, params)
    print("Public key:", pk)
    
    print("-" * 80)    
    
    # testin the encryption
    print("TESTING ENCRYPTION - SECRET KEY")
    message = generate_random_vector(n, 0, p1-1)
    print("Message:", message)
    
    # Encrypt the message
    batched = True
    c1 = encrypt_sk(sk, message, q, p1, p2, batched, params)
    print("Ciphertext:", c1)
    
    # Decrypt the message
    d1 = decrypt(sk, c1, p1, params)
    print("Decrypted message 1:", d1)
    assert message == d1

    # Encrypt the message
    print("-" * 80)    
    batched = False
    c2 = encrypt_sk(sk, message, q, p1, p2, batched, params)
    print("Ciphertext:", c2)
    
    # Decrypt the message
    d2 = decrypt(sk, c2, p1, params)
    print("Decrypted message 2:", d2)
    assert message == d2
    
    print("-" * 80)    
    
    print("TESTING ENCRYPTION - PUBLIC KEY")
    params = generate_parameters(n, q)
    # testin the encryption
    message = generate_random_vector(n, 0, p1-1)
    print("Message:", message)
    
    # Encrypt the message
    batched = False
    c3 = encrypt_pk(pk, message, q, p1, p2, batched, params)
    print("Ciphertext:", c3)
    
    # Decrypt the message
    d3 = decrypt(sk, c3, p1, params)
    print("Decrypted message 3:", d3)
    assert message == d3

    # Encrypt the message
    print("-" * 80)   
    message = generate_random_vector(n, 0, p1)
    print("Message:", message)
 
    batched = True
    c4 = encrypt_pk(pk, message, q, p1, p2, batched, params)
    print("Ciphertext:", c4)
    
    # Decrypt the message
    d4 = decrypt(sk, c4, p1, params)
    print("Decrypted message 4:", d4)
    assert message == d4
    
# ----------------------------------------------------------------------------------------

    print("-" * 80)
    print("TESTING ADDITION")
    # Generate two messages
    message1 = generate_random_vector(n, 0, p1-1)
    message2 = generate_random_vector(n, 0, p1-1)
    print("Message 1:", message1)
    print("Message 2:", message2)
    
    # Encrypt the messages
    c1 = encrypt_pk(pk, message1, q, p1, p2, False, params)
    c2 = encrypt_pk(pk, message2, q, p1, p2, False, params)
    print("Ciphertext 1:", c1)
    print("Ciphertext 2:", c2)
    
    # Add the ciphertexts
    c3 = sum_cryptograms(c1, c2)
    print("Ciphertext 3:", c3)
    
    # Decrypt the message
    d3 = decrypt(sk, c3, p1, params)
    print("Decrypted sum:", d3)
    assert add_msg(message1, message2, p1) == d3
    
# ----------------------------------------------------------------------------------------

    print("-" * 80)
    print("TESTING SUBTRACTION")
    # Generate two messages
    message1 = generate_random_vector(n, 0, p1-1)
    message2 = generate_random_vector(n, 0, p1-1)
    print("Message 1:", message1)
    print("Message 2:", message2)
    
    # Encrypt the messages
    c1 = encrypt_sk(sk, message1, q, p1, p2, False, params)
    c2 = encrypt_sk(sk, message2, q, p1, p2, False, params)
    print("Ciphertext 1:", c1)
    print("Ciphertext 2:", c2)
    
    # Add the ciphertexts
    c3 = sub_cryptograms(c1, c2)
    print("Ciphertext 3:", c3)
    
    # Decrypt the message
    d3 = decrypt(sk, c3, p1, params)
    print("Decrypted sum:", d3)
    assert sub_msg(message1, message2, p1) == d3


