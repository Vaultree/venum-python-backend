import sys
import time
from tests.yugi_test_mul import generate_formula
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
def test_scheme_mul(input):
    print("\n", "-" * 80)
    n = 4  # Dimension of the polynomial
    q = generate_modulus(2**60, 2**63, n)
    p1 = 65537
    p2 = 3
    base_decomposition = 2
    batched = False
    
    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    params = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # Generate the keys
    sk = create_sk(n, 0, 1, q)
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    # Generate the plaintexts
    m0 = [random.randint(0, p1-1) for _ in range(n)]
    m1 = [random.randint(0, p1-1) for _ in range(n)]
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params)
    c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params)
    
    d0 = decrypt(sk, c0, p1, params)
    d1 = decrypt(sk, c1, p1, params)
    
    assert m0 == d0
    assert m1 == d1
    
    # Generate the plaintexts
    m0 = [random.randint(0, 4) for _ in range(n)]
    m1 = [random.randint(0, 4) for _ in range(n)]
    
    # Encrypt the plaintexts public key
    c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params)
    c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params)
    
    d0 = decrypt(sk, c0, p1, params)
    d1 = decrypt(sk, c1, p1, params)
    
    assert m0 == d0
    assert m1 == d1
    
    for i in range(1000):
        # Generate the plaintexts
        m0 = [random.randint(0, p1-1) for _ in range(n)]
        m1 = [random.randint(0, p1-1) for _ in range(n)]
        
        # Encrypt the plaintexts public key
        # c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params)
        # c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params)
        c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params)
        c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params)
        
        # calculate the product
        prod = product(c0, c1, params)    
        c3 = relinearize(prod, rlk, c0.batched, base_decomposition, c0.q, params)
        
        dec3 = decrypt(sk, c3, p1, params)
        print("dec3: ", dec3)
        
        dec = multiply_poly_mod(m0, m1, p1)
        generate_formula(m0, m1, p1)
        print("dec: ", dec)
        
        assert dec == dec3
    
    
    
    
    
    