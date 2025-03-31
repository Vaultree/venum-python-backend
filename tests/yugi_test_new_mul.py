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
    
    # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
    params_batched = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # Generate the keys
    sk = create_sk(n, 0, 1, q)
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    # Generate the plaintexts
    m0 = [random.randint(0, p1-1) for _ in range(n)]
    m1 = [random.randint(0, p1-1) for _ in range(n)]
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
    c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params, params_batched)
    
    d0 = decrypt(sk, c0, p1, params, params_batched)
    d1 = decrypt(sk, c1, p1, params, params_batched)
    
    assert m0 == d0
    assert m1 == d1
    
    # Generate the plaintexts
    m0 = [random.randint(0, 4) for _ in range(n)]
    m1 = [random.randint(0, 4) for _ in range(n)]
    
    # Encrypt the plaintexts public key
    c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params, params_batched)
    c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params, params_batched)
    
    d0 = decrypt(sk, c0, p1, params, params_batched)
    d1 = decrypt(sk, c1, p1, params, params_batched)
    
    assert m0 == d0
    assert m1 == d1
    
    print("\n", "-" * 80)
    print(" TESTE DE STRESS - MULTIPLICAÇÃO 4096")
    print("", "-" * 80)
    n = 4096  # Dimension of the polynomial
    q = generate_modulus(2**60, 2**63, n)
    p1 = 65537
    p2 = 3
    base_decomposition = 2
    batched = False
    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    params = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
    params_batched = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # Generate the keys
    sk = create_sk(n, 0, 1, q)
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)

    tempo = 0
    total_tests = 2
    for i in range(total_tests):
        # Generate the plaintexts
        m0 = [random.randint(0, p1-1) for _ in range(n)]
        m1 = [random.randint(0, p1-1) for _ in range(n)]
        
        # Encrypt the plaintexts public key
        c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
        c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params, params_batched)
        # c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params)
        # c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params)
        
        # calculate the product
        start_time = time.perf_counter()
        prod = product(c0, c1, params)    
        c3 = relinearize(prod, rlk, c0.batched, base_decomposition, c0.q, params)
        end_time = time.perf_counter()
        tempo += end_time - start_time
        
        dec3 = decrypt(sk, c3, p1, params, params_batched)
        #print("dec3: ", dec3)
        
        dec = multiply_poly_mod(m0, m1, p1)
        #generate_formula(m0, m1, p1)
        #print("dec: ", dec)
        
        assert dec == dec3
    
    print(f"Average time: {tempo/total_tests:.6f} seconds")
    
# --------------------------------------------------------------

    n = 16  # Dimension of the polynomial
    q = generate_modulus(2**60, 2**63, n)
    p1 = 65537
    p2 = 3
    base_decomposition = 2
    batched = True
    
    print("\n", "-" * 80)
    print(" TESTANDO O MODO BATCHED - ",n, " DIMENSÕES")
    print("", "-" * 80)

    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    params = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
    params_batched = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # Generate the keys
    sk = create_sk(n, 0, 1, q)
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)

    total_tests = 1000
    for round in range(total_tests):
        # Generate the plaintexts
        m0 = [random.randint(0, p1-1) for _ in range(n)]
        m1 = [random.randint(0, p1-1) for _ in range(n)]
        
        # Encrypt the plaintexts public key
        if round % 2 == 0:
            c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
            c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params, params_batched)
        else:    
            c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params, params_batched)
            c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params, params_batched)
        
        # calculate the product
        start_time = time.perf_counter()
        prod = product(c0, c1, params)    
        c3 = relinearize(prod, rlk, c0.batched, base_decomposition, c0.q, params)
        end_time = time.perf_counter()
        tempo += end_time - start_time
        
        dec3 = decrypt(sk, c3, p1, params, params_batched)
        #print("dec3: ", dec3)
        
        dec = multiply_poly_mod_batched(m0, m1, p1)
        
        assert dec == dec3
    
    print(f"Average time: {tempo/total_tests:.6f} seconds")
    

    
    
    
    
    