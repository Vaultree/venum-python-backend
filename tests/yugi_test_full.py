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
def test_scheme_full(input):
    print("\n")
    show_title("TEST SCHEME - ENCRYPT/DECRYPT")
    dimension = [4,     8,  16,  32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]
    repeat =    [100, 100, 100, 100, 40,  40,  20,  20,    4,    4,    4,    4,     4,     4]
    
    for n in dimension:
        show_type = True
        for r in range(repeat[dimension.index(n)]):
            
            if show_type:
                print("Testing dimension: ", n)
                show_type = False
                
            q = generate_modulus(2**60, 2**63, n)
            p1 = 65537
            p2 = 3
            base_decomposition = 8
            if r%2 == 0:
                batched = False
            else:
                batched = True    
            
            # verify if p1 is equal to 1 mod 2n (necessary condition for CRT)
            if p1 % (2*n) != 1:
                print("Error: p1 must be equal to 1 mod 2n | n = ",n, "p1 mod 2n", p1 % (2*n))
                sys.exit(1)
            
            # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
            psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
            params = (psi_rev, psi_inv_rev, n_inv, bar)
            
            # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
            psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
            params_batched = (psi_rev, psi_inv_rev, n_inv, bar)
            
            # Generate the keys
            sk = create_sk(n, 0, 1, q)
            pk = generate_pk(sk, q, p1, p2, params)
            
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
            m0 = [random.randint(0, p1-1) for _ in range(n)]
            m1 = [random.randint(0, p1-1) for _ in range(n)]
            
            # Encrypt the plaintexts public key
            c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params, params_batched)
            c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params, params_batched)
            
            d0 = decrypt(sk, c0, p1, params, params_batched)
            d1 = decrypt(sk, c1, p1, params, params_batched)
            
            assert m0 == d0
            assert m1 == d1

    # **********************************************************
    show_title("TEST SCHEME - SUM/SUBTRACT")
    
    for n in dimension:
        show_type = True
        for r in range(repeat[dimension.index(n)]):
            
            if show_type:
                print("Testing dimension: ", n)
                show_type = False
                
            q = generate_modulus(2**60, 2**63, n)
            p1 = 65537
            p2 = 3
            base_decomposition = 8
            if r%2 == 0:
                batched = False
                operation = "SUM"
            elif r%2 == 1:
                batched = False  
                operation = "SUB"
            elif r%2 == 2:
                batched = True  
                operation = "SUM"
            else:
                batched = True  
                operation = "SUB"  
            
            # verify if p1 is equal to 1 mod 2n (necessary condition for CRT)
            if p1 % (2*n) != 1:
                print("Error: p1 must be equal to 1 mod 2n | n = ",n, "p1 mod 2n", p1 % (2*n))
                sys.exit(1)
            
            # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
            psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
            params = (psi_rev, psi_inv_rev, n_inv, bar)
            
            # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
            psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
            params_batched = (psi_rev, psi_inv_rev, n_inv, bar)
            
            # Generate the keys
            sk = create_sk(n, 0, 1, q)
            pk = generate_pk(sk, q, p1, p2, params)
            
            # Generate the plaintexts
            m0 = [random.randint(0, p1-1) for _ in range(n)]
            m1 = [random.randint(0, p1-1) for _ in range(n)]
            
            # Encrypt the plaintexts secret key
            c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
            c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params, params_batched)
            
            if operation == "SUM":
                sum = sum_cryptograms(c0, c1)
                dec = decrypt(sk, sum, p1, params, params_batched)
                m0m1 = add_msg(m0, m1, p1)
                assert dec == m0m1
            else:
                sub = sub_cryptograms(c0, c1)
                dec = decrypt(sk, sub, p1, params, params_batched)
                m0m1 = sub_msg(m0, m1, p1)
                assert dec == m0m1
            
            # Generate the plaintexts
            m0 = [random.randint(0, p1-1) for _ in range(n)]
            m1 = [random.randint(0, p1-1) for _ in range(n)]
            
            # Encrypt the plaintexts public key
            c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params, params_batched)
            c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params, params_batched)
            
            if operation == "SUM":
                sum = sum_cryptograms(c0, c1)
                dec = decrypt(sk, sum, p1, params, params_batched)
                m0m1 = add_msg(m0, m1, p1)
                assert dec == m0m1
            else:
                sub = sub_cryptograms(c0, c1)
                dec = decrypt(sk, sub, p1, params, params_batched)
                m0m1 = sub_msg(m0, m1, p1)
                assert dec == m0m1
    
    # **********************************************************
    show_title("TEST SCHEME - MULTIPLICATION")
    dimension = [4,     8,  16,  32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]
    repeat =    [16,   16,  16,  16,  8,   8,   4,   4,    4,    4,    2,    2,     2,     2]
    
    for n in dimension:
        show_type = True
        for r in range(repeat[dimension.index(n)]):
            
            if show_type:
                print("Testing dimension: ", n)
                show_type = False
                
            q = generate_modulus(2**62, 2**63, n)
            p1 = 65537
            p2 = 3
            base_decomposition = 8
            if r%2 == 0:
                batched = False
            else:
                batched = True    
            
            # verify if p1 is equal to 1 mod 2n (necessary condition for CRT)
            if p1 % (2*n) != 1:
                print("Error: p1 must be equal to 1 mod 2n | n = ",n, "p1 mod 2n", p1 % (2*n))
                sys.exit(1)
            
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
            
            if batched:
                mul = multiply_cryptograms(c0, c1, rlk, base_decomposition, params)
                dec = decrypt(sk, mul, p1, params, params_batched)
                m0m1 = multiply_poly_mod_batched(m0, m1, p1)
                assert dec == m0m1
            else:
                mul = multiply_cryptograms(c0, c1, rlk, base_decomposition, params)
                dec = decrypt(sk, mul, p1, params, params_batched)
                m0m1 = multiply_poly_mod(m0, m1, p1)
                assert dec == m0m1    
                        
            # Generate the plaintexts
            m0 = [random.randint(0, p1-1) for _ in range(n)]
            m1 = [random.randint(0, p1-1) for _ in range(n)]
            
            # To test the public key with more than 256 dimensions would require a larger module 
            # because the Barrett structure which is not prepared for this. Another alternative 
            # would be to use RNS (Residue Number System).
            if n > 256: # Avoiding large numbers for pk
                continue
            
            # Encrypt the plaintexts public key
            c0 = encrypt_pk(pk, m0, q, p1, p2, batched, params, params_batched)
            c1 = encrypt_pk(pk, m1, q, p1, p2, batched, params, params_batched)
            
            if batched:
                mul = multiply_cryptograms(c0, c1, rlk, base_decomposition, params)
                dec = decrypt(sk, mul, p1, params, params_batched)
                m0m1 = multiply_poly_mod_batched(m0, m1, p1)
                assert dec == m0m1
            else:
                mul = multiply_cryptograms(c0, c1, rlk, base_decomposition, params)
                dec = decrypt(sk, mul, p1, params, params_batched)
                m0m1 = multiply_poly_mod(m0, m1, p1)
                assert dec == m0m1    

# --------------------------------------------------------------
def show_title(title):
    size = 100
    print("-" * size)
    print(title)
    print("-" * size)
    
    
    
    