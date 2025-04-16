import sys
import time
import random
import pytest
from math import ceil, log

from tests.yugi_test_full import show_title
from tests.yugi_test_mul import generate_formula
from venum.ntt import *
from venum.ntt_func import *

@pytest.mark.parametrize("input", [{}])
def test_scheme_rotation(input):
    print("\n")
    show_title("TEST SCHEME - ROTATION")

    # parameters:
    n = 16
    # q = 776077649
    # p1 = 17
    p2 = 3
    base_decomposition = 2  
    batched = True
    
    q = generate_modulus(2**60,2**63,n)
    print(" Q = ", q)
    
    while True:
        p1 = generate_modulus(2**7,2**20,n)
        if p1 % 4 == 1:
            break    
    
    print("P1 = ", p1)

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
    print("Secret key: ", sk)
    
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    # Generate the plaintexts
    m0 = []
    for ct in range(n):
            m0.append(ct+1)
    print("Plaintext: ", m0)
    
    #print("Parameters: ", params)
    #print("Parameters batched: ", params_batched)
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
    print("C0: ", c0)
    
    decrypted = decrypt(sk, c0, p1, params, params_batched)
    print("Decrypted (c0): ", decrypted)
    print("-" * 80)
    
    map = []
    pot = 3
    
    primes_list = primes(n)
    primes_list = [5]
    print("Primes: ", primes)

    separator = 60
    start = False
    for pot in primes_list:
        for rot in range(n):
            # at this point we have explored all possible rotations
            if rot < (n/2):
                rotation = pot**rot  
            else:
                rotation = -pot**rot
                    
            # key rotation
            sk_rot = rotate_polynomial_coeffs(sk, rotation, n,q)
            
            print("Rotation......: ", rot)
            print("sk............: ", sk)   
            print("sk_rot........: ", view_sk_rot(sk_rot,q))   
            # ks = encrypt_sk(sk_rot, s, q, p1, p2, batched, params, params_batched)
            
            # generate new Cryptogram
            body = c0.body
            mask = c0.mask
            body = rotate_polynomial_coeffs(body, rotation, n,q)
            mask = rotate_polynomial_coeffs(mask, rotation, n,q)
            
            c1 = Cryptogram(body=body, mask=mask, batched=c0.batched, q=c0.q)
            #print("C1: ", c1)
            
            # decrypted new Cryptogram with the new secret key (rotated)
            decrypted = decrypt(sk_rot, c1, p1, params, params_batched)
            print("Decrypted.....: ", decrypted)
            print("-" * separator)
            
            # verify if the decrypted number is equal to the original plaintext
            # and if the rotation is correct
            if verify_number(decrypted, m0) == True:
                if start == True:
                    verify = 0
                    for ct in map:
                        if decrypted != ct[2]:
                            verify  += 1                     
                    
                    if verify == len(map):        
                        map.append([pot, rot, decrypted])
                else:
                    map.append([pot, rot, decrypted])
                    start = True
                
    # Show the results
    for ct in map:
        print("Power =", ct[0], "| Rotation =", ct[1], "| Decrypted = ", ct[2])
        print("-" * separator)
    print("Map size: ", len(map))

# -----------------------------------------------------------------------
# function to rotate the coefficients
def rotate_polynomial_coeffs(coeffs, d, N, q):
    tmp = []
    ret = []
    sig = []
    for i in range(N):
        # print("i =", i, d, i*d, "[",(i*d)%N,"]")
        t1 = (i*d) // N
        #print("t1 =", t1)
        signal = (-1)**(t1%2)
        sig.append(signal)
        tmp.append(((i * d) % N))
    
    #print("tmp = ", tmp)
    #print("sig = ", sig)
    for i in tmp:
        #print("I = ", coeffs[i]* sig[i])
        ret.append(mod_number(coeffs[i]* sig[i], q))
    
    return ret

# -----------------------------------------------------------------------
# show the sk_rot
def view_sk_rot(sk,q):
    n = len(sk)
    for ct in range(n):
        if sk[ct] >= (q/2):
            sk[ct] = sk[ct] - q

    return sk

# -----------------------------------------------------------------------
# verify if the number of rotations is correct
def verify_number(vet, m0):
    ret = True
    for ct in m0:
        number = ct
        if number in vet:
            ret = True
        else:
            ret = False
            break
        
    return ret

# ---------------------------------------------------------------------------
# generate prime numbers
def primes(limite):
    primos = []
    for num in range(3, limite + 1):
        eh_primo = True
        for i in range(2, int(num**0.5) + 1):
            if num % i == 0:
                eh_primo = False
                break
        if eh_primo:
            primos.append(num)
    return primos

