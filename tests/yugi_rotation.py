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

    n = 8
    q = 776077649
    p1 = 17
    p2 = 3
    base_decomposition = 2  
    batched = True
    
    # q = generate_modulus(2**28,2**30,8)
    # print(q)
    
    # test function
    assert primitive_roots(17,8) == [3, 5, 6, 7, 10, 11, 12, 14]
    
    print("Roots:", primitive_roots(p1, n))
    
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
    # sk = create_sk(n, 0, 1, q)
    sk = [0] * n
    print("Secret key: ", sk)
    
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    # Generate the plaintexts
    m0 = [1,2,3,4,5,6,7,8]
    print("Plaintext: ", m0)
    
    print("Parameters: ", params)
    print("Parameters batched: ", params_batched)
    
    # Encrypt the plaintexts secret key
    ciphertext = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
    print("Ciphertext: ", ciphertext)
    
    decrypted = decrypt(sk, ciphertext, p1, params, params_batched)
    print("Decrypted (normal): ", decrypted)
    print("*" * 60)


    roots = primitive_roots(p1, n)
    matriz(roots[0],n, p1)


# ---------------------------------------------------------------------------
# paper: https://www.arxiv.org/pdf/2503.05136 - page 145
# where J(h) is the rotation helper formula: J(h) = 5^h mod 2n, J∗(h) = −5^h mod 2n

def matriz(w, n, p1):
    
    # generate matrix of size 8x8
    # and fill it with 1's and primitive root
    mat = [ [1] * n ]
    for ct in range(n-1): 
        mat.append([w] * n)
        
    # fill the first row with 1's and func_j
    mat2 = [ [1] * n ]
    for ct in range(n-1): 
        mat2.append([func_j(1,3,n), func_j(1,2,n), func_j(1,1,n), func_j(1,0,n),
                    func_j(0,3,n), func_j(0,2,n), func_j(0,1,n), func_j(0,0,n)])
        
    pot = []
    for ct in range(n): 
        pot.append([ct] * n)
        
    print(mat)
    print(mat2)
    print(pot)  
    
    # generate matrix W (page 147)
    size = len(mat)
    print("Size: ", size)
    for ct in range(size):
        for ct2 in range(size):
            t1 = mat[ct][ct2]
            t2 = mat2[ct][ct2]
            t1 = mod_number(t1**t2, p1)
            mat[ct][ct2] = t1

    for ct in range(size):
        for ct2 in range(size):
            t1 = mat[ct][ct2]
            t2 = pot[ct][ct2]
            t1 = mod_number(t1**t2, p1)
            mat[ct][ct2] = t1

    print("Matrix W:")
    for ct in range(size):
        print(mat[ct])
# ---------------------------------------------------------------------------
def func_j(t, h, n):
    if t == 1:
        return mod_number(5**h, 2*n)
    elif t == 0:
        return mod_number(-5**h, 2*n)

# ---------------------------------------------------------------------------
def prime_factors(n):
    """Retorna um conjunto com os fatores primos de n."""
    factors = set()
    # Extrai o fator 2
    while n % 2 == 0:
        factors.add(2)
        n //= 2
    # Checa fatores ímpares a partir de 3
    f = 3
    while f * f <= n:
        while n % f == 0:
            factors.add(f)
            n //= f
        f += 2
    if n > 1:
        factors.add(n)
    return factors

# ---------------------------------------------------------------------------
def is_primitive_root(g, mod, phi, pf):
    """
    Verifica se g é raiz primitiva módulo mod.
    phi é o valor de φ(mod) e pf é o conjunto dos fatores primos de φ(mod).
    g é raiz primitiva se para cada p em pf, g^(phi/p) mod mod não é igual a 1.
    """
    for p in pf:
        if pow(g, phi // p, mod) == 1:
            return False
    return True

# ---------------------------------------------------------------------------
def primitive_roots(mod, t):
    """
    Retorna uma lista com t raízes primitivas módulo mod.
    A função inicia a verificação de candidatos a partir de 3 e avança sequencialmente.
    
    Parâmetros:
      mod: inteiro representando o módulo (normalmente um número primo)
      t: quantidade de raízes primitivas a retornar
      
    Retorno:
      Uma lista contendo t números que são raízes primitivas módulo mod.
    """
    phi = mod - 1  # Para módulo primo, φ(mod) = mod - 1
    pf = prime_factors(phi)
    roots = []
    candidate = 3
    while len(roots) < t:
        if is_primitive_root(candidate, mod, phi, pf):
            roots.append(candidate)
        candidate += 1
    return roots

