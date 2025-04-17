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
    n = 8
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
        
    # q = 65537
    # p1 = 17
    
    print("P1 = ", p1)
    print("N = ", n)

    # verify if p1 is equal to 1 mod 2n (necessary condition for CRT)
    if p1 % (2*n) != 1:
        print("Error: p1 must be equal to 1 mod 2n | n = ",n, "p1 mod 2n", p1 % (2*n))
        sys.exit(1)
    
    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    #print("Q [psi_rev] = ",psi_rev)
    params = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
    #print("p1 [psi_rev] = ",psi_rev)
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
            
    # m0 = [1, 2, 3, 4, 5, 6, 7, 8]
    print("Plaintext: ", m0)
    
    # m1 = m0.copy()
    # m2 = bit_reverse_order(m1, n)
    # m3 = bit_reverse_order(m2, n)
    # print("m1 | m3", m1, m3)
    
    #print("Parameters: ", params)
    #print("Parameters batched: ", params_batched)
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
    print("C0: ", c0)
    
    decrypted = decrypt(sk, c0, p1, params, params_batched)
    print("Decrypted (c0): ", decrypted)
    print("-" * 80)
    
    assert decrypted == m0
    
    map = []
    pot = 3

    # found the automorphism
    # 1) monta o conjunto de automorfismos    
    Z_star = [d for d in range(1, 2*n, 2) if math.gcd(d, 2*n) == 1]
    print("Automorfismos válidos (Z*_{}) = {}".format(2*n, Z_star))

    root_list = Z_star
    print("root: ", root_list)
    
    separator = 60
    start = False
    for pot in root_list:
        for rot in range(n):
            # at this point we have explored all possible rotations
            
            if rot < (n/2):
                power = pot
                rotation = power**rot  
            else:
                power = -pot
                rotation = power**rot
                    
            # key rotation
            sk_rot = rotate_polynomial_coeffs(sk, rotation, n,q)
            
            print("Rotation......: ", rot)
            print("sk............: ", sk)   
            # 2) só para exibí-la de forma humana, use uma cópia!
            signed = view_sk_rot(sk_rot.copy(), q)
            print("sk_rot (signed):", signed)
            
            # generate rotation key
            ks = encrypt_sk(sk, sk_rot, q, p1, p2, batched, params, params_batched)
            
            # generate new Cryptogram
            body = c0.body
            mask = c0.mask
            
            body = rotate_polynomial_coeffs(body, rotation, n,q)
            mask = rotate_polynomial_coeffs(mask, rotation, n,q)
            
            c1 = Cryptogram(body=body, mask=mask, batched=c0.batched, q=c0.q)
            #print("C1: ", c1)
            
            # decrypted new Cryptogram with the new secret key (rotated)
            decrypted = decrypt(sk_rot, c1, p1, params, params_batched)
            # print("Decrypted (c1).: ", decrypted)
            # print("-" * separator)
                        
            # verify if the decrypted number is equal to the original plaintext
            # and if the rotation is correct
            if verify_number(decrypted, m0) == True:
                if start == True:
                    verify = 0
                    for ct in map:
                        if decrypted != ct[3]:
                            verify  += 1                     
                    
                    if verify == len(map):        
                        map.append([power, rot, rotation, decrypted])
                else:
                    map.append([power, rot, rotation, decrypted])
                    start = True
                
    # Show the results
    print("roots: ", root_list)
    for ct in map:
        print("Power = ", pad_num(ct[0],3), "| Rot = ", pad_num(ct[1],3),"| Rotation = ", pad_num(ct[2],8) ,"| Decrypted = ", ct[3])
        # print(ct[3])
        # print("-" * separator)
    print("Map size: ", len(map))
    
    # t = verify_total_rotation(map, n)
    # print(t)
    
        # --------------------------------------------------------------------
    # # Mapeamento de rotações cíclicas reais para d ∈ Z*_2N
    # print("\n" + "="*separator)
    # print("MAPEANDO ROTAÇÕES CÍCLICAS EM FUNÇÃO DE AUTOMORFISMOS d ∈ Z*_{})".format(2*n))
    # print("="*separator)

    # Z_star = [d for d in range(1, 2*n, 2) if math.gcd(d, 2*n) == 1]
    # print("Automorfismos válidos (Z*_{}) = {}".format(2*n, Z_star))

    # rotacao_para_d = {}

    # for d in Z_star:
    #     d = d % (2*n)
    #     print(f"d = {d}")
    #     sk_rot = rotate_polynomial_coeffs(sk.copy(), d, n, q)
    #     body   = rotate_polynomial_coeffs(c0.body, d, n, q)
    #     mask   = rotate_polynomial_coeffs(c0.mask, d, n, q)

    #     c1 = Cryptogram(body=body, mask=mask, batched=c0.batched, q=c0.q)
    #     decrypted = decrypt(sk_rot, c1, p1, params, params_batched)

    #     for i in range(n):
    #         esperado = m0[i:] + m0[:i]  # rotação à esquerda de i
    #         if decrypted == esperado:
    #             rotacao_para_d[i] = d
    #             print(f"✅ Rotação de {i} posição(ões) = automorfismo d = {d}")
    #             break

    # print("\n" + "-"*separator)
    # print("TABELA: ROTACAO CÍCLICA (i) → AUTOMORFISMO d")
    # print("-"*separator)
    # for i in sorted(rotacao_para_d.keys()):
    #     print(f"Rotação {i:2d} → d = {rotacao_para_d[i]}")
    # print("-"*separator)

# -----------------------------------------------------------------------
def rotate_polynomial_coeffs(coeffs, d, N, q):
    """
    Aplica o automorfismo x -> x^d mod (x^N + 1) no polinômio dado
    pelos coeficientes `coeffs`, de comprimento N, em Z_q.
    Devolve o vetor de coeficientes rotacionado com os sinais corretos.
    """
    rotated = [0] * N
    for i in range(N):
        # índice no polinômio resultante
        j = (i * d) % N
        # quantas vezes "ultrapassou" N para determinar o sinal
        k = (i * d) // N
        sign = -1 if (k % 2) else 1
        # aplica mod q (use sua função mod_number ou %q direto)
        rotated[j] = (coeffs[i] * sign) % q
    return rotated
# -----------------------------------------------------------------------
# function to rotate the coefficients
def rotate_polynomial_coeffs2(coeffs, d, N, q):
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

# ---------------------------------------------------------------------------------------
def generate_galois_keys(sk, N, q, p1, p2, batched, params, params_batched):
    """
    Gera um dicionário de chaves Galois para todos os d em Z*_{2N}.
    Cada entrada galois_keys[d] é um ciphertext que, ao ser usado
    em key‑switching, permite aplicar o automorfismo x->x^d.
    """
    # 1) monta o conjunto de automorfismos
    Z_star = [d for d in range(1, 2*N, 2) if math.gcd(d, 2*N) == 1]

    galois_keys = {}
    for d in Z_star:
        # 2) rotaciona a sk “no anel”
        sk_rot = rotate_polynomial_coeffs(sk, d, N, q)
        # 3) encripta sk_rot sob sk (aqui sk é a chave que vai
        #    conseguir desfazer a encriptação, via decrypt)
        ks = encrypt_sk(
            sk,           # chave que poderá decriptar
            sk_rot,       # “mensagem” (o vetor de coefs da sk rotacionada)
            q, p1, p2,
            batched,
            params,
            params_batched
        )
        galois_keys[d] = ks

    return galois_keys

# ---------------------------------------------------------------------------
def poly_mul(a, b, params, q):
    """
    Multiplica dois polinômios a(x) e b(x) em R_q[x]/(x^N+1)
    usando NTT. 
    `params` = (psi_rev, psi_inv_rev, n_inv, bar).
    """
    n = len(a)
    psi_rev, psi_inv_rev, n_inv, bar = params
    
    # 1) leva para domínio NTT
    ntt_generic(a, psi_rev,q, bar)
    ntt_generic(b, psi_rev,q, bar)
    
    # 2) hadamard point‑wise
    c = [mod_number(a[i] * b[i], q) for i in range(n)]
    
    # 3) volta para domínio coeficiente
    intt_generic(c, psi_inv_rev, n_inv, q, bar)
    # 4) Barrett‑reduce (ou %q direto)
    return [mod_number(ci, q) for ci in c]

# -----------------------------------------------------------------
# function to verify the total rotation
def verify_total_rotation(map, n):
    
    vet = []
    for ct in range(n):
        vet.append(ct+1)
        
    for ct2 in range(n):
        r = []
        for ct in range(n):
            tmp = map[ct][3][ct2]
            r.append(tmp)

        for ct3 in r:
            if ct3 in vet==False:
                print("Error: ", ct3)
                return False
                
    return True

# ------------------------------------------------------------------------
def pad_num(n, width):
    """
    Retorna uma string com o número `n` alinhado à direita
    em um campo de largura `width`, preenchendo com espaços.
    """
    return f"{n:>{width}d}"