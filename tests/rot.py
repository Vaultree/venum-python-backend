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
        p1 = generate_modulus(2**7,2**10,n)
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
        m0.append((ct+1)* 1)
            
    # m0 = [2,3,1,4,7,6,8,5]
    # m0 = [1,1,1,1,1,1,1,1]
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

    print("C0: ", c0)
    
    assert decrypted == m0
    
    map = []
    pot = 3

    # found the automorphism
    # 1) monta o conjunto de automorfismos    
    Z_star = [d for d in range(1, 2*n, 2) if math.gcd(d, 2*n) == 1]
    print("Automorfismos válidos (Z*_{}) = {}".format(2*n, Z_star))

    root_list = Z_star
    # root_list = [5]
    # root_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
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
            print("P1 = ", p1, "P2 = ", p2)
            
            # generate rotation key
            ksk = generate_ks_keys(sk_rot, sk, base_decomposition, q, p1, p2, batched, params, params_batched)
            # print("KSK: ", ksk)
            # print(">>> len(ksk) esperado:", ceil_log(q, base_decomposition))
            # print(">>> len(ksk) obtido   :", len(ksk))
            # for idx, ck in enumerate(ksk):
            #     print(f"  ksk[{idx}].mask[:5] = {ck.mask[:5]}")
            #     print(f"  ksk[{idx}].body[:5] = {ck.body[:5]}")
            
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

            c2 = key_switching(ksk, base_decomposition, c1, params)
            decrypted_2 = decrypt(sk, c2, p1, params, params_batched)
            
            print("Decrypted (c1).: ", decrypted)
            print("Decrypted (c2).: ", decrypted_2)
            
            # ------------------------------------------------------
            # testando a mudança para a mesma chave:
            ksk = generate_ks_keys(sk, sk, base_decomposition, q, p1, p2, batched, params, params_batched)
            c3 = key_switching(ksk, base_decomposition, c0, params)
            decrypted_3 = decrypt(sk, c3, p1, params, params_batched)
            decrypted_0 = decrypt(sk, c0, p1, params, params_batched)
            print("Decrypted (c3).: ", decrypted_3)
            print("Decrypted (c0).: ", decrypted_0)
            print("-" * separator)
            
            # sys.exit(0)
                        
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
    
    t = 17
    raiz = 3
    pot = 5
    
    print("*" * separator)
    w = generate_w(raiz, pot, 8, t)
    
    print("*" * separator)
    w1 = generate_w1(raiz, pot, 8, t)
    
    print("*" * separator)
    d = multiply_matrices(w1, w, t)
    print("Resultado: ")
    
    mat_i = d.copy()
    
    for ct in d:
        print(ct)

    # n−1Wˆ · IRn
    modinv = modular_inverse(8//2,t)
    # modinv = 13
    
    m1 = mult_number_matrix(modinv, w, t)
    print("*" * separator)
    print("modinv: ", modinv)
    print("Resultado: ")
        
    m2 = multiply_matrices(m1, d.copy(), t)
        
    m3 = mul_matrix_vector(m2, [1,2,3,4,5,6,7,8], t)

    print("Resultado [m3]: ", m3)
    print("*" * separator)

    m4 = mul_matrix_vector(w1, m3, t)
    # m4 = mul_matrix_vector(d, m4, t)
    print("Resultado: ", m4)

    # m4 = mul_matrix_vector(w1, [13, 12, -10, -2, 18, -32, -14, 10], t)
    # print("Resultado: ", m4)
    
    # t = verify_total_rotation(map, n)
    # print(t)
    
    # m2 = multiply_matrices(w, mat_i, t)
    # print("Resultado: ")
    # for ct in m2:
    #     print(ct)
    
    
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
# generate de keys of base_decomposition
def generate_ks_keys(key_old, key_new, base_decomposition, q, p1, p2, batched, params, params_batched):

    # Validar entradas (opcional, mas recomendado)
    if not isinstance(q, (int, float)) or q <= 0:
        raise ValueError("Q deve ser um número positivo.")
    if not isinstance(base_decomposition, int) or base_decomposition < 2:
        raise ValueError("Base deve ser um inteiro maior ou igual a 2.")

    # j = math.log(q, base_decomposition)
    # j = int(math.ceil(j))
    
    j = 0
    while base_decomposition**j < q:
        j += 1
    
    keys = []
    n = len(key_old)
    if n != len(key_new):
        print("Error: The keys must have the same size")
        sys.exit(1)
        
    base_key = decompose_poly_list(key_old, base_decomposition, q)
    
    for i in range(j):
        tmp = []
        for ct in range(n):
            t1 = mod_number(key_old[ct] * base_decomposition**i, q)
            tmp.append(t1)
        
        # encrypt the old key
        # print(tmp)
        c0 = encrypt_sk(key_new, tmp, q, p1, p2, batched, params, params_batched)    
        keys.append(c0)
    
    return keys

# -----------------------------------------------------------------------
def key_switching(ksk: list, base_decomposition: int, cryptogram, params):
    """
    Aplica o Key Switching a um criptograma.

    Args:
        ksk: Lista das Key Switching Keys [RLWE(key_new, key_old*B^0), RLWE(key_new, key_old*B^1), ...].
        base_decomposition: A base B usada.
        cryptogram: O criptograma (a, b) = RLWE(key_old, m) a ser convertido.
        params: Parâmetros NTT (psi_rev, psi_inv_rev, n_inv, bar).

    Returns:
        Um novo Cryptogram (a', b') = RLWE(key_new, m).
    """
    q = cryptogram.q
    n = len(cryptogram.body)
    original_mask = cryptogram.mask # = a
    original_body = cryptogram.body # = b

    # 1. Decompor a máscara original 'a' na base B
    # a1 = [a_0, a_1, ..., a_K] onde a = Σ a_j * B^j
    a1 = decompose_poly_list(original_mask, base_decomposition, q)
    
    # print(">>> num_components esperado:", ceil_log(q, base_decomposition))
    # print(">>> len(a1) obtido       :", len(a1))
    # print(">>> cada a1[j] tem tamanho:", [len(poly) for poly in a1])

    psi_rev, psi_inv_rev, n_inv, bar = params
    num_components = len(a1) # Número de polinômios decompostos (K+1)

    # Verificação de consistência
    if len(ksk) != num_components:
         raise ValueError(f"Erro: Número de componentes decompostos ({num_components}) "
                          f"difere do número de KSKs fornecidas ({len(ksk)})")

    # --- 2. Calcular a nova máscara a' = - Σ[j=0 to K] a_j * KSa_j ---
    # Inicializa a soma Σ a_j * KSa_j
    new_mask_sum = [0] * n
    for j in range(num_components):
        a_j = a1[j]             # j-ésimo polinômio da decomposição de 'a'
        KSa_j = ksk[j].mask     # Componente 'mask' da j-ésima KSK

        # Calcular o produto a_j * KSa_j usando NTT
        product = polymul_ntt(a_j, KSa_j, q, psi_rev, psi_inv_rev, n_inv, bar)

        # Acumular na soma
        for ct in range(n):
            new_mask_sum[ct] = mod_number(new_mask_sum[ct] + product[ct], q)

    # A nova máscara a' é a negação da soma
    new_mask = [mod_number(-new_mask_sum[ct], q) for ct in range(n)]

    # --- 3. Calcular o novo corpo b' = b - Σ[j=0 to K] a_j * KSb_j ---
    # Inicializa a soma Σ a_j * KSb_j
    new_body_sum = [0] * n
    for j in range(num_components):
        a_j = a1[j]             # j-ésimo polinômio da decomposição de 'a'
        KSb_j = ksk[j].body     # Componente 'body' da j-ésima KSK

        # Calcular o produto a_j * KSb_j usando NTT
        product = polymul_ntt(a_j, KSb_j, q, psi_rev, psi_inv_rev, n_inv, bar)

        # Acumular na soma
        for ct in range(n):
            new_body_sum[ct] = mod_number(new_body_sum[ct] + product[ct], q)

    # O novo corpo b' é b - (soma)
    new_body = [mod_number(original_body[ct] - new_body_sum[ct], q) for ct in range(n)]

    # Retornar o novo criptograma (a', b')
    return Cryptogram(body=new_body, mask=new_mask, batched=cryptogram.batched, q=q)

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
        # sign = 1
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

# ------------------------------------------------------------------------
# https://www.arxiv.org/pdf/2503.05136 - page 149
def generate_w(root, j, n, t):

    mat = []
    for ct in range(n):
        mat.append([0] * n)
        
    mat[0] = [1] * n
    
    for ct in range(1,n):
        mat[ct] = [root] * n

    for ct in range(1,n):
        j3 = j**3 % (2*n)
        j2 = j**2 % (2*n)
        j1 = j**1 % (2*n)
        j0 = j**0 % (2*n)
        mat[ct][0] = (mat[ct][0]**j3)**ct % t   
        mat[ct][1] = (mat[ct][1]**j2)**ct % t   
        mat[ct][2] = (mat[ct][2]**j1)**ct % t   
        mat[ct][3] = (mat[ct][3]**j0)**ct % t   

    for ct in range(1,n):
        j3 = -j**3 % (2*n)
        j2 = -j**2 % (2*n)
        j1 = -j**1 % (2*n)
        j0 = -j**0 % (2*n)        
        mat[ct][4] = (mat[ct][4]**j3)**ct % t   
        mat[ct][5] = (mat[ct][5]**j2)**ct % t   
        mat[ct][6] = (mat[ct][6]**j1)**ct % t   
        mat[ct][7] = (mat[ct][7]**j0)**ct % t   
    
    for ct in mat:        
        print(ct)

    return mat
    
# ------------------------------------------------------------------------
# https://www.arxiv.org/pdf/2503.05136 - page 149 Matrix W*
def generate_w1(root, j, n, t):

    mat = []
    for ct in range(n):
        mat.append([1, root, root, root, root, root, root, root])
        
    for number in range(n//2):
        for ct in range(n):
            jt = j**number % (2*n)
            mat[number][ct] = (mat[number][ct]**jt)**ct % t

    for number in range(4,n):
        for ct in range(n):
            jt = -j**(number-4) % (2*n)
            mat[number][ct] = (mat[number][ct]**jt)**ct % t
    
    for ct in mat:        
        print(ct)

    return mat

# ------------------------------------------------------------------------
def mul_matrix_vector(matrix, vector, t):
    
    n = len(vector)
    row = len(matrix)
    
    if n != row:
        print("Error: The matrix and vector must have the same size")
        sys.exit(1)
    
    result = []
    for ct in range(row):
        tmp = 0
        for ct2 in range(n):
            tmp = mod_number(tmp + (matrix[ct][ct2] * vector[ct2]), t)

    # for ct in range(row):
    #     tmp = 0
    #     for ct2 in range(n):
    #         tmp = mod_number(tmp + (matrix[ct2][ct] * vector[ct]), t)
            
        result.append(tmp)
        
    return result
    
# ------------------------------------------------------------------------
def multiply_matrices(A, B, t):
    """
    Multiplica duas matrizes A e B.
    A deve ter dimensões m×p e B deve ter dimensões p×q.
    Retorna a matriz resultante C de dimensões m×q.
    """
    # Verifica se os tamanhos são compatíveis
    m = len(A)
    p = len(A[0]) if m > 0 else 0
    p2 = len(B)
    q = len(B[0]) if p2 > 0 else 0
    
    if p != p2:
        raise ValueError(f"Incompatíveis para multiplicação: A é {m}×{p}, B é {p2}×{q}")
    
    # Inicializa C com zeros
    C = [[0] * q for _ in range(m)]
    
    # Cálculo da multiplicação
    for i in range(m):
        for j in range(q):
            total = 0
            for k in range(p):
                total += ((A[i][k] * B[k][j]) % t)
            C[i][j] = total % t
    
    return C

# -----------------------------------------------------------------------
def mult_number_matrix(number, mat, t):
    
    n = len(mat)
    
    for ct in range(n):
        for ct2 in range(n):
            mat[ct][ct2] = mod_number(mat[ct][ct2] * number, t) 
            
    return mat

# -------------------------------------------------------------------------
def ceil_log(q: int, base: int) -> int:
    """Retorna o menor j tal que base**j >= q."""
    j = 0
    while base**j < q:
        j += 1
    return j
