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
        p1 = generate_modulus(2**7,2**12,n)
        if p1 % 4 == 1:
            break    
        
    # q = 65537
    # p1 = 17
    
    print("P1 = ", p1)
    print("N = ", n)

    # Gerei Q, P1, N etc.
    print(" Q = ", q)
    print("P1 = ", p1)
    print("N = ", n)

    # ... checks p1 % (2*n) etc ...

    # Generate parameters for q and p1 (batched)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    params = (psi_rev, psi_inv_rev, n_inv, bar)

    psi_rev_p1, psi_inv_rev_p1, n_inv_p1, bar_p1 = generate_parameters(n, p1)
    params_batched = (psi_rev_p1, psi_inv_rev_p1, n_inv_p1, bar_p1) # Corrected naming for clarity

    # Generate the keys
    sk = create_sk(n, 0, 1, q)
    print("Secret key: ", view_sk_rot(sk.copy(),q)) # View signed SK

    pk = generate_pk(sk, q, p1, p2, params)
    # rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params) # RLK not needed for rotation test itself, but needed for RL multiplication if you test that.

    # Generate the plaintexts
    m0 = [ct+1 for ct in range(n)]
    print("Plaintext: ", m0)

    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, True, params, params_batched) # Assuming batched=True as you used

    decrypted_initial = decrypt(sk, c0, p1, params, params_batched)
    print("Decrypted (c0): ", decrypted_initial)
    assert decrypted_initial == m0
    print("-" * 80)

    print("\nTESTANDO ROTAÇÕES COM AUTOMORFISMOS VÁLIDOS")
    print("=" * 80)

    Z_star = [d for d in range(1, 2*n, 2) if math.gcd(d, 2*n) == 1]
    print("Automorfismos válidos (d ∈ Z*_{}) = {}".format(2*n, Z_star))

    rotation_results = [] # List to store (d, decrypted_vector)

    for d in Z_star: # Iterate only through valid autom. exponents 'd'
        print("\n" + "-" * 60)
        print(f"Tentando automorfismo d = {d}")
        
        # d = d*7

        # 1) Calcular a chave secreta rotacionada s(x^d)
        sk_rot = rotate_polynomial_coeffs(sk.copy(), d, n, q)
        # print("sk rotacionada (d={}): {}".format(d, view_sk_rot(sk_rot.copy(),q)))

        # 2) Gerar Key Switching Keys de sk_rot para sk
        # ksk permite transformar um ctxt de sk_rot para sk
        # Aqui, KSK é gerada ENCRIPTANDO termos de sk_rot com sk.
        # generate_ks_keys(key_old=sk_rot, key_new=sk, ...)
        #ksk = generate_ks_keys(sk_rot, sk, base_decomposition, q, p1, p2, True, params, params_batched)
        # print(f"KSK gerada para d={d}")

        # 3) Aplicar a rotação x -> x^d DIRETAMENTE ao criptograma c0
        # O ctxt resultante c1 cifra m0(x^d) sob sk_rot
        c1 = rotate_ciphertext(c0, d, n, q)
        # print(f"Ctxt rotacionado c1 para d={d}")

        # 4) Key Switch: Converte c1 (cifrando m0(x^d) sob sk_rot)
        # para um criptograma c2 (cifrando m0(x^d) sob sk)
        #c2 = key_switching(ksk, base_decomposition, c1, params)
        # print(f"Ctxt após key switching c2 para d={d}")

        # 5) Decifrar c2 com sk original
        decrypted_rotated = decrypt(sk_rot, c1, p1, params, params_batched)
        
        # print(params_batched)

        print(f"Decrypted after rotating by d={d}: {decrypted_rotated}")

        # Store result
        rotation_results.append((d, decrypted_rotated.copy()))

    print("\n" + "="*80)
    print("Resultados das Rotações com Automorfismos d ∈ Z*_{}".format(2*n))
    print("="*80)
    for d, result_vec in rotation_results:
        print(f"d={pad_num(d,2)}: {result_vec}")

    # Find specific target rotation (optional, as it might not exist with these ops)
    target_m1 = [2,3,4,1,6,7,8,5] # Left shift 1 within halves [0,1,2,3] and [4,5,6,7]
    target_m2 = [4,1,2,3,8,5,6,7] # Right shift 1 within halves [0,1,2,3] and [4,5,6,7]
    # Add other common rotations if needed (e.g. vector reversal = [8,7,6,5,4,3,2,1])

    print("\nSearching for specific rotations...")
    found_target_m1 = None
    found_target_m2 = None
    found_reversal = None

    for d, result_vec in rotation_results:
        if result_vec == target_m1:
            found_target_m1 = d
        elif result_vec == target_m2:
            found_target_m2 = d
        elif result_vec == m0[::-1]: # Vector reversal
            found_reversal = d

    if found_target_m1 is not None:
        print(f"✅ Found target [2,3,4,1,6,7,8,5] for d = {found_target_m1}")
    else:
        print(f"❌ Target [2,3,4,1,6,7,8,5] not found with d ∈ Z*_{2*n}")

    if found_target_m2 is not None:
        print(f"✅ Found target [4,1,2,3,8,5,6,7] (R-1 within halves) for d = {found_target_m2}")
    else:
        print(f"❌ Target [4,1,2,3,8,5,6,7] (R-1 within halves) not found with d ∈ Z*_{2*n}")

    if found_reversal is not None:
        print(f"✅ Found vector reversal [8,7,6,5,4,3,2,1] for d = {found_reversal}")
    else:
        print(f"❌ Vector reversal [8,7,6,5,4,3,2,1] not found with d ∈ Z*_{2*n}")
     
# ---------------------------------------------------------------------------------------------
def rotate_ciphertext(ct, d, n, q):
    """
    Rotaciona um ciphertext no domínio de coeficientes para automorfismo x -> x^d mod (x^n+1).
    """
    new_mask = rotate_polynomial_coeffs(ct.mask, d, n, q)
    new_body = rotate_polynomial_coeffs(ct.body, d, n, q)
    return Cryptogram(body=new_body, mask=new_mask, batched=ct.batched, q=ct.q)

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
        
    # Corrigir intervalo de coeficientes
    # key_old = [mod_number(x, q) for x in key_old]
            
    for i in range(j):
        tmp = []
        expoente = mod_number(base_decomposition**i,q)
        # print("Expoente: ",expoente)
        for ct in range(n):
            t1 = mod_number(key_old[ct] * expoente, q)
            tmp.append(t1)
        
        # encrypt the old key
        # print(tmp)
        c0 = encrypt_sk(key_new, tmp, q, p1, p2, True, params, params_batched)    
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
    a1 = decompose_poly_list_rot(original_mask, base_decomposition, q)
    
    assert reconstruct_from_decomp_rot(a1, base_decomposition, q) == original_mask
    
    # Verificando a decomposicao binária:
    # posic = 0
    # decomp = [0] * n
    # for item in a1.copy():
    #     tmp = [mod_number(item[ct] * (base_decomposition ** posic), q) for ct in range(n)]
    #     decomp = [mod_number(decomp[ct] + tmp[ct], q) for ct in range(n)]
    #     posic += 1
        
    # print("decomp.....: ", decomp)
    # print("Mascara....: ", original_mask)
    # sys.exit(1)
    
    # if decomp != original_mask:
    #     print("ERRO NA DECOMPOSICAO BINÁRIA =================================")
    #     sys.exit(1)

    
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
            # if j > 0:
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
            #if j > 0:
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
        #sign = 1
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

# -------------------------------------------------------------------------
def decompose_poly_list_rot(poly, base, q):
    n = len(poly)
    max_digits = ceil_log(q, base)
    decomposition = [[] for _ in range(max_digits)]
    
    for coef in poly:
        tmp = []
        val = coef
        for _ in range(max_digits):
            tmp.append(val % base)
            val //= base
        for j in range(max_digits):
            decomposition[j].append(tmp[j])
    
    return decomposition

# ------------------------------------------------------------------------------------
def reconstruct_from_decomp_rot(a1, base, q):
    n = len(a1[0])
    res = [0] * n
    for j, poly in enumerate(a1):
        for i in range(n):
            res[i] = (res[i] + poly[i] * pow(base, j, q)) % q
    return res

# ------------------------------------------------------------------------------------
def key_switching_debug(ksk: list, base_decomposition: int, cryptogram, params):
    q = cryptogram.q
    n = len(cryptogram.body)
    a = cryptogram.mask   # = a(x)
    b = cryptogram.body   # = b(x)

    # 1) decompor a em base B
    a_decomp = decompose_poly_list_rot(a, base_decomposition, q)
    print(">>> Decomposição de a em base", base_decomposition, ":", a_decomp)

    psi_rev, psi_inv_rev, n_inv, bar = params
    num = len(a_decomp)

    new_mask_sum = [0]*n
    new_body_sum = [0]*n

    # 2) para cada componente
    for j in range(num):
        a_j = a_decomp[j]
        KSa_j = ksk[j].mask
        KSb_j = ksk[j].body

        # debug
        print(f"\n-- j = {j}")
        print(" a_j       =", a_j)
        print(" KSa_j.mask=", KSa_j)
        print(" KSa_j.body=", KSb_j)

        # produto máscara
        prod_mask = polymul_ntt(a_j, KSa_j, q, psi_rev, psi_inv_rev, n_inv, bar)
        print(" prod_mask =", prod_mask)
        # produto corpo
        prod_body = polymul_ntt(a_j, KSb_j, q, psi_rev, psi_inv_rev, n_inv, bar)
        print(" prod_body =", prod_body)

        # acumula
        for i in range(n):
            new_mask_sum[i] = mod_number(new_mask_sum[i] + prod_mask[i], q)
            new_body_sum[i] = mod_number(new_body_sum[i] + prod_body[i], q)

        print(" new_mask_sum:", new_mask_sum)
        print(" new_body_sum:", new_body_sum)

    # 3) construir saída
    new_mask = [mod_number(-x, q) for x in new_mask_sum]
    new_body = [mod_number(b[i] - new_body_sum[i], q) for i in range(n)]

    print("\n>>> new_mask final:", new_mask)
    print(">>> new_body final:", new_body)

    return Cryptogram(body=new_body, mask=new_mask, batched=cryptogram.batched, q=q)
