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
    print("ROOTS: ", roots)
    m1 = matriz(roots[0],n, p1)
    m2 = matriz_inv(roots[0],n, p1)
    
    print("=" * 60)
    print("Matrix W:")  
    show_matrix(m1)
    
    print("=" * 60)
    print("Matrix W^-1:")
    show_matrix(m2)

    print("=" * 60)    
    t3 = mul_matrix(m1, m2, p1)    
    print("Matrix W * W^-1:")
    show_matrix(t3)
    
    # n^-1 = 8^-1 mod p1
    print("=" * 60)
    p_inv = modular_inverse(n//2, p1)
    print(p_inv)
    mat = mul_number_matriz(p_inv, m1, p1)
    print("Matrix W * n^-1:")
    show_matrix(mat)
    
    # Multiply by t3 (identity matrix)
    print("=" * 60)
    mat = mul_matrix(mat, t3, p1)
    print("Matrix W * n^-1 * t3:")
    show_matrix(mat)
    
    tr = mul_matrix_vec(mat,[1,2,3,4,5,6,7,8],17)
    print("tr = ", tr)

    print("=" * 100)
    n = 8
    q = 776077649
    p1 = 257
    p2 = 3
    base_decomposition = 2  
    batched = False
    
    # q = generate_modulus(2**28,2**30,8)
    # print(q)
    
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
    m0 = [1,2,3,4,5,6,7,8]
    m1 = [0,1,0,0,0,0,0,0]
    print("Plaintext: ", m0)
    print("Plaintext: ", m1)
    
    # print("Parameters: ", params)
    # print("Parameters batched: ", params_batched)
    print("-" * 20)    
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
    c1 = encrypt_sk(sk, m1, q, p1, p2, batched, params, params_batched)
    
    c2 = multiply_cryptograms(c0, c1, rlk, base_decomposition, params)
    decrypted = decrypt(sk, c2, p1, params, params_batched)
    print("c2 = ", c2)
    print("Decrypted (c2): ", decrypted)
    print("-" * 20)
    
    decrypted = decrypt(sk, c2, p1, params, params_batched)
    print("c2 = ", c2)
    print("Decrypted (C2): ", decrypted)   
    
    # Num caso real teríamos que transformar c2 em c3 homomorficamente
    # ntt_generic(c2.body, psi_rev, q, bar) 
    # ntt_generic(c2.mask, psi_rev, q, bar)    
    # c2 = Cryptogram(body=c2.body, mask=c2.mask, batched=True, q=q)

    decrypted = decrypt(sk, c2, p1, params, params_batched)
    print("c2 = ", c2)
    print("Decrypted (C2): ", decrypted)   

    batched = True
    c3 = encrypt_sk(sk, decrypted, q, p1, p2, batched, params, params_batched)
    print("c3 = ", c3)
    print("-" * 20)
    mask = encrypt_sk(sk, [-1,1,1,1,1,1,1,1], q, p1, p2, batched, params, params_batched)
    c4 = multiply_cryptograms(c3, mask, rlk, base_decomposition, params)
    
    decrypted = decrypt(sk, c4, p1, params, params_batched)
    print("c4 = ", c4)
    print("Decrypted (C4): ", decrypted)    

    # ------------------------------------------------------
    
    print("=" * 100)
    n = 8
    q = 776077649
    p1 = 257
    p2 = 3
    base_decomposition = 2  
    batched = True
    
    # q = generate_modulus(2**28,2**30,8)
    # print(q)
    
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
    # sk = [1] * n
    print("Secret key: ", sk)
    
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    print("n = ", n)
    
    # Generate the plaintexts
    m0 = []
    for i in range(n):
        m0.append(i)
    
    print("Plaintext: ", m0)
    
    print("-" * 20)    
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)    
    decrypted = decrypt(sk, c0, p1, params, params_batched)
    print("Decrypted (C0): ", decrypted)   
    
    rot = 1
    body = rotate_polynomial_coeffs(c0.body, 5**rot, n)
    mask = rotate_polynomial_coeffs(c0.mask, 5**rot, n)
        
    s2 = rotate_polynomial_coeffs(sk, 5**rot, n)
    
    ks = encrypt_sk(sk, s2, q, p1, p2, batched, params, params_batched)   
    print("S2: ", s2)
    print("SK: ", sk)
    print("KS: ", ks)
    
    c1 = Cryptogram(body=body, mask=mask, batched=True, q=q)
    decrypted = decrypt(s2, c1, p1, params, params_batched)
    print("Decrypted (C1): ", decrypted)   

    # decomposição da chave Galois    
    kbody = decompose_poly_list(ks.body,base_decomposition,q)
    kmask = decompose_poly_list(ks.mask,base_decomposition,q)    

    print("kbody: ", kbody)
    print("kmask: ", kmask)

    # verificando a decomposição
    decomp = kbody.copy()
    lista = [0 for i in range(n)]
    i = 0
    for tmp in decomp:
        t1 = base_decomposition ** i
        t2 = [mod_number(tmp[j] * t1, q) for j in range(n)]
        lista = [mod_number(lista[j] + t2[j], q) for j in range(n)]
        i += 1
        
    if lista != ks.body:
        print("list: ", lista)
        print("body: ", ks.body)
        print("ERRO NA DECOMPOSIÇÃO - CORPO")
        sys.exit(1)

    # verificando a decomposição
    decomp = kmask.copy()
    lista = [0 for i in range(n)]
    i = 0
    for tmp in decomp:
        t1 = base_decomposition ** i
        t2 = [mod_number(tmp[j] * t1, q) for j in range(n)]
        lista = [mod_number(lista[j] + t2[j], q) for j in range(n)]
        i += 1
        
    if lista != ks.mask:
        print("list: ", lista)
        print("body: ", ks.body)
        print("ERRO NA DECOMPOSIÇÃO - MÁSCARA")
        sys.exit(1)

    a = [mod_number(ks.mask[j] * -1, q) for j in range(n)]
    b = [mod_number(c1.body[j] - ks.body[j], q) for j in range(n)]

    c2 = Cryptogram(body=b, mask=a, batched=True, q=q)
    decrypted = decrypt(sk, c2, p1, params, params_batched)
    print("Decrypted (C2): ", decrypted)
    
    # # -------------------------------------------------------
    
    # # # Parâmetros iniciais (exemplo)
    # # n = 8
    # # q = 776077649
    # # p1 = 257
    # # p2 = 3
    # # base_decomposition = 2
    # # batched = True
    # # # params e params_batched são obtidos via generate_parameters(n, q) e generate_parameters(n, p1), respectivamente.
    # # # sk é sua chave secreta.
    
    # # Generate the keys
    # sk = create_sk(n, 0, 1, q)
    # # sk = [1] * n
    # print("Secret key: ", sk)


    # # Gera um ciphertext c0 a partir de um plaintext m0 (via encrypt_sk)
    # m0 = [1, 2, 3, 4, 5, 6, 7, 8]
    # c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)

    # # Suponha que você deseja rotacionar 1 posição (para a esquerda)
    # rotations = [1]
    # # Gera o dicionário de chaves Galois para as rotações desejadas:
    # galois_keys = generate_galois_keys(sk, rotations, q, p1, p2, batched, params, params_batched)
    # # Seleciona a switching key para rotação 1:
    # ks = galois_keys[3]

    # # Calcula o fator de rotação: 5^1 mod (2*n)
    # rot_factor = pow(5, 1, 2*n)
    # # Aplica rotação ao ciphertext c0 para obter c1 (usando a função de rotação já definida)
    # c1 = Cryptogram(
    #     body=rotate_polynomial_coeffs(c0.body, rot_factor, n),
    #     mask=rotate_polynomial_coeffs(c0.mask, rot_factor, n),
    #     batched=batched,
    #     q=q
    # )

    # # Agora, aplica a key switching via a função apply_galois para obter o ciphertext final c2
    # p = (5**1) % (2*n)
    # c2 = apply_galois(c1, p, galois_keys, base_decomposition, q)

    # # c2 deverá ser decifrável com a chave secreta original
    # decrypted = decrypt(sk, c2, p1, params, params_batched)
    # print("Decrypted (after rotation and key switching):", decrypted)        
        
    # # kbody = decompose_poly_list(ks.body,base_decomposition,q)
    # # kmask = decompose_poly_list(ks.mask,base_decomposition,q)

    # # c1body = decompose_poly_list(c1.body,base_decomposition,q)
    # # c1mask = decompose_poly_list(c1.mask,base_decomposition,q)    
        
    # # print("kbody: ", kbody)
    # # print("kmask: ", kmask)
    
    # # decomp = kbody.copy()
    # # lista = [0 for i in range(n)]
    # # i = 0
    # # for tmp in decomp:
    # #     t1 = base_decomposition ** i
    # #     t2 = [mod_number(tmp[j] * t1, q) for j in range(n)]
    # #     lista = [mod_number(lista[j] + t2[j], q) for j in range(n)]
    # #     i += 1
        
    # # if lista != ks.body:
    # #     print("list: ", lista)
    # #     print("body: ", ks.body)
    # #     print("ERRO NA DECOMPOSIÇÃO")
    # #     sys.exit(1)
    # # else:
    # #     print("DECOMPOSIÇÃO OK - Body")

    # # decomp = kmask.copy()
    # # lista = [0 for i in range(n)]
    # # i = 0
    # # for tmp in decomp:
    # #     t1 = base_decomposition ** i
    # #     t2 = [mod_number(tmp[j] * t1, q) for j in range(n)]
    # #     lista = [mod_number(lista[j] + t2[j], q) for j in range(n)]
    # #     i += 1
        
    # # if lista != ks.mask:
    # #     print("list: ", lista)
    # #     print("Mask: ", ks.mask)
    # #     print("ERRO NA DECOMPOSIÇÃO")
    # #     sys.exit(1)
    # # else:
    # #     print("DECOMPOSIÇÃO OK - Mask")


# ------------------------------------------------------------------------------
def generate_galois_keys(sk, rotations, q, p1, p2, batched, params, params_batched):
    """
    Gera um dicionário de chaves Galois (switching keys) para as rotações desejadas, usando a encriptação real.
    
    Para cada rotação r na lista "rotations":
      - Se r == "col_swap": define p = M - 1, com M = 2*N (N = len(sk));
      - Se r >= 0: define p = 3^r mod (2*N);
      - Se r < 0: define p = 3^(N/2 - |r|) mod (2*N).
      
    Em seguida, computa s′ = apply_automorphism_to_poly(sk, p, q) e 
    gera a switching key ks = encrypt_sk(sk, s′, q, p1, p2, batched, params, params_batched).
    
    Retorna:
      dict: mapeia cada valor p (int) à switching key (objeto Cryptogram).
    """
    N = len(sk)
    M = 2 * N
    galois_keys = {}
    for r in rotations:
        if r == "col_swap":
            p = M - 1
        elif r >= 0:
            p = pow(3, r, M)
        else:
            p = pow(3, (N // 2 - abs(r)), M)
        # Computa s' = s(X^p)
        s_prime = apply_automorphism_to_poly(sk, p, q)
        # Realiza a encriptação real de s_prime usando a chave secreta sk
        ks = encrypt_sk(sk, s_prime, q, p1, p2, batched, params, params_batched)
        
        print("P = ", p)
        galois_keys[p] = ks
    return galois_keys

def apply_galois(ciphertext, p, galois_keys, base_decomposition, q):
    """
    Aplica a operação de key switching utilizando a chave Galois para o automorfismo X -> X^p.
    
    Procedimento conforme pseudocódigo:
      1. Aplica o automorfismo aos polinômios do ciphertext:
           c0_auto(X) = c0(X^p) mod (X^N+1)
           c1_auto(X) = c1(X^p) mod (X^N+1)
      2. Decompõe c1_auto na base w (base_decomposition), obtendo coeffs = decomposition(c1_auto, w).
      3. Para cada dígito j, utiliza o par (k0_j, k1_j) da switching key correspondente (obtida em galois_keys[p])
         para atualizar:
              new_c0 += coeffs[j] * k0_j  
              new_c1 += coeffs[j] * k1_j
      4. Retorna o novo ciphertext (new_c0, new_c1) cifrado com a chave secreta original.
    
    Args:
      ciphertext (Cryptogram): objeto contendo os polinômios (body e mask).
      p (int): elemento de Galois (automorfismo aplicado).
      galois_keys (dict): dicionário gerado por generate_galois_keys.
      base_decomposition (int): base (w) usada na decomposição.
      q (int): módulo.
      
    Retorna:
      Cryptogram: novo ciphertext resultante do key switching.
    """
    # 1. Aplicar automorfismo nos polinômios do ciphertext.
    c0_auto = apply_automorphism_to_poly(ciphertext.body, p, q)
    c1_auto = apply_automorphism_to_poly(ciphertext.mask, p, q)
    
    # 2. Obtém a switching key para p.
    evk = galois_keys.get(p)
    if evk is None:
        print("Erro: Chave Galois para p =", p, "não encontrada.")
        sys.exit(1)
    
    # 3. Decompor c1_auto na base w.
    coeffs = decompose_poly(c1_auto, base_decomposition, q)  # Retorna lista de K polinômios (cada um de tamanho N)
    
    # Inicializa new_c0 e new_c1:
    new_c0 = c0_auto[:]  # cópia de c0_auto
    N = len(c0_auto)
    new_c1 = [0] * N
    
    # Para cada dígito j da decomposição, atualiza new_c0 e new_c1:
    for j, comp in enumerate(coeffs):
        # Evk: lista de pares (k0_j, k1_j)
        k0_j, k1_j = evk[j]
        prod0 = multiply_poly_mod(comp, k0_j, q)
        prod1 = multiply_poly_mod(comp, k1_j, q)
        new_c0 = poly_add_mod(new_c0, prod0, q)
        new_c1 = poly_add_mod(new_c1, prod1, q)
    
    return Cryptogram(body=new_c0, mask=new_c1, batched=ciphertext.batched, q=q)

# ---------------------------------------------------------------------------------
def apply_automorphism_to_poly(poly, p, q):
    """
    Aplica o automorfismo que transforma X em X^p em um polinômio.
    Dados N = len(poly) e M = 2*N, para cada índice i:
       - Calcula new_index = (i*p) mod M.
       - Se new_index >= N, aplica sinal negativo e subtrai N.
    Retorna o polinômio resultante (reduzido módulo q).
    """
    N = len(poly)
    M = 2 * N
    new_poly = [0] * N
    for i in range(N):
        new_index = (i * p) % M
        sign = 1
        if new_index >= N:
            sign = -1
            new_index = new_index - N
        new_poly[new_index] = mod_number(poly[i] * sign, q)
    return new_poly

def decompose_poly(poly, base, q):
    """
    Decompõe cada coeficiente de 'poly' em uma soma:
         poly[i] = sum_{j=0}^{K-1} d_{i,j} * base^j,
    retornando K polinômios (cada um representando os dígitos para todos os coeficientes).
    K é definido como ceil(log_q / log_base).
    """
    n = len(poly)
    K = math.ceil(math.log(q, base))
    decomp = [[0] * n for _ in range(K)]
    for i in range(n):
        coeff = poly[i]
        for j in range(K):
            decomp[j][i] = coeff % base
            coeff //= base
    return decomp

def random_polynomial(n, q):
    """Gera um polinômio aleatório de grau n com coeficientes em Z_q."""
    return [random.randint(0, q-1) for _ in range(n)]

def small_gaussian_error(n, bound=10):
    """Retorna um polinômio com erro pequeno (valores aleatórios em [-bound, bound])."""
    return [random.randint(-bound, bound) for _ in range(n)]

# ---------------------------------------------------------------------------
def rotate_polynomial_coeffs(coeffs, d, N):
    tmp = []
    ret = []
    for i in range(N):
        # print("i =", i, d, i*d, "[",(i*d)%N,"]")
        t1 = d % N
        signal = -1**t1
        tmp.append(((i * d) % N) * signal)
    
    for i in tmp:
        ret.append(coeffs[i])
    
    return ret

# ---------------------------------------------------------------------------
# Funções auxiliares para operações polinomiais em R = ℤ_q[X]/(Xⁿ+1)
def poly_add_mod(p1, p2, q):
    """Soma polinomial com redução módulo q."""
    return [mod_number((a + b) , q) for a, b in zip(p1, p2)]

def poly_sub_mod(p1, p2, q):
    """Subtração polinomial com redução módulo q."""
    return [mod_number((a - b) , q) for a, b in zip(p1, p2)]

# ---------------------------------------------------------------------------
def key_switching(c1, ks, base_decomposition, q):
    """
    Realiza a operação de key switching em um criptograma RLWE conforme o Algoritmo 3 (imagem anexa).

    Dados:
      - c1: Criptograma RLWE original sob a chave z1 (objeto Cryptogram);
      - ks: Switching key – encriptação de s2 (ou de B^j · z1) sob a nova chave z2;
      - base_decomposition: base B usada na decomposição dos coeficientes;
      - q: módulo da encriptação.
      
    O algoritmo executa:
       a' = - ∑_j [B^j * (digit_j(c1.body) * digit_j(ks.body))]
       b' = c1.mask - ∑_j [B^j * (digit_j(c1.mask) * digit_j(ks.mask))]
       
    Retorna:
      Novo criptograma (objeto Cryptogram) encriptado sob a nova chave z2.
    """
    n = len(c1.body)
    
    # Decompor c1.body e c1.mask e as respectivas partes da switching key
    c1_body_decomp = decompose_poly_list(c1.body, base_decomposition, q)
    c1_mask_decomp = decompose_poly_list(c1.mask, base_decomposition, q)
    ks_body_decomp = decompose_poly_list(ks.body, base_decomposition, q)
    ks_mask_decomp = decompose_poly_list(ks.mask, base_decomposition, q)
    
    # Inicializa os acumuladores para o novo corpo (a') e para a correção da máscara (que será subtraída de c1.mask)
    new_body = [0] * n       # a' acumulado com sinal negativo
    accum_mask = [0] * n     # soma dos termos para a máscara
    
    # Para cada dígito da decomposição
    for j in range(len(c1_body_decomp)):
        # Calcula o peso para a fatia j: B^j
        weight = base_decomposition ** j
        
        # Multiplica as fatias correspondentes do criptograma e da switching key
        term_body = multiply_poly_mod(c1_body_decomp[j], ks_body_decomp[j], q)
        term_mask = multiply_poly_mod(c1_mask_decomp[j], ks_mask_decomp[j], q)
        
        # Multiplica cada termo pelo peso correspondente e reduz módulo q
        term_body = [mod_number(coef * weight, q) for coef in term_body]
        term_mask = [mod_number(coef * weight, q) for coef in term_mask]
        
        # Acumula o corpo com sinal negativo (conforme a fórmula: a' = - Σ_j (B^j * termo))
        new_body = poly_sub_mod(new_body, term_body, q)
        # Acumula os termos para a máscara
        accum_mask = poly_add_mod(accum_mask, term_mask, q)
    
    # A nova máscara é b' = c1.mask - acumulador de termos (correção)
    new_mask = poly_sub_mod(c1.mask, accum_mask, q)
    
    return Cryptogram(body=new_body, mask=new_mask, batched=c1.batched, q=q)

# -------------------------------------------------------------------------------------
def mul_number_matriz(number, m, modulus):
    """
    Multiplies a number by a matrix m.
    Returns the matrix resulting from the multiplication.
    """
    if not m:
        return []
    
    result = []
    for row in m:
        new_row = [mod_number(number * element,modulus) for element in row]
        result.append(new_row)
    
    return result
    
# ---------------------------------------------------------------------------
def show_matrix(m):
    """
    Displays the matrix m with the elements aligned by column.
    Each column uses the maximum width between its elements.    
    """
    if not m:
        return
    
    # Considerando que todas as linhas tenham o mesmo número de colunas:
    num_colunas = len(m[0])

    # Calcula a largura máxima para cada coluna
    col_widths = []
    for j in range(num_colunas):
        # Para cada coluna, converte todos os elementos para string
        # e pega o comprimento máximo
        largura = max(len(str(linha[j])) for linha in m)
        col_widths.append(largura)

    # Imprime cada linha, formatando cada elemento com o respectivo tamanho da coluna
    for linha in m:
        linha_formatada = " ".join(f"{str(valor):>{col_widths[i]}}" for i, valor in enumerate(linha))
        print(linha_formatada)
# ---------------------------------------------------------------------------
def extract_col(m,number_col):
    size = len(m)
    ret_line = []
    for ct in range(size):
        tmp = m[ct]
        ret_line.append(tmp[number_col])
    
    return ret_line            

# ---------------------------------------------------------------------------
def extract_row(m,number_row):
    return m[number_row]

# ---------------------------------------------------------------------------
def mul_col_row(col, row):
    ret = 0
    if len(col) != len(row):
        print("Error: col and row must have the same size")
        sys.exit(1)
        
    for ct in range(len(col)):
        ret += col[ct] * row[ct]
    
    return ret

# ---------------------------------------------------------------------------
def mul_matrix(m1, m2, modulus):
    """
    Multiplies two matrices m1 and m2.
    Returns the matrix resulting from the multiplication.
    """
    if len(m1) != len(m2[0]):
        print("Error: m1 and m2 must have the same size")
        sys.exit(1)
        
    result = []
    for i in range(len(m1)):
        row = []
        for j in range(len(m2[0])):
            row.append(mod_number(mul_col_row(extract_col(m1, i), extract_row(m2, j)),modulus))
        result.append(row)
    
    return result

# ---------------------------------------------------------------------------
def mul_matrix_vec(m1, m2, modulus):
    """
    Multiplies two matrices m1 and m2.
    Returns the matrix resulting from the multiplication.
    """
    if len(m1) != len(m2):
        print("Error: m1 and m2 must have the same size")
        sys.exit(1)
        
    result = []
    for i in range(len(m1)):
        for j in range(len(m2)):
            tmp = extract_row(m1, i)
            ret = 0
            for ct2 in tmp:
                ret += mod_number(ct2 * m2[j],modulus)
        result.append(mod_number(ret,modulus))
    
    return result

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
        
    # print("MAT:  ",mat)
    # print("MAT2: ",mat2)
    # print("POT: ", pot)  
    
    # generate matrix W (page 147)
    size = len(mat)
    # print("Size: ", size)
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

    # print("Matrix W:")
    # for ct in range(size):
    #     print(mat[ct])
        
    return mat

# -------------------------------------------------------------------------------
# paper: https://www.arxiv.org/pdf/2503.05136 - page 148        
def matriz_inv(w, n, p1):
    
    # generate matrix of size 8x8
    # and fill it with 1's and primitive root
    mat = []
    for ct in range(n): 
        mat.append([1,w,w,w,w,w,w,w])
        
    # fill the first row with 1's and func_j
    mat2 = [ ]
    limit = n//2
    for ct in range(limit): 
        tmp = [1,func_j(1,ct,n), func_j(1,ct,n), func_j(1,ct,n), func_j(1,ct,n),
               func_j(1,ct,n), func_j(1,ct,n), func_j(1,ct,n)]
        mat2.append(tmp)
    for ct in range(limit): 
        tmp = [1,func_j(0,ct,n), func_j(0,ct,n), func_j(0,ct,n), func_j(0,ct,n),
               func_j(0,ct,n), func_j(0,ct,n), func_j(0,ct,n)]
        mat2.append(tmp)
        
    pot = []
    for ct in range(n): 
        tmp = [0,1,2,3,4,5,6,7]
        pot.append(tmp)
        
    # print("MAT:  ",mat)
    # print("MAT2: ",mat2)
    # print("POT: ", pot)  
    
    # generate matrix W (page 147)
    size = len(mat)
    # print("Size: ", size)
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

    # print("Matrix W^-1:")
    # for ct in range(size):
    #     print(mat[ct])
        
    return mat
# ---------------------------------------------------------------------------
def func_j(t, h, n):
    if t == 1:
        return mod_number(5**h, 2*n)
    elif t == 0:
        return mod_number(-5**h, 2*n)
    else:
        sys.exit(1)

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

