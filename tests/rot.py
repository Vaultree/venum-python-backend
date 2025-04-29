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
    p1 = 17
    p2 = 3
    base_decomposition = 2  
    batched = True
    
    q = generate_modulus(2**60,2**63,n)
    print(" Q = ", q)
    
    # while True:
    #     p1 = generate_modulus(2**7,2**12,n)
    #     if p1 % 4 == 1:
    #         break    
        
    # q = 65537
    # p1 = 17
    
    print("P1 = ", p1)
    print("N = ", n)

    # verify if p1 is equal to 1 mod 2n (necessary condition for CRT)
    if p1 % (2*n) != 1:
        print("Error: p1 must be equal to 1 mod 2n | n = ",n, "p1 mod 2n", p1 % (2*n))
        sys.exit(1)

    if p1 % 4 != 1:
        print("Error: p1 must be equal to 1 mod 4!")
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
    
    # sk = [1,0,0,0,0,0,0,0]
    
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    # Generate the plaintexts
    m0 = []
    for ct in range(n):
        m0.append((ct+1)* 1)
        
    r1 = encode_input(m0, n, p1)
    r2 = decode_input(r1, p1)
    assert m0 == r2
            
    # m0 = [4, 1, 2, 3, 8, 5, 6, 7]
    # m0 = [10,20,30,40,50,60,70,80]
    # m0 = [1,1,1,1,1,1,1,1]
    print("Plaintext: ", m0)
    
    # m1 = m0.copy()
    # m2 = bit_reverse_order(m1, n)
    # m3 = bit_reverse_order(m2, n)
    # print("m1 | m3", m1, m3)
    
    #print("Parameters: ", params)
    #print("Parameters batched: ", params_batched)
    
    # Encrypt the plaintexts secret key
    c0 = encrypt_sk(sk, encode_input(m0,n,p1), q, p1, p2, batched, params, params_batched)
    print("C0: ", c0)
    
    decrypted = decrypt(sk, c0, p1, params, params_batched)
    decrypted = decode_input(decrypted, p1)            
    print("Decrypted (c0): ", decrypted)
    print("-" * 80)
    
    #sys.exit(0)

    print("C0: ", c0)
    
    assert decrypted == m0
    
    map = []
    pot = 3

    # found the automorphism
    # 1) monta o conjunto de automorfismos    
    Z_star = [d for d in range(1, 2*n, 2) if math.gcd(d, 2*n) == 1]
    print("Automorfismos válidos (Z*_{}) = {}".format(2*n, Z_star))

    root_list = Z_star
    root_list = [5]
    # root_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    print("root: ", root_list)
    
    root_list = [5]
    # start, end = -4*n, 4*n+1
    # cria uma lista de inteiros de start até end, inclusive
    # rot_vetor = list(range(start, end+1))
    
    map_rot_test = []
    map_indice = []
    
    separator = 60
    start = False
    for pot in root_list:
        for signal in [0,1]:
            for rot in range(2*n):
            # for rot in rot_vetor:
                batched = True
                
                # at this point we have explored all possible rotations
                
                # if rot < (n/2):
                #     power = pot
                #     rotation = power**rot  
                # else:
                #     power = -pot
                #     rotation = power**rot

                if signal == 0:
                    power = pot
                    rotation = power**rot  
                else:
                    power = -pot
                    rotation = power**rot
                    
                rotation = rot
                        
                # key rotation
                sk_rot = rotate_polynomial_coeffs(sk, rotation, n,q)
                
                # sk_rot = rotate_polynomial_coeffs(sk, -1, n,q)
                # sk_rot = rotate_polynomial_coeffs(sk_rot, 3, n,q)
                # sk_rot = rotate_polynomial_coeffs(sk_rot, -1, n,q)
                
                print("-" * 80)
                print("Rotation......: ", rot)
                print("sk............: ", sk)   
                # 2) só para exibí-la de forma humana, use uma cópia!
                signed = view_sk_rot(sk_rot.copy(), q)
                print("sk_rot (signed):", signed)
                print("P1 = ", p1, " | P2 = ", p2, " | Q = ", q)
                
                # generate rotation key
                ksk = generate_ks_keys(sk_rot, sk, base_decomposition, q, p1, p2, batched, params, params_batched)
                
                # para multiplicar os criptogramas com a chave sk_rot é necessário gerar as chaves
                # de relinearização.
                rlk_sk_rot = generate_rlk(sk_rot, base_decomposition, p1, p2, q, params)
                
                # print("KSK: ", ksk)
                # print(">>> len(ksk) esperado:", ceil_log(q, base_decomposition))
                # print(">>> len(ksk) obtido   :", len(ksk))
                # for idx, ck in enumerate(ksk):
                #     print(f"  ksk[{idx}].mask[:5] = {ck.mask[:5]}")
                #     print(f"  ksk[{idx}].body[:5] = {ck.body[:5]}")
                
                #generate new Cryptogram
                # body = c0.body
                # mask = c0.mask
                # body = rotate_polynomial_coeffs(body, rotation, n,q)
                # mask = rotate_polynomial_coeffs(mask, rotation, n,q)
                # c1 = Cryptogram(body=body, mask=mask, batched=c0.batched, q=c0.q)
                
                c1 = rotate_ciphertext(c0, rotation, n, q)
                
                # c1 = rotate_ciphertext(c0, -1, n, q)
                # c1 = rotate_ciphertext(c1, 3, n, q)
                # c1 = rotate_ciphertext(c1, -1, n, q)
                #print("C1: ", c1)
                
                # decrypted new Cryptogram with the new secret key (rotated)
                
                decrypted = decrypt(sk_rot, c1, p1, params, params_batched)
                decrypted = decode_input(decrypted, p1)
                print("Decrypted (c1).: ", decrypted)
                # print("-" * separator)
                
                if contains_vector(map_rot_test, decrypted) == False:
                    map_rot_test.append(decrypted)
                    map_indice.append(rotation)
                
                continue
                # print("-" * separator)
                
                c2 = key_switching(ksk, base_decomposition, c1, params)
                decrypted_2 = decrypt(sk, c2, p1, params, params_batched)
                
                decrypted_sk = decrypt(sk, c1, p1, params, params_batched)
                
                print("Decrypted with sk_rot (c1)...: ", decrypted)
                print("Decrypted with sk     (c1)...: ", decrypted_sk)
                print("Decrypted key_switching (c2).: ", decrypted_2)

                # testar mascara multiplicativa            
                c_mul = encrypt_sk(sk_rot, [1,0,0,0,0,0,0,0], q, p1, p2, batched, params, params_batched)
                decrypted_cmul = decrypt(sk_rot, c_mul, p1, params, params_batched)
                print("Decrypted (c_mul with sk_rot).: ", decrypted_cmul)
                
                c1_sk_rot = encrypt_sk(sk_rot, [-8,-7,-6,-5,-4,-3,-2,-1], q, p1, p2, batched, params, params_batched)
                c_sum = sum_cryptograms(c1, c1_sk_rot)
                c_sum_dec = decrypt(sk_rot, c_sum, p1, params, params_batched)
                print("Decrypted (c1 + c1_sk_rot = 0).: ", c_sum_dec)
                
                print("-" * 60)
                c_mul = encrypt_sk(sk_rot, [1,0,0,0,0,0,0,0], q, p1, p2, batched, params, params_batched)
                print("sk_rot: ", sk_rot)
                decrypted_c1 = decrypt(sk_rot, c1, p1, params, params_batched)
                print("Decrypted (c1 with sk_rot).: ", decrypted_c1)
                c4 = multiply_cryptograms(c1, c_mul, rlk_sk_rot, base_decomposition, params)
                decrypted_c4 = decrypt(sk_rot, c4, p1, params, params_batched)
                print("Decrypted with (c4 = c1 * c_mul).: ", decrypted_c4)
                
                print("-" * 60)
                # verificando se a multiplicação é correta
                c_mul = encrypt_sk(sk_rot, [1,0,0,0,0,0,0,0], q, p1, p2, batched, params, params_batched)
                c8_sk_rot = encrypt_sk(sk_rot, [8,7,6,5,4,3,2,1], q, p1, p2, batched, params, params_batched)
                assert [8,7,6,5,4,3,2,1] == decrypt(sk_rot, c8_sk_rot, p1, params, params_batched)
                c8 = multiply_cryptograms(c8_sk_rot, c_mul, rlk_sk_rot, base_decomposition, params)
                decrypted_c8 = decrypt(sk_rot, c8, p1, params, params_batched)
                print("Decrypted with (c8 = c8 * c_mul).: ", decrypted_c8)

                print("-" * 60)
                # verificando se a multiplicação é correta (plaintext)
                # Necessária codificação batched porque o criptograma está cifrado no modo batched
                # (parametros relativos a p1)
                (psi_rev, psi_inv_rev, n_inv, bar) = params_batched
                c_mul = [1,0,0,0,0,0,0,0]
                intt_generic(c_mul, psi_inv_rev, n_inv, p1, bar)
                
                c8_sk_rot = encrypt_sk(sk_rot, [8,7,6,5,4,3,2,1], q, p1, p2, batched, params, params_batched)
                assert [8,7,6,5,4,3,2,1] == decrypt(sk_rot, c8_sk_rot, p1, params, params_batched)
            
                # Cuidado ao passar os parametros para o criptograma (são parametros relativos a q)
                (psi_rev, psi_inv_rev, n_inv, bar) = params
                body = polymul_ntt(c8_sk_rot.body, c_mul.copy(), q, psi_rev, psi_inv_rev, n_inv, bar)
                mask = polymul_ntt(c8_sk_rot.mask, c_mul.copy(), q, psi_rev, psi_inv_rev, n_inv, bar)
                
                c9 = Cryptogram(body=body, mask=mask, batched=c8_sk_rot.batched, q=c8_sk_rot.q)        
                decrypted_c9 = decrypt(sk_rot, c9, p1, params, params_batched)
                print("Decrypted with (c9 = c9 * c_mul).: ", decrypted_c9)
                
                assert decrypted_c9 == [8,0,0,0,0,0,0,0]
                
                
                print("-" * 60)
                c_mul = encrypt_sk(sk, [1,0,0,0,0,0,0,0], q, p1, p2, batched, params, params_batched)
                decrypted_cmul = decrypt(sk, c_mul, p1, params, params_batched)
                print("Decrypted (c_mul).: ", decrypted_cmul)
                print("Decrypted c0: ", decrypt(sk, c0, p1, params, params_batched))
                c5 = multiply_cryptograms(c0, c_mul, rlk, base_decomposition, params)
                decrypted_c5 = decrypt(sk, c5, p1, params, params_batched)
                print("Decrypted with (c5 = c0 * c_mul).: ", decrypted_c5)
                
                print("=" * 90)
                
                sum4 = 0
                sum5 = 0
                for ct in range(n):
                    sum4 = mod_number(sum4 + decrypted_c4[ct], q)
                    sum5 = mod_number(sum5 + decrypted_c5[ct], q)   
                    
                if sum5 != m0[0]:
                    print("Error: ", sum5)            
                    #continue

                if sum4 != decrypted[0]:
                    print("Error: ", sum4)            
                    #continue
                
                # ------------------------------------------------------
                # Testa Key Switching identidade
                # ksk_identity = generate_ks_keys(sk, sk, base_decomposition, q, p1, p2, batched, params, params_batched)
                # c3 = key_switching(ksk_identity, base_decomposition, c0, params)
                # dec3 = decrypt(sk, c3, p1, params, params_batched)
                # dec0 = decrypt(sk, c0, p1, params, params_batched)
                # assert dec3 == dec0, "Key switching identidade falhou!"
                
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

    print("Rotações verificadas: ")
    posic = 0
    for ct in map_rot_test:
        if contains_message(ct, m0) == True:
            print(map_indice[posic],": ",ct) 
        posic +=1
                
    # Show the results
    print("roots: ", root_list)
    for ct in map:
        print("Power = ", pad_num(ct[0],3), "| Rot = ", pad_num(ct[1],3),"| Rotation = ", pad_num(ct[2],8) ,"| Decrypted = ", ct[3])
        # print(ct[3])
        # print("-" * separator)
    print("Map size: ", len(map))
    
    # t = 17
    # raiz = 3
    # pot = 5
    
    # print("*" * separator)
    # w = generate_w(raiz, pot, 8, t)
    
    # print("*" * separator)
    # w1 = generate_w1(raiz, pot, 8, t)
    
    # print("*" * separator)
    # d = multiply_matrices(w1, w, t)
    # print("Resultado: ")
    
    # mat_i = d.copy()
    
    # for ct in d:
    #     print(ct)

    # # n−1Wˆ · IRn
    # modinv = modular_inverse(8//2,t)
    # modinv = 15
    
    # m1 = mult_number_matrix(modinv, w, t)
    # print("*" * separator)
    # print("modinv: ", modinv)
    # print("Resultado: ")
        
    # m2 = multiply_matrices(m1, d.copy(), t)
        
    # m3 = mul_matrix_vector(m2, [1,2,3,4,5,6,7,8], t)

    # print("Resultado [m3]: ", m3)
    # print("*" * separator)

    # m4 = mul_matrix_vector(w1, m3, t)
    # # m4 = mul_matrix_vector(d, m4, t)
    # print("Resultado: ", m4)

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

# -------------------------------------------------------------------------------------
# Suas funções mod_number, etc.
def mod_number(a, q):
    """Computa (a mod q), garantindo resultado no intervalo [0, q-1]."""
    return (a % q + q) % q
# -------------------------------------------------------------------------------------
def modular_inverse(a, m):
    """Calcula o inverso modular de a modulo m."""
    # Implementação usando Extended Euclidean Algorithm
    # Referência: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    # g = gcd(a, m)
    # ax + my = g
    # Se g == 1, então ax = 1 mod m, e x é o inverso
    g, x, y = extended_gcd(a, m)
    if g != 1:
        # raise ValueError(f"O inverso modular não existe para {a} mod {m}")
        print(f"O inverso modular não existe para {a} mod {m}")
        sys.exit(1)
    return x % m # Retorna o inverso no intervalo [0, m-1]
# -------------------------------------------------------------------------------------
def extended_gcd(a, b):
    """Calcula gcd(a, b) e coeficientes x, y tal que ax + by = gcd(a, b)."""
    if a == 0:
        return b, 0, 1
    else:
        gcd, x, y = extended_gcd(b % a, a)
        return gcd, y - (b // a) * x, x
# -------------------------------------------------------------------------------------
# Funções de multiplicação fornecidas por você (incluindo uma corrigida para vetor)
def mul_matrix_vector(matrix, vector, t):
    n = len(vector)
    rows = len(matrix)
    cols = len(matrix[0]) if rows > 0 else 0 # Corrigir aqui para matrizes vazias
    
    if cols != n:
        print(f"Erro: As dimensões da matriz ({rows}x{cols}) e do vetor ({n}) não são compatíveis para M * v.")
        sys.exit(1)
    
    result = []
    for i in range(rows):
        tmp = 0
        for j in range(n):
            tmp = mod_number(tmp + (matrix[i][j] * vector[j]), t)
        result.append(tmp)
        
    return result
# -------------------------------------------------------------------------------------
def multiply_matrices(A, B, t):
    """Multiplica duas matrizes A e B mod t."""
    rows_A = len(A)
    cols_A = len(A[0]) if rows_A > 0 else 0
    rows_B = len(B)
    cols_B = len(B[0]) if rows_B > 0 else 0
    
    if cols_A != rows_B:
        raise ValueError(f"Incompatíveis para multiplicação: A é {rows_A}×{cols_A}, B é {rows_B}×{cols_B}")
    
    C = [[0] * cols_B for _ in range(rows_A)]
    
    for i in range(rows_A):
        for j in range(cols_B):
            total = 0
            for k in range(cols_A):
                total = mod_number(total + (A[i][k] * B[k][j]), t)
            C[i][j] = total # Total já está no módulo
    
    return C
# -------------------------------------------------------------------------------------
def mult_scalar_vector(scalar, vector, t):
    """Multiplica um escalar por um vetor mod t."""
    return [mod_number(scalar * elem, t) for elem in vector]
# -------------------------------------------------------------------------------------
def encode_input(input, n, t):
    # Calcular n^-1 mod t
    n_inv = modular_inverse(n, t)

    # Matriz W_hat da página 148/149
    W_hat = [
        [ 1,  1,  1,  1,  1,  1,  1,  1],
        [12, 14,  5,  3, 10, 11,  7,  6],
        [ 8,  9,  8,  9, 15,  2, 15,  2],
        [11,  7,  6, 10, 14,  5,  3, 12],
        [13, 13, 13, 13,  4,  4,  4,  4],
        [ 3, 12, 14,  5,  6, 10, 11,  7],
        [ 2, 15,  2, 15,  9,  8,  9,  8],
        [ 7,  6, 10, 11,  5,  3, 12, 14]
    ]

    # Matriz Identidade Reversa I_R para n=8
    IR = [
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0]
    ]


    # Passo 1: Calcular M_prod = W_hat * IR mod t
    M_prod = multiply_matrices(W_hat, IR, t)

    # Passo 2: Calcular v_temp = M_prod * v mod t
    v_temp = mul_matrix_vector(M_prod, input, t)

    # Passo 3: Calcular m = n_inv * v_temp mod t
    m_encoded = mult_scalar_vector(n_inv, v_temp, t)

    return m_encoded

# ----------------------------------------------------------------
def decode_input(input, t):
        # Matriz W_hat_star (mod 17)
    W_hat_star = [
        [ 1,  3,  9, 10, 13,  5, 15, 11],
        [ 1,  5,  8,  6, 13, 14,  2, 10],
        [ 1, 14,  9,  7, 13, 12, 15,  6],
        [ 1, 12,  8, 11, 13,  3,  2,  7],
        [ 1,  6,  2, 12,  4,  7,  8, 14],
        [ 1,  7, 15,  3,  4, 11,  9, 12],
        [ 1, 11,  2,  5,  4, 10,  8,  3],
        [ 1, 10, 15, 14,  4,  6,  9,  5]
    ]

    # Calcular v_j3 = W_hat_star * m_j3_raw mod t_mod
    v_j3_decoded = mul_matrix_vector(W_hat_star, input, t)

    return v_j3_decoded

# ---------------------------------------------------------------------------------------
def contains_vector(vectors, target):
    """
    Verifica se a lista de vetores `vectors` contém o vetor `target` como elemento.

    Argumentos:
        vectors (list of list): lista de vetores.
        target (list): vetor a buscar.

    Retorna:
        bool: True se `target` está em `vectors`, False caso contrário.

    Exemplo:
        vectors = [[1,2], [3,4]]
        contains_vector(vectors, [3,4])  # True
    """
    return target in vectors

# -----------------------------------------------------------------------------------------
def contains_message(vet, msg):    

    n = len(msg)
    for ct in range(n):
        if msg[ct] in vet:
            continue
        else:
            return False
        
    return True
