import sys
import time
import random
import pytest
from math import ceil, log, gcd

# --- Supondo que estas funções/classes existem ---
from venum.ntt import *
from venum.ntt_func import *
from tests.yugi_test_full import show_title
# from your_bfv_implementation import (
#     Cryptogram, create_sk, generate_pk, generate_rlk,
#     encrypt_sk, decrypt, generate_parameters,
#     rotate_polynomial_coeffs, view_sk_rot, mod_number,
#     generate_modulus, modInverse # modInverse precisa ser definida
# )
# --- Fim das suposições ---

# Função para calcular o inverso modular (se não existir na sua lib)
def modInverse(a, m):
    m0 = m
    y = 0
    x = 1
    if (m == 1): return 0
    while (a > 1):
        q = a // m
        t = m
        m = a % m; a = t
        t = y
        y = x - q * y
        x = t
    if (x < 0): x = x + m0
    return x if gcd(a, m0) == 1 else None # Retorna None se não houver inverso

@pytest.mark.parametrize("input", [{}])
def test_find_cyclic_rotation(input):
    print("\n")
    show_title("TEST SCHEME - FINDING CYCLIC ROTATION")

    n = 8
    p2 = 3
    base_decomposition = 2
    batched = True

    # Use módulos menores para testes mais rápidos se possível, mas mantenha as condições
    # q = generate_modulus(2**60, 2**63, n)
    q = generate_modulus(2**30, 2**32, n) # Exemplo menor
    print(" Q = ", q)

    p1 = 0
    retry_count = 0
    while True:
        # p1 = generate_modulus(2**7, 2**20, n) # Original
        p1 = generate_modulus(2**10, 2**12, n)  # Exemplo menor
        if p1 % (2 * n) == 1 and p1 > n : # Garantir p1 > n também
            break
        retry_count += 1
        if retry_count > 100:
             raise ValueError("Não foi possível encontrar p1 adequado.")

    print("P1 = ", p1)

    # Gera parâmetros
    try:
        psi_rev_q, psi_inv_rev_q, n_inv_q, bar_q = generate_parameters(n, q)
        params_q = (psi_rev_q, psi_inv_rev_q, n_inv_q, bar_q)
       
        psi_rev_p1, psi_inv_rev_p1, n_inv_p1, bar_p1 = generate_parameters(n, p1)
        params_batched = (psi_rev_p1, psi_inv_rev_p1, n_inv_p1, bar_p1)
    except Exception as e:
        print(f"Erro ao gerar parâmetros para n={n}, q={q}, p1={p1}: {e}")
        # Isso pode acontecer se p1 for muito pequeno ou não tiver as propriedades certas
        # para a NTT funcionar na biblioteca.
        pytest.skip("Falha na geração de parâmetros, pulando teste.")
        return


    sk = create_sk(n, 0, 1, q)
    # pk = generate_pk(sk, q, p1, p2, params_q) # Não estritamente necessário para este teste
    # rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params_q) # Idem

    m0 = [i + 1 for i in range(n)]
    print("Plaintext: ", m0)

    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params_q, params_batched)

    decrypted_orig = decrypt(sk, c0, p1, params_q, params_batched)
    print("Decrypted (c0): ", decrypted_orig)
    if decrypted_orig != m0:
        print("ERRO: Decriptação inicial falhou!")
        sys.exit(1)
    print("-" * 80)

    # Grupo de Galois para N=8 (2N=16)
    galois_group = [d for d in range(1, 2*n) if gcd(d, 2*n) == 1]
    print(f"(Z/{2*n}Z)^* = {galois_group}")
    print("-" * 80)

    found_left_shift = None
    found_right_shift = None

    expected_left = m0[1:] + m0[:1]
    expected_right = m0[-1:] + m0[:-1]
    print(f"Procurando por Left Shift : {expected_left}")
    print(f"Procurando por Right Shift: {expected_right}")
    print("-" * 80)


    for d in galois_group:
        if d == 1: continue # d=1 é a identidade, não rotaciona

        print(f"Testando automorfismo d = {d}")

        # Tenta rotacionar e descriptografar
        try:
            sk_rot = rotate_polynomial_coeffs(sk, d, n, q)
            body_rot = rotate_polynomial_coeffs(c0.body, d, n, q)
            mask_rot = rotate_polynomial_coeffs(c0.mask, d, n, q)
            c1 = Cryptogram(body=body_rot, mask=mask_rot, batched=c0.batched, q=c0.q)
            decrypted = decrypt(sk_rot, c1, p1, params_q, params_batched)
            print(f"  Resultado Decriptado: {decrypted}")

            if decrypted == expected_left:
                print(f"  ✅ ENCONTRADO: d={d} realiza o SHIFT LEFT por 1!")
                found_left_shift = d
            elif decrypted == expected_right:
                print(f"  ✅ ENCONTRADO: d={d} realiza o SHIFT RIGHT por 1!")
                found_right_shift = d
            else:
                 # Verificar se é alguma outra permutação interessante (opcional)
                 pass

        except Exception as e:
            print(f"  ⚠️ Erro durante a rotação/decriptação com d={d}: {e}")
            # Pode acontecer se algo der errado internamente

        print("-" * 60) # Separador para cada d testado

    print("=" * 80)
    print("Resultado da Busca:")
    if found_left_shift:
        print(f"Automorfismo para Shift Left por 1: d = {found_left_shift}")
        # Verifica se o inverso dele dá o shift right
        inv_d = modInverse(found_left_shift, 2*n)
        if inv_d == found_right_shift:
            print(f"  (Correto: {found_left_shift}^-1 mod {2*n} = {inv_d}, que foi encontrado para Shift Right)")
        elif found_right_shift:
             print(f"  (Inesperado: {found_left_shift}^-1 mod {2*n} = {inv_d}, mas {found_right_shift} foi encontrado para Shift Right)")
        else:
             print(f"  (Info: {found_left_shift}^-1 mod {2*n} = {inv_d}, mas Shift Right não foi encontrado)")

    else:
        print("Nenhum automorfismo encontrado para Shift Left por 1.")

    if found_right_shift:
        print(f"Automorfismo para Shift Right por 1: d = {found_right_shift}")
        # Verifica se o inverso dele dá o shift left
        inv_d = modInverse(found_right_shift, 2*n)
        if inv_d == found_left_shift:
             print(f"  (Correto: {found_right_shift}^-1 mod {2*n} = {inv_d}, que foi encontrado para Shift Left)")
        elif found_left_shift:
             print(f"  (Inesperado: {found_right_shift}^-1 mod {2*n} = {inv_d}, mas {found_left_shift} foi encontrado para Shift Left)")
        else:
            print(f"  (Info: {found_right_shift}^-1 mod {2*n} = {inv_d}, mas Shift Left não foi encontrado)")
    else:
        print("Nenhum automorfismo encontrado para Shift Right por 1.")

    print("=" * 80)

    # Se encontramos um, podemos gerar shifts maiores
    if found_left_shift:
        k=2
        d_k = pow(found_left_shift, k, 2*n)
        print(f"\nTestando Shift Left por {k} usando d = {found_left_shift}^{k} = {d_k} mod {2*n}")
        try:
            sk_rot_k = rotate_polynomial_coeffs(sk, d_k, n, q)
            body_rot_k = rotate_polynomial_coeffs(c0.body, d_k, n, q)
            mask_rot_k = rotate_polynomial_coeffs(c0.mask, d_k, n, q)
            c1_k = Cryptogram(body=body_rot_k, mask=mask_rot_k, batched=c0.batched, q=c0.q)
            decrypted_k = decrypt(sk_rot_k, c1_k, p1, params_q, params_batched)
            expected_k = m0[k:] + m0[:k]
            print(f"  Resultado: {decrypted_k}")
            print(f"  Esperado : {expected_k}")
            if decrypted_k == expected_k:
                print(f"  ✅ Shift Left por {k} funcionou!")
            else:
                print(f"  ❌ Shift Left por {k} falhou.")
        except Exception as e:
            print(f"  ⚠️ Erro durante a rotação/decriptação com d={d_k}: {e}")
            
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
