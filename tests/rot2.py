import sys
import time
import random
import pytest
from math import ceil, log, gcd # Adicionado gcd

# ... (resto das importações e funções auxiliares como generate_modulus, etc.)
from venum.ntt import *
from venum.ntt_func import *
from tests.yugi_test_full import show_title
# Assume que Cryptogram, create_sk, generate_pk, generate_rlk,
# encrypt_sk, decrypt, generate_parameters estão definidas em algum lugar.
# Assume que rotate_polynomial_coeffs, view_sk_rot, mod_number estão definidas.

# Função para calcular o inverso modular
def modInverse(a, m):
    m0 = m
    y = 0
    x = 1
    if (m == 1):
        return 0
    while (a > 1):
        # q é o quociente
        q = a // m
        t = m
        # m é o resto agora, processa da mesma forma que o gcd de Euclides
        m = a % m
        a = t
        t = y
        # Atualiza y e x
        y = x - q * y
        x = t
    # Garante que x seja positivo
    if (x < 0):
        x = x + m0
    return x

@pytest.mark.parametrize("input", [{}])
def test_scheme_rotation_cyclic(input):
    print("\n")
    show_title("TEST SCHEME - CYCLIC ROTATION")

    # parameters:
    n = 8 # Testando com N=8
    p2 = 3
    base_decomposition = 2
    batched = True

    q = generate_modulus(2**60,2**63,n)
    print(" Q = ", q)

    while True:
        p1 = generate_modulus(2**7,2**20,n)
        # Condição necessária para batching com x^n+1
        if p1 % (2*n) == 1:
            break

    print("P1 = ", p1)

    # Gera parâmetros para q e p1
    psi_rev_q, psi_inv_rev_q, n_inv_q, bar_q = generate_parameters(n, q)
    params_q = (psi_rev_q, psi_inv_rev_q, n_inv_q, bar_q)

    psi_rev_p1, psi_inv_rev_p1, n_inv_p1, bar_p1 = generate_parameters(n, p1)
    params_batched = (psi_rev_p1, psi_inv_rev_p1, n_inv_p1, bar_p1)

    # Gera chaves
    sk = create_sk(n, 0, 1, q)
    print("Secret key: ", sk[:8], "...") # Mostrar só uma parte

    pk = generate_pk(sk, q, p1, p2, params_q) # Usar params_q aqui
    # Rotações não precisam de rlk, mas vamos manter a geração caso precise depois
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params_q) # Usar params_q

    # Gera plaintext
    m0 = [i + 1 for i in range(n)]
    print("Plaintext: ", m0)

    # Encripta
    # Passar os parâmetros corretos é crucial
    c0 = encrypt_sk(sk, m0, q, p1, p2, batched, params_q, params_batched)
    #print("C0: ", c0)

    # Decripta para verificar
    decrypted_orig = decrypt(sk, c0, p1, params_q, params_batched) # Passar params_q e params_batched
    print("Decrypted (c0): ", decrypted_orig)
    print("-" * 80)
    
    assert decrypted_orig == m0

    if decrypted_orig != m0:
        print("ERRO: Decriptação inicial falhou!")
        # Verifique as funções generate_parameters, encrypt_sk, decrypt
        # Certifique-se que params_q e params_batched estão sendo usados corretamente
        sys.exit(1)

    # --- Testando Rotações Cíclicas ---
    print("Testando Rotação Cíclica Específica (Shift Left por 1)")

    # O automorfismo para shift left por 1 é geralmente d=3
    d_shift_left = 3
    if gcd(d_shift_left, 2*n) != 1:
         print(f"Erro: d={d_shift_left} não está em (Z/{2*n}Z)^*")
         # Tentar d = 2*n - 1 talvez? Ou verificar a teoria para N específico
         # Para N potência de 2, d=3 quase sempre funciona para um shift.
         d_shift_left = (2*n) - 1 # Tentar o outro candidato comum
         if gcd(d_shift_left, 2*n) != 1:
             print(f"Erro: Nem d=3 nem d={2*n-1} são válidos.")
             sys.exit(1)
         else:
             print(f"Usando d = {d_shift_left} para rotação.")


    # Calcula a chave secreta rotacionada
    sk_rot_left = rotate_polynomial_coeffs(sk, d_shift_left, n, q)

    # Rotaciona o ciphertext homomorficamente
    body_rot_left = rotate_polynomial_coeffs(c0.body, d_shift_left, n, q)
    mask_rot_left = rotate_polynomial_coeffs(c0.mask, d_shift_left, n, q)
    c1_left = Cryptogram(body=body_rot_left, mask=mask_rot_left, batched=c0.batched, q=c0.q)

    # Decripta com a chave rotacionada
    decrypted_left = decrypt(sk_rot_left, c1_left, p1, params_q, params_batched) # Passar params_q e params_batched
    print(f"Automorfismo d = {d_shift_left}")
    print("SK Rotacionada:", view_sk_rot(sk_rot_left[:8], q), "...")
    print("Decriptado Rotacionado (Left): ", decrypted_left)

    # Verifica se foi o shift esquerdo esperado
    expected_left = m0[1:] + m0[:1]
    print("Esperado (Left):             ", expected_left)
    if decrypted_left == expected_left:
        print("✅ Shift Left por 1 funcionou!")
    else:
        # Talvez d=3 cause shift right?
        expected_right = m0[-1:] + m0[:-1]
        if decrypted_left == expected_right:
            print("⚠️ Automorfismo d=3 resultou em Shift Right por 1!")
        else:
            print("❌ Shift Left por 1 falhou. Resultado inesperado.")

    print("-" * 80)

    # --- Testando Rotação Cíclica para a Direita ---
    print("Testando Rotação Cíclica Específica (Shift Right por 1)")

    # O automorfismo para shift right por 1 é geralmente d = 3^{-1} mod 2N
    d_shift_right = modInverse(d_shift_left, 2*n) # Calcula o inverso do d usado para left
                                                 # Se left foi d=3, right será d=11 (para N=8)
                                                 # Se left foi d=15, right será d=15 (pois 15*15 = 225 = 14*16 + 1 = 1 mod 16)

    if d_shift_right is None or gcd(d_shift_right, 2*n) != 1:
         print(f"Erro: Inverso d={d_shift_right} não é válido ou não existe.")
         sys.exit(1)

    # Calcula a chave secreta rotacionada
    sk_rot_right = rotate_polynomial_coeffs(sk, d_shift_right, n, q)

    # Rotaciona o ciphertext homomorficamente
    body_rot_right = rotate_polynomial_coeffs(c0.body, d_shift_right, n, q)
    mask_rot_right = rotate_polynomial_coeffs(c0.mask, d_shift_right, n, q)
    c1_right = Cryptogram(body=body_rot_right, mask=mask_rot_right, batched=c0.batched, q=c0.q)

    # Decripta com a chave rotacionada
    decrypted_right = decrypt(sk_rot_right, c1_right, p1, params_q, params_batched) # Passar params_q e params_batched
    print(f"Automorfismo d = {d_shift_right}")
    print("SK Rotacionada:", view_sk_rot(sk_rot_right[:8], q), "...")
    print("Decriptado Rotacionado (Right):", decrypted_right)

    # Verifica se foi o shift direito esperado
    expected_right = m0[-1:] + m0[:-1]
    print("Esperado (Right):            ", expected_right)
    if decrypted_right == expected_right:
        print("✅ Shift Right por 1 funcionou!")
    else:
         # Talvez o inverso cause shift left?
        expected_left = m0[1:] + m0[:1]
        if decrypted_right == expected_left:
            print(f"⚠️ Automorfismo d={d_shift_right} resultou em Shift Left por 1!")
        else:
            print("❌ Shift Right por 1 falhou. Resultado inesperado.")

    print("-" * 80)

    # --- Shift por k posições ---
    k = 2 # Exemplo: shift left por 2
    print(f"Testando Rotação Cíclica Específica (Shift Left por {k})")
    # Aplica o automorfismo de shift left k vezes (calculando d^k mod 2N)
    d_shift_k = pow(d_shift_left, k, 2*n) # Calcula d^k mod 2N

    sk_rot_k = rotate_polynomial_coeffs(sk, d_shift_k, n, q)
    body_rot_k = rotate_polynomial_coeffs(c0.body, d_shift_k, n, q)
    mask_rot_k = rotate_polynomial_coeffs(c0.mask, d_shift_k, n, q)
    c1_k = Cryptogram(body=body_rot_k, mask=mask_rot_k, batched=c0.batched, q=c0.q)
    decrypted_k = decrypt(sk_rot_k, c1_k, p1, params_q, params_batched)

    print(f"Automorfismo d = {d_shift_left}^{k} = {d_shift_k} (mod {2*n})")
    print(f"Decriptado Rotacionado (Left por {k}):", decrypted_k)
    expected_k = m0[k:] + m0[:k]
    print(f"Esperado (Left por {k}):            ", expected_k)
    if decrypted_k == expected_k:
       print(f"✅ Shift Left por {k} funcionou!")
    else:
       # Verificar se foi shift right por k
       expected_k_right = m0[-k:] + m0[:-k]
       if decrypted_k == expected_k_right:
           print(f"⚠️ Automorfismo d={d_shift_k} resultou em Shift Right por {k}!")
       else:
            print(f"❌ Shift Left por {k} falhou.")


# --- Funções Auxiliares (Garantir que existam e estejam corretas) ---
# Classe Cryptogram, create_sk, generate_pk, generate_rlk, encrypt_sk,
# decrypt, generate_parameters, rotate_polynomial_coeffs, view_sk_rot,
# mod_number, generate_modulus, etc.
# ... (O resto do seu código de teste, se houver)

# Certifique-se que as funções auxiliares (como generate_parameters, etc.)
# estão disponíveis e corretas. A implementação exata delas pode afetar o resultado.


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
