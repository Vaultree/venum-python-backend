import sys
import time
import random
import pytest
from math import ceil, log

from tests.yugi_test_full import show_title
from tests.yugi_test_mul import generate_formula
from venum.ntt import *
from venum.ntt_func import *

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
def random_vector(n, t):
    """
    Gera um vetor aleatório com `n` elementos em módulo `t`.
    Cada elemento será um inteiro no intervalo [0, t-1].

    Argumentos:
        n (int): tamanho do vetor.
        t (int): módulo (valor máximo é t-1).

    Retorna:
        list: vetor de inteiros aleatórios.
    """
    return [random.randrange(t) for _ in range(n)]
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("input", [{}])
def test_scheme_rotation(input):
    print("\n")
    show_title("TEST SCHEME - ROTATION")

# --- Recriação do cálculo da página 149 ---
for ct in range(100):
    
    n_dim = 8
    t_mod = 17

    # Vetor de entrada v
    v_input = [1, 2, 3, 4, 5, 6, 7, 8]
    v_input = random_vector(n_dim, t_mod)
    print("Input: ", v_input)   

    # Calcular n^-1 mod t
    n_inv = modular_inverse(n_dim, t_mod)
    print(f"N = {n_dim}, t = {t_mod}, N^-1 mod t = {n_inv}")

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

    print("-" * 60)
    print("Executando cálculo: m = n^-1 * W_hat * IR * v mod t")

    # Passo 1: Calcular M_prod = W_hat * IR mod t
    M_prod = multiply_matrices(W_hat, IR, t_mod)
    print(f"W_hat * IR (primeiras 3x3): {[row[:3] for row in M_prod][:3]}...")

    # Passo 2: Calcular v_temp = M_prod * v mod t
    v_temp = mul_matrix_vector(M_prod, v_input, t_mod)
    print(f"M_prod * v: {v_temp}")

    # Passo 3: Calcular m = n_inv * v_temp mod t
    m_encoded = mult_scalar_vector(n_inv, v_temp, t_mod)
    print(f"n^-1 * v_temp = m: {m_encoded}")

    print("ENCODE: ", m_encoded)
    print("=" * 50)

    # # Resultado esperado da página 149
    # expected_m = [13, 16, 10, 5, 9, 12, 7, 1]

    # print("-" * 60)
    # print(f"Resultado Calculado: {m_encoded}")
    # print(f"Resultado Esperado : {expected_m}")

    # assert m_encoded == expected_m, "O resultado calculado é diferente do esperado!"
    # print("✅ O resultado calculado corresponde ao esperado na página 149.")

    # ... (código anterior incluindo funcoes mod_number, mul_matrix_vector, multiply_matrices, mult_scalar_vector)

    print("\n" + "=" * 60)
    print("Executando cálculo de Decodificação Exemplo (Página 150):")
    print("v_j3 = W_hat_star * m_j3 mod t")
    print("=" * 60)

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

    # Vetor m_j3 apos divisao por Delta (cru, antes do mod 17, para input na funcao mul_matrix_vector)
    # Estes valores são dados como resultado da divisão (26, 24, -20, etc) / 2 no paper.
    # m_j3_raw = [13, 12, -10, -2, 18, -32, -14, 10]
    m_j3_raw = m_encoded

    # Modulo t
    t_mod = 17

    # Calcular v_j3 = W_hat_star * m_j3_raw mod t_mod
    v_j3_decoded = mul_matrix_vector(W_hat_star, m_j3_raw, t_mod)

    print(f"Vetor m_j3_raw (antes do mod no calculo): {m_j3_raw}")
    print(f"Resultado Calculado v_j3 : {v_j3_decoded}")

    # Resultado esperado (vetor rotacionado vr) da página 150
    # expected_v_j3 = [4, 1, 2, 3, 8, 5, 6, 7]

    expected_v_j3 = v_input

    print(f"Resultado Esperado v_j3  : {expected_v_j3}")

    assert v_j3_decoded == expected_v_j3, "O resultado decodificado calculado é diferente do esperado!"
    print("✅ O resultado da decodificação calculada corresponde ao esperado na página 150.")
    print("=" * 50)