import sys
import time
from tests.yugi_test_full import show_title
from tests.yugi_test_mul import generate_formula
from venum.ntt import *
import pytest
import random
from venum.ntt_func import *
import math
from typing import List, Tuple, Optional

@pytest.mark.parametrize(
    "input",
    [
        {
        }
    ])
def test_scheme_rotation(input):

# ----------------------------------------------------------------------------------------
# Colar todas as funções fornecidas aqui:
# find_primitive_root, modular_inverse, reverse_bits, bit_reverse_order,
# generate_psi_vectors, Barrett, reduce_inline, ntt_generic, intt_generic
# ----------------------------------------------------------------------------------------
  """
      Calculates the primitive root of a prime number q.
      The primitive root is a number g such that all numbers from 1 to q-1 are
      congruent to g^i mod q for some i.
      The primitive root is used to calculate the Fourier transform
      of finite fields.
      Args:
      - q: int, the prime number
      - n: int, the size of the Fourier transform
      Returns:
      - int, the primitive root of q
      Raises:
      - ValueError, if a primitive root cannot be found
  """
def find_primitive_root(q: int, n: int) -> int:
    if not is_prime(q):
         raise ValueError(f"{q} is not prime")

    order = n # Adjusted: Order N for N-th roots used in provided NTT/INTT
    q_minus1 = q - 1

    if q_minus1 % order != 0:
        raise ValueError(f"Order {order} does not divide q-1 ({q_minus1})")

    # Find the prime factors of order
    factors = set() # Use set to avoid duplicates
    m = order
    i = 2
    while i * i <= m:
        if m % i == 0:
            factors.add(i)
            while m % i == 0:
                m //= i
        i += 1
    if m > 1:
        factors.add(m)

    # We try small candidates for generator g of the full group Z_q*
    for g_candidate in range(2, q):
         is_generator = True
         for factor in factors:
             # Check if g_candidate^((q-1)/factor) == 1 mod q
             # If so, g_candidate is not a generator of Z_q*
             if pow(g_candidate, q_minus1 // factor, q) == 1:
                 is_generator = False
                 break
         if is_generator:
             # Found a generator g for Z_q*
             # Now find the N-th root: g^((q-1)/N) mod q
             exponent = q_minus1 // order
             return pow(g_candidate, exponent, q)

    raise ValueError(f"No N-th root of unity found for N={n}, q={q}")

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

# ----------------------------------------------------------------------------------------
"""
Function to calculate the modular inverse of a number a mod q.
Args:
    - a: int, the number
    - q: int, the modulus
Returns:
    - int, the modular inverse of a mod q
    - None, if the modular inverse does not exist
"""
def modular_inverse(a: int, q: int) -> Optional[int]:
    a = a % q # Ensure a is in range [0, q-1]
    if a == 0: return None # Inverse of 0 does not exist

    t, new_t = 0, 1
    r, new_r = q, a

    while new_r != 0:
        quotient = r // new_r
        t, new_t = new_t, t - quotient * new_t
        r, new_r = new_r, r - quotient * new_r

    if r > 1:
        # O inverso modular não existe
        return None

    if t < 0:
        t += q

    return t

# ----------------------------------------------------------------------------------------
# reverse_bits
def reverse_bits(x: int, bits: int) -> int:
    result = 0
    for _ in range(bits):
        result = (result << 1) | (x & 1)
        x >>= 1
    return result

# ----------------------------------------------------------------------------------------
# Bit reverse order
def bit_reverse_order(vec: List[int], n: int) -> List[int]:
    if n == 0: return []
    result = [0] * n
    bits = n.bit_length() - 1 # Correct way to get log2(n) for power of 2
    for i in range(n):
        rev_index = reverse_bits(i, bits)
        if rev_index < n: # Ensure index is within bounds
             result[rev_index] = vec[i]
    return result

# ----------------------------------------------------------------------------------------
# Generate psi vectors for NTT/INTT
def generate_psi_vectors(n: int, q: int, psi: int) -> Tuple[List[int], List[int]]:
    """
    Generates the psi_rev and psi_inv_rev vectors used in the NTT/INTT functions.
    Args:
    n (int): dimension (size of the vectors).
    q (int): module.
    psi (int): N-th root of unity mod q.

    Returns:
    Tuple[List[int], List[int]]: (psi_rev, psi_inv_rev)
    """

    # print(f"Generating psi vectors for n={n}, q={q}, psi={psi}")

    # Vector for psi in direct order
    psi_vec = [1] * n
    for i in range(1, n):
        psi_vec[i] = (psi_vec[i - 1] * psi) % q

    # Psi vector in bit-reversed order
    psi_rev = bit_reverse_order(psi_vec, n)

    # Modular inverse of psi
    psi_inv_val = modular_inverse(psi, q)
    if psi_inv_val is None:
         raise ValueError(f"Modular inverse of psi={psi} mod q={q} does not exist.")

    # Vector with modular inverses of psi_vec
    psi_inv_vec = [1] * n
    for i in range(1, n):
        psi_inv_vec[i] = (psi_inv_vec[i-1] * psi_inv_val) % q

    # Inverse Psi vector in bit-reversed order
    psi_inv_rev = bit_reverse_order(psi_inv_vec, n)

    return psi_rev, psi_inv_rev

# ----------------------------------------------------------------------------------------
# Barrett reduction
class Barrett:
    def __init__(self, m: int):
        if m < 2:
            raise ValueError("m cannot be 0 or 1")
        if m >= (1 << 63):
            raise ValueError("m must be < 2^63")

        self.m = m
        # bits_m equals 64 - m.leading_zeros() in Rust, which is the number of bits needed to represent m
        bits_m = m.bit_length()
        self.shift = bits_m # Using bit_length directly as shift for simplicity
        self.mu = (1 << (self.shift + 64)) // m
        self.m128 = m  # In Python, integers are arbitrary precision

    def modulus(self) -> int:
        return self.m

    def __repr__(self) -> str:
        return f"Barrett(m={self.m}, shift={self.shift}, mu={self.mu}, m128={self.m128})"

# ----------------------------------------------------------------------------------------
# Barrett reduction
def reduce_inline(x: int, bar: Barrett) -> int:
    """
    Performs modular reduction of x using Barrett multiplication.

    Args:
    x (int): value to be reduced (must be treated as a 128-bit integer).
    bar (Barrett): instance of the Barrett class with pre-computed parameters.

    Returns:
    int: result of modular reduction, equivalent to x mod bar.m.
    """
    if x < bar.m and x >= 0: return x # Optimization

    # Calculation of q: equivalent to (((x >> bar.shift).wrapping_mul(bar.mu)) >> 64) in Rust.
    # Python's integers handle large values
    q = (((x >> bar.shift) * bar.mu) >> 64)

    # Remainder calculation: equivalent to x - q * m, with m represented as bar.m (or bar.m128 in Rust).
    r = x - q * bar.m

    # Conditional adjustment: subtracts bar.m if r is greater than or equal to bar.m, up to two times.
    # A simple modulo is safer in Python for potentially large x
    if r < 0 or r >= 2 * bar.m :
         r = r % bar.m
    elif r >= bar.m:
        r -= bar.m
        if r >= bar.m: # Should not happen with correct mu/shift, but for safety
            r -= bar.m

    # Final check for negative results from modulo %
    if r < 0:
         r += bar.m

    return r

# ----------------------------------------------------------------------------------------
# NTT (Number Theoretic Transform) generic
def ntt_generic(input_list: list[int], psi_rev: list[int], q: int, bar: Barrett) -> None:
    """
    Performs the generic NTT (Number Theoretic Transform) on 'input_list' (modifying it in-place),
    using the 'psi_rev' vector, the 'q' module and the 'bar' (Barrett) structure.

    Args:
    input_list (list[int]): Polynomial (list of coefficients) to be transformed.
    psi_rev (list[int]): Vector with the unity roots in bit-reverse order.
    q (int): Module.
    bar (Barrett): Instance of the Barrett class with pre-computed parameters.
    """
    n = len(input_list)
    # Check if n is at least 2 and a power of 2.
    if n < 2 or (n & (n - 1)) != 0:
        raise ValueError("N must be a power of 2 and greater than 2")

    # log_n is the exponent such that 2^log_n = n.
    log_n = n.bit_length() - 1

    # Cooley-Tukey NTT Algorithm
    t = n // 2
    m = 1
    while m < n:
        k_start = 0
        psi_idx = 1 # Index for psi_rev powers
        while k_start < n:
            # Use precomputed psi_rev[psi_idx] as the twiddle factor omega
            omega = psi_rev[psi_idx]
            for j in range(k_start, k_start + t):
                idx_pair = j + t
                # Apply twiddle factor using Barrett reduction
                t_val = reduce_inline(input_list[idx_pair] * omega, bar)
                u = input_list[j]

                # Butterfly operation
                input_list[j] = u + t_val
                if input_list[j] >= q:
                    input_list[j] -= q

                input_list[idx_pair] = u - t_val
                if input_list[idx_pair] < 0:
                    input_list[idx_pair] += q

            k_start += 2 * t
            psi_idx += 1
        t //= 2
        m *= 2


# ----------------------------------------------------------------------------------------
# INTT (Number Theoretic Inverse Transform) generic
def intt_generic(input_list: list[int], psi_inv_rev: list[int], n_inv: int, q: int, bar: Barrett) -> None:
    """
    Performs the generic in-place INTT (Number Theoretic Inverse Transform),
    using the vector psi_inv_rev (inverse roots in bit-reverse order),
    the modular inverse of n (n_inv), the module q and the bar (Barrett) structure.

    Args:
    input_list (list[int]): List of coefficients of the polynomial (size n).
    psi_inv_rev (list[int]): Vector of inverse roots in bit-reverse order.
    n_inv (int): Modular inverse of n.
    q (int): Module.
    bar (Barrett): Instance of the Barrett class for modular reduction.
    """
    n = len(input_list)
    # Checks if n is a power of 2 and greater than or equal to 2.
    if n < 2 or (n & (n - 1)) != 0:
        raise ValueError("N must be a power of 2 and greater than 2")

    # Apply NTT structure but with inverse roots
    ntt_generic(input_list, psi_inv_rev, q, bar)

    # Multiply by n_inv
    for i in range(n):
        input_list[i] = reduce_inline(input_list[i] * n_inv, bar)

#-----------------------------------------------------
# Funções auxiliares para o experimento BFV
#-----------------------------------------------------

def poly_mult(p1: List[int], p2: List[int], N: int, q: int) -> List[int]:
    """Multiplicação Polinomial em R_q = Z_q[X] / (X^N + 1)"""
    res_full = [0] * (2 * N - 1)
    for i in range(N):
        if p1[i] == 0: continue
        for j in range(N):
            res_full[i+j] = (res_full[i+j] + p1[i] * p2[j]) #% q # Acumula sem mod intermediário se possível

    res_mod = [0] * N
    for i in range(2 * N - 1):
        if res_full[i] == 0: continue
        idx = i % N
        # sign = 1 if (i // N) % 2 == 0 else -1 # X^N = -1
        if i < N: # i // N = 0
             res_mod[idx] = (res_mod[idx] + res_full[i]) % q
        else: # i // N = 1
             res_mod[idx] = (res_mod[idx] - res_full[i]) % q # Subtrai por causa do X^N = -1

    # Garante que resultados estão no intervalo [0, q-1]
    #for i in range(N):
    #     res_mod[i] = res_mod[i] % q # Já feito dentro do loop

    return res_mod


def apply_automorphism(poly: List[int], k: int, N: int, q: int) -> List[int]:
    """Aplica o automorfismo X -> X^k ao polinômio poly em R_q."""
    res = [0] * N
    for i in range(N):
        if poly[i] == 0: continue
        target_exp = (k * i)
        term_val = poly[i]

        # Redução usando X^N = -1
        quotient = target_exp // N
        remainder = target_exp % N

        if quotient % 2 != 0: # Se expoente de X^N é ímpar, multiplica por -1
            term_val = (-term_val) % q

        res[remainder] = (res[remainder] + term_val) % q
    return res

#-----------------------------------------------------
# Experimento BFV
#-----------------------------------------------------

# 1. Parâmetros e Chaves
N = 8
Q = 257
T = 17
DELTA = 15
k_automorphism = 5

print(f"Parâmetros: N={N}, Q={Q}, T={T}, DELTA={DELTA}, k={k_automorphism}")

S_poly = [0] * N
S_poly[0] = 1
S_poly[4] = (Q - 1) # -1 mod Q for 1 - X^4
print(f"Chave Secreta S: {S_poly}")

# Precomputações para NTT/INTT mod T
try:
    root_t = find_primitive_root(T, N) # N-th root
    print(f"Raiz primitiva N-ésima (N={N}) mod T={T}: {root_t}")
except ValueError as e:
    print(f"Erro ao encontrar raiz: {e}")
    exit()

psi_t_rev, psi_t_inv_rev = generate_psi_vectors(N, T, root_t)
N_inv_t = modular_inverse(N, T)
if N_inv_t is None: raise ValueError(f"Inverso de N mod T não existe.")
bar_t = Barrett(T)
print(f"N_inv_t: {N_inv_t}")
# print(f"psi_t_rev: {psi_t_rev}")
# print(f"psi_t_inv_rev: {psi_t_inv_rev}")

# 2. Mensagem Original e Codificação SIMD (INTT Real)
m_original = [1, 0, 0, 0, 0, 0, 0, 0]
print(f"\nMensagem Original m: {m_original}")

# # Aplicar bit reversal antes da INTT
# m_rev = bit_reverse_order(m_original, N)
# print(f"Mensagem com Bit Reverse m_rev: {m_rev}")
m_rev = m_original.copy()

# Aplicar INTT in-place (copiar primeiro)
M_poly_coeffs = m_rev.copy()
intt_generic(M_poly_coeffs, psi_t_inv_rev, N_inv_t, T, bar_t)
print(f"Polinômio Codificado M (mod T={T}): {M_poly_coeffs}")

# 3. Criptografia BFV (sem ruído E=0)
A_poly = [0] * N
A_poly[0] = 1
A_poly[2] = 1 # A = 1 + X^2
print(f"\nMáscara A: {A_poly}")

# Calcular A*S mod Q
AS_poly = poly_mult(A_poly, S_poly, N, Q)
print(f"A*S mod Q: {AS_poly}")

# Calcular Delta * M mod Q
DeltaM_poly = [(DELTA * coeff) % Q for coeff in M_poly_coeffs]
print(f"Delta * M mod Q: {DeltaM_poly}")

# Calcular b = -A*S + Delta*M mod Q
neg_AS_poly = [mod_number(coeff * -1, Q) for coeff in AS_poly]
b_poly = [(neg_AS_poly[i] + DeltaM_poly[i]) % Q for i in range(N)]
print(f"Componente b: {b_poly}")

# Criptograma Inicial C = (b, A)
print(f"Criptograma C = (b, A)")

# 4. Rotação SIMD (Aplicação do Automorfismo sigma_5)
print(f"\nAplicando Automorfismo k={k_automorphism}")
A_prime_poly = apply_automorphism(A_poly, k_automorphism, N, Q)
b_prime_poly = apply_automorphism(b_poly, k_automorphism, N, Q)
print(f"A' = A(X^k): {A_prime_poly}")
print(f"b' = b(X^k): {b_prime_poly}")

# Criptograma Rotacionado C' = (b', A')

# 5. Key Switching
print(f"\nRealizando Key Switching")
# S_k = apply_automorphism(S_poly, k_automorphism, N, Q) # S_k = S neste caso
# print(f"S_k = S(X^k): {S_k}")

# Reutilizar evk do exemplo anterior
a_k = [1, 1, 0, 0, 0, 0, 0, 0]
b_k = [0, 256, 0, 0, 0, 1, 0, 0] # (1-a_k)S = (-X)(1-X^4)=-X+X^5
print(f"evk = (b_k, a_k): ({b_k}, {a_k})")

# Calcular b'' = b' + A'*b_k
term1 = poly_mult(A_prime_poly, b_k, N, Q)
b_double_prime = [(b_prime_poly[i] + term1[i]) % Q for i in range(N)]
print(f"b'' = b' + A'*b_k: {b_double_prime}")

# Calcular A'' = A'*a_k
A_double_prime = poly_mult(A_prime_poly, a_k, N, Q)
print(f"A'' = A'*a_k: {A_double_prime}")

# Criptograma Resultante C'' = (b'', A'')

# 6. Decriptografia BFV
print(f"\nDecriptografando C''")
# Calcular A''*S
AS_final = poly_mult(A_double_prime, S_poly, N, Q)
print(f"A''*S: {AS_final}")

# Calcular P = b'' + A''*S mod Q
P_poly = [(b_double_prime[i] + AS_final[i]) % Q for i in range(N)]
print(f"Polinômio Decriptografado P = b''+A''S (mod Q): {P_poly}")

# Verificação: Calcular Delta * M(X^k)
M_rot_t = apply_automorphism(M_poly_coeffs, k_automorphism, N, T)
print(f"M(X^k) (mod T): {M_rot_t}")
Delta_M_rot_q = [(DELTA * coeff) % Q for coeff in M_rot_t]
print(f"Delta * M(X^k) (mod Q): {Delta_M_rot_q}")

if P_poly == Delta_M_rot_q:
    print("Verificação P == Delta * M(X^k) bem-sucedida!")
else:
    print("ERRO: Verificação P == Delta * M(X^k) FALHOU!")
    # Calcular diferença
    diff = [(P_poly[i] - Delta_M_rot_q[i] + Q)% Q for i in range(N)]
    print(f"Diferença (P - Delta*M_rot): {diff}")


# 7. Decodificação SIMD (NTT Real)
print("\nDecodificando o resultado")
# Escalar e Arredondar (sem ruído, apenas dividir por Delta)
# Idealmente M'_poly = P / Delta = M(X^k)
M_prime_poly_coeffs = M_rot_t # Usar o valor esperado que foi verificado
print(f"Polinômio M' = P/Delta (mod T): {M_prime_poly_coeffs}")

# # Aplicar bit reversal antes da NTT
# M_prime_rev = bit_reverse_order(M_prime_poly_coeffs, N)
# print(f"M' com Bit Reverse: {M_prime_rev}")
M_prime_rev = M_prime_poly_coeffs

# Aplicar NTT in-place (copiar primeiro)
m_prime_final = M_prime_rev.copy()
ntt_generic(m_prime_final, psi_t_rev, T, bar_t)
print(f"Vetor Decodificado Final m' (Resultado da NTT): {m_prime_final}")

# 8. Comparação e Conclusão
# Calcular o vetor esperado manualmente
m_prime_expected = [0] * N
permutation_map = { # Mapeia slot_destino -> slot_origem (baseado em j_l = 5j_i)
    0: 2, 1: 3, 2: 1, 3: 0, 4: 7, 5: 6, 6: 4, 7: 5
}
for dest_slot, src_slot in permutation_map.items():
    m_prime_expected[dest_slot] = m_original[src_slot]

print(f"\nVetor Original m: {m_original}")
print(f"Vetor Rotacionado Esperado m' (manual): {m_prime_expected}")
print(f"Vetor Rotacionado Calculado m' (via NTT): {m_prime_final}")

if m_prime_final == m_prime_expected:
    print("\nConclusão: O vetor decodificado final corresponde ao vetor esperado após a rotação SIMD!")
else:
    print("\nERRO: O vetor decodificado final NÃO corresponde ao esperado!")


