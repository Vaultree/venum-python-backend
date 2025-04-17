from typing import Optional
import math
from typing import List, Tuple

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
    order = 2 * n
    q_minus1 = q - 1

    # Find the prime factors of order
    factors = []
    m = order
    i = 2
    while i * i <= m:
        if m % i == 0:
            factors.append(i)
            while m % i == 0:
                m //= i
        i += 1
    if m > 1:
        factors.append(m)

    # We try to use g = 3 as a candidate
    g = 3
    is_primitive_root = True
    for factor in factors:
        exponent = q_minus1 // factor
        if pow(g, exponent, q) == 1:
            is_primitive_root = False
            break

    if not is_primitive_root:
        # If g = 3 is not a primitive root, we test other small candidates
        for g_candidate in range(5, 1000):
            is_primitive_root = True
            for factor in factors:
                exponent = q_minus1 // factor
                if pow(g_candidate, exponent, q) == 1:
                    is_primitive_root = False
                    break
            if is_primitive_root:
                exp = q_minus1 // order  # exp = (q - 1) / (2n)
                return pow(g_candidate, exp, q)
        raise ValueError("No primitive roots found!")
    else:
        exp = q_minus1 // order  # exp = (q - 1) / (2n)
        return pow(g, exp, q)

# ----------------------------------------------------------------------------------------
def get_prime_factors(num):
    """Helper function to get distinct prime factors of num."""
    factors = set()
    d = 2
    temp_num = num
    while d * d <= temp_num:
        if temp_num % d == 0:
            factors.add(d)
            while temp_num % d == 0:
                temp_num //= d
        # Otimização: só testar 2 e depois ímpares
        d = d + 1 if d == 2 else d + 2
    if temp_num > 1:
        factors.add(temp_num)
    return list(factors)

# --------------------------------------------------------------------------------
def find_primitive_root2(q: int, n: int) -> int:
    """
    Finds a primitive 2n-th root of unity modulo q.
    Requires q to be prime and q = 1 (mod 2n).
    """
    if q % (2 * n) != 1:
        raise ValueError(f"q={q} must be 1 mod 2n={2*n}")

    order = 2 * n
    exponent = (q - 1) // order # Exponent to potentially lower the order from q-1 to 2n

    # Get distinct prime factors of the target order 2n
    prime_factors_of_order = get_prime_factors(order)

    # Iterate through candidates for the base g
    g = 2
    while g < q:
        # Calculate the candidate 2n-th root of unity
        omega = pow(g, exponent, q)

        # omega == 1 cannot be primitive if order > 1
        if omega == 1:
            g += 1
            continue

        # Check if omega has order *exactly* 2n
        is_primitive = True
        for p in prime_factors_of_order:
            test_exponent = order // p
            if pow(omega, test_exponent, q) == 1:
                is_primitive = False
                break # Failed the test for this prime factor p

        if is_primitive:
            return omega # Found a primitive 2n-th root

        g += 1 # Try the next base

    raise ValueError(f"No primitive {order}-th root of unity found for q={q}")
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
    result = [0] * n
    bits = int(math.log2(n))
    for i in range(n):
        rev_index = reverse_bits(i, bits)
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
    psi (int): basis for generating the elements.

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

    # Vector with modular inverses of psi_vec
    psi_inv = [modular_inverse(x, q) for x in psi_vec]
    psi_inv_rev = bit_reverse_order(psi_inv, n)

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
        self.shift = bits_m
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
    # Calculation of q: equivalent to (((x >> bar.shift).wrapping_mul(bar.mu)) >> 64) in Rust.
    q = (((x >> bar.shift) * bar.mu) >> 64)
    # Remainder calculation: equivalent to x - q * m, with m represented as bar.m (or bar.m128 in Rust).
    r = x - q * bar.m

    # Conditional adjustment: subtracts bar.m if r is greater than or equal to bar.m, up to two times.
    if r >= bar.m:
        r -= bar.m
        if r >= bar.m:
            r -= bar.m
            if r >= bar.m:
                r %= bar.m

    return r

# ----------------------------------------------------------------------------------------
# Multiplication of two polynomials modulo q and modulo (X^n + 1)
def multiply_poly_mod(a: list[int], b: list[int], q: int) -> list[int]:
    """
    Multiplies two polynomials modulo q and modulo (X^n + 1).
    Each polynomial is represented as a list of integers, where the index corresponds to the degree 
    of the monomial.
    Multiplication is done conventionally and, for terms of degree >= n, the sign is inverted
    (due to the congruence X^n ≡ -1). Then, each coefficient is reduced modulo q.

    Args:
    a (list[int]): Coefficients of the first polynomial.
    b (list[int]): Coefficients of the second polynomial.
    q (int): Module for reducing the coefficients.

    Returns:
    list[int]: List with the coefficients of the resulting polynomial.    
    """
    
    n = len(a)
    if len(b) != n:
        raise ValueError("Polynomials must be the same size.")
    
    result = [0] * n

    for i in range(n):
        a_i = a[i]
        for j in range(n):
            b_j = b[j]
            idx = (i + j) % n
            # If i+j >= n, it means that there was a "turn" in the polynomial, so we multiply it by -1.
            sign = -1 if (i + j) >= n else 1
            product = sign * a_i * b_j
            result[idx] = (result[idx] + product) % q

    return result

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

    # For each iteration i: m = 2^i and t = n / 2^(i+1)
    for i in range(log_n):
        m = 1 << i
        t = n // (1 << (i + 1))
        # Iterate over blocks of size 2*t
        for block_index, block_start in enumerate(range(0, n, 2 * t)):
            # Get the multiplicative factor for the current block
            s = psi_rev[m + block_index]
            # Apply "butterfly" to each pair in the block
            for k in range(t):
                idx_left = block_start + k
                idx_right = block_start + t + k
                u = input_list[idx_left]
                b_val = input_list[idx_right]
                # Multiply b_val by s and apply modular reduction via Barrett
                v = reduce_inline(b_val * s, bar)
                summ = u + v

                # Update the left half element with the modular sum
                input_list[idx_left] = summ - q if summ >= q else summ

                # Update the right half element with the modular difference
                if u >= v:
                    input_list[idx_right] = u - v
                else:
                    input_list[idx_right] = q - (v - u)

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
    
    # log_n is the exponent such that 2^log_n = n
    log_n = int(math.log2(n))
    
    # First part: iterations with pairs (m, t)
    # For i from 0 to log_n - 2 (limit = log_n - 1)
    limit = log_n - 1
    for i in range(limit):
        m = n >> (i + 1)   # Equivalent to n / 2^(i+1)
        t = 1 << i         # Equivalent to 2^i
        
        # Process blocks of 2*t elements.
        for j, block_start in enumerate(range(0, n, 2 * t)):
            # Get the psi factor for the current block.
            s = psi_inv_rev[m + j]
            # Process each pair within the block.
            for k in range(t):
                left_index = block_start + k
                right_index = block_start + t + k
                u = input_list[left_index]
                v = input_list[right_index]
                soma = u + v
                # Update the left half element
                input_list[left_index] = soma - q if soma >= q else soma
                # Update the right half element:
                # Compute (u + q - v) * s, and apply modular reduction.
                input_list[right_index] = reduce_inline((u + q - v) * s, bar)
    
    # Second part: last round (when m == 1)
    # Pre-compute s_n_inv to avoid repeated multiplications.
    s_n_inv = reduce_inline(psi_inv_rev[1] * n_inv, bar)
    half = n // 2
    # Separates the vector into two halves.
    for i in range(half):
        u = input_list[i]
        v = input_list[i + half]
        soma = u + v
        input_list[i] = reduce_inline(soma * n_inv, bar)
        input_list[i + half] = reduce_inline((u + q - v) * s_n_inv, bar)
        
# ----------------------------------------------------------------------------------------
# NTT (Number Theoretic Transform) - Cooley-Tukey Algorithm
# Espera entrada em ordem NATURAL, produz saída em ordem BIT-REVERSA.
def ntt_generic2(input_list: list[int], psi_rev: list[int], q: int, bar: Barrett) -> None:
#def ntt_cooley_tukey(input_list: list[int], psi_rev: list[int], q: int, bar: Barrett) -> None:
    """
    Performs NTT using the Cooley-Tukey algorithm (decimation-in-time).
    Modifies 'input_list' in-place.
    Assumes natural order input and produces bit-reversed order output.

    Args:
    input_list (list[int]): Polynomial coefficients (natural order).
    psi_rev (list[int]): Roots of unity in bit-reversed order.
    q (int): Modulus.
    bar (Barrett): Barrett reduction object for modulus q.
    """
    n = len(input_list)
    if n < 1 or (n & (n - 1)) != 0:
        raise ValueError("N must be a power of 2 and >= 1")
    if n == 1:
        return # NTT of size 1 is identity

    t = n
    m = 1
    while m < n:
        t //= 2
        if t == 0: # Proteção caso n não seja potência de 2 (já verificado, mas seguro)
             break
        for i in range(m):
            # w_m^i, mas indexado usando a ordem bit-reverse de psi_rev
            # O índice m+i em psi_rev corresponde à raiz necessária nesta etapa
            s = psi_rev[m + i]
            # Itera pelos pares dentro dos blocos
            for j in range(i, n, 2 * m): # Stride é 2*m, começando em i
                 k = j + m # Índice do par
                 # Butterfly:
                 # u = input_list[j]
                 # v = reduce_inline(input_list[k] * s, bar) # v = a[k] * w
                 # input_list[j] = (u + v) % q -> u + v if u+v < q else u+v-q
                 # input_list[k] = (u - v + q) % q -> u - v if u >= v else u - v + q

                 u = input_list[j]
                 # Calcula v = input_list[k] * s mod q usando Barrett
                 v = reduce_inline(input_list[k] * s, bar)

                 sum_val = u + v
                 input_list[j] = sum_val - q if sum_val >= q else sum_val

                 # diff_val = u - v mod q
                 if u >= v:
                     input_list[k] = u - v
                 else:
                     input_list[k] = u + q - v
        m *= 2

# ----------------------------------------------------------------------------------------
# INTT (Inverse Number Theoretic Transform) - Gentleman-Sande Algorithm
# Espera entrada em ordem BIT-REVERSA, produz saída em ordem NATURAL.
def intt_generic2(input_list: list[int], psi_inv_rev: list[int], n_inv: int, q: int, bar: Barrett) -> None:
# def intt_gentleman_sande(input_list: list[int], psi_inv_rev: list[int], n_inv: int, q: int, bar: Barrett) -> None:
    """
    Performs INTT using the Gentleman-Sande algorithm (decimation-in-frequency).
    Modifies 'input_list' in-place.
    Assumes bit-reversed order input and produces natural order output.

    Args:
    input_list (list[int]): Polynomial coefficients in NTT domain (bit-reversed order).
    psi_inv_rev (list[int]): Inverse roots of unity in bit-reversed order.
    n_inv (int): Modular inverse of n modulo q.
    q (int): Modulus.
    bar (Barrett): Barrett reduction object for modulus q.
    """
    n = len(input_list)
    if n < 1 or (n & (n - 1)) != 0:
        raise ValueError("N must be a power of 2 and >= 1")
    if n == 1:
        # Apenas aplica a escala n_inv para n=1
        if n_inv != 1: # Otimização: não multiplica se n_inv for 1
             input_list[0] = reduce_inline(input_list[0] * n_inv, bar)
        return

    t = 1
    m = n
    while m > 1:
        h = m // 2 # h = m/2
        # Itera pelos diferentes 's' necessários para esta etapa
        for i in range(h):
            # Raiz s = (psi^-1)^(bit_rev(i)) para a etapa m/2
            # O índice h+i em psi_inv_rev corresponde a esta raiz
            s = psi_inv_rev[h + i]
            # Itera pelos pares dentro dos blocos
            for j in range(i, n, m): # Stride é m, começando em i
                k = j + h # Índice do par
                # Butterfly:
                # u = input_list[j]
                # v = input_list[k]
                # input_list[j] = (u + v) % q
                # input_list[k] = ((u - v + q) * s) % q

                u = input_list[j]
                v = input_list[k]

                sum_val = u + v
                input_list[j] = sum_val - q if sum_val >= q else sum_val

                # diff_val = (u - v + q) % q
                if u >= v:
                    diff_val = u - v
                else:
                    diff_val = u + q - v

                # input_list[k] = (diff_val * s) % q usando Barrett
                input_list[k] = reduce_inline(diff_val * s, bar)
        # t *= 2 # Não usado no cálculo do índice j diretamente, mas representa o tamanho do bloco atual
        m //= 2

    # Multiplica todos os coeficientes pelo inverso de n (mod q) no final
    if n_inv != 1: # Otimização: não multiplica se n_inv for 1
        for i in range(n):
            input_list[i] = reduce_inline(input_list[i] * n_inv, bar)
# ----------------------------------------------------------------------------------------
# Pointwise multiplication of two polynomials
def pointwise_mul_inplace2(a: list[int], b: list[int], bar: Barrett) -> None:
    """
    Performs the in-place point-by-point multiplication of the elements of list 'a' by the 
    corresponding elements of list 'b', applying modular reduction via the reduce_inline function.

    Args:
    a (list[int]): List of integers that will be modified in-place.
    b (list[int]): List of integers used as a multiplier.
    bar (Barrett): Instance of the Barrett class to perform modular reduction.
    """
    if len(a) != len(b):
        raise ValueError("Lists must be the same size.")
    
    for i in range(len(a)):
        a[i] = reduce_inline(a[i] * b[i], bar)

# ----------------------------------------------------------------------------------------
# Pointwise multiplication of two polynomials
def polymul_ntt(a: list[int], b: list[int], q: int,
                psi_rev: list[int], psi_inv_rev: list[int],
                n_inv: int, bar: Barrett) -> list[int]:
    """
    Multiplies two polynomials using NTT and INTT.

    Each polynomial is represented as a list of coefficients of size n,
    where n is a power of 2.

    Args:
    a (list[int]): Coefficients of the first polynomial.
    b (list[int]): Coefficients of the second polynomial.
    q (int): Modulus.
    psi_rev (list[int]): Vector with the roots of unity for NTT (in bit-reverse order).
    psi_inv_rev (list[int]): Vector with the inverse roots for INTT (in bit-reverse order).
    n_inv (int): Modular inverse of n.
    bar (Barrett): Instance of the Barrett class for modular reduction.

    Returns:
    list[int]: Polynomial resulting from multiplication.   
    """
    if len(a) != len(b):
        raise ValueError("Polynomials must be the same size.")
    
    n = len(a)
    # Make copies so as not to modify the original polynomials
    A = a.copy()
    B = b.copy()

    # Apply NTT (direct transform) to A and B
    ntt_generic(A, psi_rev, q, bar)
    ntt_generic(B, psi_rev, q, bar)

    # Point-to-point multiplication in the NTT domain
    pointwise_mul_inplace2(A, B, bar)

    # Apply INTT (inverse transform) to return to the normal domain
    intt_generic(A, psi_inv_rev, n_inv, q, bar)

    return A

# ----------------------------------------------------------------------------------------
# Generate parameters for NTT and INT
def generate_parameters(n: int, q: int) -> tuple[list[int], list[int], int, Barrett]:
    """
    Generates the necessary parameters for NTT and INTT.

    Args:
    n (int): Dimension (size of the polynomials), must be a power of 2.
    q (int): Module, must satisfy q ≡ 1 (mod 2n).

    Returns:
    tuple: A tuple containing:
    - psi_rev (list[int]): Vector of psi in bit-reversed order.
    - psi_inv_rev (list[int]): Vector of the inverse of psi in bit-reversed order.
    - n_inv (int): Modular inverse of n modulo q.
    - bar (Barrett): Instance of the Barrett class for the modulo q.    
    """
    
    if q % (2 * n) != 1:
        raise ValueError("q must be 1 mod 2n")
    
    # Calculate the primitive root psi
    psi = find_primitive_root(q, n)
    #print(f" * * * * * * Primitive root psi = {psi}", " for q =", q)
    
    # Compute the modular inverse of n (n_inv)
    n_inv = modular_inverse(n, q)
    if n_inv is None:
        raise ValueError("Modular inverse of n does not exist for the given q")
    
    # Generates the psi_rev and psi_inv_rev vectors
    # psi_n_th = pow(psi, 2, q)
    # psi_rev, psi_inv_rev = generate_psi_vectors(n, q, psi_n_th)
    
    psi_rev, psi_inv_rev = generate_psi_vectors(n, q, psi)
    
    # Instantiate the Barrett structure for the q module
    bar = Barrett(q)
    
    return psi_rev, psi_inv_rev, n_inv, bar

# ----------------------------------------------------------------------------------------
# Calculate the modular of a number a mod q (use this function to avoid negative results)
def mod_number(a: int, q: int) -> int:
    """
    Returns the modulus of 'a' with respect to 'q', ensuring a result
    in the range [0, q-1], even if 'a' is negative.

    Args:
    a (int): Integer (can be negative or positive).
    q (int): Modulus (must be positive).

    Returns:
    int: The result of the operation a mod q.
    """
    return ((a % q) + q) % q

# ----------------------------------------------------------------------------------------
# Multiply two polynomials modulo q (batched)
def multiply_poly_mod_batched(a: list[int], b: list[int], q: int) -> list[int]:
    s = len(a)
    if s != len(b):
        raise ValueError("Polynomials must be the same size.")
    
    for ct in range(s):
        a[ct] = mod_number(a[ct] * b[ct], q)

    return a
    