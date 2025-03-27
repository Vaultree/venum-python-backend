import random
from symtable import Symbol
from venum.ntt import *
from sympy import Poly, symbols
from sympy.polys.domains import GF

class Cryptogram:
    def __init__(self, mask: list[int], body: list[int]) -> None:
        self.mask = mask   # Atributo público do tipo list[int]
        self.body = body   # Atributo público do tipo list[int]

    def __repr__(self) -> str:
        return f"CustomObject(mask={self.mask}, body={self.body})"
# --------------------------------------------------------------
# function encode CRT
def encode_crt(crtnumber: list[int], base: list[int]) -> int:
    """
    Implements the Chinese Remainder Theorem algorithm to find the number
    that satisfies the system of congruences given by the remainders in 'crtnumber'
    and moduli in 'base'.

    Parameters:
    crtnumber (list[int]): Vector with the remainders of the congruences.
    base (list[int]): Vector with the moduli (must contain coprime numbers).

    Returns:
    int: Unique solution of the system of congruences modulo the product of the elements of 'base'.
    """
    
    # Calculate the product of all modules
    prod = 1
    for b in base:
        prod *= b

    size = len(base)
    # Intermediate vectors
    mult = []  # prod // base[i]
    vet1 = []  # mod_number(mult[i], base[i])
    vet2 = []  # modular_inverse vet1[i] modulus base[i]

    for i in range(size):
        tmp = prod // base[i]
        mult.append(tmp)
        # Calculate vet1[i] = mod128(tmp, base[i])
        tmp_mod = mod_number(tmp, base[i])
        vet1.append(tmp_mod)
        # Calculate vet2[i] = modinv(vet1[i], base[i])
        tmp_inv = modular_inverse(tmp_mod, base[i])
        vet2.append(tmp_inv)

    # Total sum: ∑ (crtnumber[i] * mult[i] * vet2[i])
    sum_tot = 0
    base_tot = 1
    for i in range(size):
        sum_tot += crtnumber[i] * mult[i] * vet2[i]
        base_tot *= base[i]

    # Returns the final result, guaranteed to be in the range [0, base_tot - 1]
    ret = mod_number(sum_tot, base_tot)
    return ret

# --------------------------------------------------------------
# function to convert object Poly to list
def poly_to_vector(poly, dimension: int) -> list:
    """
    Converts a Poly object into a vector (list) of coefficients with length 'dimension'.
    The coefficients are extracted and, if the polynomial has a degree lower than 'dimension - 1',
    zeros are added to fill up to the desired dimension, preserving leading and trailing zeros.

    Args:
    poly: Object representing a polynomial. This can be, for example, a sympy.Poly.
    dimension (int): Desired size of the vector of coefficients.

    Returns:
    list: List of coefficients of the polynomial with length equal to 'dimension'.
    """
    # Try using the all_coeffs() method (e.g. sympy.Poly)
    if hasattr(poly, "all_coeffs") and callable(poly.all_coeffs):
    # all_coeffs() returns the coefficients in descending order (highest degree first),
    # then we reverse them to get ascending order (degree 0 first)
        coeffs = list(reversed(poly.all_coeffs()))
    elif hasattr(poly, "coeffs"):
        coeffs = poly.coeffs
    else:
        raise ValueError("The Poly object does not have a known method or attribute to extract the coefficients.")
    
    # If the extracted vector has fewer elements than the desired dimension,
    # add zeros to the end.
    if len(coeffs) < dimension:
        coeffs.extend([0] * (dimension - len(coeffs)))
    # If it has more elements than the dimension, truncate the vector.
    elif len(coeffs) > dimension:
        coeffs = coeffs[:dimension]

    return coeffs

# --------------------------------------------------------------
# function to convert list to object Poly
def vector_to_poly(vec: list[int], q: int, variable='x') -> Poly:
    
    # Se 'variable' for uma string, converte-a para um Symbol usando a função symbols.
    if isinstance(variable, str):
        variable = symbols(variable)
    
    poly = Poly(reversed(vec), variable, domain=GF(q))
    return poly

# --------------------------------------------------------------
# create SK function
def create_sk(n: int, q: int) -> list[int]:
    """
    Generates a secret key (SK) for the Ring-LWE cryptosystem.
    The secret key is a vector of 'n' elements with values between 0 and 'q-1'.

    Parameters:
    n (int): Dimension of the secret key.
    q (int): Modulus of the Ring-LWE cryptosystem.

    Returns:
    list[int]: Secret key vector with 'n' elements.
    """
    return generate_random_vector(n, 0, q-1)

# --------------------------------------------------------------
# generate a random vector
def generate_random_vector(n, min_val, max_val):
    """
    Gera um vetor aleatório com n elementos.
    
    Parâmetros:
    - n (int): número de elementos do vetor.
    - min_val (float): valor mínimo do intervalo.
    - max_val (float): valor máximo do intervalo.
    
    Retorna:
    - list: vetor com n números aleatórios entre min_val e max_val.
    """
    return [random.randint(min_val, max_val) for _ in range(n)]

# --------------------------------------------------------------
# generate a mask vector
def generate_mask_vector(n, q):
    """
    Generates a random mask vector with 'n' elements in the range [0, q-1].
    
    Parameters:
    n (int): Dimension of the mask vector.
    q (int): Modulus of the Ring-LWE cryptosystem.
    
    Returns:
    list[int]: Mask vector with 'n' elements.
    """
    return generate_random_vector(n, 0, q-1)

# --------------------------------------------------------------
# generate a noise CRT vector
def generate_noise_crt(n: int, p1: int, p2: int) -> list[int]:
    
    noise_vector = generate_random_vector(n, 0, p2)
    crt_noise = []
    for ct in range(n):
        tmp = encode_crt([0, noise_vector[ct]], [p1, p2])
        crt_noise.append(tmp)
    
    return crt_noise

# --------------------------------------------------------------
# generate public key function
def generate_pk(sk: list[int], n: int, q: int, p1: int, p2: int, params: tuple[list[int], list[int], int, Barrett]) -> Cryptogram:
    """
    Generates a public key (PK) for the Ring-LWE cryptosystem.
    The public key is a vector of 'n' elements with values between 0 and 'q-1'.

    Parameters:
    sk (list[int]): Secret key vector.
    n (int): Dimension of the secret key and public key.
    q (int): Modulus of the Ring-LWE cryptosystem.
    p1 (int): First prime number for CRT encoding.
    p2 (int): Second prime number for CRT encoding.

    Returns:
    list[int]: Public key vector with 'n' elements.
    """
    # Generate a mask vector
    mask = generate_mask_vector(n, q)
    
    # Generate a noise vector
    noise = generate_noise_crt(n, p1, p2)
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params
    
    # multiply mask * secret key
    mask_key = polymul_ntt(mask, sk, q, psi_rev, psi_inv_rev, n_inv, bar)
    
    # add noise to the result
    body = [mod_number(mask_key[i] + noise[i], q) for i in range(n)]
    
    # convert the mask for negative multiply to -1
    mask = [mod_number(-1 * mask[i]) for i in range(n)]
    
    # Calculate the public key
    pk = Cryptogram(body, mask)
    
    return pk