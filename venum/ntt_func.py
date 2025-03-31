import random
from symtable import Symbol
import sys
from venum.ntt import *
from sympy import Poly, symbols
from sympy.polys.domains import GF
from dataclasses import dataclass, field


@dataclass
class PublicKey:
    body: List[int] = field(default_factory=list)
    mask: List[int] = field(default_factory=list)
    q: int = 0
@dataclass
class Cryptogram:
    body: List[int] = field(default_factory=list)
    mask: List[int] = field(default_factory=list)
    batched: bool = False
    q: int = 0
    
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
def create_sk(n: int, minimo: 0, maximo: int, q: int) -> list[int]:
    """
    Generates a secret key (SK) for the Ring-LWE cryptosystem.
    The secret key is a vector of 'n' elements with values between 0 and 'q-1'.

    Parameters:
    n (int): Dimension of the secret key.
    q (int): Modulus of the Ring-LWE cryptosystem.

    Returns:
    list[int]: Secret key vector with 'n' elements.
    """
    
    key = generate_random_vector(n, minimo, maximo)
    for ct in range(n):
        key[ct] = mod_number(key[ct], q)
    
    return key

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
    
    noise_vector = generate_random_vector(n, 0, p2-1)
    crt_noise = []
    for ct in range(n):
        tmp = encode_crt([0, noise_vector[ct]], [p1, p2])
        crt_noise.append(tmp)
    
    return crt_noise

# --------------------------------------------------------------
# generate a noise CRT vector
def encode_msg_crt(msg: List[int], n: int, p1: int, p2: int) -> list[int]:
    
    noise_vector = generate_random_vector(n, 0, p2-1)
    message = []
    for ct in range(n):
        tmp = encode_crt([msg[ct], noise_vector[ct]], [p1, p2])
        message.append(tmp)
    
    return message

# --------------------------------------------------------------
# generate public key function
def generate_pk(sk: List[int], q: int, p1: int, p2: int, params: tuple[list[int], list[int], int, Barrett]) -> PublicKey:
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
    n = len(sk)
    
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
    mask = [mod_number(-1 * mask[i], q) for i in range(n)]
    
    # Calculate the public key
    pk = PublicKey(body=body, mask=mask, q=q)
    
    return pk

# -----------------------------------------------------
# Gerar número primo que atenda à condição específica
def is_prime(n, k=16):
    """
    Teste de primalidade probabilístico de Miller-Rabin.
    
    Parâmetros:
      n: inteiro a ser testado.
      k: número de iterações (quanto maior, maior a precisão do teste).
    
    Retorna:
      True se n é provavelmente primo, False se é composto.
    """
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    
    # Escreve n-1 como d * 2^s
    s = 0
    d = n - 1
    while d % 2 == 0:
        s += 1
        d //= 2
    
    # Executa k iterações do teste
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

# --------------------------------------------------------------------------------
def generate_modulus(inicio, fim, n):
    """
    Procura no intervalo [inicio, fim] o primeiro número primo que satisfaça:
       primo % (2 * n) == 1.
    Se o primo inicial não atender à condição, soma 2 ao candidato e testa novamente,
    até encontrar um primo válido ou ultrapassar o limite 'fim'.
    
    Parâmetros:
      inicio: início do intervalo de busca.
      fim: fim do intervalo de busca.
      n: inteiro usado na condição (primo % (2*n) == 1).
    
    Retorna:
      O número primo que atende à condição ou None se nenhum for encontrado.
    """
    candidato = None
    
    # aleatoriza o início da busca
    inicio = random.randint(0, fim + 1)
    
    # Encontra o primeiro primo no intervalo
    for i in range(inicio, fim + 1):
        if is_prime(i):
            candidato = i
            break

    if candidato is None:
        print("Nenhum número primo encontrado no intervalo.")
        return None

    # Testa a condição para o primo encontrado, incrementando de 2 se necessário
    while candidato <= fim:
        if candidato % (2 * n) == 1:
            return candidato
        candidato += 2
        # Garante que o novo candidato seja primo
        while candidato <= fim and not is_prime(candidato):
            candidato += 2

    print("Nenhum primo que atenda à condição foi encontrado no intervalo.")
    return None

# -----------------------------------------------------
def encrypt_sk(sk: List[int], msg: List[int], q: int, p1: int, p2: int, batched: bool, params: tuple[list[int], list[int], int, Barrett], params_batched: tuple[list[int], list[int], int, Barrett]) -> Cryptogram:
    """
    Encrypts a message using the secret key and a mask.
    
    Parameters:
    sk (list[int]): Secret key vector.
    msg (list[int]): Message vector.
    n (int): Dimension of the vectors.
    params: Parameters for NTT/INTT.
    q (int): Modulus of the Ring-LWE cryptosystem.
    p1 (int): First prime number for CRT encoding.      
    p2 (int): Second prime number for CRT encoding.
    batched (bool): If True, the encryption is done in batch mode.
    params_batched: Parameters for NTT/INTT in batched mode.
    
    Returns:
    Cryptogram: Encrypted message.
    """
    
    #print("Start message: ", msg)
    # verify the cryptography mode batched or not
    if batched:
        msg = msg.copy()
        psi_rev, psi_inv_rev, n_inv, bar = params_batched
        intt_generic(msg, psi_inv_rev, n_inv, p1, bar)  
        #print("Batched message: ", msg)
        
    n = len(sk)
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params
    
    # Generate a mask vector
    mask = generate_mask_vector(n, q)
    
    # Generate a noise vector
    lmessage = encode_msg_crt(msg, n, p1, p2)
    
    # multiply mask * secret key
    mask_key = polymul_ntt(mask, sk, q, psi_rev, psi_inv_rev, n_inv, bar)
    
    # add noise to the result
    body = [mod_number(mask_key[i] + lmessage[i], q) for i in range(n)]
    
    # convert the mask for negative multiply to -1
    mask = [mod_number(-1 * mask[i], q) for i in range(n)]
    
    # Calculate the public key
    crypto = Cryptogram(body=body, mask=mask, batched=batched,q=q)
    
    return crypto

# ------------------------------------------------------------------------------
def decrypt(sk: list[int], crypto: Cryptogram, p1: int, params: tuple[list[int], list[int], int, Barrett], params_batched: tuple[list[int], list[int], int, Barrett]) -> list[int]:
    """
    Decrypts a message using the secret key and a mask.
    Parameters:
    sk (list[int]): Secret key vector.
    crypto (Cryptogram): Encrypted message.
    p1 (int): First prime number for CRT encoding.
    params (tuple): Parameters for NTT/INTT.
    Returns:
    list[int]: Decrypted message.
    """
    n = len(sk)
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params
    
    # calculate AS: mask * secret key
    mask_key = polymul_ntt(crypto.mask, sk, crypto.q, psi_rev, psi_inv_rev, n_inv, bar)
    
    # Subtract the result from the body of the cryptogram
    result = [mod_number(crypto.body[i] + mask_key[i], crypto.q) for i in range(n)]
    
    # Decode the result using CRT (verify negative numbers)
    ret = []
    limit = crypto.q // 2
    for ct in result:
        tmp = ct
        if tmp > limit:
            tmp -= crypto.q
            
        tmp = mod_number(tmp, p1)
        # if tmp < 0 or tmp >= p1:
        #     print("Error: Decryption failed.", "ct", tmp, "| p1 = ", p1)
        ret.append(tmp % p1) # here decrypt the message CRT
        
    if crypto.batched:
        #print("Batched message antes - decifragem: ", ret)
        psi_rev, psi_inv_rev, n_inv, bar = params_batched
        ntt_generic(ret, psi_inv_rev, p1, bar) 
        ret = ret[::-1]
        #print("Batched message depois - decifragem: ", ret)
        return ret

    return ret
# ------------------------------------------------------------------------------
# function encrypt pk
def encrypt_pk(pk: PublicKey, msg: List[int], q: int, p1: int, p2: int, batched: bool, params: tuple[list[int], list[int], int, Barrett], params_batched: tuple[list[int], list[int], int, Barrett]) -> Cryptogram:
    """
    Encrypts a message using the public key.
    
    Parameters:
    pk (PublicKey): Public key.
    msg (list[int]): Message vector.
    q (int): Modulus of the Ring-LWE cryptosystem.
    p1 (int): First prime number for CRT encoding.
    p2 (int): Second prime number for CRT encoding.
    batched (bool): If True, the encryption is done in batch mode.
    params (tuple): Parameters for NTT/INTT.
    params_batched (tuple): Parameters for NTT/INTT in batch
    
    Returns:
    Cryptogram: Encrypted message.
    """
    n = len(pk.body)
    
    # verify the cryptography mode batched or not
    if batched:
        msg = msg.copy()
        psi_rev, psi_inv_rev, n_inv, bar = params_batched
        intt_generic(msg, psi_inv_rev, n_inv, p1, bar)  
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params
    
    # Generate a noise vector
    message = encode_msg_crt(msg, n, p1, p2)

    # generate vector u
    u = generate_random_vector(n, 1, 2)
    
    # generate noise for mask
    noise_mask = generate_noise_crt(n, p1, p2)
    noise_mask = [mod_number(0, q) for i in range(n)]
    
    # generate noise for body
    noise_body = generate_noise_crt(n, p1, p2)
    noise_body = [mod_number(0, q) for i in range(n)]
    
    # calculate AS: mask * secret key
    pk0 = polymul_ntt(pk.body, u, q, psi_rev, psi_inv_rev, n_inv, bar)
    pk1 = polymul_ntt(pk.mask, u, q, psi_rev, psi_inv_rev, n_inv, bar)
    
    # add to noise
    pk0 = [mod_number(pk0[i] + noise_body[i], q) for i in range(n)]
    pk1 = [mod_number(pk1[i] + noise_mask[i], q) for i in range(n)]
    
    # add to message
    body = [mod_number(pk0[i] + message[i], q) for i in range(n)]
        
    # Calculate the public key
    crypto = Cryptogram(body=body, mask=pk1, batched=batched,q=q)
    
    return crypto

# ----------------------------------------------------------------------------------
# sum of cryptograms
def sum_cryptograms(c1: Cryptogram, c2: Cryptogram):
    
    n = len(c1.body)
    if n != len(c2.body):
        raise ValueError("The two ciphertexts must have the same length")
    
    if c1.q != c2.q:
        raise ValueError("The two ciphertexts must have the same modulus")  
    
    if c1.batched != c2.batched:
        raise ValueError("The two ciphertexts must have the same batched mode")
    
    body = []
    mask = []
    for i in range(n):
        body.append(mod_number(c1.body[i] + c2.body[i] , c1.q))
        mask.append(mod_number(c1.mask[i] + c2.mask[i] , c1.q))
    
    return Cryptogram(body=body, mask=mask, batched=c1.batched, q=c1.q)

# ----------------------------------------------------------------------------------
# function add messages
def add_msg(m1: List[int], m2: List[int], q):
    
    n = len(m1)
    if n != len(m2):
        raise ValueError("The two messages must have the same length")
    
    return [(mod_number(m1[i] + m2[i], q)) for i in range(n)]

# ----------------------------------------------------------------------------------
# sum of cryptograms
def sub_cryptograms(c1: Cryptogram, c2: Cryptogram):
    
    n = len(c1.body)
    if n != len(c2.body):
        raise ValueError("The two ciphertexts must have the same length")
    
    if c1.q != c2.q:
        raise ValueError("The two ciphertexts must have the same modulus")  
    
    if c1.batched != c2.batched:
        raise ValueError("The two ciphertexts must have the same batched mode")
    
    body = []
    mask = []
    for i in range(n):
        body.append(mod_number(c1.body[i] - c2.body[i] , c1.q))
        mask.append(mod_number(c1.mask[i] - c2.mask[i] , c1.q))
    
    return Cryptogram(body=body, mask=mask, batched=c1.batched, q=c1.q)

# ----------------------------------------------------------------------------------
# function add messages
def sub_msg(m1: List[int], m2: List[int], q):
    
    n = len(m1)
    if n != len(m2):
        raise ValueError("The two messages must have the same length")
    
    return [(mod_number(m1[i] - m2[i], q)) for i in range(n)]

# ----------------------------------------------------------------------------------
# generate de relinearization key
def generate_rlk(sk: List[int], base: int, p1: int, p2: int, q: int, params: tuple[list[int], list[int], int, Barrett]) -> list[Cryptogram]:
    
    digit_count = math.log(q, base)
    digit_count = math.ceil(digit_count)
    
    n = len(sk)
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params

    rlk = []
    
    # calculate S^2
    sk2 = polymul_ntt(sk, sk, q, psi_rev, psi_inv_rev, n_inv, bar)
    
    for i in range(digit_count):
        
        mask = generate_mask_vector(n, q)
        mask2 = mask.copy()
        noise = generate_noise_crt(n, p1, p2)
        
        # calculate AS
        mask_secret = polymul_ntt(mask, sk, q, psi_rev, psi_inv_rev, n_inv, bar)
        
        # add noise to the result
        mask_secret = [mod_number(mask_secret[i] + noise[i], q) for i in range(n)]
        
        message = sk2.copy()
        
        # multiply the message by the base
        tmp = base ** i
        message = [mod_number(tmp * message[i], q) for i in range(n)]
            
        body = [mod_number(mask_secret[i] + message[i], q) for i in range(n)]
        mask = [mod_number(-1 * mask2[i], q) for i in range(n)]
        
        rlk.append(Cryptogram(body=body, mask=mask, batched=False, q=q))
    
    return rlk

# ----------------------------------------------------------------------------------
# calculate product of two cryptograms
def product(cripto0: Cryptogram, cripto1: Cryptogram, params: tuple[list[int], list[int], int, Barrett]) -> (List[int], List[int],List[int]): # type: ignore
    
    n = len(cripto0.body)
    if n != len(cripto1.body):
        raise ValueError("The two ciphertexts must have the same length")
    
    if cripto0.q != cripto1.q:
        raise ValueError("The two ciphertexts must have the same modulus")
    
    if cripto0.batched != cripto1.batched:
        raise ValueError("The two ciphertexts must have the same batched mode")
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params
    
    # body * body:
    c0 = polymul_ntt(cripto0.body, cripto1.body, cripto0.q, psi_rev, psi_inv_rev, n_inv, bar)
    
    # (body * mask) + (mask * body):
    t1 = polymul_ntt(cripto0.body, cripto1.mask, cripto0.q, psi_rev, psi_inv_rev, n_inv, bar)
    t2 = polymul_ntt(cripto0.mask, cripto1.body, cripto0.q, psi_rev, psi_inv_rev, n_inv, bar)
    c1 = [mod_number(t1[i] + t2[i], cripto0.q) for i in range(n)]
    
    # mask * mask:
    c2 = polymul_ntt(cripto0.mask, cripto1.mask, cripto0.q, psi_rev, psi_inv_rev, n_inv, bar)
    
    return (c0, c1, c2)
   
# ----------------------------------------------------------------------------------
def decompose_poly_list(poly: list[int], base: int, q: int) -> list[list[int]]:
    """
    Decompõe um polinômio (representado como uma lista de inteiros) em uma soma de polinômios,
    onde os coeficientes de cada polinômio componente são os dígitos da representação
    dos coeficientes do polinômio original na base 'base'. A operação trabalha em GF(modulo).

    Cada polinômio é representado como uma lista de inteiros, onde o i-ésimo elemento é o coeficiente de x^i.
    Assim, para cada i temos:
    
        poly[i] = components[0][i] + components[1][i]*base + components[2][i]*base^2 + ... + components[num_components-1][i]*base^(num_components-1)
    
    Parâmetros:
      poly          : list[int]
                      Polinômio de entrada (coeficiente de x^i é poly[i]).
      base          : int
                      Base usada para decomposição (por exemplo, 2 para binário ou 10 para decimal).
      num_components: int
                      Número de dígitos/componentes a serem extraídos.
      modulo        : int
                      Módulo usado para trabalhar em GF(modulo). Os coeficientes são reduzidos módulo 'modulo'.
    
    Retorna:
      Uma lista de listas de inteiros, onde o j-ésimo elemento (0 ≤ j < num_components) é o polinômio
      componente correspondente ao dígito extraído para a potência base^j.
    
    Lança:
      ValueError: Se algum coeficiente do polinômio requer mais dígitos do que 'num_components'.
    """
    # Reduz os coeficientes do polinômio no corpo GF(modulo)
    poly_mod = [c % q for c in poly]
    n = len(poly_mod)
    
    # Inicializa num_components listas, cada uma com n coeficientes (um para cada termo do polinômio)
    num_components = math.log(q, base)
    num_components = math.ceil(num_components)
    
    components = [[0] * n for _ in range(num_components)]
    
    # print("componentes: ", num_components)
    
    # Para cada coeficiente, extrai os dígitos na base 'base'
    for i, coeff in enumerate(poly_mod):
        a = coeff
        for j in range(num_components):
            digit = a % base
            # Em GF(modulo), garantimos que o dígito também esteja reduzido, embora geralmente base < modulo
            components[j][i] = digit % q
            a //= base
            # print("a = ", a)
            
        if a != 0:
            raise ValueError(f"O coeficiente do termo x^{i} requer mais de {num_components} dígitos na base {base}.")
    
    return components

# ----------------------------------------------------------------------------------
# relinearization
def relinearize(prod: (List[int], List[int], List[int]), rlk: list[Cryptogram], batched: bool, base_decomposition: int, q: int, params: tuple[list[int], list[int], int, Barrett]) -> Cryptogram: # type: ignore
    
    (c0, c1, c2) = prod
    n = len(c0)
    
    # decompose the c2
    quad_decomposed = decompose_poly_list(c2, base_decomposition, q)
    
    # decomp = quad_decomposed.copy()
    # lista = [0 for i in range(n)]
    # i = 0
    # for tmp in decomp:
    #     t1 = base_decomposition ** i
    #     t2 = [mod_number(tmp[j] * t1, q) for j in range(n)]
    #     lista = [mod_number(lista[j] + t2[j], q) for j in range(n)]
    #     i += 1
        
    # if lista != c2:
    #     print("list: ", lista)
    #     print("c2: ", c2)
    #     print("ERRO NA DECOMPOSIÇÃO")
    #     sys.exit(1)
    # else:
    #     print("DECOMPOSIÇÃO OK")
        
    
    # collect the parameters
    psi_rev, psi_inv_rev, n_inv, bar = params

    # start mask and body
    mask = [0 for i in range(n)]
    body = [0 for i in range(n)]
    
    for aux_key, component in zip(rlk, quad_decomposed):   
        # print("component: ", component, "tamanho: ", len(component))
          
        tmp = polymul_ntt(aux_key.mask, component, q, psi_rev, psi_inv_rev, n_inv, bar)
        mask = [mod_number(mask[i] + tmp[i], q) for i in range(n)]
        tmp = polymul_ntt(aux_key.body, component, q, psi_rev, psi_inv_rev, n_inv, bar)
        body = [mod_number(body[i] + tmp[i], q) for i in range(n)]
        
    mask = [mod_number(mask[i] + c1[i], q) for i in range(n)]
    body = [mod_number(body[i] + c0[i], q) for i in range(n)]   
    
    return Cryptogram(body=body, mask=mask, batched=batched, q=q)

# -----------------------------------------------------------------------------------------
# function to multiply two polynomials
def multiply_cryptograms(cripto0: Cryptogram, cripto1: Cryptogram, rlk: list[Cryptogram], base_decomposition, params: tuple[list[int], list[int], int, Barrett]) -> Cryptogram:
    
    # calculate the product
    prod = product(cripto0, cripto1, params)
    
    # calculate the relinearization
    c3 = relinearize(prod, rlk, cripto0.batched, base_decomposition, cripto0.q, params)
    
    return c3   
 