import random

from .glwe import GlweDistribution, GlweSample
from .logging import logger

# from sympy import Poly
from sympy import Poly, GF
from sympy.abc import x
from sympy.polys.specialpolys import random_poly

import math
from typing import Iterable, Tuple


class SecretKey:
    """
    A secret key for the scheme.

    Attributes:
    - dist: The GLWE distribution used to generate the secret key.
    - secret_poly: a secret polynomial.
    """

    def __init__(self, dist: GlweDistribution, secret_poly: Poly):
        self._dist = dist
        self.secret_poly = secret_poly

    def __repr__(self):
        return f'SecretKey({self.secret_poly})'

    @classmethod
    def rand(cls, dist: GlweDistribution, modulus=None) -> 'SecretKey':
        """
        Generates a random secret key.

        Args:
        - dist: A GLWE distribution to use.
        - modulus: The modulus to use for the secret key. If None, the modulus
          of the distribution is used.

        Returns:
        - A secret key.
        """

        secret = dist.sample_polynomial(modulus)
        logger.debug(f"Generating random secret key from: {secret}")
        return cls(dist, secret)

    @property
    def dist(self):
        return self._dist


class PublicKey:
    """
    A public key for the scheme.

    Attributes:
    - glwe_sample: A GLWE sample representing the public key.
    """

    def __init__(self, glwe_sample: GlweSample):
        self.glwe_sample = glwe_sample

    def __repr__(self):
        return f'PublicKey(mask={self.glwe_sample.mask}, '
        f'body={self.glwe_sample.body})'

    @classmethod
    def from_secret_key(cls, secret_key: SecretKey) -> 'PublicKey':
        """
        Generates a public key from a secret key.

        Args:
        - secret_key: The secret key to derive the public key from.

        Returns:
        - A public key.
        """

        sample = secret_key.dist.sample_zero_encryption(secret_key.secret_poly)
        logger.debug(f"Generating public key with sample: {sample}")
        return cls(sample)


def gen_key_pair(dist: GlweDistribution,
                 modulus=None) -> Tuple[SecretKey, PublicKey]:
    """
    Generates a key pair.

    Args:
    - dist: The GLWE distribution to use.
    - modulus: The modulus to use for the secret key. If None, the modulus of
      the distribution is used.

    Returns:
    - A tuple (sk, pk) where sk is the secret key and pk is the public key.
    """

    sk = SecretKey.rand(dist, 2)
    pk = PublicKey.from_secret_key(sk)
    return sk, pk


class RelinKey:
    """
    A relinearization key for the scheme. Used to normalize ciphertexts.

    Attributes:
    - aux_keys: A list of auxiliary keys as described in the reference paper.
    - base: The base/radix used to generate the auxiliary keys.
    """

    def __init__(self, aux_keys: Iterable[GlweSample], base: int):
        self.aux_keys = aux_keys
        self.base = base

    @staticmethod
    def _compute_aux_keys(sk: SecretKey, base: int) -> Iterable[GlweSample]:
        digit_count = math.log(sk.dist.params.ciphertext_modulus, base)
        digit_count = math.ceil(digit_count)
        aux_keys = []
        sk2 = sk.secret_poly ** 2
        for i in range(digit_count):
            mask = sk.dist.sample_mask()
            crt_noise = (sk.dist.sample_crt_noise()
                         .set_domain(sk.dist.cipher_ring))
            masked_secret = mask * sk.secret_poly
            masked_secret = masked_secret % sk.dist.poly_modulus
            noisy_secret = masked_secret + crt_noise
            noisy_secret = noisy_secret % sk.dist.poly_modulus
            message = base ** i * sk2
            message = message % sk.dist.poly_modulus
            body = (noisy_secret + message) % sk.dist.poly_modulus
            aux_keys.append(GlweSample(mask=-mask, body=body))
        return aux_keys

    @staticmethod
    def _compute_aux_keys_crt(sk: SecretKey, base: int) -> Iterable[GlweSample]:
        digit_count = math.log(sk.dist.params.ciphertext_modulus, base)
        digit_count = math.ceil(digit_count)
        aux_keys = []
        sk2 = sk.secret_poly ** 2
        #sk2 = sk.secret_poly * sk.secret_poly
        
        sk2 = sk2 % sk.dist.poly_modulus
        
        # print("=" * 80)
        # print("GERANDO A CHAVE DE RELINEARIZAÇÃO")
        
        # print("SK =", sk.secret_poly );
        # print("SK2 =", sk2 );
        
        # print("modulus: ", sk.dist.poly_modulus)
        
        for i in range(digit_count):
            # mask = sk.dist.sample_mask() 
            grau = sk.dist.params.dimension  # dimensao
            # print("size: ", grau)

            # gerando uma mascara randomica
            tmp = [random.randint(1, sk.dist.params.ciphertext_modulus) for _ in range(grau)]
            mask = Poly(reversed(tmp), x, domain=sk.dist.cipher_ring)

            # mask = Poly(reversed([1,0,0,0]), x, domain=sk.dist.params.ciphertext_modulus)
            
            # print("mask: AUX KEY ", mask)
                        
            # crt_noise = (sk.dist.sample_crt_noise()
            #              .set_domain(sk.dist.cipher_ring))
            masked_secret = mask * sk.secret_poly
            masked_secret = masked_secret % sk.dist.poly_modulus
            
            # print("noise: ", crt_noise)
            
            #noisy_secret = masked_secret + crt_noise
            
            # ruido simples para chave de relinearização:
            # ruido = [random.randint(0, 1) for _ in range(4)]
            # ruido = Poly(reversed(ruido), x, domain=sk.dist.cipher_ring)
            # print("ruido da chave de relinearização: ", ruido)
            # noisy_secret = masked_secret + ruido
            
            # aqui geramos o ruído da chave de relinearização
            ruido = gerar_ruido(sk.dist.params.dimension, 3)
            
            print("Ruido simples - Key: ", ruido)
            size = sk.dist.params.dimension
            ruido_key = []
            for ct in range(size):
                tmp = encode_crt([0, ruido[ct]], [sk.dist.params.plaintext_modulus, 3])
                ruido_key.append(tmp)
            
            # ruido_key = Poly(reversed(ruido), x, domain=sk.dist.cipher_ring)
            ruido_key = Poly(reversed(ruido_key), x, domain=sk.dist.cipher_ring)
            print("Ruido CRT - Key: ", ruido_key)
            
            # AS + (0,e)
            noisy_secret = masked_secret + ruido_key
            noisy_secret = noisy_secret % sk.dist.poly_modulus
            
            # print("noisy_secret (antes)..: ", masked_secret)
            # print("noisy_secret (depois).: ", noisy_secret)
            
            
            # print("base ** i: ", base ** i)
            # print("base: ", base)
            # print("i: ", i)
            # print("sk2: ", sk2)
            
            # w^i * sk^2
            message = (base ** i) * sk2
            message = message % sk.dist.poly_modulus
            
            # print("message after: ", message)
            # print("-----------------------------------------------------")             
            body = (noisy_secret + message) % sk.dist.poly_modulus
            
            # mask = mask * -1
            # print("mask: AUX KEY ", mask)
            # print("mask: AUX KEY (NEG)", -mask)
            
            # print("i = ", i)
            # print("mask: AUX KEY ", mask)
            # print("body: AUX KEY ", body)
            # print("-----------------------------------------------------")
            
            aux_keys.append(GlweSample(mask=-mask , body=body))
        return aux_keys
    
    @classmethod
    def from_secret_key(cls, secret_key: SecretKey,
                        base: int = 2) -> 'RelinKey':
        """
        Generates a relinearization key from a secret key.

        Args:
        - secret_key: The secret key to derive the relinearization key from.
        - base: The base to use for the relinearization key. Defaults to 2.

        Returns:
        - A relinearization key.
        """

        aux_keys = cls._compute_aux_keys(secret_key, base)
        return cls(aux_keys, base)

    @classmethod
    def from_secret_key_crt(cls, secret_key: SecretKey,
                        base: int = 2) -> 'RelinKey':
        """
        Generates a relinearization key from a secret key.

        Args:
        - secret_key: The secret key to derive the relinearization key from.
        - base: The base to use for the relinearization key. Defaults to 2.

        Returns:
        - A relinearization key.
        """

        aux_keys = cls._compute_aux_keys_crt(secret_key, base)
        return cls(aux_keys, base)

    def digit_count(self):
        """
        The number of parts in the relinearization key decomposition.
        """

        return len(self.aux_keys)

# Faz o módulo de número positivo ou negativo
def mod128(a: int, b: int) -> int:
    """
    Calcula a operação módulo de forma a garantir um resultado não negativo,
    equivalente à função Rust apresentada.
    
    Parâmetros:
        a (int): Número inteiro (pode ser negativo).
        b (int): Número inteiro positivo (equivalente a u128 em Rust).
    
    Retorna:
        int: O resultado de ((a % b) + b) % b, sempre não negativo.
    """
    return ((a % b) + b) % b

# retorna o inverso multiplicativo de um numero
def modinv(a: int, q: int) -> int:
    """
    Calcula o inverso multiplicativo de 'a' no módulo 'm', ou seja, encontra um inteiro x tal que (a * x) % m == 1.
    
    Parâmetros:
        a (int): Número inteiro.
        m (int): Módulo (inteiro positivo).
    
    Retorna:
        int: O inverso multiplicativo de 'a' módulo 'm'.
    
    Levanta:
        ValueError: Se o inverso não existir (quando gcd(a, m) != 1).
    """
    # Ajusta 'a' para o intervalo [0, m-1]
    a = a % q
    if a == 0:
        raise ValueError("Não existe inverso multiplicativo para 0 no módulo dado.")
    
    # Inicializa os coeficientes para o algoritmo estendido
    t, new_t = 0, 1
    r, new_r = q, a

    # Algoritmo estendido de Euclides
    while new_r != 0:
        quotient = r // new_r
        t, new_t = new_t, t - quotient * new_t
        r, new_r = new_r, r - quotient * new_r

    # Se r > 1, 'a' e 'm' não são coprimos, logo o inverso não existe
    if r > 1:
        raise ValueError("O inverso multiplicativo não existe pois 'a' e 'm' não são coprimos.")
    
    # Ajusta t para ser positivo
    if t < 0:
        t = t + q
    
    return t

def encode_crt(crtnumber: list[int], base: list[int]) -> int:
    """
    Implementa o algoritmo do Teorema Chinês do Resto para encontrar o número
    que satisfaz o sistema de congruências dado pelos restos em 'crtnumber'
    e modulos em 'base'.
    
    Parâmetros:
        crtnumber (list[int]): Vetor com os restos das congruências.
        base (list[int]): Vetor com os módulos (deve conter números coprimos entre si).
    
    Retorna:
        int: Solução única do sistema de congruências módulo o produto dos elementos de 'base'.
    """
    # Calcula o produto de todos os módulos
    prod = 1
    for b in base:
        prod *= b

    size = len(base)
    # Vetores intermediários
    mult = []  # Cada elemento: prod // base[i]
    vet1 = []  # mod128(mult[i], base[i])
    vet2 = []  # inverso multiplicativo de vet1[i] módulo base[i]

    for i in range(size):
        tmp = prod // base[i]
        mult.append(tmp)
        # Calcula vet1[i] = mod128(tmp, base[i])
        tmp_mod = mod128(tmp, base[i])
        vet1.append(tmp_mod)
        # Calcula vet2[i] = modinv(vet1[i], base[i])
        tmp_inv = modinv(tmp_mod, base[i])
        vet2.append(tmp_inv)

    # Soma total: ∑ (crtnumber[i] * mult[i] * vet2[i])
    sum_tot = 0
    base_tot = 1
    for i in range(size):
        sum_tot += crtnumber[i] * mult[i] * vet2[i]
        base_tot *= base[i]

    # Retorna o resultado final, garantido no intervalo [0, base_tot - 1]
    ret = mod128(sum_tot, base_tot)
    return ret

def gerar_ruido(n: int, limite: int) -> list[int]:
    """
    Gera um vetor com n números inteiros aleatórios no intervalo [0, limite).
    Ou seja, cada posição do vetor terá um valor entre 0 e limite-1.
    
    Parâmetros:
        n (int): Número de elementos no vetor.
        limite (int): Limite superior (não incluso) para os valores.
        
    Retorna:
        list[int]: Uma lista contendo n inteiros no intervalo [0, limite).
    """
    return [random.randrange(limite) for _ in range(n)]