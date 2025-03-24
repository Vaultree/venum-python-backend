
from .rns import RnsBasis
from .crt import CrtEncoder
from .logging import logger

from sympy import Poly, GF
from sympy.abc import x
from sympy.polys.specialpolys import random_poly

import random
from dataclasses import dataclass


@dataclass
class EncryptionParameters:
    """
    Parameters for the encryption scheme.
    """

    dimension: int
    ciphertext_modulus: int
    plaintext_modulus: int
    noise_modulus: int
    seed: int = None

    def __post_init__(self):
        if (self.plaintext_modulus * self.noise_modulus
                >= self.ciphertext_modulus):
            raise ValueError(
                'Invalid parameters: plaintext_modulus * noise_modulus '
                '>= ciphertext_modulus')


class GlweSample:
    """
    A sample from the GLWE distribution.

    The sample is a pair of polynomials (mask, body).
    """

    def __init__(self, mask, body):
        self.mask = mask
        self.body = body

    @staticmethod
    def _compute_mask(mask, noise, u, poly_modulus):
        new_mask = mask * u + noise
        new_mask = new_mask % poly_modulus
        return new_mask

    @staticmethod
    def _compute_body(body, noise, message, u, poly_modulus):
        new_body = body * u + message + noise
        new_body = new_body % poly_modulus
        return new_body

    @classmethod
    def compute_sample(cls, mask, mask_noise, body, body_noise,
                       message, u, poly_modulus):
        """
        Compute a new sample from the given parameters. The new sample
        corresponds to an encryption of the given message.
        """

        new_mask = cls._compute_mask(mask, mask_noise, u,
                                     poly_modulus)
        new_body = cls._compute_body(body, body_noise, message, u,
                                     poly_modulus)
        return cls(mask=new_mask, body=new_body)

    @classmethod
    def _compute_zero_sample(
            cls, mask: Poly, secret: Poly,
            crt_noise: Poly, poly_modulus: Poly):
        body = (mask * secret +
                crt_noise.set_domain(poly_modulus.domain)) % poly_modulus
        return cls(mask=-mask, body=body)

    def __repr__(self):
        return f'GlweSample(mask={self.mask}, body={self.body})'


class GlweDistribution:
    def __init__(self, params: EncryptionParameters):
        """
        Initialize the GLWE distribution with the given parameters.

        Args:
        - params: the encryption parameters.
        """

        if params.seed is not None:
            random.seed(params.seed)
            logger.warning(f"Setting random seed to {params.seed}")
        self.params = params
        self.plaintext_ring = GF(params.plaintext_modulus, symmetric=False)
        self.cipher_ring = GF(params.ciphertext_modulus, symmetric=False)
        self.poly_modulus = Poly(
            x ** params.dimension + 1, x, domain=self.cipher_ring)
        crt_basis = RnsBasis(
            [self.params.plaintext_modulus, self.params.noise_modulus])
        self.crt_encoder = CrtEncoder(crt_basis, self.plaintext_ring)

    def sample_polynomial(self, modulus=None):
        """
        Sample a polynomial with coefficients in the given modulus.

        Args:
        - modulus: the modulus for the coefficients. If None, the
            ciphertext modulus is used.

        Returns:
        - a polynomial with coefficients in the given modulus.
        """

        modulus = modulus or self.params.ciphertext_modulus
        degree = self.params.dimension - 1
        return random_poly(
            x,
            n=degree,
            inf=0,
            sup=modulus - 1
        ).as_poly(domain=self.cipher_ring)

    def sample_mask(self):
        """
        Sample a mask polynomial.

        Returns:

        - a polynomial with coefficients in the ciphertext modulus.
        """

        return self.sample_polynomial()

    def sample_noise(self):
        """
        Sample a noise polynomial.

        Returns:

        - a polynomial with coefficients in the noise modulus.
        """

        # FIX: for security reasons, noise should not be sampled from
        # a uniform distribution.
        return self.sample_polynomial(self.params.noise_modulus)

    def sample_crt_noise(self):
        """
        Sample a CRT noise polynomial.

        Returns:

        - a CRT-encoded noise polynomial.
        """

        noise = self.sample_noise()
        logger.debug(f"Sampled CRT noise: {noise}")
        crt_noise = self.crt_encoder.encode_pure_noise(noise)
        logger.debug(f"CRT noise: {crt_noise}")
        return crt_noise

    def sample_zero_encryption(self, secret: Poly):
        """
        Produces a random GLWE sample corresponding to an encryption of
        zero message.

        Args:
        - secret: the secret polynomial to use for the encryption.

        Returns:
        - a GLWE sample corresponding to an encryption of zero.
        """

        # mask = self.sample_mask()
        # size = len(mask)
        # for ct in range(size):
        #     mask[ct] = 2
        
        #crt_noise = self.sample_crt_noise()
        
        # mask = Poly(reversed([100,200,300,400]), x, domain=self.cipher_ring)
        tmp = [random.randint(1, self.params.ciphertext_modulus-1) for _ in range(self.params.dimension)]
        mask = Poly(reversed(tmp), x, domain=self.cipher_ring)
        
        # crt_noise = Poly(reversed([0,0,0,0]), x, domain=self.cipher_ring)
        
        print("mask (sample zero encryption)......: ", mask)
        
        ruido = gerar_ruido(self.params.dimension, 3)
        print("ruido (noise geração da chave publica): ", ruido)
        size = self.params.dimension
        crt_noise = []
        for ct in range(size):
            tmp = encode_crt([0, ruido[ct]], [self.params.plaintext_modulus, 3])
            crt_noise.append(tmp)
            
        crt_noise = Poly(reversed(crt_noise), x, domain=self.cipher_ring)
        print("crt_noise (sample zero encryption).: ", crt_noise)
        
        # print("mask (sample zero encryption)......: ", mask)
        # print("crt_noise (sample zero encryption).: ", crt_noise)
        
        return GlweSample._compute_zero_sample(
            mask, secret, crt_noise, self.poly_modulus)

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