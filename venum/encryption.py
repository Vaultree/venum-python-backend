import sys
from .logging import logger
from .glwe import GlweSample, GlweDistribution
from .key import SecretKey, PublicKey, RelinKey
from .numeric import radix_decompose_poly

# from sympy import Poly
from sympy.abc import x

from typing import Iterable
from sympy import Poly, GF



class Cipher:
    """
    A class representing a GLWE ciphertext.

    Attributes:
    - glwe_sample: A GlweSample object representing the ciphertext.
    """

    def __init__(self, glwe_sample):
        self.glwe_sample = glwe_sample

    def __repr__(self):
        return f'Cipher(mask={self.glwe_sample.mask}, '
        f'body={self.glwe_sample.body})'


class Encryptor:
    """
    A class for handling encryption and decryption of messages.
    """

    def __init__(self, dist: GlweDistribution, plaintext_encoder):
        """
        Initializes an Encryptor object.

        Args:
        - dist: A GlweDistribution object representing the distribution
          used for encryption.
        - plaintext_encoder: An object that encodes and decodes messages
          according to the `venum.plaintext_encoding.Encoder` interface.
        """

        self._dist = dist
        self.plaintext_encoder = plaintext_encoder

    @property
    def dist(self):
        return self._dist

    def _encrypt_pk(self, pk: PublicKey, message: Iterable[int],
                    plaintext_encoder=None) -> Cipher:
        """
        Encrypts a message.

        Args:
        - pk: A PublicKey object representing the public key.
        - message: An iterable of integers representing the message.
        - plaintext_encoder: An object that encodes and decodes messages
          according to the `venum.plaintext_encoding.Encoder` interface.
          If None, the default encoder is used.

        Returns:
        - A Cipher object representing the encrypted message.
        """

        logger.debug(f'Encrypting message: {message}')

        plaintext_encoder = plaintext_encoder or self.plaintext_encoder

        message = plaintext_encoder.encode(message)
        logger.debug(f'encoded message: {message}')

        crt_message = self.dist.crt_encoder.encode_pure_message(
            message).set_domain(self.dist.cipher_ring)
        logger.debug(f'crt_message: {crt_message}')

        crt_noise1 = (self.dist.sample_crt_noise()
                      .set_domain(self.dist.cipher_ring))
        crt_noise2 = (self.dist.sample_crt_noise()
                      .set_domain(self.dist.cipher_ring))

        u = self.dist.sample_polynomial(modulus=2)
        logger.debug(f'sampled u: {u}')

        logger.debug(f"using public key: {pk.glwe_sample}")
        sample = GlweSample.compute_sample(
            mask=pk.glwe_sample.mask,
            mask_noise=crt_noise2,
            body=pk.glwe_sample.body,
            body_noise=crt_noise1,
            message=crt_message,
            u=u,
            poly_modulus=self.dist.poly_modulus,
        )
        return Cipher(sample)

    def _encrypt_sk(self, sk: SecretKey, message: Iterable[int],
                    plaintext_encoder=None) -> Cipher:
        logger.debug(f'Encrypting message: {message}')

        plaintext_encoder = plaintext_encoder or self.plaintext_encoder

        message = plaintext_encoder.encode(message)
        logger.debug(f'encoded message: {message}')
        
        # forcando o ruido == 0
        ruido = plaintext_encoder.encode([0,0,0,0])
        
        crt_message = self.dist.crt_encoder.encode(
            message,
            # self.dist.sample_noise()
            ruido
        ).set_domain(self.dist.cipher_ring)
        
        print("message cifragem.: ", message)
        print("message CRT......: ", message)
        
        if crt_message != message:
            print("ERRO NA ENCODING =================================")
        
        # TODO: extract this into an easily testable function
        zero_sample = sk.dist.sample_zero_encryption(sk.secret_poly)
        
        print("zero_sample mask: ", zero_sample.mask)
        print("zero_sample body: ", zero_sample.body)
        
        body = zero_sample.body + crt_message
        mask = zero_sample.mask
        sample = GlweSample(body=body, mask=mask)
        return Cipher(sample)

    def encrypt(self, key: SecretKey | PublicKey, message: Iterable[int],
                plaintext_encoder=None) -> Cipher:
        if isinstance(key, SecretKey):
            logger.debug("Secret key encryption")
            return self._encrypt_sk(key, message, plaintext_encoder=plaintext_encoder)
        else:
            logger.debug("Public key encryption")
            return self._encrypt_pk(key, message, plaintext_encoder=plaintext_encoder)

    def decrypt(self, sk: SecretKey, cipher: Cipher) -> Iterable[int]:
        """
        Decrypt_s a ciphertext.

        Args:
        - sk: A SecretKey object representing the secret key.
        - cipher: A Cipher object representing the ciphertext.

        Returns:
        - An iterable of integers representing the decrypted message.
        """

        logger.debug(f"{cipher}")
        cipher_mask = cipher.glwe_sample.mask
        cipher_body = cipher.glwe_sample.body

        crt_message = (cipher_body + cipher_mask * sk.secret_poly)
        crt_message = crt_message % self.dist.poly_modulus
        logger.debug(f"{crt_message}")
        
        print("crt_message decifragem: ", crt_message)
        
        noisy_message = self.dist.crt_encoder.decode(crt_message)
        
        print("noise_message: ", noisy_message)
        
        logger.debug(f"{noisy_message}")

        noiseless_coefs = [rns[0] for rns in noisy_message]

        # correct for dimension size
        len_diff = self.dist.params.dimension - len(noiseless_coefs)
        if len_diff > 0:
            noiseless_coefs.extend([0] * len_diff)

        noiseless_message = Poly.from_list(noiseless_coefs,
                                           x, domain=self.dist.plaintext_ring)

        logger.debug(f"{noiseless_message}")
        cleartext = self.plaintext_encoder.decode(noiseless_message)
        return cleartext


class Rank2Cipher:
    """
    A class representing a rank-2 ciphertext, representing a non-normalized
    Cipher containing squared terms. Usually produced as an intermediate step
    during homomorphic multiplication.

    Attributes:
    - constant: A Poly object representing the constant term over the secret.
    - linear: A Poly object representing the linear term over the secret.
    - quadratic: A Poly object representing the quadratic term over the secret.
    """

    def __init__(self, constant: Poly, linear: Poly, quadratic: Poly):
        self.constant = constant
        self.linear = linear
        self.quadratic = quadratic

    def relinearize(self, relin_key: RelinKey, poly_modulus) -> Cipher:
        """
        Relinearizes the rank-2 ciphertext into a normalized Cipher.

        Args:
        - relin_key: A RelinKey object representing the relinearization key.
        - poly_modulus: A Poly object representing the polynomial modulus.

        Returns:
        - A Cipher object representing the relinearized ciphertext.
        """

        cipher_ring = poly_modulus.domain
        quad_decomposed = radix_decompose_poly(
            poly=self.quadratic,
            radix=relin_key.base,
            num_components=relin_key.digit_count(),
            domain=cipher_ring
        )
        
        print("radix: ", relin_key.base);
        print("digit_count: ", relin_key.digit_count());
        print("Domain: ", cipher_ring);
        
        # decomposicao binaria
        # quad_decomposed = decompose_poly(self.quadratic, relin_key.base, relin_key.digit_count(), cipher_ring);
        
        # print("passou aqui");
        
        print("quadratico: ", self.quadratic);
        quad_decomposed = [x for x in quad_decomposed];
        print("Decomposicão: ", sum(1 for _ in quad_decomposed));
        # sys.exit(1)
        
        # if sum(1 for _ in quad_decomposed) != relin_key.digit_count():
        #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        #     print(len(quad_decomposed), relin_key.digit_count())
        #     print("ERRO NA DECOMPOSICAO BINÁRIA (numero de elementos diferentes )=================================")
        #     sys.exit(1)
        
        decomp = Poly([0], x, domain=cipher_ring)
        posic = 0
        for ct in quad_decomposed:
            #print("ct: ", ct)
            if posic == 0:
                decomp = ct * (2 ** posic)
            else:
                decomp = decomp + (ct * (2 ** posic))
                decomp = decomp % poly_modulus
            posic += 1
            
        print("decomp.....: ", decomp)
        print("quadratico.: ", self.quadratic)
        
        if decomp != self.quadratic:
            print("ERRO NA DECOMPOSICAO BINÁRIA =================================")
            sys.exit(1)
        
        mask = Poly([0], x, domain=cipher_ring)
        body = Poly([0], x, domain=cipher_ring)
        for aux_key, component in zip(relin_key.aux_keys,
                                      quad_decomposed):
            mask += aux_key.mask * component
            mask = mask % poly_modulus
            body += aux_key.body * component
            body = body % poly_modulus
            
        mask += self.linear
        body += self.constant
        
        print("Poly modulus >>>>>>>>>>>>>>>>: ", poly_modulus)

        mask = mask % poly_modulus
        body = body % poly_modulus
        return Cipher(GlweSample(mask=mask, body=body))
        # return Cipher(GlweSample(mask=body, body=mask))


def decompose_poly(poly: Poly, base: int, num_components: int, modulo: int):
    """
    Decompõe um polinômio em uma soma de polinômios cujos coeficientes 
    correspondem aos dígitos da representação dos coeficientes do polinômio
    original na base 'base', trabalhando em GF(modulo).

    Parâmetros:
      poly          : Poly
                      Polinômio de entrada (univariado)
      base          : int
                      Base usada para decomposição (ex: 2 para binário)
      num_components: int
                      Número de componentes (dígitos a serem extraídos)
      modulo        : int
                      Módulo usado para reduzir os coeficientes (ex: 281474972188673)

    Retorna:
      Uma lista (Iterable) de objetos Poly, onde o j-ésimo polinômio corresponde
      à parte dos coeficientes associada a base^j.
    """
    
    # Recria o polinômio no domínio do corpo finito GF(modulo)
    #poly_mod = Poly(poly.as_expr(), poly.gens, domain=GF(modulo))
    
    poly_mod = poly
    
    # Obtém um dicionário dos termos: chaves são tuplas de expoentes; valores são os coeficientes
    coeff_dict = poly_mod.as_dict()
    
    # Inicializa um dicionário para cada componente (cada "casa" na base)
    comp_coeffs = {j: {} for j in range(num_components)}
        
    # Para cada termo do polinômio, decompor o coeficiente na base fornecida
    for monom, coeff in coeff_dict.items():
        # Para polinômios univariados, o monômio é representado como uma tupla (expoente,)
        exp = monom[0]
        a = coeff  # coeficiente já em GF(modulo)
        # Extração dos dígitos
        for j in range(num_components):
            digit = a % base
            a = a // base
            # Armazena o dígito no dicionário do componente j para o monômio com expoente 'exp'
            comp_coeffs[j][(exp,)] = digit
        # Se ainda sobrar valor, significa que a decomposição exige mais dígitos do que 'num_components'
        if a != 0:
            raise ValueError(f"O coeficiente do termo x^{exp} exige mais de {num_components} dígitos na base {base}.")
    
    # Constrói os polinômios correspondentes a cada componente, mantendo o mesmo conjunto de variáveis
    polys = []
    for j in range(num_components):
        p = Poly.from_dict(comp_coeffs[j], poly.gens, domain=modulo)
        polys.append(p)
    
    return polys