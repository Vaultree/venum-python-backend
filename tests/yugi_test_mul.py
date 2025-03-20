import sys
from tests.yugi_test_enc_sk import vetor_aleatorio
from venum.glwe import EncryptionParameters, GlweDistribution
from venum.encryption import Encryptor
from venum.plaintext_encoding import PolynomialEncoder
from venum.key import gen_key_pair, RelinKey
from venum.evaluation import Evaluator

from sympy import Poly
from sympy.abc import x
import pytest
import random


# @pytest.mark.skip(reason="multiplication needs fixing")
@pytest.mark.parametrize(
    "input",
    [
        {
            "params": EncryptionParameters(
                dimension=4,
                # ciphertext_modulus=1400472361734830353,
                ciphertext_modulus=281474972188673,
                plaintext_modulus=65537,
                noise_modulus=3,
                #seed=1,
            ),
            "lhs": [0, 0, 0, 0],
            "rhs": [0, 0, 0, 0],
            # [xˆ0, x^1, x^2, x^3]
        }
    ])
def test_mul(input):
    params, lhs, rhs = input["params"], input["lhs"], input["rhs"]

    dist = GlweDistribution(params)

    sk, pk = gen_key_pair(dist)
    encryptor = Encryptor(dist, PolynomialEncoder(dist))
    
    print("SK...: ", sk)
    
    slot1 = []
    slot2 = []
    slot3 = []
    slot4 = []
    quant1 = 0
    quant2 = 0
    quant3 = 0
    quant4 = 0
    total = 100
    quantidade_decifragens_corretas = 0
    
    for i in range(total):
        # print("ITERACAO: ", i)
            
        # lhs = vetor_aleatorio(4, 1, 10)
        rhs = vetor_aleatorio(4, 1, params.plaintext_modulus-1)
        lhs = vetor_aleatorio(4, 1, params.plaintext_modulus-1)

        rhs = vetor_aleatorio(4, 1, 64)
        lhs = vetor_aleatorio(4, 1, 64)
         
        number = random.randint(1, params.plaintext_modulus-1)
        
        # lhs = [number,0,0,0]    
        # lhs = [0,number,0,0]  
        # lhs = [0,0,number,0] 
        # lhs = [0,0,0,number]        
        # lhs = [10,20,30,40]
        
        # verificando... 
        expected = (Poly(reversed(lhs), x, domain=dist.plaintext_ring) *
                Poly(reversed(rhs), x, domain=dist.plaintext_ring) %
                dist.poly_modulus.set_domain(dist.plaintext_ring))
        expected = list(reversed(expected.all_coeffs()))
        
        print("VALORES..: ", lhs, rhs)
        print("ESPERADO.: ", expected)
        
        generate_formula(lhs, rhs, params.plaintext_modulus)
        #sys.exit(1)
        
        lhs_cipher = encryptor.encrypt(sk, lhs)
        rhs_cipher = encryptor.encrypt(sk, rhs)
        
        d1 = encryptor.decrypt(sk, lhs_cipher)
        d2 = encryptor.decrypt(sk, rhs_cipher)
        
        # verificando a decifragem
        print("VALORES ENTRADA.....: ", lhs, rhs)
        print("VALORES DECIFRADOS..: ", d1, d2)
        print("Produto esperado....: ", expected)
        generate_formula(lhs, rhs, params.plaintext_modulus)
        assert d1 == lhs
        assert d2 == rhs
        
        print("c1 mask: ", lhs_cipher.glwe_sample.mask)
        print("c1 body: ", lhs_cipher.glwe_sample.body)

        print("c2 mask: ", rhs_cipher.glwe_sample.mask)
        print("c2 body: ", rhs_cipher.glwe_sample.body)

        relin_key = RelinKey.from_secret_key_bfv(sk)
        eval = Evaluator(dist, relin_key)
        cipher_result = eval.mul(lhs_cipher, rhs_cipher)
        
        print("c3 mask: ", cipher_result.glwe_sample.mask)
        print("c3 body: ", cipher_result.glwe_sample.body)
        
        decrypted = encryptor.decrypt(sk, cipher_result)
        
        generate_formula(lhs, rhs, params.plaintext_modulus)
        print("RESULTADO: ", decrypted)
        print("ESPERADO.: ", expected)
        
        t = vector_difference(decrypted, expected)
        print("DIFERENCA: ", t)
        
        if t[0] > 0:
            quant1 += 1
        if t[1] > 0:
            quant2 += 1
        if t[2] > 0:            
            quant3 += 1
        if t[3] > 0:
            quant4 += 1
        
        if t[0] not in slot1 and t[0] > 0:
            slot1.append(t[0])
        if t[1] not in slot2 and t[1] > 0:
            slot2.append(t[1])
        if t[2] not in slot3 and t[2] > 0:
            slot3.append(t[2])
        if t[3] not in slot4 and t[3] > 0:
            slot4.append(t[3])
        
        if decrypted == expected:
            quantidade_decifragens_corretas += 1
            #sys.exit(1)
            
        #assert decrypted == expected
        print("--------------------------------------------------")
    
    print("slot1 [Diferenças distintas verificadas]: ", slot1)
    print("slot2 [Diferenças distintas verificadas]: ", slot2)
    print("slot3 [Diferenças distintas verificadas]: ", slot3)
    print("slot4 [Diferenças distintas verificadas]: ", slot4)
    print("Quantidade de erros no slot1: ", quant1)
    print("Quantidade de erros no slot2: ", quant2)
    print("Quantidade de erros no slot3: ", quant3)
    print("Quantidade de erros no slot4: ", quant4)
    print("Total de testes: ", total)
    print("Quantidade de decifragens corretas: ", quantidade_decifragens_corretas
          , " de ", total)
    print("--------------------------------------------------")
    sys.exit(1)


def generate_formula(a, b, q):
    # String original da fórmula
    formula = "PolynomialMod[PolynomialMod[(A*x^0 + B*x^1 + C*x^2 + D*x^3) * (E*x^0 + F*x^1 + G*x^2 + H*x^3), x^4+1], q]"
    
    # Substitui os coeficientes do primeiro polinômio (A, B, C, D)
    for letter, pos in zip(["A", "B", "C", "D"], range(4)):
        formula = formula.replace(letter, str(a[pos]))
    
    # Substitui os coeficientes do segundo polinômio (E, F, G, H)
    for letter, pos in zip(["E", "F", "G", "H"], range(4)):
        formula = formula.replace(letter, str(b[pos]))
    
    # Substitui o módulo 'q'
    formula = formula.replace("q", str(q))
    
    print("Formula Wolfram Alpha:", formula)

def vector_difference(v1, v2):
    """
    Recebe dois vetores (listas) de números e retorna um novo vetor contendo
    a diferença elemento a elemento: v1[i] - v2[i].

    Parâmetros:
      v1: list
          Primeiro vetor.
      v2: list
          Segundo vetor, que deve ter o mesmo tamanho de v1.

    Retorna:
      list: vetor resultante da subtração de cada elemento de v1 pelo elemento correspondente em v2.

    Lança:
      ValueError se os vetores tiverem tamanhos diferentes.
    """
    if len(v1) != len(v2):
        print(v1)
        print(v2)
        raise ValueError("Os vetores devem ter o mesmo tamanho.")
    
    return [a - b for a, b in zip(v1, v2)]