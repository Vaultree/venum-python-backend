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
                dimension=4096,
                ciphertext_modulus=1400472361734830353,
                # ciphertext_modulus=281474972188673,
                # plaintext_modulus=65537,
                plaintext_modulus=12289,
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
    
    # envenenando os parametros
    params.dimension = 16
    params.ciphertext_modulus=1400472361734830353
    #params.ciphertext_modulus = gerar_primo(2**60, 2**61, params.dimension)
    params.plaintext_modulus = 65537
    total = 1000  # TOTAL DE TESTES
    chave_publica = True

    dist = GlweDistribution(params)

    sk, pk = gen_key_pair(dist)
    encryptor = Encryptor(dist, PolynomialEncoder(dist))
    
    print("SK...: ", sk)
    print("-" *  100)
    
    slot1 = []
    slot2 = []
    slot3 = []
    slot4 = []
    quant1 = 0
    quant2 = 0
    quant3 = 0
    quant4 = 0 
    quantidade_decifragens_corretas = 0
    
    slots = [0] * params.dimension
    quant = [0] * params.dimension
    
    for i in range(total):
        # print("ITERACAO: ", i)
            
        # rhs = vetor_aleatorio(params.dimension, 1, 16)
        # lhs = vetor_aleatorio(params.dimension, 1, 16)
        
        rhs = vetor_aleatorio(params.dimension, 1, params.plaintext_modulus-1)
        lhs = vetor_aleatorio(params.dimension, 1, params.plaintext_modulus-1)        
        
        # lhs = gerar_vetor_especial(params.dimension, i % params.dimension, params.plaintext_modulus-1)
        
        n1 = random.randint(0, params.dimension-1)
        n2 = random.randint(0, params.dimension-1)
        n3 = random.randint(0, params.dimension-1)
        n4 = random.randint(0, params.dimension-1)
        
        # lhs = gerar_vetor_especial_multi(params.dimension, [n1, n2, n3, n4], params.plaintext_modulus-1)

        #rhs = vetor_aleatorio(4, 1, 64)
        #lhs = vetor_aleatorio(4, 1, 64)
         
        number = random.randint(1, params.plaintext_modulus-1)
        
        # lhs = [number,0,number,0]    
        # lhs = [0,number,0,0]  
        # lhs = [0,0,number,0] 
        # lhs = [0,0,0,number]        
        # lhs = [10,20,30,40]
        
        # checking...
        expected = (Poly(reversed(lhs), x, domain=dist.plaintext_ring) *
                Poly(reversed(rhs), x, domain=dist.plaintext_ring) %
                dist.poly_modulus.set_domain(dist.plaintext_ring))
        expected = list(reversed(expected.all_coeffs()))
        
        # filling with zeros
        while True:
            if len(expected) < params.dimension:
                expected.append(0)
            else:
                break
            
        # print("VALORES..: ", lhs, rhs)
        # print("ESPERADO.: ", expected)
        
        #generate_formula(lhs, rhs, params.plaintext_modulus)
        #sys.exit(1)
        
        # lhs_cipher = encryptor.encrypt(sk, lhs)
        # rhs_cipher = encryptor.encrypt(sk, rhs)

        if chave_publica:
            print("ATENÇÃO: UTILIZANDO CHAVE PÚBLICA ********************************************")
            lhs_cipher = encryptor.encrypt(pk, lhs)
            rhs_cipher = encryptor.encrypt(pk, rhs)
        else:
            print("ATENÇÃO: UTILIZANDO CHAVE PRIVADA ********************************************")
            lhs_cipher = encryptor.encrypt(sk, lhs)
            rhs_cipher = encryptor.encrypt(sk, rhs)
        
        d1 = encryptor.decrypt(sk, lhs_cipher)
        d2 = encryptor.decrypt(sk, rhs_cipher)
        
        # verificando a decifragem
        # print("VALORES ENTRADA.....: ", lhs, rhs)
        # print("VALORES DECIFRADOS..: ", d1, d2)
        # print("Produto esperado....: ", expected)
        #generate_formula(lhs, rhs, params.plaintext_modulus)
        assert d1 == lhs
        assert d2 == rhs
        
        # print("c1 mask: ", lhs_cipher.glwe_sample.mask)
        # print("c1 body: ", lhs_cipher.glwe_sample.body)

        # print("c2 mask: ", rhs_cipher.glwe_sample.mask)
        # print("c2 body: ", rhs_cipher.glwe_sample.body)

        relin_key = RelinKey.from_secret_key_crt(sk)
        eval = Evaluator(dist, relin_key)
        cipher_result = eval.mul(lhs_cipher, rhs_cipher)
        
        # print("c3 mask: ", cipher_result.glwe_sample.mask)
        # print("c3 body: ", cipher_result.glwe_sample.body)
        
        decrypted = encryptor.decrypt(sk, cipher_result)
        
        print("RESULTADO.: ", decrypted)
        print("ESPERADO..: ", expected)
        print("LHS.......: ", lhs)
        print("RHS.......: ", rhs)
        print("FÓRMULA WOLFRAM ALPHA: (VÁLIDO SOMENTE PARA 4 DIMENSÕES) ")
        generate_formula(lhs, rhs, params.plaintext_modulus)
        
        if decrypted == expected:
            quantidade_decifragens_corretas += 1
            
        # if decrypted != expected:
        #     sys.exit(1)

        # verificando diferenças
        t = vector_difference(decrypted, expected)
        if t != [0] * params.dimension:
            print("DIFERENCA: ", t)
        else:
            print("Nenhuma diferença encontrada.")  
        
        if params.dimension == 4:
            if t[0] != 0:
                quant1 += 1
            if t[1] != 0:
                quant2 += 1
            if t[2] != 0:            
                quant3 += 1
            if t[3] != 0:
                quant4 += 1
                
            for i in range(4):
                while t[i] < 0:
                    t[i] += params.plaintext_modulus
            
            if t[0] not in slot1 and t[0] != 0:
                slot1.append(t[0])
            if t[1] not in slot2 and t[1] != 0:
                slot2.append(t[1])
            if t[2] not in slot3 and t[2] != 0:
                slot3.append(t[2])
            if t[3] not in slot4 and t[3] != 0:
                slot4.append(t[3])
        else:
            limit = params.dimension
            for i in range(limit):
                if t[i] != 0:
                    quant[i] += 1
                    
                for i in range(limit):
                    while t[i] < 0:
                        t[i] += params.plaintext_modulus
                        
                if t[i] not in slots and t[i] != 0:
                    slots.append(t[i])
            #sys.exit(1)
            
        #assert decrypted == expected
        # print("--------------------------------------------------")
        print("-" * 50)
    
    if params.dimension == 4:
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
    else:
        if slots != [0] * params.dimension:
            print("slots [Diferenças distintas verificadas]: ", slots)
        else:
            print("Nenhuma diferença no slot foi encontrado.")
        
        if quant != [0] * params.dimension:
            print("Quantidade de erros: ", quant)
        else:
            print("Nenhum erro encontrado.")    
        print("Total de testes: ", total)
        print("Quantidade de decifragens corretas: ", quantidade_decifragens_corretas
            , " de ", total)    
    print("=" * 100)
    
    print("DIMENSÃO............: ", params.dimension)
    print("MÓDULO CIPHERTEXT...: ", params.ciphertext_modulus)
    print("MÓDULO DO CLEARTEXT.: ", params.plaintext_modulus)
    # sys.exit(1)


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
    
    print(formula)

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

import random

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

def gerar_primo(inicio, fim, n):
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

def gerar_vetor_especial(n: int, pos: int, l: int) -> list[int]:
    """
    Gera um vetor (lista) de dimensão n com todos os valores iguais a 0,
    exceto um único valor na posição especificada por 'pos'. Nesse elemento,
    o valor é gerado aleatoriamente no intervalo [0, l).

    Parâmetros:
        n (int): Dimensão do vetor.
        pos (int): Posição em que o valor não-zero será inserido.
        l (int): Limite superior para o valor gerado (valor máximo possível é l-1).

    Retorna:
        list[int]: Vetor gerado.

    Levanta:
        ValueError: Se a posição 'pos' não estiver entre 0 e n-1.
    """
    if pos < 0 or pos >= n:
        raise ValueError("A posição especificada está fora do intervalo do vetor.")

    vetor = [0] * n
    vetor[pos] = random.randrange(l)
    return vetor

def gerar_vetor_especial_multi(n: int, posicoes: list[int], l: int) -> list[int]:
    """
    Gera um vetor (lista) de dimensão n com todos os valores iguais a 0,
    exceto nas posições indicadas em 'posicoes'. Em cada posição especificada,
    o valor é gerado aleatoriamente no intervalo [0, l).
    
    Parâmetros:
        n (int): Dimensão do vetor.
        posicoes (list[int]): Lista com as posições onde o vetor terá valores não nulos.
        l (int): Limite superior para os valores gerados (valor máximo possível é l-1).
    
    Retorna:
        list[int]: Vetor gerado.
    
    Levanta:
        ValueError: Se alguma posição em 'posicoes' estiver fora do intervalo [0, n-1].
    """
    vetor = [0] * n
    
    for pos in posicoes:
        if pos < 0 or pos >= n:
            raise ValueError(f"A posição {pos} está fora do intervalo do vetor de dimensão {n}.")
        vetor[pos] = random.randrange(l)
    
    return vetor