import sys
import time
from tests.yugi_test_full import show_title
from tests.yugi_test_mul import generate_formula
from venum.ntt import *
import pytest
import random
from venum.ntt_func import *

@pytest.mark.parametrize(
    "input",
    [
        {
        }
    ])
def test_scheme_rotation(input):
    print("\n")
    show_title("TEST SCHEME - ROTATION")

    n = 4
    # q = generate_modulus(2**30, 2**32, n)
    q=585316169
    p1 = 257
    p2 = 3
    base_decomposition = 2
    batched = True
    
    # verify if p1 is equal to 1 mod 2n (necessary condition for CRT)
    if p1 % (2*n) != 1:
        print("Error: p1 must be equal to 1 mod 2n | n = ",n, "p1 mod 2n", p1 % (2*n))
        sys.exit(1)
    
    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    params = (psi_rev, psi_inv_rev, n_inv, bar)
    
    # To use batched mode it is necessary to calculate psi_rev, psi_inv_rev, n_inv and bar for module p1
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, p1)
    params_batched = (psi_rev, psi_inv_rev, n_inv, bar)
    
    pot = []
    for ct in range(0, 256):
        tmp = 64**ct % p1
        if pot.__contains__(tmp) == False:
            pot.append(tmp)
    pot.sort()
    
    print("pot: ", pot)
    print("PSI_REV:", psi_rev)
    print("PSI_INV_REV:", psi_inv_rev)
    
    
    # Generate the keys
    # sk = create_sk(n, 0, 1, q)
    sk = [1,0,1,0]
    print("Secret key: ", sk)
    
    pk = generate_pk(sk, q, p1, p2, params)
    rlk = generate_rlk(sk, base_decomposition, p1, p2, q, params)
    
    # Generate the plaintexts
    m0 = [1,2,3,4]
    # m0 = [1,0,0,0]
    print("Plaintext: ", m0)
    
    print("Parameters: ", params)
    print("Parameters batched: ", params_batched)
    
    # Encrypt the plaintexts secret key
    ciphertext = encrypt_sk(sk, m0, q, p1, p2, batched, params, params_batched)
    print("Ciphertext: ", ciphertext)
    
    # Rotate the ciphertext
    r = 1
    ciphertext_rot = rotate_cyphertext(ciphertext, r)
    print("Ciphertext rotated: ", ciphertext_rot)
    
    sk_rot = rotate_left(sk, r)
    decrypted = decrypt(sk_rot, ciphertext_rot, p1, params, params_batched)
    print("Decrypted rotated (sk_rot)..: ", decrypted)

    decrypted = decrypt(sk, ciphertext_rot, p1, params, params_batched)
    print("Decrypted rotated (sk)......: ", decrypted)
    
    decrypted = decrypt(sk, ciphertext, p1, params, params_batched)
    print("Decrypted (normal): ", decrypted)
    print("*" * 60)
    
    # ********************************

# ----------------------------------------------------------------------------------
def rotate_cyphertext(ciphertext: Cryptogram, n):
    """
    Rotaciona o vetor ciphertext para a esquerda por n posições.
    
    Parâmetros:
      ciphertext (list): Lista a ser rotacionada.
      n (int): Número de posições para rotacionar.
      
    Retorna:
      list: O vetor rotacionado.
    """
    a = ciphertext.mask
    b = ciphertext.body
    
    # Rotaciona o vetor ciphertext para a esquerda por n posições
    a_rot = rotate_left(a, n)
    b_rot = rotate_left(b, n)
    
    # Cria um novo objeto Cryptogram com os vetores rotacionados
    ciphertext_rot = Cryptogram(body=b_rot, mask=a_rot, batched=ciphertext.batched, q=ciphertext.q)
    return ciphertext_rot
    
# ----------------------------------------------------------------------------------
def rotate_left(vector, n):
    """
    Rotaciona o vetor para a esquerda por n posições.
    
    Parâmetros:
      vector (list): Lista a ser rotacionada.
      n (int): Número de posições para rotacionar.
      
    Retorna:
      list: O vetor rotacionado.
    """
    if not vector:  # Verifica se a lista está vazia
        return vector
    n = n % len(vector)  # Garante que n está dentro do tamanho da lista
    
    # return [vector[1], vector[0], vector[3], vector[2]]
    
    return vector[n:] + vector[:n]
        

