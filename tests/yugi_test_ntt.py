import sys
import time
from tests.yugi_test_mul import generate_formula, gerar_primo
from venum.ntt import find_primitive_root, generate_parameters, modular_inverse, multiply_poly_mod, polymul_ntt
import pytest
import random

@pytest.mark.parametrize(
    "input",
    [
        {
            "modulus": 1400472361734830353,
            "dimension": 1024,
            "number": 162577
        }
    ])
def test_ntt(input):
    modulus, dimension = input["modulus"], input["dimension"]
    
    # calculate the primitive root of the modulus
    root = find_primitive_root(modulus, dimension)
    print("\nPrimitive root: ", root)
    
    # check if the modular inverse exists
    number = input["number"]
    inverse = modular_inverse(number, modulus)
    print("\nModular inverse: ", inverse)
    assert number * inverse % modulus == 1

    # check if the modular inverse exists
    for ct in range(10):
        number = random.randint(0, modulus - 1)
        inverse = modular_inverse(number, modulus)
        print("Modular inverse: ", inverse)
        assert number * inverse % modulus == 1
        
    # --------------------------------------------------------------
    # Test multiply_poly_mod
    a = [11, 12, 13, 14]
    b = [47, 36, 25, 1027]
    # Módulo q
    q = 65537
    resultado = multiply_poly_mod(a, b, q)
    print("Result of multiplying polynomials: ", resultado)
    generate_formula(a, b, q)
    
    # --------------------------------------------------------------
    # Testing the code to multiply using NTT/INTT
    
    print("-" * 80)
    n = 4  # Dimension of the polynomial
    q = 12289
    
    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    
    # Example of n-dimensional polynomials
    # Polynomial A(x) = 1 + 2x + 3x^2 + 4x^3 ...
    # Polynomial B(x) = 8 + 7x + 6x^2 + 5x^3 ...
    a = generate_random_vector(n, 0, q - 1)
    b = generate_random_vector(n, 0, q - 1) 

    print("Polynomial A:", a)
    print("Polynomial B:", b)
    print("q...........:", q)
    generate_formula(a, b, q)
    
    # Multiplication of polynomials using NTT/INTT
    result = polymul_ntt(a, b, q, psi_rev, psi_inv_rev, n_inv, bar)
    print("Result of multiplication: ", result)
    
    # check if the result is correct
    assert multiply_poly_mod(a, b, q) == result
    
    # --------------------------------------------------------------
    # Test the multiplication of polynomials using NTT/INTT
    
    n = 4
    while True:
        limit = 10
        if n > 2048:
            limit = 1
            
        for _ in range(limit):
            print("-" * 80)
            q = gerar_primo(2**62, 2**63, n)
            
            # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
            psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
            
            # Example of n-dimensional polynomials
            # Polynomial A(x) = 1 + 2x + 3x^2 + 4x^3 ...
            # Polynomial B(x) = 8 + 7x + 6x^2 + 5x^3 ...
            
            # max_point = 2**16
            max_point = q-1
            a = generate_random_vector(n, 0, max_point)
            b = generate_random_vector(n, 0, max_point) 

            #print("Polynomial A:", a)
            #print("Polynomial B:", b)
            print("q...........:", q)
            print("Dimensions..:", n)
            
            # Multiplication of polynomials using NTT/INTT
            result = polymul_ntt(a, b, q, psi_rev, psi_inv_rev, n_inv, bar)
            # print("Result of multiplication: ", result)
            
            # check if the result is correct
            assert multiply_poly_mod(a, b, q) == result
            
        n = n * 2  # Dimension of the polynomial
        if n > 16384:
            break
        
    # Test Timing:
    print("-" * 80)
    n = 1024  # Dimension of the polynomial
    q = gerar_primo(2**62, 2**63, n)
    
    # Generate the required parameters (vectors psi_rev, psi_inv_rev, n_inv and the Barrett structure)
    psi_rev, psi_inv_rev, n_inv, bar = generate_parameters(n, q)
    
    max_point = q-1
    a = generate_random_vector(n, 0, max_point)
    b = generate_random_vector(n, 0, max_point) 

    # Multiplication of polynomials using NTT/INTT
    total_rounds= 1000
    start_time = time.perf_counter()
    for _ in range(total_rounds):
        result = polymul_ntt(a, b, q, psi_rev, psi_inv_rev, n_inv, bar)
        
    end_time = time.perf_counter()
    elapsed_time = (end_time - start_time) / total_rounds

    # print("Polinômio resultante:", result)
    print("Execution time: {:.6f} seconds".format(elapsed_time))

# --------------------------------------------------------------
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
