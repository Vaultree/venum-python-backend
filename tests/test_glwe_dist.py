from venum.glwe import EncryptionParameters, GlweDistribution, GlweSample

import pytest
from sympy import Poly
from sympy.abc import x


@pytest.fixture
def params():
    return EncryptionParameters(
        dimension=4,
        ciphertext_modulus=383,
        plaintext_modulus=127,
        noise_modulus=3,
        seed=0
    )


def test_poly_coeffs_within_bounds(params):
    dist = GlweDistribution(params)
    poly = dist.sample_polynomial()
    assert all(-params.ciphertext_modulus <= x <=
               params.ciphertext_modulus for x in poly.all_coeffs())


def test_poly_dimension_within_bounds(params):
    dist = GlweDistribution(params)
    poly = dist.sample_polynomial()
    assert poly.degree() == params.dimension - 1


@pytest.fixture
def zero_sample_polys():
    return {
        "mask": Poly(x**3 + 1, x),
        "secret": Poly(x**2 + 1, x),
        "crt_noise": Poly(x + 1, x),
        "expected_body": Poly(x**3 + x**2 + 2),
    }


def test_zero_sample(params, zero_sample_polys):
    dist = GlweDistribution(params)

    mask = zero_sample_polys["mask"].set_domain(dist.cipher_ring)
    secret = zero_sample_polys["secret"].set_domain(dist.cipher_ring)
    crt_noise = zero_sample_polys["crt_noise"].set_domain(dist.cipher_ring)
    expected_body = zero_sample_polys["expected_body"].set_domain(
        dist.cipher_ring)

    sample = GlweSample._compute_zero_sample(
        mask, secret, crt_noise, dist.poly_modulus
    )
    assert sample.mask == -mask
    assert sample.body == expected_body


@pytest.fixture
def sample_polys():
    """Fixture to provide test polynomials."""
    poly_modulus = Poly(x**3 + 1, x)  # Example modulus polynomial
    mask = Poly(x + 1, x)
    noise = Poly(2*x + 3, x)
    body = Poly(x**2 + 2, x)
    mask_noise = Poly(3*x + 2, x)
    body_noise = Poly(x + 4, x)
    message = Poly(2*x**2 + x + 1, x)
    u = Poly(x + 1, x)
    return {
        "poly_modulus": poly_modulus,
        "mask": mask,
        "noise": noise,
        "body": body,
        "mask_noise": mask_noise,
        "body_noise": body_noise,
        "message": message,
        "u": u
    }


def test_compute_mask(sample_polys):
    mask = sample_polys["mask"]
    noise = sample_polys["noise"]
    u = sample_polys["u"]
    poly_modulus = sample_polys["poly_modulus"]

    result = GlweSample._compute_mask(
        mask, noise,
        u, poly_modulus)
    expected = (mask * u + noise) % poly_modulus
    assert result == expected


def test_compute_body(sample_polys):
    body = sample_polys["body"]
    noise = sample_polys["body_noise"]
    message = sample_polys["message"]
    u = sample_polys["u"]
    poly_modulus = sample_polys["poly_modulus"]

    result = GlweSample._compute_body(
        body, noise,
        message, u, poly_modulus)
    expected = (body * u + message + noise) % poly_modulus
    assert result == expected


def test_compute_sample(sample_polys):
    mask = sample_polys["mask"]
    mask_noise = sample_polys["mask_noise"]
    body = sample_polys["body"]
    body_noise = sample_polys["body_noise"]
    message = sample_polys["message"]
    u = sample_polys["u"]
    poly_modulus = sample_polys["poly_modulus"]

    sample = GlweSample.compute_sample(
        mask, mask_noise, body, body_noise,
        message, u, poly_modulus
    )
    expected_mask = (mask * u + mask_noise) % poly_modulus
    expected_body = (body * u + message + body_noise) % poly_modulus
    assert sample.mask == expected_mask
    assert sample.body == expected_body
