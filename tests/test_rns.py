import pytest
from venum.rns import RnsBasis


@pytest.fixture
def rns_basis():
    return RnsBasis([3, 5, 7])


def test_initialization(rns_basis):
    rns_number = rns_basis.to_rns(10)
    expected = [10 % m for m in rns_basis.moduli]
    assert rns_number.residues == expected, "RNS initialization failed."


def test_addition(rns_basis):
    x = 10
    y = 6
    a = rns_basis.to_rns(x)
    b = rns_basis.to_rns(y)
    result = a + b
    expected = [(x + y) % m for m in rns_basis.moduli]
    assert result.residues == expected, \
        f"Expected addition result {expected}, got {result.residues}"


def test_subtraction(rns_basis):
    x = 10
    y = 6
    a = rns_basis.to_rns(x)
    b = rns_basis.to_rns(y)
    result = a - b
    expected = [(x - y) % m for m in rns_basis.moduli]
    assert result.residues == expected, \
        f"Expected subtraction result {expected}, got {result.residues}"


def test_multiplication(rns_basis):
    x = 10
    y = 6
    a = rns_basis.to_rns(x)
    b = rns_basis.to_rns(y)
    result = a * b
    expected = [(x * y) % m for m in rns_basis.moduli]
    assert result.residues == expected, \
        f"Expected multiplication result {expected}, got {result.residues}"


def test_to_int(rns_basis):
    x = 10
    a = rns_basis.to_rns(x)
    assert a.to_int() == x, f"Expected int value {x}, got {int(a)}"
