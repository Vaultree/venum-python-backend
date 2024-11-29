from venum import numeric

import pytest


@pytest.mark.parametrize("number, radix, expected", [
    # radix 2
    (0, 2, [0]),
    (1, 2, [1]),
    (10, 2, [1, 0, 1, 0]),
    (100, 2, [1, 1, 0, 0, 1, 0, 0]),
    (10000, 2, [1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0]),
    (1234512345, 2, [1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1,
     0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1]),

    # radix 10
    (0, 10, [0]),
    (1, 10, [1]),
    (10, 10, [1, 0]),
    (100, 10, [1, 0, 0]),
    (10000, 10, [1, 0, 0, 0, 0]),
    (1234512345, 10, [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]),
])
def test_radix_decomposition(number, radix, expected):
    for i in range(len(expected)):
        actual = numeric.nth_digit(number, radix, i)
        assert actual == expected[-1 - i]
