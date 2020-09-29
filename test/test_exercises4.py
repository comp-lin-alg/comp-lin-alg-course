'''Tests for the fourth exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m, n', [(20, 7), (40, 13), (9, 87)])
def test_operator_2_norm(m, n):
    random.seed(8473*m + 9283*n)
    A = random.randn(m, n)

    norm1 = cla_utils.operator_2_norm(A)
    u, s, v = np.linalg.svd(A)
    norm2 = s[0]

    assert(np.abs(norm1 - norm2) < 1.0e-6)

@pytest.mark.parametrize('m', [20, 40, 233])
def test_cond(m):
    random.seed(3168*m)
    A = random.randn(m, m)

    cond1 = cla_utils.cond(A)
    cond2 = np.linalg.cond(A)

    assert(np.abs(cond1 - cond2) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
