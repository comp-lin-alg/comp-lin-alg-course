'''Tests for the eighth exercise set.'''
import pytest
from cla_utils import pow_it
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_pow_it(m):
    random.seed(1302*m)
    A = random.randn(m, m)
    x0 = random.randn(m)
    xi = pow_it(A, x0, tol=1.0e-4, maxit=10000)
    r = np.dot(A, xi)
    r /= np.linalg.norm(r)
    r -= xi/np.linalg.norm(xi)
    assert(np.linalg.norm(r) < 1.0e-4)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
