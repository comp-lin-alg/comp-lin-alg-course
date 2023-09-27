'''Tests for the tenth exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m, k', [(20, 4), (40, 20), (70, 13)])
def test_arnoldi(m, k):
    A = random.randn(m, m) + 1j*random.randn(m, m)
    b = random.randn(m) + 1j*random.randn(m)

    Q, H = cla_utils.arnoldi(A, b, k)
    assert(Q.shape == (m, k+1))
    assert(H.shape == (k+1, k))
    assert(cla_utils.norm((Q.conj().T)@Q - np.eye(k+1)) < 1.0e-6)
    assert(cla_utils.norm(A@Q[:,:-1] - Q@H) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_GMRES(m):
    A = random.randn(m, m)
    b = random.randn(m)

    x, _ = cla_utils.GMRES(A, b, maxit=1000, tol=1.0e-3)
    assert(cla_utils.norm(np.dot(A, x) - b) < 1.0e-3)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
