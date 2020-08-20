'''Tests for the eighth exercise set.'''
import pytest
from cla_utils import arnoldi
from numpy import random
import numpy as np


@pytest.mark.parametrize('m, k', [(20, 4), (40, 20), (70, 13)])
def test_arnoldi(m, k):
    A = random.randn(m, m) + 1j*random.randn(m, m)
    b = random.randn(m) + 1j*random.randn(m)

    Q, H = arnoldi(A, b, k)
    assert(Q.shape == (m, k))
    assert(H.shape == (k+1, k))
    assert(np.abs(np.dot(Q, np.conj(Q).T) - np.eye(m)) < 1.0e-6)
    assert(np.linalg.norm(np.dot(A, Q) - np.dot(Q[:,:-1], H)) < 1.0e-6)

    
@pytest.mark.parametrize('m', [20, 204, 18])
def test_GMRES(m):
    A = random.randn(m, m)
    b = random.randn(m)

    x, _ = arnoldi(A, b, maxit=1000, tol=1.0e-3)
    assert(np.linalg.norm(np.dot(A, x) - b) < 1.0e-3)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
