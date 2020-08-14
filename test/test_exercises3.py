'''Tests for the third exercise set.'''
import pytest
from cla_utils import householder, householder_solve, householder_solve_qr
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 40, 87])
def test_householder(m):
    random.seed(1878*m)
    A = random.randn(m, m)
    A0 = 1.0*A  # make a deep copy
    R = householder(A0)
    assert(np.allclose(mat, np.triu(R))  # check R is upper triangular
    assert(np.linalg.norm(np.dot(R.T, R) - np.dot(A.T, A)) < 1.0e-6)


@pytest.mark.parametrize('m, k', [(20, 7), (40, 13), (87, 9)])
def test_householder_solve(m):
    random.seed(2432*m + 7438*k)
    A = random.randn(m, m)
    x0 = random.randn(m, k)
    b = np.dot(A, x0)
    x = householder_solve(A, b)
    assert(np.linalg.norm(x - x0) < 1.0e-6)


@pytest.mark.parametrize('m, n', [(20, 7), (40, 13), (87, 9)])
def test_householder_qr(m, n):
    random.seed(4732*m + 1238*k)
    A = random.randn(m, n)
    A0 = 1*A
    Q, R = householder_solve_qr(A0)

    # check orthogonality
    assert(linalg.norm(np.dot(Q.T, Q) - np.eye(m)) < 1.0e-6)
    # check upper triangular
    assert(np.allclose(mat, np.triu(R))
    # check QR factorisation
    assert(linalg.norm(np.dot(Q, R) - A) < 1.0e-6)


@pytest.mark.parametrize('m, n', [(20, 7), (40, 13), (87, 9)])
def test_householder_ls(m, n):
    random.seed(8473*m + 9283*k)
    A = random.randn(m, n)
    b = random.randn(m)
           
    x = householder_ls(A, b)

    #check normal equation residual
    assert( np.linalg.norm(np.dot(A.T, np.dot(A, x) - b) < 1.0e-6))


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
