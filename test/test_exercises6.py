'''Tests for the fifth exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np
from cla_utils import LU_inplace, get_Lk, solve_L, solve_U


@pytest.mark.parametrize('m, k', [(20, 4), (204, 100), (18, 7)])
def test_get_Lk(m, k):
    random.seed(9752*m)
    lk = random.randn(k)
    Lk = cla_utils.get_Lk(m, lk)
    assert(np.count_nonzero(Lk) == m + k - 1)

    b = random.randn(m)
    x = np.dot(Lk, b)
    assert(np.abs(x[0:k-1]-b[0:k-1]) < 1.0e-6)
    for i in range(k,m):
        assert(np.abs(x[i] - b[i] - np.dot(b[i+1:], lk)) < 1.0e-6)


@pytest.mark.parametrize('m', [3, 20, 204, 18])
def test_LU_inplace(m):
    random.seed(8564*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    cla_utils.LU_inplace(A)
    L = np.eye(m)
    i1 = np.tril_indices(m, k=-1)
    L[i1] = A[i1]
    #print(L)
    U = np.triu(A)
    A1 = np.dot(L, U)
    print(A1)
    err = (A1 - A0)
    assert(np.linalg.norm(err) < 1.0e-6)


@pytest.mark.parametrize('m, k', [(20, 4), (204, 100), (18, 7)])
def test_solve_L(m, k):
    random.seed(1002*m + 2987*k)
    b = random.randn(m, k)
    L = np.tril(random.randn(m, m))
    x = cla_utils.solve_L(L, b)
    err1 = b - np.dot(L, x)
    assert(np.linalg.norm(err1) < 1.0e-6)
    A = random.randn(m, m)
    x = solve_L(A, b)
    err2 = b - np.dot(A, x)
    assert(np.linalg.norm(err2) < 1.0e-6)


@pytest.mark.parametrize('m, k', [(20, 4), (204, 100), (18, 7)])
def test_solve_U(m, k):
    random.seed(1002*m + 2987*k)
    b = random.randn(m, k)
    U = np.triu(random.randn(m, m))
    x = cla_utils.solve_U(U, b)
    err1 = b - np.dot(U, x)
    assert(np.linalg.norm(err1) < 1.0e-6)
    A = random.randn(m, m)
    err2 = b - np.dot(A, x)
    assert(np.linalg.norm(err2) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_inverse_LU(m):
    random.seed(5422*m)
    A = random.randn(m, m)
    A0 = 1.0*A

    Ainv = cla_utils.inverse_LU(A0)
    err = np.dot(Ainv, A) - np.eye(m)
    assert(np.linalg.norm(err) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
