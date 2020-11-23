'''Tests for the seventh exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_perm(m):
    random.seed(4533*m)
    p = random.permutation(np.arange(m))
    for k in range(20):
        i, j = random.choice(np.arange(m), 2, replace=False)
        p0 = 1.0*p
        x = random.randn(m)
        cla_utils.perm(p0, i, j)
        assert(np.abs(x[int(p0[i])]-x[int(p[j])]) < 1.0e-6)
        assert(np.abs(x[int(p0[j])]-x[int(p[i])]) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_LUP_inplace(m):
    random.seed(8364*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    p = cla_utils.LUP_inplace(A)
    L = np.eye(m)
    i1 = np.tril_indices(m, k=-1)
    L[i1] = A[i1]
    U = np.triu(A)
    A1 = np.dot(L, U)
    A0 = A0[p, :]
    err = A1 - A0
    assert(np.linalg.norm(err) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_solve_LUP(m):
    random.seed(8364*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    b = random.randn(m)
    x = cla_utils.solve_LUP(A, b)
    assert(np.linalg.norm(b - np.dot(A0, x)) < 1.0e-6)


@pytest.mark.parametrize('m', [3, 9, 18])
def test_det_LUP(m):
    random.seed(1477*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    detA = cla_utils.det_LUP(A)
    _, s, _ = np.linalg.svd(A0)
    assert(np.linalg.norm(detA - np.prod(s)) < 1.0e-6)

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
