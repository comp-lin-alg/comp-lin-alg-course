'''Tests for the first exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


# Test the basic matvec
@pytest.mark.parametrize('m, n', [(20, 20), (40, 20), (20, 45)])
def test_basic_matvec(m, n):
    random.seed(1878*m + 1950*n)
    A = random.randn(m, n)
    x = random.randn(n)

    b0 = A@x
    b = cla_utils.basic_matvec(A, x)

    assert(cla_utils.norm(b-b0) < 1.0e-6)


# Test the basic matvec
@pytest.mark.parametrize('m, n', [(20, 20), (40, 20), (20, 45)])
def test_column_matvec(m, n):
    random.seed(1878*m + 1950*n)
    A = random.randn(m, n)
    x = random.randn(n)

    b0 = A@x
    b = cla_utils.column_matvec(A, x)

    assert(cla_utils.norm(b-b0) < 1.0e-6)


@pytest.mark.parametrize('m, n', [(20, 20), (40, 20), (20, 45)])
def test_rank2_matrix(m, n):
    random.seed(1451*m + 1901*n)
    u1 = 1/np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))
    u2 = 1/np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))
    v1 = 1/np.sqrt(2)*(random.randn(n) + 1j*random.randn(n))
    v2 = 1/np.sqrt(2)*(random.randn(n) + 1j*random.randn(n))

    A = cla_utils.rank2(u1, u2, v1, v2)
    a1 = random.randn(m)
    a2 = random.randn(n)

    n1 = np.vdot(a1, A@a2)
    n2 = np.vdot(a1, u1)*np.vdot(v1, a2) + \
        np.vdot(a1, u2)*np.vdot(v2, a2)

    assert(np.abs(n1-n2) < 1.0e-7)


@pytest.mark.parametrize('m', [10, 20, 200])
def test_rank1pert_inv(m):
    random.seed(1001*m)
    u = 1/np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))
    v = 1/np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))

    A = np.eye(m) + np.outer(u, v.conj())
    Ainv = cla_utils.rank1pert_inv(u, v)

    x = 1/np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))

    y = A@x
    err = x - Ainv@y

    assert(cla_utils.norm(err)<1.0e-7)


@pytest.mark.parametrize('m', [3, 7, 20, 43])
def test_ABiC(m):
    random.seed(1348*m)
    B = random.randn(m, m)
    B_sym = B + B.T
    C = random.randn(m, m)
    C_ssym = C - C.T
    A = B_sym + 1j*C_ssym
    xr = random.randn(m)
    xi = random.randn(m)

    Ahat = B_sym
    Ahat[np.triu_indices(m, 1)] = C_ssym[np.triu_indices(m, 1)]
    zr, zi = cla_utils.ABiC(Ahat, xr, xi)

    z = zr + 1j*zi
    x = xr + 1j*xi
    err = z - A@x

    assert(cla_utils.norm(err) < 1.0e-7)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
