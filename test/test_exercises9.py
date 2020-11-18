'''Tests for the nine exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_pow_it(m):
    random.seed(1302*m)
    A = random.randn(m, m)
    A = 0.5*(A + A.T)
    x0 = random.randn(m)
    xi, lambda0 = cla_utils.pow_it(A, x0, tol=1.0e-6, maxit=10000)
    assert(np.linalg.norm(A@xi-lambda0*xi) < 1.0e-3)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_inverse_it(m):
    random.seed(1302*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A = 0.5*(A + np.conj(A).T)
    e, _ = np.linalg.eig(A)
    x0 = random.randn(m)
    mu = e[m//2] + random.randn() + 1j*random.randn()
    xi, li = cla_utils.inverse_it(A, x0, mu, tol=1.0e-8, maxit=10000)
    es = np.abs(e - mu)
    i1 = np.argsort(es)
    ll = e[i1[0]]
    assert(np.abs(ll - li) < 1.0e-6)
    r = np.dot(A, xi)
    assert(np.linalg.norm(r - li*xi) < 1.0e-4)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_rq_it(m):
    random.seed(1302*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A = 0.5*(A + np.conj(A).T)
    e, _ = np.linalg.eig(A)
    x0 = random.randn(m)
    mu = e[m//2] + random.randn() + 1j*random.randn()
    xi, li = cla_utils.rq_it(A, x0, tol=1.0e-8, maxit=10000)
    r = np.dot(A, xi)
    assert(np.linalg.norm(r - li*xi) < 1.0e-4)


@pytest.mark.parametrize('m', [20, 30, 18])
def test_pure_QR(m):
    random.seed(1302*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A = 0.5*(A + A.conj().T)
    A0 = 1.0*A
    A2 = cla_utils.pure_QR(A0, maxit=10000, tol=1.0e-5)
    #check it is still Hermitian
    assert(np.linalg.norm(A2 - np.conj(A2).T) < 1.0e-4)
    #check for upper triangular
    assert(np.linalg.norm(A2[np.tril_indices(m, -1)])/m**2 < 1.0e-5)
    #check for conservation of trace
    assert(np.abs(np.trace(A0) - np.trace(A2)) < 1.0e-6)

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
