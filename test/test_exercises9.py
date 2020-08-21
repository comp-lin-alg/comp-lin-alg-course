'''Tests for the eighth exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_pow_it(m):
    random.seed(1302*m)
    A = random.randn(m, m)
    x0 = random.randn(m)
    xi = cla_utils.pow_it(A, x0, tol=1.0e-4, maxit=10000)
    r = np.dot(A, xi)
    r /= np.linalg.norm(r)
    r -= xi/np.linalg.norm(xi)
    assert(np.linalg.norm(r) < 1.0e-4)


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
    r /= np.linalg.norm(r)
    r -= xi/np.linalg.norm(xi)
    assert(np.linalg.norm(r) < 1.0e-4)


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
    r /= np.linalg.norm(r)
    r -= xi/np.linalg.norm(xi)
    assert(np.linalg.norm(r) < 1.0e-4)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_pure_QR(m):
    random.seed(1302*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A0 = 1.0*A
    A2 = cla_utils.pure_QR(QR, maxit=10, tol=1.0e-100)
    #check it is still Hermitian
    assert(np.linalg.norm(A2 - np.conj(A2).T) < 1.0e-4)
    #check for orthogonality
    x0 = random.randn(m)
    assert(np.linalg.norm(np.dot(A2, x0))-np.linalg(np.dot(A0, x0)) < 1.0e-6)
    #check for conservation of trace
    assert(np.abs(np.tr(A0) - np.tr(A2)) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
