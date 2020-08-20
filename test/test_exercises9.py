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


@pytest.mark.parametrize('m', [20, 204, 18])
def test_inverse_it(m):
    random.seed(1302*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A = 0.5*(A + np.conj(A).T)
    e, _ = np.linalg.eig(A)
    x0 = random.randn(m)
    mu = e[m//2] + random.randn() + 1j*random.randn()
    xi, li = inverse_it(A, x0, mu, tol=1.0e-8, maxit=10000)
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
    xi, li = rq_it(A, x0, tol=1.0e-8, maxit=10000)
    r = np.dot(A, xi)
    r /= np.linalg.norm(r)
    r -= xi/np.linalg.norm(xi)
    assert(np.linalg.norm(r) < 1.0e-4)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
