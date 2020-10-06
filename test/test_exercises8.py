'''Tests for the eighth exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_Q1AQ1s(m):
    random.seed(4373*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    Ah = cla_utils.Q1AQ1s(A)
    assert(np.abs(np.trace(A0) - np.trace(Ah)) < 1.0e-6)
    b = random.randn(m)
    x0 = np.dot(A0, b)
    xh = np.dot(Ah, b)
    assert(np.abs(np.linalg.norm(x0) - np.linalg.norm(xh)) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_hessenberg(m):
    random.seed(4373*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    cla_utils.hessenberg(A)
    # check preserves trace
    assert(np.abs(np.trace(A0) - np.trace(A)) < 1.0e-6)
    b = random.randn(m)
    x0 = np.dot(A0, b)
    xh = np.dot(A, b)
    # check transformation was via unitary transformations
    assert(np.abs(np.linalg.norm(x0) - np.linalg.norm(xh)) < 1.0e-6)
    # check Hessenberg structure
    assert(np.linalg.norm(A[np.tril_indices(m, -2)]) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_hessenbergQ(m):
    random.seed(3213*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    Q = cla_utils.hessenbergQ(A)
    assert(np.abs(np.trace(A0) - np.trace(A)) < 1.0e-6)
    b = random.randn(m)
    x0 = np.dot(A0, b)
    xh = np.dot(A, b)
    assert(np.abs(np.linalg.norm(x0) - np.linalg.norm(xh)) < 1.0e-6)
    # check Hessenberg structure
    assert(np.linalg.norm(A[np.tril_indices(m, -2)]) < 1.0e-6)
    # check the Schur factorisation
    assert(np.linalg.norm(A - np.dot(Q, np.dot(A, np.conj(Q).T))))


@pytest.mark.parametrize('m', [20, 204, 18])
def test_ev(m):
    random.seed(3213*m)
    A = random.randn(m, m)
    A = 0.5*(A + np.conj(A).T)
    A0 = 1.0*A
    ee, V = cla_utils.ev(A)
    for i in range(m):
        ee1 = np.dot(V[:, i], np.dot(A, V[:, i]))
        assert(np.abs(ee[i] - ee1) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
