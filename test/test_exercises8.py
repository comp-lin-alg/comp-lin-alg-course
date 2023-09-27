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
    # check preserves trace (similarity transform)x
    assert(np.abs(np.trace(A0) - np.trace(Ah)) < 1.0e-6)
    # check transformation was via unitary transformations
    assert(np.abs(cla_utils.norm(Ah) - cla_utils.norm(A0)) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_hessenberg(m):
    random.seed(4373*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    cla_utils.hessenberg(A)
    # check preserves trace
    assert(np.abs(np.trace(A0) - np.trace(A)) < 1.0e-6)
    # check transformation was via unitary transformations
    assert(np.abs(cla_utils.norm(A) - cla_utils.norm(A0)) < 1.0e-6)
    # check Hessenberg structure
    assert(cla_utils.norm(A[np.tril_indices(m, -2)]) < 1.0e-6)


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
    assert(np.abs(cla_utils.norm(A) - cla_utils.norm(A0)) < 1.0e-6)
    # check Hessenberg structure
    assert(cla_utils.norm(A[np.tril_indices(m, -2)]) < 1.0e-6)
    # check orthogonality
    assert(cla_utils.norm(Q.T@Q - np.eye(m)) < 1.0e-6)
    # check the Schur factorisation
    assert(cla_utils.norm(A - np.dot((Q.conj()).T, np.dot(A0, Q))) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_ev(m):
    random.seed(3213*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    V = cla_utils.ev(A)
    #check that V and AV are aligned
    norm = cla_utils.norm
    for i in range(m):
        v = V[:, i]
        Av = A0@v
        v /= v[0]
        Av /= Av[0]
        assert(norm(Av - v) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
