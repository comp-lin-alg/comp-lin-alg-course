'''Test integration using the quadrature module.'''
import pytest
from cla_utils import basic_matvec, column_matvec, rank2
from numpy import random
import numpy as np


# Test the basic matvec
@pytest.mark.parametrize('m, n', [(20, 20), (40, 20), (20, 45)])
def test_basic_matvec(m, n):
    random.seed(1878*m + 1950*n)    
    A = random.randn(m, n)
    x = random.randn(n)

    b0 = A.dot(x)
    b = basic_matvec(A, x)

    assert(np.linalg.norm(b-b0)<1.0e-6)

# Test the basic matvec
@pytest.mark.parametrize('m, n', [(20, 20), (40, 20), (20, 45)])
def test_rank2_matrix(m, n):
    random.seed(1451*m + 1901*n)
    u1 = np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))
    u2 = np.sqrt(2)*(random.randn(m) + 1j*random.randn(m))
    v1 = np.sqrt(2)*(random.randn(n) + 1j*random.randn(n))
    v2 = np.sqrt(2)*(random.randn(n) + 1j*random.randn(n))

    A = rank2(u1, u2, v1, v2)
    a1 = random.randn(m)
    a2 = random.randn(n)

    n1 = np.vdot(a1, A.dot(a2))
    n2 = np.vdot(a1, u1)*np.vdot(v1, a2) + \
        np.vdot(a1, u2)*np.vdot(v2, a2)

    assert(np.abs(n1-n2)<1.0e-7)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
