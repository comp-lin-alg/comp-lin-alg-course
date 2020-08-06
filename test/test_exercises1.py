'''Test integration using the quadrature module.'''
import pytest
from cla_utils import basic_matvec, column_matvec
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


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
