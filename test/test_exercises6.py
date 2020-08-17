'''Tests for the fifth exercise set.'''
import pytest
from cla_utils import LU_inplace
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_LU_inplace(m):
    random.seed(8564*m)
    A = random.randn(m, m)
    A0 = 1.0*A
    LU_inplace(A)
    L = np.eye(m)
    i1 = np.tril_indices(m, k=-1)]
    L[i1] = A[i1]
    U = np.triu(A)
    A1 = np.dot(L, U)
    assert(np.abs(A1 - A0) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
