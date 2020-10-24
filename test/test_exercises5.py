'''Tests for the fifth exercise set.'''
import pytest
import cla_utils
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_solve_R(m):
    random.seed(8323*m)
    A = random.randn(m, m)
    Q, R = np.linalg.qr(A)
    b = random.randn(m)
    x0 = cla_utils.solve_R(R, b)
    x1 = np.linalg.solve(R, b)
    err = x0 - x1

    assert(np.linalg.norm(err) < 1.0e-6)

    #does this test everything we want to exercise 5?


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
