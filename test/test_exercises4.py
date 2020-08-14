'''Tests for the fourth exercise set.'''
import pytest
from cla_utils import operator_2_norm
from numpy import random
import numpy as np


@pytest.mark.parametrize('m, n', [(20, 7), (40, 13), (9, 87)])
def test_householder_ls(m, n):
    random.seed(8473*m + 9283*k)
    A = random.randn(m, n)

    norm1 = operator_2_norm(A)
    u, s, v = linalg.svd(A)
    norm2 = s[0]

    assert(np.abs(norm1 - norm2) < 1.0e-6)


if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
