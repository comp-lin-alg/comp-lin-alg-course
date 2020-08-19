'''Tests for the eighth exercise set.'''
import pytest
from cla_utils import Q1AQ1s, hessenberg
from numpy import random
import numpy as np


@pytest.mark.parametrize('m', [20, 204, 18])
def test_Q1AQ1s(m):
    random.seed(4373*m)
    A = random.randn(m, m)
    A0 = 1.0*A0
    Ah = Q1AQ1s(A)    
    assert(np.abs(np.trace(A0) - np.trace(Ah)) < 1.0e-6)
    b = random.randn(m, m)
    x0 = np.dot(A0, b)
    xh = np.dot(Ah, b)
    assert(np.abs(np.linalg.norm(x0) - np.linalg.norm(xh)) < 1.0e-6)


@pytest.mark.parametrize('m', [20, 204, 18])
def test_hessenberg(m):
    random.seed(4373*m)
    A = random.randn(m, m)
    A0 = 1.0*A0
    H = hessenberg(A)
    assert(np.abs(np.trace(A0) - np.trace(H)) < 1.0e-6)
    

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)
