import numpy as np


def randomQ(m):
    """
    Produce a random orthogonal mxm matrix.

    Inputs 

    :param m: the matrix dimension parameter.
    
    Outputs

    :return Q: the mxm numpy array containing the orthogonal matrix.
    """
    Q, R = linalg.qr(random.randn(m, m))
    return Q


def randomR(m):
    """
    Produce a random upper triangular mxm matrix.

    Inputs 

    :param m: the matrix dimension parameter.
    
    Outputs

    :return R: the mxm numpy array containing the upper triangular matrix.
    """
    
    A = random.randn(m, m)
    return numpy.triu(A)


def backward_stability_householder(m):
    """
    Verify backward stability for QR factorisation using Householder for
    real mxm matrices.

    Inputs

    :param m: the matrix dimension parameter.
    """
    # repeat the experiment a few times to capture typical behaviour
    for k in range(20):
        Q1 = randomQ(m)
        R1 = randomR(m)

        raise NotImplementedError


def solve_R(R, b):
    """
    Solve the system Rx=b where R is an mxm upper triangular matrix 
    and b is an m dimensional vector.

    Inputs

    :param A: an mxm-dimensional numpy array
    :param b: an m-dimensional numpy array

    Outputs

    :param x: an m-dimensional numpy array
    """
                     
    raise NotImplementedError


def back_stab_solve_R(m):
    """
    Verify backward stability for back substitution for
    real mxm matrices.

    Inputs

    :param m: the matrix dimension parameter.
    """
    # repeat the experiment a few times to capture typical behaviour
    for k in range(20):
        A = random.randn(m, m)
        R = np.triu(A)

        raise NotImplementedError


def back_stab_householder_solve(m):
    """
    Verify backward stability for the householder algorithm
    for solving Ax=b for an m dimensional square system.

    Inputs

    :param m: the matrix dimension parameter.
    """
    raise NotImplementedError
