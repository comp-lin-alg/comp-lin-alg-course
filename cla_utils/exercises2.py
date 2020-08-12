import numpy as np


def orthog_cpts(v, Q):
    """
    Given a vector v and an orthonormal set of vectors q_1,...q_n,
    compute v = r + u_1q_1 + u_2q_2 + ... + u_nq_n
    for scalar coefficients u_1, u_2, ..., u_n and
    residual vector r

    Inputs

    :param v: an m-dimensional numpy array
    :param Q: an mxn-dimensional numpy array whose columns are the
    orthonormal vectors

    Outputs

    :param r: an m-dimensional numpy array containing the residual
    :param u: an n-dimensional numpy array containing the coefficients
    """

    m, n = Q.shape

    r = np.zeros(m, dtype=np.float64)
    u = np.zeros(n, dtype=np.float64)
    return r, u


def solveQ(Q, b):
    """
    Given a unitary mxm matrix Q and a vector b, solve Qx=b for x.

    Inputs

    :param Q: an mxm dimensional numpy array containing the unitary matrix
    :param b: the m dimensional array for the RHS

    Outputs

    :param x" m dimensional array containing the solution.
    """

    raise NotImplementedError

    return x
