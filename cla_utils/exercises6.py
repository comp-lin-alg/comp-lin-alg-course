import numpy as np


def get_Lk(m, lvec):
    """Compute the lower triangular row operation mxm matrix L_k 
    which has ones on the diagonal, and below diagonal entries
    in column k given by lvec (k is inferred from the size of lvec).

    :param m: integer giving the dimensions of L.
    :param lvec: a k-1 dimensional numpy array.

    :return Lk: an mxm dimensional numpy array.

    """
                     
    raise NotImplementedError


def LU_inplace(A):
    """Compute the LU factorisation of A, using the in-place scheme so
    that the strictly lower triangular components of the array contain
    the strictly lower triangular components of L, and the upper
    triangular components of the array contain the upper triangular
    components of U.

    :param A: an mxm-dimensional numpy array

    """
                     
    raise NotImplementedError


def solve_L(L, b):
    """
    Solve systems Lx_i=b_i for x_i with L lower triangular, i=1,2,\ldots,k

    :param L: an mxm-dimensional numpy array, assumed lower triangular
    :param b: an mxk-dimensional numpy array, with ith column containing \
    b_i

    :return x: an mxk-dimensional numpy array, with ith column containing \
    the solution x_i

    """
                     
    raise NotImplementedError


def solve_U(U, b):
    """
    Solve systems Ux_i=b_i for x_i with U upper triangular, i=1,2,\ldots,k

    :param U: an mxm-dimensional numpy array, assumed upper triangular
    :param b: an mxk-dimensional numpy array, with ith column containing \
    b_i

    :return x: an mxk-dimensional numpy array, with ith column containing \
    the solution x_i

    """
                     
    raise NotImplementedError


def inverse_LU(A):
    """
    Form the inverse of A via LU factorisation.

    :param A: an mxm-dimensional numpy array.

    :return Ainv: an mxm-dimensional numpy array.

    """
                     
    raise NotImplementedError
