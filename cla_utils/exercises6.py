import numpy as np


def getLk(m, lvec):
    """Compute the lower triangular row operation mxm matrix L_k 
    which has ones on the diagonal, and below diagonal entries
    in column k given by lvec (k is inferred from the size of lvec).

    Inputs

    :param m: integer giving the dimensions of L.
    :param lvec: a k-1 dimensional numpy array.

    Outputs

    :param Lk: an mxm dimensional numpy array.

    """
                     
    raise NotImplementedError

    return Lk


def LU_inplace(A):
    """Compute the LU factorisation of A, using the in-place scheme so
    that the strictly lower triangular components of the array contain
    the strictly lower triangular components of L, and the upper
    triangular components of the array contain the upper triangular
    components of U.

    Inputs

    :param A: an mxm-dimensional numpy array

    Outputs

    No outputs as this is an in-place implementation.

    """
                     
    raise NotImplementedError
