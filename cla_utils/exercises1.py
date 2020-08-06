import numpy as np


def basic_matvec(A, x):
    """
    Elementary matrix-vector multiplication.

    :param A: an mxn-dimensional numpy array
    :param x: an n-dimensional numpy array

    returns an m-dimensional numpy array which is the product of A with x

    This should be implemented using a double loop over the entries of A
    """

    pass


def column_matvec(A, x):
    """
    Matrix-vector multiplication using the representation of the product
    Ax as linear combinations of the columns of A, using the entries in 
    x as coefficients.


    :param A: an mxn-dimensional numpy array
    :param x: an n-dimensional numpy array

    returns an m-dimensional numpy array which is the product of A with x

    This should be implemented using a single loop over the entries of x
    """
