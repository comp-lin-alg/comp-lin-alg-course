import numpy as np
import timeit
import numpy.random as random

#pre-construct a matrix in the namespace to use in tests
random.seed(1651)
A0 = random.randn(500, 500)
x0 = random.randn(500)

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

    pass


def timeable_basic_matvec():
    """
    Doing a matvec example with the basic_matvec that we can
    pass to timeit.
    """
    
    b = basic_matvec(A0, x0)


def timeable_column_matvec():
    """
    Doing a matvec example with the column_matvec that we can
    pass to timeit.
    """
    
    b = column_matvec(A0, x0)


def timeable_numpy_matvec():
    """
    Doing a matvec example with the builtin numpy matvec so that
    we can pass to timeit.
    """

    b = A0.dot(x0)

    
def time_matvecs():
    """
    Get some timings for matvecs.
    """

    print("Timing for basic_matvec")
    print(timeit.Timer(timeable_basic_matvec).timeit(number=1))
    print("Timing for column_matvec")
    print(timeit.Timer(timeable_column_matvec).timeit(number=1))
    print("Timing for numpy matvec")
    print(timeit.Timer(timeable_numpy_matvec).timeit(number=1))
