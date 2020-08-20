import numpy as np


def get_A3():
    """
    Return A3 matrix for investigating power iteration.
    
    :return A3: a 3x3 numpy array.
    """

    return array([[ 0.68557183+0.46550108j,  0.12934765-0.1622676j ,
                    0.24409518+0.25335939j],
                  [ 0.1531015 +0.66678983j,  0.45112492+0.18206976j,
                    -0.02633966+0.43477693j],
                  [-0.10817164-1.16879196j, -0.18446849+0.03755672j,
                   0.06430325-0.44757084j]])


def get_B3():
    """
    Return B3 matrix for investigating power iteration.

    :return B3: a 3x3 numpy array.
    """
    return array([[ 0.46870499+0.37541453j,  0.19115959-0.39233203j,
                    0.12830659+0.12102382j],
                  [ 0.90249603-0.09446345j,  0.51584055+0.84326503j,
                    -0.02582305+0.23259079j],
                  [ 0.75419973-0.52470311j, -0.59173739+0.48075322j,
                    0.51545446-0.21867957j]])


def pow_it(A, x0, tol, maxit, store_iterations = False):
    """
    For a matrix A, apply the power iteration algorithm with initial
    guess x0, until either 

    ||r|| < tol where

    r = Ax/||x|| - x/||x||,

    or the number of iterations exceeds maxit.

    :param A: an mxm numpy array
    :param x0: the starting vector for the power iteration
    :param tol: a positive float, the tolerance
    :param maxit: integer, max number of iterations
    :param store_iterations: if True, then return the entire sequence \
    of power iterates, instead of just the final iteration. Default is \
    False.

    :return x0: an m dimensional numpy array, or \
    if store_iterations, an mxmaxit dimensional numpy array.
    """

    xi = 1.0*x0
    m = x0.size
    if store_iterations:
        xt = np.zeros((m, maxit+1), dtype=np.float64)
        xt[:, 0] = x0
    for i in range(maxit):
        xi = np.dot(A, xi)
        xi /= np.linalg.norm(xi)
        r = np.dot(A, xi)
        r /= np.linalg.norm(r)
        r -= xi
        rmag = np.linalg.norm(r)
        if store_iterations:
            xt[:, i+1] = xi
        if rmag < tol:
            break
    if store_iterations:
        return xt
    else:
        return xi
