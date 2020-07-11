.. default-role:: math

LU decomposition
================

In this section we look at the some other algorithms for solving the
equation `Ax=b` when `A` is invertible. On the one hand the `QR`
factorisation has great stability properties. On the other, it can be
beaten by other methods for speed when there is particular structure
to exploit (such as lots of zeros in the matrix). In this section, we
explore the the family of methods that go right back to the technique
of Gaussian elimination, that you will have been familiar with since
secondary school.

The computational way to view Gaussian elimination is through the LU
decomposition of an invertible matrix, `A=LU`, where `L` is lower
triangular (`l_{ij}=0` for `j<i`) and `U` is upper triangular
(`u_{ij}=0` for `j>i`). Here we use the symbol `U` instead of `R` to
emphasise that we are looking as square matrices.  The process of
obtaining the `LU` decomposition is very similar to the Householder
algorithm, in that we repeatedly left multiply `A` by matrices to
transform below-diagonal entries in each column to zero, working from
the first to the last column. The difference is that whilst the
Householder algorithm left multiplies with unitary matrices, here,
we left multiply with lower triangular matrices.

The first step puts zeros below the first entry in the first column.

   .. math::

      A_1 = L_1A = \begin{pmatrix}
      u_1 & v_2^1 & v_2^1 & \ldots & v_n^1 \\
      \end{pmatrix},

      \,
      u_1 = \begin{pmatrix} u_{11} \\ 0 \\ \ldots \\ 0\end{pmatrix}.

Then, the next step puts zeros  below the second entry in the second
column.

   .. math::

      A_2 = L_2L_1A = \begin{pmatrix}
      u_1 & u_2 & v_2^2 & \ldots & v_n^2 \\
      \end{pmatrix},

      \,
      u_2 = \begin{pmatrix} u_{12} \\ u_{22} \\ 0 \\ \ldots \\ 0 \\
      \end{pmatrix}.

After repeated left multiplications we have

   .. math::

      A_n = \underbrace{L_n\ldots L_2L_1}A = U.

If we assume (we will show this later) that all these lower triangular
matrices are invertible, we can define

   .. math::

      L = (L_n\ldots L_2L_1)^{-1} = L_1^{-1}L_2^{-1}\ldots L_n^{-1},

      \mbox{ so that }

      L^{-1} = L_n\ldots L_2L_1.

Then we have `L^{-1}A = U`, i.e. `A=LU`.

So, we need to find lower triangular matrices `L_k` that do not change
the first `k-1` rows, and 
