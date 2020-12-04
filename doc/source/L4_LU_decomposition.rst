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

An algorithm for LU decomposition
---------------------------------

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454095315>`_.

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

      A_n = {L_n\ldots L_2L_1}A = U.

This process of transforming `A` to `U` is called Gaussian elimination.

If we assume (we will show this later) that all these lower triangular
matrices are invertible, we can define

   .. math::

      L = (L_n\ldots L_2L_1)^{-1} = L_1^{-1}L_2^{-1}\ldots L_n^{-1},

      \mbox{ so that }

      L^{-1} = L_n\ldots L_2L_1.

Then we have `L^{-1}A = U`, i.e. `A=LU`.

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454096015>`_.

So, what's the advantage of writing `A=LU`? Well, we can define
`y=Ux`.  Then, we can solve `Ax=b` in two steps, first solving `Ly=b`
for `y`, and then solving `Ux=y` for `x`. The latter equation is an
upper triangular system that can be solved by the back
substitution algorithm we introduced for QR factorisation. The former
equation can be solved by forward substitution, derived in an analogous
way, written in pseudo-code as follows.

* `x_1 \gets b_1/L_{11}`
* FOR `i= 2` TO `m`

  * `x_i \gets (b_i - \sum_{k=1}^iL_{ik}x_k)/L_{ii}`

Forward substitution has an operation count that is identical to back
substitution, by symmetry, i.e. `\mathcal{O}(m^2)`. In contrast, we
shall see shortly that the Gaussian elimination process has an
operation count `\mathcal{O}(m^3)`. Hence, it is much cheaper to solve
a linear system with a given `LU` factorisation than it is to form `L`
and `U` in the first place. We can take advantage of this in the
situation where we have to solve a whole sequence of linear systems
`Ax=b_i`, `i=1,2,\ldots,K`, with the same matrix `A` but different
right hand side vectors. In this case we can pay the cost of forming
`LU` once, and then use forward and back substitution to cheaply solve
each system. This is particularly useful when we need to repeatedly
solve systems as part of larger iterative algorithms, such as time
integration methods or Monte Carlo methods.

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454096580>`_.

So, we need to find lower triangular matrices `L_k` that do not change
the first `k-1` rows, and transforms the `k`-th column `x_k` of `A_k`
as follows.

   .. math::

      Lx_k = L\begin{pmatrix}
      x_{1k}\\
      \vdots\\
      x_{kk}\\
      x_{k+1,k}\\
      \vdots\\
      x_{m,k}\\
      \end{pmatrix}
      = \begin{pmatrix}
      x_{1k}\\
      \vdots\\
      x_{kk}\\
      0 \\
      \vdots\\
      0 \\
      \end{pmatrix}.

As before with the Householder method, we see that we need the top-left
`k\times k` submatrix of `L` to be the identity (so that it doesn't change
the first `k` rows). We claim that the following matrix transforms
`x_k` to the required form.

   .. math::

      L_k = \begin{pmatrix}
      1 & 0 & 0 & \ldots & 0 & \ldots & \ldots & \ldots & 0 \\
      0 & 1 & 0 & \ldots & 0 & \ldots & \ldots& \vdots & 0 \\
      0 & 0 & 1 & \ldots & 0 & \ldots & \ldots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & \vdots & \vdots & \vdots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & 1 & 0 & \ldots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & -l_{k+1,k} & 1 & \ldots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & -l_{k+2,k} & 0 & \ddots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & \vdots & 0 & \ldots & \ddots & 0 \\
      \vdots & \ddots & \ddots & \ddots & -l_{m,k} & 0 & \ldots & \ldots &1 \\
      \end{pmatrix},

      \quad

      l_k = \begin{pmatrix}
      0 \\
      0 \\
      0 \\
      \vdots \\
      0 \\
      l_{k+1,k}=x_{k+1,k}/x_{kk} \\
      l_{k+2,k}= x_{k+2,k}/x_{kk} \\
      \vdots\\
      l_{m,k} = x_{m,k}/x_{kk} \\
      \end{pmatrix}.

This has the identity block as required, and we can verify that `L_k`
puts zeros in the entries of `x_k` below the diagonal by first writing
`L_k = I - l_ke_k^*`. Then,

   .. math::

      L_kx_k = I - l_ke_k^* = x_k - l_k\underbrace{(e_k^*x_k)}_{=x_{kk}},

which subtracts off the below diagonal entries as required. Indeed,
multiplication by `L_k` implements the row operations that are performed
to transform below diagonal elements of `A_k` to zero during Gaussian
elimination.

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454097320>`_.

The determinant of a lower triangular matrix is equal to the trace
(product of diagonal entries), so `\det(L_k)=1`, and consequently
`L_k` is invertible, enabling us to define `L^{-1}` as above.
To form `L` we need to multiply the inverses of all the `L_k` matrices
together, also as above. To do this, we first note that `l_k^*e_k=0`
(because `l_k` is zero in the only entry that `e_k` is nonzero). Then
we claim that `L_k^{-1}=I + l_ke_k^*`, which we verify as follows.

   .. math::

      (I + l_ke_k^*)L_k =       (I + l_ke_k^*)(I - l_ke_k^*)
      = I + l_ke_k^* - l_ke_k^* + (l_ke_k^*)(l_ke_k*)

      = I + \underbrace{l_k(e_k^*l_k)e_k*}_{=0} = I,

as required. Similarly if we multiply the inverse lower triangular
matrices from two consecutive iterations, we get

   .. math::

      L_k^{-1}L_{k+1}^{-1} = (I + l_ke_k^*)(I + l_{k+1}e_{k+1}^*)
      = I + l_ke_k^* + l_{k+1}e_{k+1}^* + l_k\underbrace{(e_k^*l_{k+1})}_{=0}e_{k+1}^*

      = I + l_ke_k^* + l_{k+1}e_{k+1}^*,

since `e_k^*l_{k+1}=0` too, as `l_{k+1}` is zero in the only place
where `e_k` is nonzero. If we iterate this argument, we get

   .. math::

      L = I + \sum_{i=1}^{m-1}l_ie_i^*.

Hence, the `k`th column of `L` is the same as the `k`th column of `L_k^{-1}`,
i.e.,

   .. math::

      L = \begin{pmatrix}
      1 & 0 & 0 & \ldots & 0 & \ldots & \ldots & \ldots & 0 \\
      l_{21} & 1 & 0 & \ldots & 0 & \ldots & \ldots& \vdots & 0 \\
      l_{31} & l_{32} & 1 & \ldots & 0 & \ldots & \ldots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & \vdots & \vdots & \vdots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & 1 & 0 & \ldots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & l_{k+1,k} & 1 & \ldots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & l_{k+2,k} & l_{k+2,k+1} & \ddots & \vdots & 0 \\
      \vdots & \ddots & \ddots & \ddots & \vdots & l_{m-1,k+1} & \ldots & \ddots & 0 \\
      \vdots & \ddots & \ddots & \ddots & l_{m,k} & l_{m,k+1} & \ldots & \ldots &1 \\
      \end{pmatrix}.

In summary, we can compute entries of `L` during the Gaussian elimination
process of transforming `A` to `U`. Note that the matrices `L_1,L_2,\ldots`
should not be explicitly formed during the elimination process, they are just
a mathematical concept to translate from the row operations into the final
`L` matrix.

.. proof:exercise::

   Having said that, let's take a moment to compute some examples
   using the `L_1,L_2,\ldots` matrices (to help with understanding).
   The :func:`cla_utils.exercises6.get_Lk` function has been left
   unimplemented. It should return one of these matrices given the
   `l_k` entries.  The test script ``test_exercises6.py`` in the
   ``test`` directory will test this function.

   Once it passes the tests, experiment with the inverse and
   multiplication properties above, to verify that they work.

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454098164>`_.

The Gaussian elimination algorithm is written in pseudo-code as
follows. We start by copying `A` into `U`, and setting `L` to
an identity matrix, and then work "in-place" i.e. replacing values
of `U` and `L` until they are completed. In a computer implementation,
this memory should be preallocated and then written to instead of
making copies (which carries overheads).

* `U \gets A`
* `L \gets I`
* FOR `k=1` TO `m-1`

  * for `j=k+1` TO `m`

    * `l_{jk} \gets u_{jk}/u_{kk}`
    * `u_{j,k:m} \gets u_{j,k:m} - l_{jk}u_{k,k:m}`
  * END FOR
* END FOR

To do an operation count for this algorithm, we note that the
dominating operation is the update of `U` inside the `j` loop. This
requires `m-k+1` multiplications and subtractions, and is iterated
`m-k` times in the `j` loop, and this whole thing is iterated from
`j=k+1` to `m`. Hence the asymptotic operation count is

   .. math::

      N_{\mbox{FLOPs}} = \sum_{k=1}^{m-1}\sum_{j=k+1}^m 2(m-k+1),

      = \sum_{k=1}^{m-1}2(m-k+1)\underbrace{\sum_{j={k+1}}^m 1}_{=m-k}

      = \sum_{k=1}^{m-1}2m^2 - 4mk + 2k^2

      \sim 2m^3 -4\frac{m^3}{2} + \frac{2m^3}{3} = \frac{2m^3}{3}.

.. proof:exercise::

   Since the diagonal entries of `L` are all ones, the total amount of
   combined memory required to store `L` and `U` is the same as the
   amount of memory required to store `A`. Further, each iteration of
   the LU factorisation algorithm computes one column of `L` and one
   rows of `U`, and the corresponding column an row of `A` are not
   needed for the rest of the algorithm. This creates the opportunity
   for a memory-efficient 'in-place' algorithm in which the matrix `A`
   is modified until it contains the values for `L` and `U`.

   The :func:`cla_utils.exercises6.LU_inplace` function has been left
   unimplemented. It should implement this in-place low-storage
   procedure, applying the changes to the provided matrix `A`.  The
   test script ``test_exercises6.py`` in the ``test`` directory will
   test this function.

.. proof:exercise::

   The LU factorisation requires 3 loops (this is why it has a cubic
   FLOP count). In the algorithm above, there are two explicit loops
   and one explicit one (in the slice notation). It is possible to
   rewrite this in a single loop, using an outer product. Identify
   this outer product, and update
   :func:`cla_utils.exercises6.LU_inplace` to make use of this
   reformulation (using :func:`numpy.outer`). Do you notice any
   improvement in speed?

.. proof:exercise::

   The functions :func:`cla_utils.exercises6.solve_L` and
   :func:`cla_utils.exercises6.solve_U` have been left unimplemented.
   They should use forward and backward substitution to solve lower
   and upper triangular systems respectively. The interfaces are set
   so that multiple right hand sides can be provided and solved at the
   same time. The functions should only use one loop over the columns
   of `L` (or `U`), to efficiently solve the multiple problems. The
   test script ``test_exercises6.py`` in the ``test`` directory will
   test these functions.

.. proof:exercise::

   Propose an algorithm to use the LU factorisation to compute the
   inverse of a matrix.  The functions
   :func:`cla_utils.exercises6.inverse_LU` has been left unimplemented.
   Complete it using your algorithm, using functions developed in the
   previous exercises where possible. The test script
   ``test_exercises6.py`` in the ``test`` directory will test these
   functions.

Pivoting
--------

.. hint::

   Video recordings for this material is available `here
   <https://player.vimeo.com/video/454098919>`_, and then
   `here
   <https://player.vimeo.com/video/454108809>`_.

Gaussian elimination will fail if a zero appears on the diagonal,
i.e. we get `x_{kk}=0` (since then we can't divide by it). Similarly,
Gaussian elimination will amplify rounding errors if `x_{kk}` is very
small, because a small error becomes large after dividing by `x_{kk}`.
The solution is to reorder the rows in `A_k` so that that `x_{kk}` has
maximum magnitude. This would seem to mess up the `LU` factorisation
procedure. However, it is not as bad as it looks, as we will now
see.

The main tool is the permutation matrix.

.. proof:definition:: Permutation matrix

   An `m\times m` permutation matrix has precisely one entry equal to
   1 in every row and column, and zero elsewhere.

A compact way to store a permutation matrix `P` as a size `m` vector
`p`, where `p_i` is equal to the number of the column containing the 1
entry in row `i` of `P`.  Multiplying a vector `x` by a permutation
matrix `P` simply rearranges the entries in `x`, with `(Px)_i =
x_{p_i}`.

During Gaussian elimination, say that we are at stage `k`, and
`(A_k)_{kk}` is not the largest magnitude entry in the `k`th column of
`A_k`. We reorder the rows to fix this, and this is what we call
*pivoting*. Mathematically this reordering is equivalent to
multiplication by a permutation matrix `P_k`. Then we continue the
Gaussian elimination procedure by left multiplying by `L_k`, placing
zeros below the diagonal in column `k` of `P_kA_k`.

In fact, `P_k` is a very specific type of permutation matrix, that only
swaps two rows. Therefore, `P_k^{-1}=P_k`, even though this is not
true for general permutation matrices.

We can pivot at every stage of the procedure, producing a permutation
matrix `P_k`, `k=1,\ldots, {m-1}` (if no pivoting is necessary at a given
stage, then we just take the identity matrix as the pivoting matrix
for that stage). Then, we end up with the result of Gaussian elimination
with pivoting,

   .. math::

      L_{m-1}P_{m-1}\ldots L_2P_2L_1P_1 = U.

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454109227>`_.

This looks like it has totally messed up the LU factorisation, because
`LP` is not lower triangular for general lower triangular matrix `L`
and permutation matrix `P`. However, we can save the situation, by
trying to swap all the permutation matrices to the right of all of the
`L` matrices. This does change the `L` matrices, because matrix-matrix
multiplication is not commutative. However, we shall see that it does
preserve the lower triangular matrix structure.

To see how this is done, we focus on how things look after two stages
of Gaussian elimination. We have

   .. math::

      A_2 = L_2P_2L_1P_1 = L_2\underbrace{P_2L_1P_2}_{=L_1^{(2)}}P_2P_1
      = L_2L_1^{(2)}P_2P_1,

having used `P_2^{-1}=P_2`. Left multiplication with `P_2` exchanges
row 2 with some other row `j` with `j>2`. Hence, right multiplication
with `P_2` does the same thing but with columns instead of rows.
Therefore, `L_1P_2` is the same as `L_1` but with column 2 exchanged
with column `j`. Column 2 is just `e_2` and column `j` is just `e_j`,
so now column 2 has the 1 in row `j` and column `j` has the 1 in
row 2. Then, `P_2L_1P_2` exchanges row 2 of `L_1P_2` with row `j` of
`L_1P_2`. This just exchanges `l_{12}` with `l_{1j}`, and swaps the
1s in columns 2 and `j` back to the diagonal. In summary, `P_2L_1P_2`
is the same as `L_1` but with `l_{12}` exchanged with `l_{1j}`.

Moving on to the next stage, and we have

   .. math::

      A_3 = L_3P_3L_2L_1P_2P_1 = L_3\underbrace{P_3L_2P_3}_{=L_2^{(3)}}
      \underbrace{P_3L_1P_3}_{=L_1^{(3)}}P_3P_2P_1.

By similar arguments we see that `L_2^{(3)}` is the same as `L_2` but
with `l_{23}` exchanged with `l_{2j}` for some (different) `j`, and
`L_2^{(3)}` is the same as `L_2^{(2)}` with `l_{13}` exchanged with
`l_{1j}`. After iterating this argument, we can obtain

   .. math::

      \underbrace{L_{m-1}^{(m-1)}\ldots L_2^{(m-1)}L_1^{(m-1)}}_{L^{-1}}
      \underbrace{P_{m-1}\ldots P_2P_1}_P = U,

where we just need to keep track of the permutations in the `L`
matrices as we go through the Gaussian elimination stages. These `L`
matrices have the same structure as the basic LU factorisation, and hence
we obtain

   .. math::

      L^{-1}PA = U \implies PA = LU.

This is equivalent to permuting the rows of `A` using `P` and then
finding the LU factorisation using the basic algorithm (except we
can't implement it like that because we only decide how to build `P`
during the Gaussian elimination process).

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454109660>`_.

The LU factorisation with pivoting can be expressed in the following
pseudo-code.

* `U\gets A`
* `L\gets I`
* `P\gets I`
* FOR `k=1` TO `m-1`

  * Choose `i\geq k` to maximise `|u_{ik}|`
  * `u_{k,k:m} \leftrightarrow u_{i,k:m}` (row swaps)
  * `l_{k,1:k-1} \leftrightarrow l_{i,1:k-1}` (row swaps)
  * `p_{k,1:m} \leftrightarrow p_{i,1:m}`
  * FOR `j=k+1` TO `m`

    * `l_{jk} \gets u_{jk}/u_{kk}`
    * `u_{j,k:m} \gets u_{j,k:m} - l_{jk}u_{k,k:m}`
  * END FOR
* END FOR

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454110324>`_.

To solve a system `Ax=b` given the a pivoted LU factorisation `PA=LU`,
we left multiply the equation by `P` and use the factorisation get
`LUx=Pb`. The procedure is then as before, but `b` must be permuted to
`Pb` before doing the forwards and back substitutions.

We call this strategy *partial pivoting*. In contrast, *complete
pivoting* additionally employs permutations `Q_k` on the right that
swap columns of `A_k` as well as the rows swapped by the permutations
`P_k`. By similar arguments, one can obtain the LU factorisation with
complete pivoting, `PAQ=LU`.

.. proof:exercise::

   The function :func:`cla_utils.exercises7.perm` has been left
   unimplemented. It should take an `m\times m` permutation matrix
   `P`, stored as a vector of indices `p\in\mathbb{N}^m` so that
   `(Px)_i = x_{p_i}`, `i=1,2,\ldots, m`, and replace it with the
   matrix `P_{i,j}P` (also stored as a vector of indices) where
   `P_{i,j}` is the permutation matrix that exchanges the entries `i`
   and `j`. The test script ``test_exercises7.py`` in the ``test``
   directory will test this function.


.. proof:exercise::

   The function :func:`cla_utils.exercises7.LUP_inplace` has been left
   unimplemented. It should extend the in-place algorithm for LU
   factorisation (with the outer-product formulation, if you managed
   it) to the LUP factorisation. As well as computing L and U "in
   place" in the array where the input A is stored, it will compute a
   permutation matrix, which can and should be constructed using
   :func:`cla_utils.exercises7.perm`.The test script
   ``test_exercises7.py`` in the ``test`` directory will test this
   function.


.. proof:exercise::

   The function :func:`cla_utils.exercises7.solve_LUP` has been left
   unimplemented. It should use the LUP code that you have written to
   solve the equation `Ax=b` for `x` given inputs `A` and `b`.  The
   test script ``test_exercises7.py`` in the ``test`` directory will
   test this function.

.. proof:exercise::

   Show how to compute the determinant of `A` from the LUP
   factorisation in `\mathcal{O}(m)` time (having already constructed
   the LUP factorisation which costs `\mathcal{O}(m^3)`). Complete the
   function :func:`cla_utils.exercises7.det_LUP` to implement this
   computation. The test script ``test_exercises7.py`` in the ``test``
   directory will test this function.

Stability of LU factorisation
-----------------------------

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454110810>`_.

To characterise the stability of LU factorisation, we quote the following
result.

.. proof:theorem::

   Let `\tilde{L}` and `\tilde{U}` be the result of the Gaussian
   elimination algorithm implemented in a floating point number system
   satisfying axioms I and II. If no zero pivots are encountered, then

      .. math::

	 \tilde{L}\tilde{U} = A + \delta A

   where

      .. math::

	 \frac{\|\delta A\|}{\|L\|\|U\|} = \mathcal{O}(\varepsilon),

   for some perturbation `\delta A`.

The algorithm is backward stable if `\|L\|\|U\|=\mathcal{O}(\|A\|)`,
but there will be problems if `|L\|\|U\|\gg \|A\|`. For a proof of this
result, see the textbook by Golub and van Loan.

A similar result exists for pivoted LU. The main extra issue is that
small changes could potentially lead to a different pivoting matrix
`\tilde{P}` which is then `O(1)` different from `P`. This is characterised
in the following result (which we also do not prove).

.. proof:theorem::

   Let `\tilde{P}`, `\tilde{L}` and `\tilde{U}` be the result of the
   partial pivoted Gaussian elimination algorithm implemented in a
   floating point number system satisfying axioms I and II. If no zero
   pivots are encountered, then

      .. math::

	 \tilde{L}\tilde{U} = A + \delta A

   where

      .. math::

	 \frac{\|\delta A\|}{\|A\|} = \mathcal{O}(\rho\varepsilon),

   for some perturbation `\delta A`, and where `\rho` is the growth
   factor,

      .. math::

	 \rho = \frac{\max_{ij}|u_{ij}|}{|a_{ij}|}.

Thus, partial pivoting (and complete pivoting turns out not to help
much extra) can keep the entries in `L` under control, but there can
still be pathological cases where entries in `U` can get large,
leading to large `\rho` and unstable computations.

Taking advantage of matrix structure
------------------------------------

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454111577>`_.

The cost of the standard Gaussian elimination algorithm to form `L`
and `U` is `\mathcal{O}(m^3)`, which grows rather quickly as `m`
increases. If there is structure in the matrix, then we can often
exploit this to reduce the cost. Understanding when and how to exploit
structure is a central theme in computational linear algebra.
Here we will discuss some examples of structure to be exploited.

When `A` is a lower or upper triangular matrix then we can use
forwards or back substitution, with `\mathcal{O}(m^2)` operation count
as previously discussed.

When `A` is a diagonal matrix, i.e. `A_{ij}=0` for `i\ne j`, it only
has `m` nonzero entries, that can be stored as a vector,
`(A_{11},A_{22},\ldots,A_{mm})`. In this case, `Ax=b` can be solved in
`m` operations, just by setting `x_i=b_i/A_{ii}`, for
`i=1,2,\ldots,m`.

Similarly, if `A \in \mathcal{C}^{dm\times dm}` is block diagonal,
i.e.

   .. math::

      A = \begin{pmatrix}
      B_{1} & 0 & \ldots & 0 \\
      0 & B_{2} & \ldots & 0 \\
      \vdots & \vdots & \ddots & 0 \\
      0 & 0 & \ldots & B_{m}
      \end{pmatrix},

where `B_{i}\in\mathcal{C}^{d\times d}` for `i=1,2,\ldots,m`. The inverse
of `A` is

   .. math::

      A = \begin{pmatrix}
      B_{1}^{-1} & 0 & \ldots & 0 \\
      0 & B_{2}^{-1} & \ldots & 0 \\
      \vdots & \vdots & \ddots & 0 \\
      0 & 0 & \ldots & B_{m}^{-1}
      \end{pmatrix}.

A generalisation of a diagonal matrix is a banded matrix, where
`A_{ij}=0` for `i>j+p` and for `i<j-q`. We call `p` the upper
bandwidth of `A`; `q` is the lower bandwidth. When the matrix is
banded, there are already zeros below the diagonal of `A`, so we know
that the corresponding entries in the `L_k` matrices will be zero.
Further, because there are zeros above the diagonal of `A`, these do
not need to be updated when applying the row operations to those
zeros.

.. proof:exercise::

   Construct the `100\times 100` matrix `A` as follows: take `A=3I`,
   then set `A_{1,i}=1`, for `i=1,\ldots,100`. Then set `A_{i,1}=i` for
   `i=1,\ldots,100`.  Using your own LU factorisation, compute the LU
   factorisation of `A`. What
   do you observe about the number of non-zero entries in `L` and `U`?
   Explain this using what you have just learned about banded
   matrices. Can the situation be improved by pivoting? (Just think about
   it, don't need to implement it.)

The Gaussian elimination algorithm (without pivoting) for a banded
matrix is given as pseudo-code below.

* `U \gets A`
* `L \gets I`
* FOR `k=j+1` TO `\min(k+p,m)`

  * `l_{jk} \gets u_{jk}/u_{kk}`
  * `n \gets \min(k+q, m)`
  * `u_{j,k:n} \gets u_{j,k:n}- l_{jk}u_{k,k:n}`
  * END FOR
* END FOR

The operation count for this banded matrix algorithm is
`\mathcal{O}(mpq)`, which is linear in `m` instead of cubic!
Further, the resulting matrix `L` has lower bandwidth `p`
and `U` has upper bandwidth `q`. This means that we can also
exploit this structure in the forward and back substitution
algorithms as well. For example, the forward substitution algorithm
is given as pseudo-code below.

* `x_1 \gets b_1/L_{11}`
* FOR `k=2` TO `m`

  * `j \gets \max(1, k-p)`
  * `x_k \gets \frac{b_k -L_{k,j:k-1}x_{j:k-1}}{L_{kk}}`
* END FOR

This has an operation count `\mathcal{O}(mp)`. The story is
very similar for the back substitution.

.. hint::

   A video recording for this material is available `here
   <https://player.vimeo.com/video/454112153>`_.

Another example that we have already encountered is unitary matrices
`Q`. Since `Q^{-1}=Q^*`, solving the system `Qx=b` is just the cost of
applying `Q^*`, with operation count `\mathcal{O}(m^2)`.

An important matrix that we shall encounter later is an upper
Hessenberg matrix, that has a lower bandwidth of 1, but no particular
zero structure above the diagonal. In this case, the `L` matrix is
still banded (with lower bandwidth 1) but the `U` matrix is not.  This
means that there are still savings due to the zeros in `L`, but work
has to be done on the entire column of `U` above the diagonal, and so
solving an upper Hessenberg system has operation count
`\mathcal{O}(m^2)`.

Cholesky factorisation
----------------------

An example of extra structure which we shall discuss in a bit more
detail is the case of Hermitian positive definite matrices. Recall
that a Hermitian matrix satisfies `A^*=A`, whilst positive definite
means that

   .. math::

      x^*Ax > 0, \, \forall \|x\|>0.

When `A` is Hermitian positive definite, it is possible to find an
upper triangular matrix `R` such that `A=R^*R`, which is called the
Cholesky factorisation. To show that it is possible to compute
the Cholesky factorisation, we start by assuming that `A` has
a 1 in the top-left hand corner, so that

   .. math::

      A = \begin{pmatrix}
      1 & w^* \\
      w & K \\
      \end{pmatrix}

where `w` is a `m-1` vector containing the rest of the first column
of `A`, and `K` is an `(m-1)\times(m-1)` Hermitian positive
definite matrix. (Exercise: show that `K` is Hermitian positive
definite.)

After one stage of Gaussian elimination, we have

   .. math::

      \underbrace{\begin{pmatrix}
      1 & 0 \\
      -w & I \\
      \end{pmatrix}}_{L_1^{-1}}
      \underbrace{
      \begin{pmatrix}
      1 & w^* \\
      w & K \\
      \end{pmatrix}}_{A}
      =
      \begin{pmatrix}
      1 & w^* \\
      0 & K - ww^* \\
      \end{pmatrix}.

Further,

   .. math::

      \begin{pmatrix}
      1 & w^* \\
      0 & K - ww^* \\
      \end{pmatrix}=
      \underbrace{
      \begin{pmatrix}
      1 & 0 \\
      0 & K - ww^* \\
      \end{pmatrix}}_{A_1}
      \underbrace{
      \begin{pmatrix}
      1 & w^T \\
      0 & I \\
      \end{pmatrix}}_{(L_1^{-1})^*=L_1^{-*}},

so that `A = L_1^{-1}A_1L_1^{-*}`. If `a_{11} \neq 1`, we at least
know that `a_{11}= e_1^*Ae_1>0`, and the factorisation becomes

   .. math::

      A =
      \underbrace{\begin{pmatrix} \alpha & 0 \\
      w/\alpha & I \\
      \end{pmatrix}}_{R_1^T}
      \underbrace{
      \begin{pmatrix}
      1 & 0 \\
      0 & K - \frac{ww^*}{a_{11}} \\
      \end{pmatrix}}_{A_1}
      \underbrace{
      \begin{pmatrix}
      \alpha & w/\alpha \\
      0 & I \\
      \end{pmatrix}}_{R_1},

where `\alpha=\sqrt{a_{11}}`. We can check that `A_1` is positive
definite, since

   .. math::

      x^*A_1x = x^*R_1^{-*}AR_1x = (R_1^{-1}x)^*AR_1x = y^*Ay > 0, \mbox{ where }
      y = R_1x.

Hence, `K-{ww^*}/{a_{11}}` is positive definite, since

   .. math::

      r^*\left(K-\frac{ww^*}{a_{11}}\right)r = \begin{pmatrix} 0 \\ r \\ \end{pmatrix}^*
      A_1 \begin{pmatrix} 0 \\ r \\ \end{pmatrix} > 0,

and hence we can now perform the same procedure all over again to `K -
{ww^*}/a_{11}`. By induction we can always continue until we have the
required Cholesky factorisation, which is unique (since there were no
choices to be made at any step).

We can then present the Cholesky factorisation as pseudo-code.

* `R\gets A`
* FOR `k=1` TO `m`

  * FOR `j=k+1` to `m`

    * `R_{j,j:m} \gets R_{j,j:m} - R_{k,j:m}\bar{R}_{kj}/R_{kk}`
  * `R_{k,k:m} \gets R_{k,k:m}/\sqrt{R_{k:k}}`

The operation count of the Cholesky factorisation is dominated
by the operation inside the `j` loop, which has one division,
`m-j+1` multiplications, and `m-j+1` subtractions, giving
`\sim 2(m-j)` FLOPs. The total operation count is then

   .. math::

      N_{\mbox{FLOPs}} = \sum_{k=1}^m\sum_{j=k+1}^m
      \sim \frac{1}{3}m^3.
