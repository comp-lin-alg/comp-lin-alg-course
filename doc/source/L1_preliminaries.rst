.. default-role:: math

===============
 Preliminaries
===============

In this preliminary section, we revise a few key linear algebra
concepts that will be used in the rest of the course, emphasising
the column space of matrices. We will quote some standard results
that should be found in an undergraduate linear algebra course.

.. hint::

   Before you attempt any exercises, you need to make sure that you
   have everything you need set up on your computer. See the checklist
   in :ref:`the exercises section <comp_exercises>`.

Matrices, vectors and matrix-vector multiplication
==================================================

.. details:: Supplementary video

   .. vimeo:: 450145459

We will consider the multiplication of a vector

   .. math::

      x = \begin{pmatrix} x_1 \\
      x_2 \\
      \vdots \\
      x_n \\
      \end{pmatrix}, \quad x_i \in \mathbb{C}, \, i=1,2,\ldots,n,
      \mbox{ i.e. } x \in \mathbb{C}^n,

by a matrix

   .. math::

      A = \begin{pmatrix}
      a_{11} & a_{12} & \ldots & a_{1n} \\
      a_{21} & a_{22} & \ldots & a_{2n} \\
      \vdots & \vdots & \ddots & \vdots \\
      a_{m1} & a_{m2} & \ldots & a_{mn} \\
      \end{pmatrix},

i.e. `A\in \mathbb{C}^{m\times n}`. `A` has `m` rows and `n` columns
so that the product

   .. math::

      b = Ax

produces `b \in \mathbb{C}^m`, defined by

   .. math::
      :label: matvec

      b_i = \sum_{j=1}^n a_{ij}x_j, \, i=1,2,\ldots,m.

In this course it is important to
consider the general case where `m \neq n`, which has many applications
in data analysis, curve fitting etc. We will usually state generalities
in this course for vectors over the field `\mathbb{C}`, noting where things
specialise to `\mathbb{R}`.

.. details:: Supplementary video

   .. vimeo:: 450156255

We can quickly check that the map `x \to Ax` given by matrix
multiplication is a linear map from `\mathbb{C}^n \to \mathbb{C}^m`, since
it is straightforward to check from the definition that

   .. math::

      A(\alpha x + y) = \alpha Ax + Ay,

for all `x,y \in \mathbb{C}^n` and `\alpha\in \mathbb{C}`. (Exercise:
show this for yourself.)

.. details:: Supplementary video

   .. vimeo:: 450157385

It is very useful to interpret matrix-vector multiplication as a linear
combination of the columns of `A` with coefficients taken from the entries
of `x`. If we write `A` in terms of the columns,

   .. math::

      A = \begin{pmatrix}
      a_1 & a_2 & \ldots & a_n \\
      \end{pmatrix},

where

   .. math::

      a_i \in \mathbb{C}^m, \, i=1,2,\ldots,n,

then

   .. math::

      b = \sum_{j=1}^n x_j a_j,

i.e. a linear combination of the columns of `A` as described above.



   
.. details:: Supplementary video

   .. vimeo:: 450161699

We can extend this idea to matrix-matrix multiplication. Taking
`A\in \mathbb{C}^{m\times l}`, `C\in \mathbb{C}^{l\times n}`,
`B\in \mathbb{C}^{m\times n}`, with `B=AC`, then the components of
`B` are given by

   .. math::

      b_{ij} = \sum_{k=1}^l a_{ik}c_{kj}, \quad 1\leq i \leq m, \,
      1\leq j \leq n.

Writing `b_j \in \mathbb{C}^m` as the jth column of `B`, for `1\leq j \leq n`,
and `c_j` as the jth column of `C`,
we see that

   .. math::

      b_j = Ac_j.

This means that the jth column of `B` is the matrix-vector product of
`A` with the jth column of `C`. This kind of "column thinking" is very
useful in understanding computational linear algebra algorithms.



   
.. details:: Supplementary video

   .. vimeo:: 450162431

An important example is the outer product of two vectors, `u \in
\mathbb{C}^m` and `v \in \mathbb{C}^n`. Here it is useful to see these
vectors as matrices with one column, i.e. `u \in \mathbb{C}^{m\times
1}` and `v \in \mathbb{C}^{n\times 1}`. The outer product is `u v^T
\in \mathbb{C}^{m\times n}`. The columns of `v^T` are just single numbers
(i.e. vectors of length 1), so viewing this as a matrix multiplication
we see

   .. math::

      uv^T = \begin{pmatrix}
      uv_1 & uv_2 & \ldots & uv_n
      \end{pmatrix},

which means that all the columns of `uv^T` are multiples of `u`. We will
see in the next section that this matrix has rank 1. In the complex
number case, the transpose $^T$ is replaced by the adjoint $^*$ which
is the complex conjugate of the transpose. There will be more about this
later.

.. _ex-basic-matvec:

.. proof:exercise::

   The :func:`cla_utils.exercises1.basic_matvec` function has been left
   unimplemented. To finish the function, add code so that it
   computes the matrix-vector product `b=Ax` from inputs `A` and `x`.
   In this first implementation, you should simply implement
   :eq:`matvec` with a double nested for loop (one for the sum over `j`,
   and one for the `i` elements of `b`). Run this script to test your code
   (and all the exercises from this exercise set)::

      py.test test/test_exercises1.py

   from the Bash command line. Make sure you commit your modifications
   and push them to your fork of the course repository.

.. hint::

  Don't forget to activate the virtual environment before running the tests to make sure that you have access to all the necessary packages

.. _ex-column-matvec:

.. proof:exercise::

   The :func:`cla_utils.exercises1.column_matvec` function has been
   left unimplemented.  To finish the function, add code so that it
   computes the matrix-vector product `b=Ax` from inputs `A` and `x`.
   This second implementation should use the column-space formulation
   of matrix-vector multiplication, i.e., `b` is a weighted sum of the
   columns of `A` with coefficients given by the entries in `x`.  This
   should be implemented with a single for loop over the entries of
   `x`. The test script ``test_exercises1.py`` will also test
   this function.

.. hint::

   It will be useful to use the Python "slice" notation, for
   example::

     A[:, 3]

   will return the 4th (since Python numbers from zero) column of `A`.
   For more information, see the `Numpy documentation on slicing.
   <https://numpy.org/doc/stable/reference/arrays.indexing.html>`_

.. proof:exercise::

   The :func:`cla_utils.exercises1.time_matvecs` function computes
   the execution time for these two implementations for some example
   matrices and compares them with the built-in Numpy matrix-vector
   product. Run this function and examine the output. You should
   observe that the basic implementation is much slower than the
   built-in implementation. This is because built-in Numpy operations
   use compiled C code that is wrapped in Python, which avoids the
   overheads of run-time interpretation of the Python code and
   manipulation of Python objects. Numpy is really useful for
   computational linear algebra programming because it preserves the
   readability and flexibility of Python (writing code that looks much
   more like maths, access to object-oriented programming models)
   whilst giving near-C speed if used appropriately.  You can read
   more about the advantages of using Numpy `here
   <https://numpy.org/devdocs/user/whatisnumpy.html>`_.  You
   should also observe that the column implementation is somewhere
   between the speed of the basic implementation and the built-in
   implementation. This is because (if you did it correctly), each
   iteration of the for loop involves adding an entire array (a
   scaling of one of the columns of `A`) to another array (where `b`
   is being calculated). This will also use compiled C code through
   Numpy, removing some (but not all) of the Python overheads in the
   basic implementation.

   In this course, we will present algorithms in the notes that generally
   do not express the way that Numpy should be used to implement them.
   In these exercises you should consider the best way to make use of
   Numpy built-in operations (which will often make the code more maths-like
   and readable, as well as potentially faster).

Range, nullspace and rank
=========================



   
.. details:: Supplementary video

   .. vimeo:: 450162984


In this section we'll quickly rattle through some definitions and results.

.. proof:definition:: Range

   The range of `A`, `\mbox{range}(A)`, is the set of vectors that can
   be expressed as `Ax` for some `x`.

The next theorem follows as a result of the column space
interpretation of matrix-vector multiplication.

.. proof:theorem::

   `\mbox{range}(A)` is the vector space spanned by the columns of `A`.

.. proof:definition:: Nullspace

   The nullspace `\mbox{null}(A)` of `A` (or kernel) is the set of
   vectors `x` satisfying `Ax=0`, i.e.

   .. math::

      \mbox{null}(A) = \{x \in \mathbb{C}^n: Ax=0\}.



   
.. details:: Supplementary video

   .. vimeo:: 450166119

.. proof:definition:: Rank

   The column rank `\mbox{rank}(A)` of `A` is the dimension of the
   column space of `A`.  The row rank `\mbox{rank}(A)` of `A` is the
   dimension of the row space of `A`. It can be shown that the column
   rank and row rank of a matrix are equal, so we shall just refer
   to the rank.

If

   .. math::

      A = \begin{pmatrix}
      a_1 & a_2 & \ldots & a_n \\
      \end{pmatrix},

the column space of `A` is `\mbox{span}(a_1,a_2,\ldots,a_n)`.

.. proof:definition::

   An `m\times n` matrix `A` is full rank if it has maximum possible rank
   i.e. rank equal to `\min(m, n)`.

If `m\geq n` then `A` must have `n` linearly independent columns to be
full rank. The next theorem is then a consequence of the column space
interpretation of matrix-vector multiplication.

.. proof:theorem::

   An `m\times n` matrix `A` is full rank if and only if it maps no two
   distinct vectors to the same vector.

.. proof:definition::

   A matrix `A` is called nonsingular, or invertible, if it is a square
   matrix (`m=n`) of full rank.

.. proof:exercise::

   The :func:`cla_utils.exercises1.rank2` function has been left
   unimplemented.  To finish the function, add code so that it
   computes the rank-2 matrix `A = u_1v_1^* + u_2v_2^*` from
   `u_1,u_2\in \mathbb{C}^m` and `v_1,v_2 \in \mathbb{C}^n`. As you
   can see, the function needs to implement this rank-2 matrix by
   first forming two matrices `B` and `C` from the inputs,
   and then forming `A` as the product of `B` and `C`. The
   test script ``test_exercises1.py`` in the ``test`` directory will also test this function.

   To measure the rank of `A`, we can use the built-in rank
   function::

     r = numpy.linalg.matrix_rank(A)

   and we should find that the rank is equal to 2. Can you explain why
   this should be the case (use the column space interpretation of
   matrix-matrix multiplication)?

Invertibility and inverses
==========================



   
.. details:: Supplementary video

   .. vimeo:: 450171203

This means that an invertible matrix has columns that form a basis for
`\mathbb{C}^m`. Given the canonical basis vectors defined by

   .. math::

      e_j = \begin{pmatrix}
      0 \\
      \ldots \\
      0 \\
      1 \\
      0 \\
      \ldots \\
      0 \\
      \end{pmatrix},

i.e. `e_j` has all entries zero except for the jth entry which is 1, we can
write

   .. math::

      e_j = \sum_{k=1}^m z_{jk} a_k, \quad 1\leq j \leq m.

In other words,

   .. math::

      I =
      \begin{pmatrix}
      e_1 & e_2 & \ldots & e_m
      \end{pmatrix}

      = ZA.

We call `Z` a (left) inverse of `A`.

.. exercise::

   Prove that `Z` is the unique left inverse of `A`, and prove that `Z`
   is also the unique right inverse of `A`, satisfying `I = AZ`.

We write `Z=A^{-1}`.

The first four parts of the next theorem are a consequence of what
we have so far, and we shall quote the fifth and sixth (see a linear algebra
course).

.. proof:theorem::

   Let `A \in \mathbb{C}^{m\times m}`. Then the following are equivalent.

   #. `A` has an inverse.
   #. `\mbox{rank}(A)=m`.
   #. `\mbox{range}(A)=\mathbb{C}^m`.
   #. `\mbox{null}(A)=\{0\}`.
   #. 0 is not an eigenvalue of `A`.
   #. The determinant `\det(A)\neq 0`.



   
.. details:: Supplementary video

   .. vimeo:: 450172407

Finding the inverse of a matrix can be seen as a change of basis. Considering
the equation `Ax= b`, we have `x = A^{-1}b` for invertible `A`. We have
seen already that `b` can be written as

   .. math::

      b = \sum_{j=1}^m x_j a_j.

Since the columns of `A` span `\mathbb{C}^m`, the entries of `x` thus
provide the unique expansion of `b` in the columns of `A` which form a
basis.  Hence, whilst the entries of `b` give basis coefficients for
`b` in the canonical basis `(e_1,e_2,\ldots,e_m)`, the entries of `x`
give basis coefficients for `b` in the basis given by the columns of `A`.

.. proof:exercise::

   For matrices of the form, `A = I + uv^*`, where `I` is the `m\times
   m` identity matrix, and `u,v \in \mathbb{C}^m`, show that whenever
   `A` is invertible, the inverse is of the form `A^{-1} = I + \alpha uv^*`
   where `\alpha \in \mathbb{C}`, and calculate the form of `\alpha`.

   The :func:`cla_utils.exercises1.rank1pert_inv` function has been
   left unimplemented.  To finish the function, add code so that it
   computes `A^{-1}` using your formula (and not any built-in matrix
   inversion routines). The test script ``test_exercises1.py`` in the
   ``test`` directory will also test this function.

   Add a function to :mod:`cla_utils.exercises1` that measures the
   time to compute the inverse of `A` for an input matrix of size 400,
   and compare with the time to compute the inverse of `A` using the built-in
   inverse::

     numpy.linalg.inv(A)

   What do you observe? Why do you think this is? We will examine the
   cost of general purpose matrix inversion algorithms later.


Adjoints and Hermitian matrices
===============================



   
.. details:: Supplementary video

   .. vimeo:: 450173092

.. proof:definition:: Adjoint

   The adjoint (or Hermitian conjugate) of `A\in \mathbb{C}^{m\times n}`
   is a matrix `A^* \in \mathbb{C}^{n\times m}` (sometimes written
   `A^\dagger` or `A'`), with

   .. math::

      a^*_{ij} = \bar{a_{ji}},

   where the bar denotes the complex conjugate of a complex number. If
   `A^* = A` then we say that `A` is Hermitian.

   For real matrices, `A^*=A^T`. If `A=A^T`, then we say that the matrix
   is symmetric.

The following identity is very important when dealing with adjoints.

.. proof:theorem::

   For matrices `A`, `B` with compatible dimensions (so that they can
   be multiplied),

   .. math::

      (AB)^* = B^*A^*.

.. proof:exercise::

   (This is an advanced exercise if the other exercises are complete.
   If you are behind on the exercises please skip this one.)
   
   Consider a matrix `A=B + iC` where `B,C\in\mathbb{R}^{m\times m}`
   and `A` is Hermitian. Show that `B=B^T` and `C=-C^T`. To save
   memory, instead of storing values of `A` (`m\times m` complex
   numbers to store), consider equivalently storing a real-valued
   `m\times m` array `\hat{A}` with `\hat{A}_{ij}=B_{ij}` for `i\geq j`
   and `\hat{A}_{ij}=C_{ij}` for `i<j`.

   The :func:`cla_utils.exercises1.ABiC` function has been left
   unimplemented. It should implement matrix vector multiplication
   `z=Ax`, returning the real and imaginary parts of `z`, given the
   real and imaginary parts of `x` as inputs, and given the real array
   `\hat{A}` as above. You should implement the multiplication using
   real arithmetic only, with just one loop over the entries of `x`,
   using the column space interpretation of matrix-vector
   multiplication. The test script ``test_exercises1.py`` in the
   ``test`` directory will also test this function.

.. hint::

   You can use the Python "slice" notation, to assign into a slice
   of an array, for example::

     x[3:5] = y[3:5]

   will copy the 4th and 5th entries of `y` (Python numbers from zero,
   and the upper limit of the slice is the first index value not to
   use.  For more information, see the `Numpy documentation on
   slicing.
   <https://numpy.org/doc/stable/reference/arrays.indexing.html>`_


Inner products and orthogonality
================================



   
.. details:: Supplementary video

   .. vimeo:: 450172520

The inner product is a critical tool in computational linear algebra.

.. proof:definition:: Inner product

   Let `x,y\in \mathbb{C}^m`. Then the inner product of `x` and `y` is

   .. math::

      x^*y = \sum_{i=1}^m \bar{x}_iy_i.

(Exercise: check that the inner product is bilinear, i.e. linear in
both of the arguments.)

We will frequently use the natural norm derived from the inner product
to define size of vectors.

.. proof:definition:: 2-Norm

   Let `x\in \mathbb{C}^m`. Then the 2-norm of `x` is

   .. math::

      \|x\| = \sqrt{\sum_{i=1}^m |x_i|^2} = \sqrt{x^*x}.

Orthogonality will emerge as an early key concept in this course.

.. proof:definition:: Orthogonal vectors

   Let `x,y\in \mathbb{C}^m`. The two vectors are orthogonal if
   `x^*y=0`.

   Similarly, let `X`, `Y` be two sets of vectors. The two sets
   are orthogonal if

   .. math::

      x^*y = 0,\quad \forall x\in X, \, y\in Y.

   A set `S` of vectors is itself orthogonal if

   .. math::

      x^*y = 0,\quad\forall x,y \in S.

   We say that `S` is orthonormal if we also have `\|x\|=1`
   for all `x\in S`.

Orthogonal components of a vector
=================================



   
.. details:: Supplementary video

   .. vimeo:: 450184086

Let `S=\{q_1,q_2,\ldots,q_n\}` be an orthonormal set of vectors in
`\mathbb{C}^m`, and take another arbitrary vector `v\in \mathbb{C}^m`.
Now take

   .. math::

      r = v - (q_1^*v)q_1 - (q_2^*v)q_2 - \ldots - (q_n^*v)q_n.

Then, we can check that `r` is orthogonal to `S`, by calculating
for each `1\leq i \leq n`,

   .. math::

      q^*_ir = q_i^*v - (q_1^*v)(q_i^*q_1) - \ldots - (q_n^*v)(q_i^*q_n)

      = q_i^*v - q_i^*v = 0,

since `q_i^*q_j=0` if `i\neq j`, and 1 if `i=j`.
Thus,

   .. math::

      v = r + \sum_{i=1}^n (q_i^*v)q_i
      = r + \sum_{i=1}^n \underbrace{(q_i q_i^*)}_{\mbox{rank-1 matrix}}v.

If `S` is a basis for `\mathbb{C}^m`, then `n=m` and `r=0`, and we have

   .. math::

      v = \sum_{i=1}^m (q_i q_i^*)v.

.. proof:exercise::

   The :func:`cla_utils.exercises2.orthog_cpts` function has been left
   unimplemented. It should implement the above computation, returning
   `r` and the coefficients of the component of `v` in each
   orthonormal direction. The test script ``test_exercises2.py`` in
   the ``test`` directory will test this function.

Unitary matrices
================



   
.. details:: Supplementary video

   .. vimeo:: 450184373

.. proof:definition:: Unitary matrices

   A matrix `Q\in \mathbb{C}^{m\times m}` is unitary if `Q^* =Q^{-1}`.

   For real matrices, a matrix `Q`  is orthogonal if `Q^T=Q^{-1}`.

.. proof:theorem::

   The columns of a unitary matrix `Q` are orthonormal.

.. proof:proof::

   We have `I = Q^*Q`. Then using the column space interpretation
   of matrix-matrix multiplication,

   .. math::

      e_j = Q^*q_j,

   where `q_j` is the jth column of `Q`. Taking row i of `e_j`, we have

   .. math::

      \delta_{ij} = q_i^*q_j, \mbox{ where }
      \delta_{ij} = \left\{
      \begin{array}{ccc}
      1 & \mbox{if} & i=j, \\
      0 & \mbox{otherwise} & \\
      \end{array}\right. .

Extending a theme from earlier, we can interpret `Q^*=Q^{-1}` as
representing a change of orthogonal basis. If `Qx = b`, then
`x=Q^*b` contains the coefficients of `b` expanded in the basis
given by the orthonormal columns of `Q`.

.. proof:exercise::

   The :func:`cla_utils.exercises2.solveQ` function has been left
   unimplemented. Given a square unitary matrix `Q` and a vector `b`
   it should solve `Qx=b` using information above (it is not expected
   to work when `Q` is not unitary or square). The test script
   ``test_exercises2.py`` in the ``test`` directory will test this
   function.

   Add a function to :mod:`cla_utils.exercises2` that measures the
   time to solve `Qx=b` using ``solveQ`` for an input matrix of sizes 100,
   200, 400,
   and compare with the times to solve the equation using the general purpose
   solve (which uses LU factorisation, which we will discuss later)::

     x = numpy.linalg.solve(Q, b)

   What did you expect and was it observed?

   A quick way to get an orthogonal matrix is to take a general matrix $A$
   and find the QR factorisation, which we will cover in the next section.

     Q, R = numpy.linalg.qr(A)

   returns two matrices, of which `Q` is orthogonal.

Vector norms
============



   
.. details:: Supplementary video

   .. vimeo:: 450184674

Various vector norms are useful to measure the size of a vector.
In computational linear algebra we need them for quantifying errors
etc.

.. proof:definition:: Norms

   A norm is a function `\|\cdot\|:\mathbb{C}^m \to \mathbb{R}`, such that

   #. `\|x\|\geq 0`, and `\|x\|=0\implies x =0.`
   #. `\|x+y\| \leq \|x\| + \|y\|` (triangle inequality).
   #. `\|\alpha x\| = |\alpha|\|x\|` for all `x \in \mathbb{C}^m`
      and `\alpha \in \mathbb{C}`.

We have already seen the 2-norm, or Euclidean norm, which is part of a
larger class of norms called p-norms, with

   .. math::

      \|x\|_p = \left(\sum_{i=1}^m |x_i|^p\right)^{1/p}, \quad

for real `p>0`. We will also consider weighted norms

   .. math::

      \|x\|_{W,p} = \|Wx \|_p,

where `W` is a matrix.

Projectors and projections
==========================



   
.. details:: Supplementary video

   .. vimeo:: 450185110

.. proof:definition:: Projector

   A projector `P` is a square matrix that satisfies `P^2=P`.

If `v \in \mbox{range}(P)`, then there exists `x` such that
`Px = v`. Then,

   .. math::

      Pv = P(Px) = P^2x = Px = v,

and hence multiplying by `P` does not change `v`.

Now suppose that `Pv \neq v` (so that `v\notin \mbox{range}(P)`).
Then,

   .. math::

      P(Pv - v) = P^2v - Pv = Pv - Pv = 0,

which means that `Pv-v` is the nullspace of `P`. We have

   .. math::

      Pv -v = -(I-P)v.

.. proof:definition:: Complementary projector

   Let `P` be a projector. Then we call `I-P` the complementary projector.

To see that `I-P` is also a projector, we just calculate,

   .. math::

      (I-P)^2 = I^2 - 2P + P^2 = I - 2P + P = I - P.

If `Pu=0`, then `(I-P)u = u`.

In other words, the nullspace
of `P` is contained in the range of `I-P`.

On the other hand, if `v` is in the range of `I-P`,  then
there exists some `w` such that

   .. math::

      v = (I-P)w = w - Pw.

We have

   .. math::

      Pv = P(w-Pw) = Pw - P^2w = Pw - Pw = 0.

Hence, the range of `I-P` is contained in the nullspace of `P`.
Combining these two results we see that the range of `I-P`
is equal to the nullspace of `P`. Since `P` is the complementary
projector to `I-P`, we can repeat the same argument to show
that the range of `P` is equal to the nullspace of `I-P`.

We see that a projector `P` separates `\mathbb{C}^m` into two
subspaces, the nullspace of `P` and the range of `P`. In fact the
converse is also true: given two subspaces `S_1` and `S_2`
of `\mathbb{C}^m` with `S_1 \cap S_2 = \{0\}`, then there
exists a projector `P` whose range is `S_1` and whose nullspace
is `S_2`.



   
.. details:: Supplementary video

   .. vimeo:: 450185494

Now we introduce orthogonality into the concept of projectors.

.. proof:definition:: Orthogonal projector

   `P` is an orthogonal projector if

   .. math::

      (Pv)^*(Pv-v) = 0, \, \forall v \in \mathbb{C}^m.

In this case, `P` separates the space into two orthogonal subspaces.

Constructing orthogonal projectors from sets of orthonormal vectors
===================================================================

Let `\{q_1,\ldots,q_n\}` be an orthonormal set of vectors in
`\mathbb{C}^m`. We write

.. math::

   \hat{Q} = \begin{pmatrix}
   q_1 & q_2 & \ldots & q_n \\
   \end{pmatrix}.

Previously we showed that for any `v\in \mathbb{C}^m`, we have

.. math::

   v = \underbrace{r}_{\mbox{Orthogonal to column space of }\hat{Q}} +
   \underbrace{\sum_{i=1}^n (q_iq^*_i)v}_{\mbox{in the column space of }\hat{Q}}.

Hence, the map

.. math::

   v \mapsto Pv = \underbrace{\sum_{i=1}^n (q_iq^*_i)}_{=P}v,

is an orthogonal projector. In fact, `P` has very simple form.

.. _orthogonal_projector:

.. proof:theorem::

   The orthogonal projector `P` takes the form

   .. math::

      P = \hat{Q}\hat{Q}^*.

.. proof:proof::

   From the change of basis interpretation of multiplication by
   `\hat{Q}^*`, the entries in `\hat{Q}^*v` gives coefficients of the
   projection of `v` onto the column space of `\hat{Q}` when expanded
   using the columns as a basis. Then, multiplication by `\hat{Q}`
   gives the projection of `v` expanded again in the canonical basis.
   Hence, multiplication by `\hat{Q}\hat{Q}^*` gives exactly the same
   result as multiplication by the formula for `P` above.

This means that `\hat{Q}\hat{Q}^*` is an orthogonal projection onto
the range of `\hat{Q}`. The complementary projector is `P_{\perp} =
I - \hat{Q}\hat{Q}^*` is an orthogonal projection onto the nullspace
of `\hat{Q}`.

An important special case is when `\hat{Q}` has just one column,
and then

.. math::

   P = q_1q_1^*, \, P_{\perp}=I - q_1q_1^*.

We notice that `P^* = (\hat{Q}\hat{Q}^*) = \hat{Q}\hat{Q}^* = P`.
In fact the following is true.

.. proof:theorem::

   `P=P^*` if and only if `Q` is orthogonal.

.. proof:exercise::

   The :func:`cla_utils.exercises2.orthog_proj` function has been left
   unimplemented. Given an orthonormal set `q_1,q_2,\ldots,q_n`, it
   should provide the orthogonal projector `P`. The test script
   ``test_exercises2.py`` in the ``test`` directory will also test
   this function.
