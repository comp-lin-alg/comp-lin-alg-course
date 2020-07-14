.. default-role:: math

Finding eigenvalues of matrices
===============================

We start with some preliminary terminology.  A vector `x\in
\mathbb{C}^m` is an *eigenvector* of a square matrix `A\in
\mathbb{C}^{m\times m}` with *eigenvalue* `\lambda` if `Ax=\lambda
x`. An eigenspace is the subspace `E_{\lambda}\subset\mathbb{C}^m`
containing all eigenvectors of `A` with eigenvalue `\lambda`.

There are a few reasons why we are interested in computing
eigenvectors and eigenvalues of a matrix `A`.

#. Eigenvalues and eigenvectors encode information about `A`.
#. Eigenvalues play an important role in stability calculations
   in physics and engineering.
#. We can use eigenvectors to underpin the solution of linear systems
   involving `A`.
#. ...

How to find eigenvalues?
------------------------

The method that we first encounter in our mathematical education is to
find solutions of `(A-\lambda I)x = 0`, which implies that
`\det(A-\lambda I)=0`. This gives a degree `m` polynomial to solve for
`\lambda`, called the *characteristic polynomial*. Unfortunately,
there is no general solution for polynomials of degree 5 or greater
(from Galois theory). Further, the problem of finding roots of
polynomials is numerically unstable. All of this means that we should
avoid using polynomials finding eigenvalues. Instead, we should try to
apply transformations to the matrix `A` to a form that means that the
eigenvalues can be directly extracted.

The eigenvalue decomposition of a matrix `A` finds a nonsingular matrix
`X` and a diagonal matrix `\Lambda` such that

   .. math::

      A = X\Lambda X^{-1}.

The diagonal entries of `\Lambda` are the eigenvalues of `A`. Hence,
if we could find the eigenvalue decomposition of `A`, we could just
read off the eigenvalues of `A`; the eigenvalue decomposition is
"eigenvalue revealing". Unfortunately, it is not always easy or even
possible to transform to an eigenvalue decomposition. Hence we shall
look into some other eigenvalue revealing decompositions of `A`.

We quote the following result that explains when an eigenvalue
decomposition can be found.

.. proof:theorem::

   An `m\times m` matrix `A` has an eigenvalue decomposition if and
   only if it is non-defective, meaning that the geometric
   multiplicity of each eigenvalue (the dimension of the eigenspace
   for that eigenvalue) is equal to the algebraic multiplicity (the
   number of times that the eigenvalue is repeated as a root in the
   characteristic polynomial `\det(I\lambda - A)=0`.

If the algebraic multiplicity is greater than the geometric
multiplicity for any eigenvalue of `A`, then the matrix is defective,
the eigenspaces do not span `\mathbb{C}^m`, and an eigenvalue
decomposition is not possible.

This all motivates the search for other eigenvalue revealing
decompositions of `A`.

.. proof:definition:: Similarity transformations

   For `X\in \mathbb{C}^{m\times m}` a nonsingular matrix, the map
   `A\mapsto X^{-1}AX` is called a similarity transformation of `A`.
   Two matrices `A` and `B` are *similar* if `B=X^{-1}AX`.

The eigenvalue decomposition shows that (when it exists), `A` is similar
to `\Lambda`. The following result shows that it may be useful to examine
other similarity transformations.

.. proof:theorem::

   Two similar matrices `A` and `B` have the same characteristic polynomial,
   eigenvalues, and geometric multiplicities.

.. proof:proof::

   See a linear algebra textbook.

The goal is to find a similarity transformation such that `A` is
transformed to a matrix `B` that has some simpler structure where the
eigenvalues can be easily computed (with the diagonal matrix of the
eigenvalue decomposition being one example).

One such transformation comes from the Schur factorisation.

.. proof:definition:: Schur factorisation

   A Schur factorisation of a square matrix `A` takes the form `A =
   QTQ^*`, where `Q` is unitary (and hence `Q^*=Q^{-1}`) and `T` is
   upper triangular.

It turns out that, unlike the situation for the eigenvalue
decomposition, the following is true.

.. proof:theorem::

   Every square matrix has a Schur factorisation.

This is useful, because the characteristic polynomial of an upper
triangular matrix is just `\prod_{i=1}^m (\lambda-T_{ii})`, i.e.  the
eigenvalues of `T` are the diagonal entries
`(T_{11},T_{22},\ldots,T_{mm})`. So, if we can compute the Schur
factorisation of `A`, we can just read the eigenvalues from the diagonal
matrices of `A`.

There is a special case of the Schur factorisation, called the unitary
diagonalisation

.. proof:definition:: Unitary diagonalisation

   A unitary diagonalisation of a square matrix `A` takes the form `A =
   Q\Lambda Q^*`, where `Q` is unitary (and hence `Q^*=Q^{-1}`) and `\Lambda`
   is diagonal.

A unitary diagonalisation is a Schur factorisation *and* an eigenvalue
decomposition.
   
.. proof:theorem::

   A Hermitian matrix is unitary diagonalisable, with real `\Lambda`.

Hence, if we have a Hermitian matrix, we can follow a Schur
factorisation strategy (such as we shall develop in this section), and
obtain an eigenvalue decomposition as a bonus.

Transformations to Schur factorisation
--------------------------------------

Just as for the QR factorisations, we will compute the Schur
factorisation successively, with multiplication by a sequence of
unitary matrices `Q_1,Q_2,\ldots`. There are two differences for the
Schur factorisation. First, the matrices must be multiplied not just
on the left but also on the right with the inverse, i.e.

   .. math::

      A \mapsto \underbrace{Q_1^*AQ_1}_{A_1} \mapsto \underbrace{Q_2^*Q_1^*AQ_2Q_1}_{A_2}, \ldots

At each stage, we have a similarity transformation,

   .. math::

      A = \underbrace{Q_1Q_2\ldots Q_k}_{=Q}A_k\underbrace{Q_k^*\ldots Q_2^*Q_1^*}_{=Q^*},

\emph{i.e.} `A` is similar to `A_k`. Second, the successive sequence is
infinite, i.e. we will develop an iterative method that converges in
the limit `k\to\infty`.  We should terminate the iterative method
when `A_k` is sufficiently close to being upper triangular (which

We should not be surprised by this news, since if the successive
sequence were finite, we would have derived an explicit formula for
computing the eigenvalues of the characteristic polynomial of `A`
which is explicit in general. 

In fact, there are two stages to this process. The first stage, which
is finite (takes `m-1` steps) is to use similarity transformations to
upper Hessenberg form (`H_{ij}=0` for `i>j+1`). If `A` is Hermitian,
then `H` will be tridiagonal. This stage is not essential but it makes
the second, iterative, stage much faster.

Similarity transformation to upper Hessenberg form
--------------------------------------------------

We already know how to use a unitary matrix to set all entries to zero
below the diagonal in the first column of `A` by left multiplication
`Q^*_1A`, because this is the Householder algorithm. The problem is
that we then have to right multiply by `Q_1` to make it a similarity
transformation, and this puts non-zero entries back in the column
again. The easiest way to see this is to write
`Q_1^*AQ_1=(Q_1^*(Q_1^*A)^*)^*`. `(Q_1^*A)^*` has zeros in the first
row to the right of the first entry. Then, `Q_1^*(Q_1^*A)` creates
linear combinations of the first column with the other columns,
filling the zeros in with non-zero values again. Then finally taking
the adjoint doesn't help with these non-zero values. Again, we
shouldn't be surprised that this is impossible, because if it was,
then we would be able to build a finite procedure for computing
eigenvalues of the characteristic polynomial, which is impossible in
general.

A slight modification of this idea (and the reason that we can
transform to upper Hessenberg form) is to use a Householder rotation
`Q_1^*` to set all entries to zero below the *second* entry in the
first column. This matrix leaves the first row unchanged, and hence
right multiplication by `Q_1` leaves the first column unchanged. We
can create zeros using `Q_1^*` and `Q_1` will not destroy them. This
procedure is then repeated with multiplication by `Q_2^*`, which
leaves the first two rows unchanged and puts zeros below the third
entry in the second column, which are not spoiled by right
multiplication by `Q_2`. Hence, we can transform `A` to a similar
upper Hessenberg matrix `H` in `m-2` iterations.

This reduction to Hessenberg form can be expressed in the following
pseudo-code.

* FOR `k=1` TO `m-2`

  * `x\gets A_{k+1:m,k}`
  * `v_k\gets \mbox{sign}(x_1)\|x\|_2e_1 + x`
  * `v_k\gets v_k/\|v\|_2`
  * `A_{k+1:m,k:m} \gets A_{k+1:m,k:m}- 2v_k(v_k^*A_{k+1:m,k:m})`
  * `A_{k:m,k+1:m} \gets A_{k:m,k+1:m}- 2(A_{k:m,k+1:m,k:m}v_k)v_k^*`
* END FOR

Note the similarities and differences with the Householder algorithm
for computing the QR factorisation.

To calculate the operation count, we see that the algorithm is
dominated by the two updates to `A`, the first of which applies a
Householder reflection to rows from the left, and the second applies
the same reflections to columns from the right.

The left multiplication applies a Householder
reflection to the last `m-k` rows, requiring 4 FLOPs per
entry. However, these rows are zero in the first `k-1` columns,
so we can skip these and just work on the last `m-k+1` entries
of each of these rows.

Then, the total operation count for the left multiplication is

   .. math::

      4 \times \sum_{k=1}^{m-1} (m-k)(m-k+1) \sim \frac{4}{3}m^3.
  
