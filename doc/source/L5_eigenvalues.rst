.. default-role:: math

Finding eigenvalues of matrices
===============================

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454117340>`_.

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

.. proof:example::

   Consider the `m\times m` diagonal matrix

      .. math::

	 A = \begin{pmatrix}
	 1 & 0 & \ldots & \ldots & 0 \\
	 0 & 2 & \ldots & \ldots &0 \\
	 0 & 0 & \ddots & \ddots & \vdots \\
	 \vdots & \vdots & \ldots & \ddots & \vdots \\
	 0 & 0 & \ldots & \ldots & m \\
	 \end{pmatrix}.
	 
   The characteristic polynomial of `A` is

      .. math::

	 (\lambda-1)(\lambda-2)\ldots(\lambda-m),

   and the eigenvalues are clearly `1,2,3,\ldots,m`. This is called
   the Wilkinson Polynomial. Numpy has some tools for manipulating
   polynomials which are useful here. When an `m\times m` array is
   passed in to :func:`numpy.poly`, it returns an array of
   coefficients of the polynomial. For `m=20`, obtain this array and
   then perturb the coefficients `a_k \to \tilde{a}_k =
   a_k(1+10^{-10}r_k)` where `r_k` are randomly sampled normally
   distributed numbers with mean 0 and variance 1. :func:`numpy.roots`
   will compute the roots of this perturbed polynomial. Plot these
   roots as points in the complex plane. Repeat this 100 times,
   superposing the root plots on the same graph. What do you observe?
   What does it tell you about this problem and what should we
   conclude about the wisdom of finding eigenvalues using
   characteristic polynomials?

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454118485>`_.

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

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454122744>`_.

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

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454122918>`_.

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

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454123177>`_.

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

i.e. `A` is similar to `A_k`. Second, the successive sequence is
infinite, i.e. we will develop an iterative method that converges in
the limit `k\to\infty`.  We should terminate the iterative method
when `A_k` is sufficiently close to being upper triangular (which
can be measured by checking some norm on the lower triangular part
of `A` and stopping when it is below a tolerance).

We should not be surprised by the news that the method needs to be
iterative, since if the successive sequence were finite, we would have
derived an explicit formula for computing the eigenvalues of the
characteristic polynomial of `A` which is explicit in general.

In fact, there are two stages to this process. The first stage, which
is finite (takes `m-1` steps) is to use similarity transformations to
upper Hessenberg form (`H_{ij}=0` for `i>j+1`). If `A` is Hermitian,
then `H` will be tridiagonal. This stage is not essential but it makes
the second, iterative, stage much faster.

Similarity transformation to upper Hessenberg form
--------------------------------------------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454123306>`_.

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

.. proof:exercise::

   The :func:`cla_utils.exercises8.Q1AQ1s` function has been left
   uncompleted. It should apply the Householder transformation `Q_1`
   to the input `A` (without forming `Q_1` of course) that transforms
   the first column of `A` to have zeros below the diagonal, and then
   apply a transformation equivalent to right multiplication by
   `Q_1^*` (again without forming `Q_1`).  The test script
   ``test_exercises8.py`` in the ``test`` directory will test this
   function.

   Experiment with the output of this function. What happens to the
   first column?
      
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

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454123643>`_.

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

.. proof:exercise::

   The :func:`cla_utils.exercises8.hessenberg` function has been left
   unimplemented. It should implement the algorithm above, using only
   one loop over `k`. It should work "in-place", changing the input
   matrix. At the left multiplication, your implementation should
   exploit the fact that zeros do not need to be recomputed where
   there are already expected to be zeros. The test script
   ``test_exercises8.py`` in the ``test`` directory will test this
   function.

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
  
The right multiplication does the same operations but now there are no
zeros to take advantage of, so all `m` entries in the each of the last
`m-k` columns need to be manipulated. With 4 FLOPs per entry, this becomes

   .. math::

      4\times \sum_{k=1}^{m-1} m(m-k) \sim \frac{10}{3}m^3 FLOPs.

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454123926>`_.
      
In the Hermitian case, the Hessenberg matrix becomes tridiagonal, and
these extra zeros can be exploited, leading to an operation count
`\sim 4m^3/3`.

It can be shown that this transformation to a Hessenberg matrix is
backwards stable, i.e. in a floating point implementation, it gives
`\tilde{Q},\tilde{H}` such that

   .. math::

      \tilde{Q}\tilde{H}\tilde{Q}^* = A + \delta A, \,
      \frac{\|\delta A\|}{\|A\|}=\mathcal{O}(\varepsilon),

for some `\delta A`.

.. proof:exercise::

   The :func:`cla_utils.exercises8.hessenbergQ` function has been left
   unimplemented. It should implement the Hessenberg algorithm again
   (you can just copy paste the code from the previous exercise) but
   it should also return the matrix `Q` such that `QHQ^*=A`. You need
   to work out how to alter the algorithm to construct this. The test
   script ``test_exercises8.py`` in the ``test`` directory will test
   this function.

.. proof:exercise::

   The :func:`cla_utils.exercises8.ev` function has been left
   unimplemented. It should return the eigenvectors of
   `A` by first reducing to Hessenberg form, using the functions you
   have already created, and then calling
   :func:`cla_utils.exercises8.hessenberg_ev`, which computes the
   eigenvectors of upper Hessenberg matrices (do not
   edit this function!). The test script ``test_exercises8.py`` in the
   ``test`` directory will test this function.

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454124279>`_.

In the next few sections we develop the iterative part of the
transformation to the upper triangular matrix `T`. This algorithm
works for a broad class of matrices, but the explanation is much
easier for the case of real symmetric matrices, which have real
eigenvalues and orthogonals eigenvectors (which we shall normalise to
`\|q_i\|=1`, `i=1,2,\ldots,m`). The idea is that we will have already
transformed to Hessenberg form, which will be tridiagonal in this
case. Before describing the iterative transformation, we will discuss
a few key tools in explaining how it works.

Rayleigh quotient
-----------------

The first tool that we shall consider is the Rayleigh quotient. If
`A\in \mathbb{C}^{m\times m}` is a real symmetric matrix, then the
Rayleigh quotient of a vector `x \in \mathbb{C}^{m}` is defined as

   .. math::

      r(x) = \frac{x^TAx}{x^Tx}.

If `x` is an eigenvector of `A`, then

   .. math::

      r(x) = \frac{x^T\lambda x}{x^Tx} = \lambda,

i.e. the Rayleigh quotient gives the corresponding eigenvalue.

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454124455>`_.

If `x` is not exactly an eigenvector of `A`, but is just close to one,
we might hope that `r(x)` is close to being an eigenvalue. To
investigate this we will consider the Taylor series expansion of
`r(x)` about an eigenvector `q_J` of `A`. We have

   .. math::

      \nabla r(x) = \frac{2}{x^Tx}\left(Ax-r(x)x\right),

which is zero when `x=q_J`, because then `r(q_J)=\lambda_J`:
eigenvectors of `A` are stationary points of `r(x)`! Hence, the Taylor
series has vanishing first order term,

   .. math::

      r(x) = r(q_J) + (x-q_J)^T\underbrace{\nabla r(q_J)}_{=0}
      + \mathcal{O}(\|x-q_J\|^2),

i.e.

   .. math::

      r(x) - r(q_J) = \mathcal{O}(\|x-q_J\|^2), \quad \mbox{as }
      x \to q_J.

The Rayleigh quotient gives a quadratically accurate estimate to the
eigenvalues of `A`.

.. proof:exercise::

   Add a function to :mod:`cla_utils.exercises8` that investigates
   this property by:

   #. Forming a Hermitian matrix `A`,
   #. Finding an eigenvector `v` of `A` with eigenvalue `\lambda` (you can use :func:`numpy.linalg.eig` for this),
   #. Choosing a perturbation vector `r`, and perturbation parameter `\epsilon>0`,
   #. Comparing the Rayleigh quotient of `v + \epsilon r` with `\lambda`,
   #. Plotting (on a log-log graph, use :func:`matplotlib.pyplot.loglog`) the error in estimating the eigenvalue as a function of `\epsilon`.

   The best way to do this is to plot the computed data values as points,
   and then superpose a line plot of `a\epsilon^k` for appropriate
   value of `k` and `a` chosen so that the line appears not to far away
   from the points on the same scale. This means that we can check
   by eye if the error is scaling with `\epsilon` at the expected rate.
   
Power iteration
---------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454124710>`_.

Power iteration is a simple method for finding the eigenvalue of
`A` with largest eigenvalue (in magnitude). It is based on the following
idea. We expand a vector `v` in eigenvectors of `A`,

   .. math::

      v = a_1q_1 + a_2q_2 + \ldots a_mq_m,

where we have ordered the eigenvalues so that `|\lambda_1|\geq |\lambda_2|
\geq |\lambda 3| \geq \ldots \geq |\lambda_m`.

Then,

   .. math::

      Av = a_1\lambda_1q_1 + a_2\lambda_2q_2 + \ldots a_m\lambda_m q_m,

and hence, repeated applications of `A` gives

   .. math::

      A^kv = \underbrace{AA\ldots A}_{k\mbox{ times}}v

      = a_1\lambda^k_1q_1 + a_2\lambda^k_2q_2 + \ldots a_m\lambda^k_m q_m.

If `|\lambda_1|>|\lambda_2|`, then provided that `a_1=q_1^Tv\neq 0`,
the first term `a_1\lambda^k_1q_1` rapidly becomes larger than all of
the others, and so `A^kv \approx a_1\lambda^k_1 q_1`, and we can
normalise to get `q_1 \approx A^kv/\|A^kv\|`. To keep the magnitude of
the estimate from getting too large or small (depending on the size of
`\lambda_1` relative to 1), we can alternately apply `A` and normalise,
which gives the power iteration. Along the way, we can use the Rayleigh
quotient to see how our approximation of the eigenvalue is going.

* Set `v_0` to some initial vector (hoping that `\|q_1^Tv_0\|>0`).
* FOR `k=1,2,\ldots`

  * `w\gets Av^{k-1}`,
  * `v^k\gets w/\|w\|`,
  * `\lambda^{(k)} \gets (v^k)^TAv^k`.

Here we have used the fact that `\|v^k\|=1`, so there is no need to
divide by it in the Rayleigh quotient. We terminate the power
iteration when we decide that the changes in `\lambda` indicate
that the error is small. This is guided by the following result.

.. _power_iteration:

.. proof:theorem::

   If `|\lambda_1|> |\lambda_2|` and `\|q_1^Tv_0\|>0`, then after
   `k` iterations of power iteration, we have

      .. math::

	 \|v^k - \pm q_1\| = \mathcal{O}\left(
	 \left|\frac{\lambda_2}{\lambda_1}\right|^k\right),
	 \quad |\lambda^{(k)} - \lambda_1|=
	 \mathcal{O}\left(\left|\frac{\lambda_2}{\lambda_1}\right|^{2k}\right),

   as `k\to\infty`. At each step `\pm` we mean that the result holds
   for either `+` or `-`.

.. proof:proof::

   We have already shown the first equation using the Taylor series, and
   the second equation comes by combining the Taylor series error with
   the Rayleigh quotient error.

The `\pm` feature is a bit annoying, and relates to the fact that the
normalisation does not select `v^k` to have the direction as `q_1`.

.. proof:exercise::

   The :func:`cla_utils.exercises9.pow_it` function has been left
   unimplemented. It should apply power iteration to a given matrix
   and initial vector, according to the docstring. The test script
   ``test_exercises9.py`` in the ``test`` directory will test this
   function.

.. proof:exercise::

   The functions :func:`cla_utils.exercises9.A3` and
   :func:`cla_utils.exercises9.B3` each return a 3x3 matrix, `A_3` and
   `B_3` respectively. Apply :func:`cla_utils.exercises9.pow_it` to
   each of these functions. What differences in behaviour do you
   observe? What is it about `A_3` and `B_3` that causes this?

   
Inverse iteration
-----------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454124799>`_.

Inverse iteration is a modification of power iteration so that we can
find eigenvalues other than `\lambda_1`. To do this, we use the fact
that eigenvectors `q_j` of `A` are also eigenvectors of `(A - \mu
I)^{-1}` for any `\mu\in \mathbb{R}` not an eigenvalue of `A`
(otherwise `A-\mu I` is singular). To show this, we write

   .. math::

      (A - \mu I)q_j = (\lambda_j - \mu)q_j
      \implies (A - \mu I)^{-1}q_j = \frac{1}{\lambda_j - \mu}q_j.

Thus `q_j` is an eigenvalue of `(A - \mu I)^{-1}` with eigenvalue
`1/(\lambda_j - \mu)`. We can then apply power iteration to `(A-\mu
I)^{-1}` (which requires a matrix solve per iteration), which
converges to an eigenvector `q` for which `1/|\lambda-\mu|` is
smallest, where `\lambda` is the corresponding eigenvalue. In other
words, we will find the eigenvector of `A` whose eigenvalue is closest
to `\mu`.

This algorithm is called inverse iteration, which we express in
pseudo-code below.

* `v^{0}\gets` some initial vector with `\|v^0\|=1`.

* FOR `k=1,2,\ldots`

  * SOLVE `(A-\mu I)w = v^{k-1}` for `w`
  * `v^k\gets w/\|w\|`
  * `\lambda^{(k)} \gets (v^k)^TAv^k`

We can then directly extend :numref:`Theorem
{number}<power_iteration>` to the inverse iteration algorithm.
We conclude that the convergence rate is not improved relative
to power iteration, but now we can "dial in" to different
eigenvalues by choosing `\mu`.

.. proof:exercise::

   The :func:`cla_utils.exercises9.inverse_it` function has been left
   unimplemented. It should apply inverse iteration to a given matrix
   and initial vector, according to the docstring. The test script
   ``test_exercises9.py`` in the ``test`` directory will test this
   function.

.. proof:exercise::

   Using the `A_3` and `B_3` matrices, explore the inverse iteration
   using different values of `\mu`. What do you observe?

Rayleigh quotient iteration
---------------------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454303115>`_.

Since we can use the Rayleigh quotient to find an approximation of an
eigenvalue, and we can use an approximation of an eigenvalue to find
the nearest eigenvalue using inverse iteration, we can combine them
together. The idea is to start with a vector, compute the Rayleigh
quotient, use the Rayleigh quotient for `\mu`, then do one step of
inverse iteration to give an updated vector which should now be closer
to an eigenvector. Then we iterate this whole process. This is called
the Rayleigh quotient iteration, which we express in pseudo-code
below.

   * `v^{0}` some initial vector with `\|v^0\|=1`.
   * `\lambda^{(0)} \gets (v^0)^TAv^0`
   * FOR `k=1,2,\ldots`
  
     * SOLVE `(A-\lambda^{(k-1)} I)w = v^{k-1}` for `w`
     * `v^k\gets w/\|w\|`
     * `\lambda^{(k)} \gets (v^k)^TAv^k`

This dramatically improves the convergence since if
`\|v^k-q_J\|=\mathcal{O}(\delta)` for some small `\delta`, then the
Rayleigh quotient gives `|\lambda^{(k)}-q_J|=\mathcal{O}(\delta^2)`.
Then, inverse iteration gives an estimate

.. math::

   \|v^{k+1}-\pm q_J\| = \mathcal{O}(|\lambda^{(k)}-\lambda_J|
   \|v^k-q_J\|) = \mathcal{O}(\delta^3).

Thus we have cubic convergence, which is super fast!

.. proof:exercise::

   The :func:`cla_utils.exercises9.rq_it` function has been left
   unimplemented. It should apply inverse iteration to a given matrix
   and initial vector, according to the docstring. The test script
   ``test_exercises9.py`` in the ``test`` directory will test this
   function.

.. proof:exercise::

   The interfaces to :func:`cla_utils.exercises9.inverse_it` and
   :func:`cla_utils.exercises9.rq_it` have been designed to optionally
   provide the iterated values of the eigenvector and eigenvalue.  For
   a given initial condition (and choice of `\mu` in the case of
   inverse iteration), compare the convergence speeds of the
   eigenvectors and eigenvalues, using some example matrices of
   different sizes (don't forget to make them Hermitian).

The pure QR algorithm
---------------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454124953>`_.

We now describe the QR algorithm, which will turn out to be an
iterative algorithm that converges to the diagonal matrix (upper
triangular matrix for the general nonsymmetric case) that `A` is
similar to. Why this works is not at all obvious at first, and
we shall explain this later. For now, here is the algorithm
written as pseudo-code.

* `A^{(0)} \gets A`
* FOR `k=1,2,\ldots`

  * FIND `Q^{(k)},R^{(k)}` such that `Q^{(k)}R^{(k)}=A^{(k-1)}` (USING QR FACTORISATION)
  * `A^{(k)} = R^{(k)}Q^{(k)}`

.. proof:exercise::

   The :func:`cla_utils.exercises9.pure_QR` function has been left
   unimplemented. It should implement the pure QR algorithm as above,
   using your previous code for finding the QR factorisation using
   Householder transformations. You should think about avoiding
   unecessary allocation of new numpy arrays inside the loop. The
   method of testing for convergence has been left as well. Have a
   think about how to do this and document your implementation. The
   test script ``test_exercises9.py`` in the ``test`` directory will
   test this function.

.. proof:exercise::

   Investigate the behaviour of the pure QR algorithm applied to the
   functions provided by :func:`cla_utils.exercises9.get_A100`,
   :func:`cla_utils.exercises9.get_B100`,
   :func:`cla_utils.exercises9.get_C100`, and 
   :func:`cla_utils.exercises9.get_D100`. You can use
   :func:`matplotlib.pyplot.pcolor` to visualise the entries,
   or compute norms of the components of the matrices below the diagonal,
   for example. What do you observe? How does this relate to the structure
   of the four matrices?
    
The algorithm simply finds the QR factorisation of `A`, swaps Q and R,
and repeats. We call this algorithm the "pure" QR algorithm, since it
can be accelerated with some modifications that comprise the
"practical" QR algorithm that is used in practice.

We can at least see that this is computing similarity transformations since

   .. math::

      A^{(k)} = R^{(k)}Q^{(k)} = (Q^{(k)})^*Q^{(k)}R^{(k)}Q^{(k)} = (Q^{(k)})^*A^{(k-1)}Q^{(k)},

so that `A^{(k)}` is similar to `A^{(k-1)}` and hence to `A^{(k-2)}` and all
the way back to `A`. But why does `A^{(k)}` converge to a diagonal matrix?
To see this, we have to show that the QR algorithm is equivalent to
another algorithm called simultaneous iteration.

Simultaneous iteration
----------------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454125180>`_.

One problem with power iteration is that it only finds one
eigenvector/eigenvalue pair at a time. Simultaneous iteration is a
solution to this. The starting idea is simple: instead of working on
just one vector `v`, we pick a set of linearly independent vectors
`v_1^{0},v_2^{0},\ldots,v_n^{0}` and repeatedly apply `A` to each of
these vectors. After a large number applications and normalisations in
the manner of the power iteration, we end up with a linear independent
set `v_1^{k},v_2^{k},\ldots,v_n^{k}`, `n\leq m`. All of the vectors in this set
will be very close to `q_1`, the eigenvector with largest magnitude of
corresponding eigenvalue. We can choose `v_1^{k}` as our approximation
of `q_1`, and project this approximation of `q_1` from the rest of the
vectors `v_2^{k},v_3^{k},\ldots v_m^{k}`.  All the remaining vectors
will be close to `q_2`, the eigenvector with the next largest
magnitude of corresponding eigenvalue. Similarly we can choose the
first one of the remaining projected vectors as an approximation of
`q_2` and project it again from the rest.

We can translate this idea to matrices by defining `V^{(0)}` to be the
matrix with columns given by the set of initial `v`s. Then after `k`
applications of `A`, we have `V^{(k)}=A^{k} V^{(0)}`. By the column space
interpretation of matrix-matrix multiplication, each column of `V^{(k)}`
is `A^{k}` multiplied by the corresponding column of `V^{(0)}`. To make the
normalisation and projection process above, we could just apply the
Gram-Schmidt algorithm, sequentially forming an orthonormal spanning
set for the columns of `V^{(k)}` working from left to right.  However, we
know that an equivalent way to do this is to form the (reduced) QR
factorisation of `V^{(k)}`, `\hat{Q}^{(k)}\hat{R}^{(k)}=V^{(k)}`; the columns of
`\hat{Q}^{(k)}` give the same orthonormal spanning set.  Hence, the
columns of `\hat{Q}^{(k)}` will converge to eigenvectors of `A`, provided
that:

#. The first `n` eigenvalues of `A` are distinct in absolute value:
   `|\lambda_1| > |\lambda_2| > \ldots > |\lambda_n|`. If we want to find
   all of the eigenvalues `n=m`, then all the absolute values of the
   eigenvalues must be distinct.
#. The `v` vectors can be expressed as a linear sum of the first `n`
   eigenvectors `q_1,\ldots,q_n` in a non-degenerate way. This turns
   out  to be equivalent (we won't show it here) to the condition that
   `\hat{Q}^TV^{(0)}` has an LU factorisation (where `\hat{Q}` is the
   matrix whose columns are the first `n` eigenvectors of `A`).

One problem with this idea is that it is not numerically stable.  The
columns of `V^{(k)}` rapidly become a very ill-conditioned basis for the
spanning space of the original independent set, and the values of
eigenvectors will be quickly engulfed in rounding errors. There is a
simple solution to this though, which is to orthogonalise after
each application of `A`. This is the simultaneous iteration algorithm,
which we express in the following pseudo-code.

* TAKE A UNITARY MATRIX `\hat{Q}^{(0)}`
* FOR `k=1,2,\ldots`

  * `Z\gets A\hat{Q}^{(k-1)}`
  * FIND `Q^{(k)},R^{(k)}` such that `Q^{(k)}R^{(k)}=Z` (USING QR FACTORISATION)

This is mathematically equivalent to the process we described above,
and so it converges under the same two conditions listed above.
    
We can already see that this looks rather close to the QR algorithm.
The following section confirms that they are in fact equivalent.

The pure QR algorithm and simultaneous iteration are equivalent
---------------------------------------------------------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454125393>`_.

To be precise, we will show that the pure QR algorithm is equivalent
to simultaneous iteration when the initial independent set is the
canonical basis `I`, i.e. `Q^{(0)}=I`. From the above, we see that
that algorithm converges provided that `Q^T` has an LU decomposition,
where `Q` is the limiting unitary matrix that simultaneous iteration
is converging to.  To show that the two algorithms are equivalent, we
append them with some auxiliary variables, which are not needed for
the algorithms but are needed for the comparison.

To simultaneous iteration we append a running similarity transformation
of `A`, and a running product of all of the `R` matrices.

* `{Q'}^{(0)} \gets I`
* FOR `k=1,2,\ldots`

  * `Z\gets A{Q'}^{(k-1)}`
  * FIND `{Q'}^{(k)},R^{(k)}` such that `{Q'}^{(k)}R^{(k)}=Z` (USING QR FACTORISATION)
  * `A^{(k)} = ({Q'}^{(k)})^TA{Q'}^{(k)}`
  * `{R'}^{(k)} = R^{(k)}R^{(k-1)}\ldots R^{(1)}`

To the pure QR factorisation we append a running product of the `Q^{k}`
matrices, and a running product of all of the `R` matrices (again).

* `A^{(0)} \gets A`
* FOR `k=1,2,\ldots`

  * FIND `Q^{(k)},R^{(k)}` such that `Q^{(k)}R^{(k)}=A^{(k-1)}` (USING QR FACTORISATION)
  * `A^{(k)} = R^{(k)}Q^{(k)}`
  * `{Q'}^{(k)} = Q^{(1)}Q^{(2)}\ldots Q^{(k)}`
  * `{R'}^{(k)} = R^{(k)}R^{(k-1)}\ldots R^{(1)}`

.. proof:theorem:: pure QR and simultaneous iteration with `I` are equivalent

   The two processes above generate identical sequences of matrices
   `{R'}^{(k)}`, `{Q'}^{(k)}` and `A^{(k)}`, which are related by
   `A^{k} = {Q'}^{(k)}{R'}^{(k)}` (the `k`-th power of `A`, not
   `A^{(k)}`!), and `A^{(k)}=({Q'}^{(k)})^TA{Q'}^{(k)}`.


.. proof:proof::

   We prove by induction. At `k=0`, `A_k={R'}^{(k)}={Q'}^{(k)}=0`. Now
   we assume that the inductive hypothesis is true for `k`, and aim to
   deduce that it is true for `k+1`.

   For simultaneous iteration, we immediately have the simularity
   formula for `A^{(k)}` by definition, and we just need to verify the QR
   factorisation of `A^k`. From the inductive hypothesis,

      .. math::

	 A^k = AA^{k-1} = A{Q'}^{(k-1)}{R'}^{(k-1)}
	 = Z{R'}^{(k-1)} = {Q'}^{(k)}\underbrace{R^{(k)}{R'}^{(k-1)}}_{={R'}^{(k)}}
	 = {Q'}^{(k)}{R'}^{(k)},

   as required (using the definition of `Z` and then the definition of
   `{R'}^{(k)}`).

   For the QR algorithm, we again use the inductive hypothesis on the
   QR factorization of `A^k` followed by the inductive hypothesis on
   the similarity transform to get

      .. math::

	 A^k = AA^{k-1} =A{Q'}^{(k-1)}{R'}^{(k-1)}
	 {Q'}^{(k-1)}A^{(k-1)}{R'}^{(k-1)} =
	 {Q'}^{(k-1)}Q^{(k)}R^{(k)}{R'}^{(k-1)}
	 = {Q'}^{(k)}{R'}^{(k)},

   where we used the algorithm definitions in the third equality and
   then the definitions of `{Q'}^{(k)}` and `{R'}^{(k)}`. To verify
   the similarity transform at iteration `k` we use the algorithm definitions
   to write

      .. math::

	 A^{(k)} = R^{(k)}Q^{(k)} = (Q^{(k)})^TQ^{(k)}R^{(k)}Q^{(k)}
	 = ({Q'}^{(k)})^TA({Q'})^{(k)},

   as required.

This theorem tells us that the QR algorithm will converge under the
conditions that simultaneous iteration converges. It also tells us
that the QR algorithm finds an orthonormal basis (the columns of
`{Q'}^{(k)}`) from the columns of each power of `A^k`; this is how
it relates to power iteration.

The practical QR algorithm
--------------------------

.. hint::
   
   A video recording for this material is available `here
   <https://player.vimeo.com/video/454125822>`_.

The practical QR algorithm for real symmetric matrices has a number of
extra elements that make it fast. First, recall that we start by
transforming to tridiagonal (symmetric Hessenberg) form. This cuts
down the numerical cost of the steps of the QR algorithm. Second, the
Rayleigh quotient algorithm idea is incorporated by applying shifts
`A^{(k)}-\mu^{(k)}I`, where `\mu^{(k)}` is some eigenvalue
estimate. Third, when an eigenvalue is found (i.e. an eigenvalue
appears accurately on the diagonal of `A^{(k)}`) the off-diagonal
components are very small, and the matrix decouples into a block
diagonal matrix where the QR algorithm can be independently applied to
the blocks (which is cheaper than doing them all together). This final
idea is called deflation.

A sketch of the practical QR algorithm is as follows.

* `A^{(0)} \gets` TRIDIAGONAL MATRIX
* FOR `k=1,2,\ldots`

  * PICK A SHIFT `\mu^{(k)}` (discussed below)
  * `Q^{(k)}R^{(k)} = A^{(k-1)} - \mu^{(k)}I` (from QR factorisation)
  * `A^{(k)} = R^{(k)}Q^{(k)} + \mu^{(k)}I`
  * IF `A^{(k)}_{j,j+1}\approx 0` FOR SOME `j`

    * `A_{j,j+1}\gets 0`
    * `A_{j+1,j}\gets 0`

    * continue by applying the practical QR algorithm to
      the diagonal blocks `A_1` and `A_2` of
      
      .. math::

	 A_k =
	 \begin{pmatrix}
	 A_1 & 0 \\
	 0 & A_2 \\
	 \end{pmatrix}

One possible way to select the shift `\mu^{(k)}` is to calculate a
Rayleigh quotient with `A` using the last column `q_m^{(k)}` of `{Q'}^{(k)}`,
which then gives cubic convergence for this eigenvector and
eigenvalue. In fact, this is just `A_{mm}^k`,

   .. math::

      A_{mm}^{(k)} = e_m^TA^{(k)}e_m = e_m^T({Q'}^{(k)})^TA{Q'}^{(k)}e_m
      = (q_m^{(k)})^TAq_{m} = \mu^{(k)}.

This is very cheap, we just read off the bottom right-hand corner
from `A^{(k)}`! This is called the Rayleigh quotient shift.

It turns out that the Rayleigh quotient shift is not guaranteed to
work in all cases, so there is an alernative approach called the
Wilkinson shift, but we won't discuss that here.
