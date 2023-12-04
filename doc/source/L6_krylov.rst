.. default-role:: math

Iterative Krylov methods for `Ax=b`
===================================

.. details:: Supplementary video

   .. vimeo:: 454126320

In the previous section we saw how iterative methods are necessary
(but can also be fast) for eigenvalue problems `Ax=\lambda x`.
Iterative methods can also be useful for solving linear systems
`Ax=b`, generating a sequence of vectors `x^k` that converge to the
solution. We shall examine Krylov subspace methods, where each
iteration mainly involves a matrix-vector multiplication. For dense
matrices, matrix-vector multiplication costs `\mathcal{O}(m^2)`, but
often (e.g. numerical solution of PDEs, graph problems, etc.)  `A` is
sparse (i.e. has zeros almost everywhere) and the matrix-vector
multiplication costs `\mathcal{O}(m)`. The goal is then to find a
method where only a few iterations are necessary before the error is
very small, so that the solver has total cost `\mathcal{O}(mN)` where
`N` is the total number of iterations, hopefully small.

Since we only need the result of a matrix-vector multiplication, it is
even possible to solve a linear system without storing `A`
explicitly. Instead one can just provide a subroutine that implements
matrix-vector multiplication in some way; this is called a
"matrix-free" implementation.

Krylov subspace methods
-----------------------

.. details:: Supplementary video

   .. vimeo:: 454126582

In this section we will introduce Krylov subspace methods for solving
`Ax=b` (we will not specialise to real or symmetric matrices
here). The idea is to approximate the solution using the basis

.. math::

   (b, Ab, A^2b, A^3b, \ldots, A^kb)

whose span is called a Krylov subspace. After each iteration the
Krylov subspace grows by one dimension. As we have already seen from
studying power iteration, the later elements in this sequence will get
very parallel (they will all be approximating the eigenvector with
largest magnitude of eigenvalue). Hence, we once again need to resort
to orthogonalising the basis. We could simply take the QR
factorisation of this basis, but the Arnoldi iteration coming up
next also provides a neat way to solve the equation when projected
onto the Krylov subspace.

Arnoldi iteration
-----------------

The key to Krylov subspace methods turns out to be the transformation
of `A` to an upper Hessenberg matrix by orthogonal similarity
transforms, so that `A=QHQ^*`. We have already looked at using
Householder transformations to do this in the previous section.
The Householder technique is not so suitable for large dimensional
problems, so instead we look at a way of proceeding column by
column, just like the Gram-Schmidt method does for finding
`QR` factorisations.

We do this by rewriting `AQ=QH`. The idea is that at iteration `n` we
only look at the first `n` columns of `Q`, which we
call `\hat{Q}_n`. When `m` is large, this is a significant saving:
`mn \ll m^2`.
To execute the iteration, it turns out that
we should look at the `(n+1)\times n` upper left-hand section of `H`,
i.e.

   .. math::

      \tilde{H}_n = \begin{pmatrix}
      h_{11} & \ldots & h_{1n} \\
      h_{21} & \ddots & \vdots \\
      0 & \ddots & h_{nn} \\
      0 & 0 & h_{n+1,n} \\
      \end{pmatrix}

Then, `A\hat{Q}_n = \hat{Q}_{n+1}\tilde{H}_n`.

.. details:: Supplementary video

   .. vimeo:: 454127181

Using the column space interpretation of matrix-matrix multiplication,
we see that the `n`-th column is

   .. math::

      Aq_n = h_{1n}q_1 + h_{2n}q_2 + \ldots + h_{n,n}q_n + h_{n+1,n}q_{n+1}.

This formula shows us how to construct the non-zero entries of the
nth column of `H`; this defines the Arnoldi algorithm which we
provide as pseudo-code below.

* `q_1\gets b/\|b\|`
* FOR `n=1,2,\ldots`

  * `v\gets Aq_n`
  * FOR `j=1` TO `n`

    * `h_{jn}\gets q_j^*v`
    * `v \gets v - h_{jn}q_j`
  * END FOR
  * `h_{n+1,n} \gets \|v\|`
  * `q_{n+1} \gets v/\|v\|`
* END FOR

.. proof:exercise::

   `(\ddagger)` The :func:`cla_utils.exercises10.arnoldi` function has
   been left unimplemented. It should implement the Arnoldi algorithm
   using Numpy array operations where possible, and return the `Q` and
   `H` matrices after the requested number of iterations is complete.
   What is the minimal number of Python for loops possible?

   The test
   script ``test_exercises10.py`` in the ``test`` directory will test
   this function.
  
If we were to form the QR factorisation of the `m\times n` Krylov
matrix

   .. math::

      K_n = \begin{pmatrix}
      b & Ab & \ldots & A^{n_1}b \\
      \end{pmatrix}

then we would get `Q=Q_n`. Importantly, in the Arnoldi iteration, we
never form `K_n` or `R_n` explicitly, since these are very
ill-conditioned and not useful numerically.

.. details:: Supplementary video

   .. vimeo:: 454136990
   
But what is the use of the `\tilde{H}_n` matrix? Applying
`\hat{Q}_n^*` to `A\hat{Q}_n = \hat{Q}_{n+1}\tilde{H}_n` gives

   .. math::

      \hat{Q}_n^*A\hat{Q}_n = \hat{Q}_n^*\hat{Q}_{n+1}\tilde{H}_n,

      = \begin{pmatrix}
      1 & 0 & \ldots & 0 & 0 \\
      0 & \ddots & \ddots & \vdots & \vdots \\
      \vdots & \ddots & \ddots & \vdots & \vdots \\
      0 & \ldots & \ldots & 1 & 0 \\ 
      \end{pmatrix}
      \tilde{H}_n
      = H_n,

where `H_n` is the `n\times n` top left-hand corner of `H`.

.. details:: Supplementary video

   .. vimeo:: 454171516

The intrepretation of this is that `H_n` is the orthogonal projection
of `A` onto the Krylov subspace `\mathrm{span}(K_n)`. To see this, take any vector `v`,
and project `Av` onto the the Krylov subspace `\mathrm{span}(K_n)`.

   .. math::

      PAv = \hat{Q}_n\hat{Q}_n^*v.

Then, changing basis to the orthogonal basis gives

   .. math::

      \hat{Q}_n^*(\hat{Q}_n\hat{Q}_n^*A)\hat{Q}_n = \hat{Q}_n^*A\hat{Q}_n
      = H_n.

GMRES
-----

.. details:: Supplementary video

   .. vimeo:: 454171559

The Generalised Minimum Residual method (GMRES), due to Saad (1986),
exploits these properties of the Arnoldi iteration. The idea is
that we build up the orthogonal basis for the Krylov subspaces
one by one, and at each iteration we solve the projection of
`Ax=b` onto the Krylov basis as a least squares problem, until
the residual `\|Ax-b\|` is below some desired tolerance.

To avoid the numerical instabilities that would come from using the
basis `(b,Ab,A^2b,\ldots)`, we use the Arnoldi iteration to build an
orthonormal basis, and seek approximate solutions of the form `x_n =
\hat{Q}_ny` for `y\in\mathbb{C}^n`. We then seek the value of `y`
that minimises the residual

   .. math::

      \mathcal{R}_n = \|AQ_ny - b\|.

This explains the Minimum Residual part of the name. We also see from
this definition that the residual cannot increase with iterations,
because it only increases the subspace where we seek a solution.

This problem can be simplified further by using `AQ_n = Q_{n+1}\tilde{H}_n`,
so

   .. math::

      \mathcal{R}_n = \|\hat{Q}_{n+1}\tilde{H}_n y - b\|.

Remembering that `b=\|b\|q_1`, we see that the entire residual is in
the column space of `\hat{Q}_{n+1}`, and hence we can invert
on the column space using `\hat{Q}_{n+1}^*` which does not change
the norm of the residual due to the orthonormality.

   .. math::

      \mathcal{R}_n = \|\tilde{H}_n y - \hat{Q}_{n+1}^*b\|
      = \|\tilde{H}_n y - e_1\|b\|\|.

Finding `y` to minimise `\mathcal{R}_n` requires the solution of a
least squares problem, which can be computed via QR factorisation
as we saw much earlier in the course.

.. details:: Supplementary video

   .. vimeo:: 454171921

We are now in position to present the GMRES algorithm as pseudo-code.

* `q_1 \gets b/\|b\|`
* FOR `n=1,2,\dots`

  * APPLY STEP `n` OF ARNOLDI
  * FIND `y` TO MINIMIZE `\|\tilde{H}_ny - \|b\|e_1\|`
  * `x_n \gets \hat{Q}_ny`
  * CHECK IF `\mathcal{R}_n <` TOL
* END FOR

.. proof:exercise::

   `(\ddagger)` The :func:`cla_utils.exercises10.GMRES` function has
   been left unimplemented. It should implement the basic GMRES
   algorithm above, using one loop over the iteration count. You can
   paste code from your :func:`cla_utils.exercises10.arnoldi`
   implementation, and you should use your least squares code to solve
   the least squares problem.  The test script ``test_exercises10.py``
   in the ``test`` directory will test this function.

.. proof:exercise::

   The least squares problem in GMRES requires the QR
   factorisation of `H_k`. It is wasteful to rebuild this from scratch
   given that we just computed the QR factorisation of
   `H_{k-1}`. Modify your code so that it recycles the QR
   factorisation, applying just one extra Householder rotation per
   GMRES iteration. Don't forget to check that it still passes the
   test.

.. hint::

   Don't get confused by the two Q matrices involved in GMRES! There
   is the Q matrix for the Arnoldi iteration, and the Q matrix for
   the least squares problem. These are not the same.

Convergence of GMRES
--------------------

.. details:: Supplementary video

   .. vimeo:: 454198706

The algorithm presented as pseudocode is the way to implement GMRES
efficiently. However, we can make an alternative formulation
of GMRES using matrix polynomials.

We know that `x_n\in \mathrm{span}(K_n)`, so we can write

   .. math::

      x_n = c_0b + c_1Ab + c_2A^2b + \ldots + c_{n-1}A^{n-1}b
      = p'(A)b,

where

   .. math::

      p'(z) = c_0 + c_1z + c_2z^2 + \ldots + \ldots c_{n-1}z^{n-1} \implies
      p'(A) = c_0I + c_1A + c_2A^2 + \ldots + c_{n-1}A^{n-1}.

Here we have introduced the idea of a matrix polynomial, where the
kth power of `z` is replaced by the kth power of `A`.

The residual becomes

   .. math::

      r_n = b - Ax_n = b - Ap'(A)b = (I - Ap'(A))b
      = p(A)b,

where `p(z) = 1 - zp'(z)`. Thus, the residual is a matrix polynomial
`p` of `A` applied to `b`, where `p\in \mathcal{P}_n`, and

   .. math::

      \mathcal{P}_n = \{\mbox{degree }\leq n\mbox{ polynomials with }
      p(0)=1\}.

Hence, we can recast iteration `n` of GMRES as a polynomial
optimisation problem: find `p_n\in \mathcal{P}_n` such that
`\|p_n(A)b\|` is minimised.  We have

   .. math::

      \|r_n\|  = \|p_n(A)b\| \leq \|p_n(A)\|\|b\|
      \implies \frac{\|r_n\|}{\|b\|} \leq \|p_n(A)\|.

Assuming that `A` is diagonalisable, `A=V\Lambda V^{-1}`, then
`A^s=V\Lambda^sV^{-1}`, and

   .. math::

      \|p_n(A)\| = \|Vp_n(\Lambda^s)V^{-1}\|
      \leq \underbrace{\|V\|\|V^{-1}\|}_{=\kappa(V)}
      \|p_n\|_{\Lambda(A)},

where `\Lambda(A)` is the set of eigenvalues of `A`, and

   .. math::

      \|p\|_{\Lambda(A)} = \sup_{x\in \Lambda(A)}\|p(x)\|.

Hence we see that GMRES will converge quickly if `V` is
well-conditioned, and `p(x)` is small for all `x\in \Lambda(A)`.  This
latter condition is not trivial due to the `p(0)=1` requirement. One
way it can happen is if `A` has all eigenvalues clustered in a small
number of groups, away from $0$. Then we can find a low degree
polynomial that passes through 1 at `x=0`, and 0 near each of the
clusters. Then GMRES will essentially converge in a small number of
iterations (equal to the degree of the polynomial). There are problems
if the eigenvalues are scattered over a wide region of the complex
plane: we need a very high degree polynomial to make `p(x)` small at
all the eigenvalues and hence we need a very large number of
iterations. Similarly there are problems if eigenvalues are very close
to zero.

.. proof:exercise::

   Investigate the convergence of the matrices provided by the
   functions :func:`cla_utils.exercises10.get_AA100`,
   :func:`cla_utils.exercises10.get_BB100`, and
   :func:`cla_utils.exercises10.get_CC100`, by looking at the residual
   after each iteration (which should be provided by
   :func:`cla_utils.exercises10.GMRES` as specified in the docstring).
   What do you observe? What is it about the three matrices that
   causes this different behaviour?

Preconditioned GMRES
--------------------

.. details:: Supplementary video

   .. vimeo:: 454218547

This final topic has been a strong focus of computational linear algebra
over the last 30 years. Typically, the matrices that we want to solve
do not have eigenvalues clustered in a small number of groups, and so
GMRES is slow. The solution (and the challenge) is to find a matrix
`\hat{A}` such that `\hat{A}x = y` is cheap to solve (diagonal, or triangular, or
something else) and such that `\hat{A}^{-1}A` *does* have eigenvalues clustered
in a small number of groups (e.g. `\hat{A}` is a good approximation of `A`, so
that `\hat{A}^{-1}A\approx I` which has eigenvalues all equal to 1). We call
`\hat{A}` the preconditioning matrix, and the idea is to apply GMRES to
the (left) preconditioned system

   .. math::

      \hat{A}^{-1}Ax = \hat{A}^{-1}b.

GMRES on this preconditioned system is equivalent to the following algorithm,
called preconditioned GMRES.

* SOLVE `\hat{A}\tilde{b}=b`.
* `q_1 \gets \tilde{b}/\|\tilde{b}\|`
* FOR `n=1,2,\dots`

  * SOLVE `\hat{A}v = Aq_n`
  * FOR `j=1` TO `n`

    * `h_{jn}=q_j^*v`
    * `v = v - h_{jn}q_j`
  * END FOR
  * `h_{n+1,n} \gets \|v\|`
  * `q_{n+1}`\gets v/\|v\|`
  * FIND `y` TO MINIMIZE `\|\tilde{H}_ny - \|\tilde{b}\|e_1\|`
  * `x_n \gets \hat{Q}_ny`
  * CHECK IF `\mathcal{R}_n <` TOL
* END FOR

.. proof:exercise::

   Show that this algorithm is equivalent to GMRES applied to the
   preconditioned system.
  
The art and science of finding preconditioning matrices `\hat{A}` (or
matrix-free procedures for solving `\hat{A}x=y`) for specific problems
arising in data science, engineering, physics, biology etc. can
involve ideas from linear algebra, functional analysis, asymptotics,
physics, etc., and represents a major activity in scientific computing
today.

Knowing when to stop
--------------------

We should stop an iterative method when the error is sufficiently small.
But, we don't have access the the exact solution, so we can't compute the
error. Two things we can look at are:

* The residual `{r}^k=A{x}^k-{b}`, or
* The pseudo-residual `{s}^k = {x}^{k+1}-{x}^k`,
  which both tend to zero as `{x}^k\to{x}^*` provided that `A` is
  invertible.
  
How do their sizes relate to the size of `{e}^k={x}^k-{x}^*`?

.. math::
   {e}^k  = & {x}^* - {x}^k \\
    = & A^{-1}(A{x}^*-A{x}^k) \\
    = & A^{-1}({b}-A{x}^k) \\
    = & A^{-1}{r}^k,

so `\|e^k\| \leq \|A^{-1}\|\|r^k\|`.

The relative error `\|e^k\|/\|x^*\|` satisfies

.. math::

   \frac{\|e^k\|}{\|x^*\|}
   = \frac{\|A^{-1}r^k\|}{\|x\|}
   \leq \|A^{-1}\|\frac{\|r^k\|}{\|x\|}
   \leq \|A^{-1}\|\|A\|\frac{\|r^k\|}{\|b\|},

so the relative error is bounded from above by
the condition number of `\|A\|` multiplied by
the relative residual `\|r^k\|/\|b\|`. If the condition
number is large, it is possible to have a small residual
but still have a large condition number.

Similar results hold for the pseudoresidual, but are
more complicated to show for the case of GMRES.
