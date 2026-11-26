.. default-role:: math

==================
Errata for 2026/27
==================

Here we'll list all the updates to the exercises and coursework Classroom
repositories for Autumn 2026. An up to date version is in the
`public repository <https://github.com/comp-lin-alg/comp-lin-alg-course>`_.

Chapter 3: QR factorisation
---------------------------

* In :doc:`L2_QR_factorisation`, the slice `A[0:3, 3:6]` is a
  `3\times 3` submatrix, rather than `4\times 4`.

Chapter 4: Analysing algorithms
-------------------------------

* In :doc:`L3_analysing_algorithms`, in the operation count for
  modified Gram-Schmidt, the inner product is `r_{ij}\gets q_i^*v_j`,
  rather than `q_i^*v_i`.
* In :doc:`L3_analysing_algorithms`, modified Gram-Schmidt has `n-i`
  inner iterations over `j` for each `i`, rather than `n-i-1`.
* In :doc:`L3_analysing_algorithms`, the Frobenius norm is
  `\|A\|_F=\sqrt{\sum_{i,j}|A_{ij}|^2}`; the absolute value was
  missing, which matters for complex matrices.
* In :doc:`L3_analysing_algorithms`, the absolute and relative
  condition numbers are defined with a limit over small
  perturbations, `\lim_{\delta\to 0}\sup_{0<\|\delta x\|\leq\delta}`,
  rather than a supremum over all `\delta x\neq 0`, so that they
  measure the conditioning of `f` at `x`.
* In :doc:`L3_analysing_algorithms`, in the example on the roots of
  `x^2-2x+c`, the perturbed coefficient should be `c=1-\delta c`
  (rather than `1+\delta c`) for the roots to be `1\pm\sqrt{\delta c}`.
* In :doc:`L3_analysing_algorithms`, machine epsilon
  `\varepsilon=2^{-53}` is half the distance between 1 and the next
  floating point number (which is `2^{-52}`), not the distance itself.
* In :doc:`L3_analysing_algorithms`, in the backward stability of
  the Householder algorithm, Householder transformations are
  reflections, not rotations.
* In :doc:`L3_analysing_algorithms`, the exercise on the stability of
  the Householder algorithm has been corrected: `R_1` is a random
  upper triangular matrix, and the exercise now states that
  `A=Q_1R_1` is formed and factorised as `A=Q_2R_2` before comparing.
* In :doc:`L3_analysing_algorithms`, each step of back substitution
  requires `m-i` multiplications and subtractions, rather than
  `m-i-1`.
* In :doc:`L3_analysing_algorithms`, backward stability of computing
  `\tilde{Q}^*b` means that `(\tilde{Q}+\delta Q)\tilde{y}=b` for some
  small `\delta Q`; this does not follow from
  `\tilde{y}=(\tilde{Q}+\delta Q)^*b`, as was stated.
* In :doc:`L3_analysing_algorithms`, in the proof that solving
  `Ax=b` via QR is backward stable, the result should read
  `b=(A+\Delta A)\tilde{x}`; a bracket was misplaced.

Chapter 5: Finding eigenvalues of matrices
------------------------------------------

* In :doc:`L5_eigenvalues`, after the Schur factorisation theorem,
  the eigenvalues are read from the diagonal of `T`, not from the
  diagonal of `A`.
* In :doc:`L5_eigenvalues`, the stopping criterion for the iterative
  transformation to Schur form should measure the strictly lower
  triangular part of `A_k`, rather than the lower triangular part
  of `A`.
* In :doc:`L5_eigenvalues`, the argument for why the transformation
  to Schur form must be iterative should conclude that an explicit
  formula for the roots of the characteristic polynomial is
  impossible in general (not "explicit in general").
* In :doc:`L5_eigenvalues`, the reduction to upper Hessenberg form
  takes `m-2` steps, rather than `m-1`.
* In :doc:`L5_eigenvalues`, in the explanation of why a single
  similarity transformation cannot triangularise `A`, the product
  should be `Q_1^*(Q_1^*A)^*` rather than `Q_1^*(Q_1^*A)`, and it
  forms linear combinations of rows, not columns.
* In :doc:`L5_eigenvalues`, in the reduction to Hessenberg form,
  Householder transformations are reflections, not rotations.
* In :doc:`L5_eigenvalues`, in the Hessenberg pseudo-code, the
  normalisation should be `v_k\gets v_k/\|v_k\|_2`, rather than
  `v_k/\|v\|_2`.
* In :doc:`L5_eigenvalues`, in the operation count for the reduction
  to Hessenberg form, the right multiplication costs `\sim 2m^3`
  FLOPs; `\frac{10}{3}m^3` is the total for the left and right
  multiplications together.
* In :doc:`L5_eigenvalues`, in the proof of convergence of power
  iteration, the eigenvector estimate follows from the eigenvector
  expansion of `A^kv^0`, not from the Taylor series; the eigenvalue
  estimate then follows from the quadratic accuracy of the Rayleigh
  quotient.
* In :doc:`L5_eigenvalues`, inverse iteration converges to the
  eigenvector for which `1/|\lambda-\mu|` is largest (not smallest),
  i.e. whose eigenvalue is closest to `\mu`.
* In :doc:`L5_eigenvalues`, in the convergence of Rayleigh quotient
  iteration, the eigenvalue error is `|\lambda^{(k)}-\lambda_J|`,
  rather than `|\lambda^{(k)}-q_J|`.
* In :doc:`L5_eigenvalues`, the exercise on
  :func:`cla_utils.exercises9.rq_it` should ask for Rayleigh quotient
  iteration, rather than inverse iteration.
* In :doc:`L5_eigenvalues`, in the hint on complex QR factorisation,
  `\mbox{sign}(x_1)` should be replaced by `e^{i\arg(x_1)}`, rather
  than `-e^{i\arg(x_1)}`; the minus sign reintroduces the
  cancellation that the choice of sign is meant to avoid.
* In :doc:`L5_eigenvalues`, the convergence condition for
  simultaneous iteration with `n<m` vectors also requires
  `|\lambda_n|>|\lambda_{n+1}|`.
* In :doc:`L5_eigenvalues`, in the simultaneous iteration
  pseudo-code, `\hat{Q}^{(0)}` is an `m\times n` matrix with
  orthonormal columns (unitary only when `n=m`), and the QR
  factorisation computes `\hat{Q}^{(k)},\hat{R}^{(k)}`.
* In :doc:`L5_eigenvalues`, in the connection between the QR
  algorithm and inverse iteration, the indices should be `m` rather
  than `n` (`e_m`, `\tilde{r}^{(k)}_{mm}`), and the last column of
  `{Q'}^{(k)}` is the result of `k` steps of inverse iteration with
  `\mu=0` applied to `A^T`, starting from `e_m`.
* In :doc:`L5_eigenvalues`, in the description of shifts in the
  practical QR algorithm, it is the last diagonal entry of `A^{(k)}`
  (not of `R^{(k)}`) that converges, via inverse iteration, to the
  eigenvalue of `A` of smallest magnitude.
* In :doc:`L5_eigenvalues`, the Rayleigh quotient shift computed
  from `A^{(k)}` is `A^{(k)}_{mm}=(q_m^{(k)})^TAq_m^{(k)}`, and it is
  used as the shift `\mu^{(k+1)}` for the next iteration (rather than
  `\mu^{(k)}`).

Chapter 6: Krylov methods
-------------------------

* In :doc:`L6_krylov`, the last column of the Krylov matrix `K_n`
  should be `A^{n-1}b`, rather than `A^{n_1}b`.
* In :doc:`L6_krylov`, the orthogonal projection of `Av` onto the
  Krylov subspace is `PAv=\hat{Q}_n\hat{Q}_n^*Av`; the `A` was
  missing on the right-hand side.
* In :doc:`L6_krylov`, in the exercise on recycling the QR
  factorisation in GMRES, the least squares matrix is `\tilde{H}_k`
  rather than `H_k` (and `\tilde{H}_{k-1}` rather than `H_{k-1}`).
  Householder transformations are reflections, not rotations.
* In :doc:`L6_krylov`, in the GMRES convergence bound,
  `p_n(A)=Vp_n(\Lambda)V^{-1}`, rather than `Vp_n(\Lambda^s)V^{-1}`.
* In :doc:`L6_krylov`, in "Knowing when to stop", invertibility of
  `A` is not needed for the residual and pseudo-residual to tend to
  zero as `x^k\to x^*`; it is needed for the converse, that a small
  residual implies a small error.
* In :doc:`L6_krylov`, in the bound on the relative error, the
  denominators should be `\|x^*\|` rather than `\|x\|`.
* In :doc:`L6_krylov`, in "Knowing when to stop", the relative
  error is bounded by the condition number of `A` (not of `\|A\|`)
  times the relative residual; if the condition number is large, a
  small residual can coexist with a large relative error (not a
  large condition number).

Chapter 7: Preconditioning Krylov methods
-----------------------------------------

* In :doc:`L7_preconditioning`, in the proof that a stationary method
  does not converge when `\rho(C)\geq 1`, the eigenvalue condition
  `|\lambda|>1` should be `|\lambda|\geq 1`, to include the case
  `\rho(C)=1`.
* In :doc:`L7_preconditioning`, in the proof that `\rho(C)<1` implies
  convergence, `\lambda^k` in the eigenvector expansion should be
  `\lambda_i^k`.
* In :doc:`L7_preconditioning`, in the proof that the eigenvalues of
  the symmetrically preconditioned matrix lie in `[1-c,1+c]`, the
  square root `S` should satisfy `SS=M_s^{-1}`, rather than `SS=M_s`.
* In :doc:`L7_preconditioning`, the SOR convergence condition
  `0<\omega 2` should be `0<\omega<2`.
* In :doc:`L7_preconditioning`, in the finite difference example,
  `f_k=f(k/n)/n^2` should be `f_k=f(k/(n+1))/(n+1)^2`, consistent
  with the grid spacing `h=1/(n+1)`.
* In :doc:`L7_preconditioning`, the example "Jacobi iteration for the
  example matrix" contained swapped spatial and eigenvector indices
  and used sine instead of cosine when combining exponentials.
  The corrected eigenvalues of `D^{-1}A` are
  `\lambda_l=1-\cos(l\pi/(n+1))`, for `l=1,\ldots,n`.
  Hence the maximum is `1+\cos(\pi/(n+1))`, attained at `l=n`,
  and the upper bound for the scaled Jacobi parameter is
  `\omega<2/(1+\cos(\pi/(n+1)))`.
* In :doc:`L7_preconditioning`, in the introduction to Chebyshev
  acceleration, the coefficients `c_{jk}` should be sought for
  `j=0,\ldots,k`, including the coefficient of the initial guess `x^0`.
  The eigenvalues to be made small are those of the matrix `p_k(C)`,
  rather than of the polynomial `p_k` itself.
  After changing variables to the interval `[\alpha,\beta]`, the
  leading coefficient is `2^{2k-1}/(\beta-\alpha)^k`, rather than
  `2^{2k-1}/(\beta-\alpha)`.
  The accelerated iterate is `y^k=\sum_{j=0}^kc_{jk}x^j`, rather
  than a sum with `x^k` in every term.
* In :doc:`L7_preconditioning`, in the commented-out "Accelerated
  Jacobi's method" block, the iteration matrix is
  `-M^{-1}N=I-D^{-1}A`; the minus sign was missing.
  Also, if `\lambda` is an eigenvalue of `I-D^{-1}A`, then
  `1-\lambda` is an eigenvalue of `D^{-1}A`, rather than `D^{-1}`.
* In :doc:`L7_preconditioning`, the symmetric convergence analysis
  should use `\lambda_{\max}=\|A\|_2`, rather than `\|A\|_2^2`.
  The positivity argument for `\hat{\mu}_{\min}` uses the square of
  the smallest singular value of `M^{-1}A^{1/2}`, rather than its
  minimum eigenvalue.
* In :doc:`L7_preconditioning`, the CG preconditioner must be symmetric
  positive definite. The theorem placing the preconditioned eigenvalues
  in `[1-c,1+c]` also assumes that `A` is symmetric.
* In :doc:`L7_preconditioning`, the eigenvalue clustering proof uses
  the relation `Cv=(1-\lambda)v`, without assuming an eigendecomposition.
  Eigenvalue clustering alone does not guarantee fast GMRES convergence;
  for diagonalizable matrices, eigenvector conditioning also enters
  the convergence bound. The simplified proof of stationary convergence
  when `\rho(C)<1` covers the diagonalizable case.
