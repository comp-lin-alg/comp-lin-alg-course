.. default-role:: math

Preconditioning Krylov methods
==============================

In this section we will discuss some preconditioners and how to
analyse them. The most important question is how quickly the
preconditioned GMRES algorithm converges for a given matrix and
preconditioner. We will focus on the link between stationary
iterative methods and preconditioners.

Stationary iterative methods
---------------------------

As we have already discussed, given a matrix equation `Ax=b`,
iterative methods provide a way of obtaining a (hopefully) better
approximate solution `{x}^{k+1}` from a previous approximate
`{x}^k`. Stationary iterative methods are defined from splittings
as follows.

.. proof:definition:: Stationary iterative methods

   A stationary iterative method is constructed from matrices `M` and
   `N`with `A=M+N`. Then the iterative method is defined by

   .. math::

      M{x}^{k+1}=-N{x}^k+{b}.

The word "stationary" refers to the fact that exactly the same thing
is done at each iteration. This contrasts with Krylov methods such as
GMRES, where the sequence of operations depends on the previous
iterations (e.g. a different size least square system is solve in each
GMRES iteration).

In this section we will introduce/recall some classic stationary
methods.

..
   \frame{
   \frametitle{2D Poisson equation}
   Consider the following problem:
   \begin{eqnarray*}
     -\nabla^2 p(x,y) &=& f(x,y), \quad 0\leq x,y \leq 1, \\
     p(x,y) & = & 0, \quad x=0\,\mathrm{or}\,x=1\,\mathrm{or}\,y=0\,\mathrm{or}\,y=1.
   \end{eqnarray*}
   Finite difference approximation is
   \begin{eqnarray*}
    4p_{i,j} - p_{i-1,j} - p_{i+1,j} - p_{i,j+1} - p_{i,j-1} &=& \Delta x^2 f_{i,j}, 
    \quad 0<i,j<N, \\
    p_{i,j} =  0, && \quad i=0\,\mathrm{or}\,N\,\mathrm{or}\,j=0\,\mathrm{or}\,N.
    \\
   \end{eqnarray*}
   }

   \frame{
   \begin{eqnarray*}
     4p_{i,j} - p_{i-1,j} - p_{i+1,j} - p_{i,j+1} - p_{i,j-1} &=& \Delta x^2 f_{i,j}, 
     \quad 0<i,j<N, \\
     p_{i,j} =  0, && \quad i=0\,\mathrm{or}\,N\,\mathrm{or}\,j=0\,\mathrm{or}\,N.
     \\
   \end{eqnarray*}.
   If we write 
   \[
   {p}=\left(p_{1,1},\ldots,p_{N-1,1},
   p_{1,2},\ldots,p_{N-1,2},
   \ldots,
   p_{1,N-1},\ldots,p_{N-1,N-1}\right),
   \]
   then equation becomes
   \[
   \begin{pmatrix}
   D+2I & -I & \ldots & 0 \\
   -I & D+2I & \ddots & \vdots \\
   \vdots & \ddots & \ddots & -I \\
   0 & \ldots & -I & D + 2I \\
   \end{pmatrix}
   {p} = 
   {b},
   \,
   D = 
   \begin{pmatrix}
   2 & -1 & \ldots & 0 \\
   -1 & 2 & \ddots & \vdots \\
   \vdots & \ddots & \ddots & -1 \\
   0 & \ldots & -1 & 2 \\
   \end{pmatrix}.
   \]
   }

.. proof:definition:: Richardson iteration

   For a chosen parameter `\omega>0`, take `M=I/\omega`. This
   defines the iterative method given by

   .. math::
   
      {x}^{k+1} = {x}^k+\omega\left({b}-A{x}^k\right).

Richardson, L.F. (1910). *The approximate arithmetical solution
by finite differences of physical problems involving differential
equations, with an application to the stresses in a masonry
dam*. Philos. Trans. Roy. Soc. London Ser. A 210: 307-357.
      
This approach is convenient for parallel computing, because each entry in
`x^{k+1}` can be updated independently, once `Ax^k` has been evaluated.

..
   \frame{
   \frametitle{Richardson's method for the Poisson problem}
   Finite difference approximation is
   \begin{eqnarray*}
   4p_{i,j} - p_{i-1,j} - p_{i+1,j} - p_{i,j+1} - p_{i,j-1} &=& \Delta x^2 f_{i,j}, 
   \quad 0<i,j<N, \\
   p_{i,j} =  0, && \quad i=0\,\mathrm{or}\,N\,\mathrm{or}\,j=0\,\mathrm{or}\,N.
   \\
   \end{eqnarray*}
   Richardson iteration gives:
   \[
   p_{i,j}^{k+1} =  p_{i,j}^k  -\omega\left(4p_{i,j}^k - p_{i-1,j}^k - p_{i+1,j}^k - p_{i,j+1}^k -
   p_{i,j-1}^k + \Delta x^2 f_{i,j}\right).
   \]
   }


.. proof:definition:: Jacobi's method
   
   Split `A=L+D+U` with `L` strictly lower triangular, `D`
   diagonal and `U` strictly upper triangular, i.e.

   .. math::
      
      L_{ij}=0, \, j\geq i, \quad D_{ij}=0, \, i\neq j, \quad U_{ij}, \, i\geq j.

   Then, Jacobi's method is

   .. math::

     D{x}^{k+1} = {b}-(L+U){x}^k.
   
`D` is very cheap to invert because it is diagonal; entries in
`x^{k+1}` can be updated independently once `(L+U)x^k` has been evaluated.

Jacobi, C.G.J. (1845). *Ueber eine neue Aufloesungsart der bei der
Methode der kleinsten Quadrate vorkommenden linearen Gleichungen*,
Astronomische Nachrichten, 22, 297-306.

.. \frame{
   \frametitle{Jacobi's method for the Poisson problem}
   Finite difference approximation is
   \begin{eqnarray*}
   4p_{i,j} - p_{i-1,j} - p_{i+1,j} - p_{i,j+1} - p_{i,j-1} &=& \Delta x^2 f_{i,j}, 
   \quad 0<i,j<N, \\
   p_{i,j} =  0, && \quad i=0\,\mathrm{or}\,N\,\mathrm{or}\,j=0\,\mathrm{or}\,N.
   \\
   \end{eqnarray*}
   Jacobi iteration gives:
   \[
   4p_{i,j}^{k+1} =  p_{i-1,j}^k + p_{i+1,j}^k + p_{i,j+1}^k +
   p_{i,j-1}^k + \Delta x^2 f_{i,j}.
   \]
   }

.. proof:definition:: Gauss-Seidel Method

   Split `A=L+D+U` with `L` strictly lower triangular, `D`
   diagonal and `U` strictly upper triangular.
   The Gauss-Seidel method (forward or backwards) is
   
   .. math::

      (L+D){x}^{k+1} = {b}-U{x}^k, \quad\mbox{or},\quad
      (U+D){x}^{k+1} = {b}-L{x}^k.

Each Gauss-Seidel iteration requires the solution of a triangular
system by forward/backward substitution.

.. proof:exercise::

   Show that forward Gauss-Seidel is a modification Jacobi's method
   but using new values as soon as possible.

.. \frame{
   \frametitle{Gauss-Seidel method for the Poisson problem}
   Finite difference approximation is
   \begin{eqnarray*}
   4p_{i,j} - p_{i-1,j} - p_{i+1,j} - p_{i,j+1} - p_{i,j-1} &=& \Delta x^2 f_{i,j}, 
   \quad 0<i,j<N, \\
   p_{i,j} =  0, && \quad i=0\,\mathrm{or}\,N\,\mathrm{or}\,j=0\,\mathrm{or}\,N.
   \\
   \end{eqnarray*}
   For our ordering, Gauss-Seidel iteration gives:
   \[
   4p_{i,j}^{k+1} =  p_{i-1,j}^{k+1} + p_{i+1,j}^k + p_{i,j+1}^k +
   p_{i,j-1}^{k+1} + \Delta x^2 f_{i,j}.
   \]
   }


.. proof:definition:: Scaled Gauss-Seidel method

   We introduce a scaling/relaxation parameter `\omega>0` and
   take `M=D/\omega+L`, so that

   .. math::
      \left(\frac{1}{\omega}D+L\right){x}^{k+1} 
      = {b}+\left(\left(\frac{1}{\omega}-1\right)D-U\right){x}^k.

For `\omega=1`, we recover Gauss-Seidel. For `1<\omega<2`, we often
obtain faster convergence. This is called Successive Over-Relaxation
(SOR).  The optimal value of `\omega` is known for some problems.
This was state of the art for numerical solution of PDEs in the 50s
and 60s.

* Richardson and Jacobi are *simultaneous displacement
  methods*: updates can be done simultaneously (e.g. on a GPU).
  Changing variables by a permutation does not alter the algorithm.

* Gauss-Seidel and SOR are *successive displacement methods*:
  we can only overwrite the old vector with the new one element by element.
  Successive displacement methods usually converge faster, and changing
  variables by a permutation does alter the algorithm.

..
   \frame{
   \frametitle{Red-black ordering for Poisson}
   For the Poisson equation, we can use an alternative ordering. 
   \begin{enumerate}
   \item {red} points with `i+j` even.
   \item {black} points with `i+j` odd.
   \end{enumerate}
   Order the vector with all the red points first, then all the black
   points.\\
   Since the red points depend on the black points and vice-versa,
   this is completely parallel within the red points and within the black points.
   }

Using splitting methods as preconditioners
------------------------------------------

A (non-symmetric) preconditioner `\hat{A}` can be built from a
splitting method by applying one iteration with initial guess
`{z}^0={0}`. Then

.. math::   
   \hat{A}^{-1}A{x} = {z},

where

.. math::
   M{z} = -N{z}^0 + A{x} = A{x},
   
i.e. `\hat{A}=M`. Later we shall see how to relate convergence
properties of splitting methods to the convergence of
preconditioned CG using `\hat{A}=M`.

Symmetric iterative methods
---------------------------

Consider a symmetric matrix `A=A^T`.xs
If we can build iterative methods from the splitting
`A=M+N`, then we can also build iterative methods from the splitting
`A=A^T=M^T+N^T`. We can then combine them together.

.. proof:definition:: Symmetric iterative method

   Given a splitting `A=M+N`, a symmetric method performs one
   stationary iteration using `M+N`, followed by one stationary
   iteration using `M^T+N^T`, i.e.

   .. math::

      M{x}^{k+\frac{1}{2}}=-N{x}^k + {b}, \quad 
      M^T{x}^{k+1}=-N^T{x}^{k+\frac{1}{2}} + {b}.
   
.. proof:example:: Symmetric Successive Over-Relaxation (SSOR).

   For a symmetric matrix `A=L+D+U`, `L=U^T`. The symmetric version
   of SOR is then
   
   .. math::
      
      (L+\frac{1}{\omega}D){x}^{k+\frac{1}{2}}&=\left(\left(\frac{1}{\omega}-1\right)
      D-U\right){x}^k + {b}, \\
      (U+\frac{1}{\omega}D){x}^{k+1}&=\left(\left(\frac{1}{\omega}-1\right)
      D-L\right){x}^{k+\frac{1}{2}} + {b}.

Some Krylov methods, notably the Conjugate Gradient method, require
the preconditioner `\hat{A}` to be symmetric.
We can build symmetric preconditioners from symmetric splitting methods.
Write the symmetric iteration as a single step with
`{x}^0={0}`.

.. math::

   M^T{x}^{1}&=(M-A)^T{x}^{\frac{1}{2}} + {b}, \\
   &= (M-A)^TM^{-1}{b} + {b}, \\
   & = (M^T + M-A)M^{-1}{b},

so that

.. math::

   x^1 = M^{-T}(M^T + M - A)M^{-1}b,
   
i.e. `\hat{A}=M^{-T}(M^T + M-A)M^{-1}`.

.. proof:example:: Symmetric Gauss-Seidel preconditioner

   `\hat{A} = (L+D)^{-T}D(L+D)^{-1}`.

.. A = L+D+L^T, M=L+D, (L+D)^T + (L+D) - L+D+L^T = D.
   
Convergence criteria for stationary methods
-------------------------------------------

In this section we will look at the convergence of stationary
methods. This is relevant because it relates directly to the
convergence properties of the corresponding preconditioned Krylov
method when the stationary method is used as a preconditioner.

For a splitting `A=M+N`, recall that the iterative method is

.. math::

   M{x}^{k+1} = -N{x}^k + {b}.

On the other hand, the solution `{x}^*` of `A{x}={b}` satisfies

.. math::

   M{x}^* = -N{x}^* + {b}.

Subtracting these two equations gives

.. math::

   M{e}^{k+1} = -N{e}^k, \quad {e}^k = {x}^*-{x}^k,

so

.. math::

   {e}^{k+1}=C{e}^k \implies {e}^k = C^k{e}^0, \quad
   C:=-M^{-1}N = -M^{-1}(A-M) = I - M^{-1}A.

`C` is called the *iteration matrix*.

For a symmetric iterative method,

.. math::
   M{x}^{k+\frac{1}{2}}=-N{x}^k + {b}, \quad
   M^T{x}^{k+1}=-N^T{x}^{k+\frac{1}{2}} + {b},

we subtract `Ax^*=b` from both equations to get

.. math::
   M{e}^{k+\frac{1}{2}}=-N{e}^k, \quad
   M^T{e}^{k+1}=-N^T{e}^{k+\frac{1}{2}}.

Then eliminating `e^{k+1/2}` gives

.. math::

   M^Te^{k+1} = N^TM^{-1}Ne^k,

i.e. the iteration matrix is 
   
.. math::

   C = M^{-T}N^TM^{-1}N

.. proof:exercise::

   Show that
   
.. math::

   C = I-\left(M_s\right)^{-1}A,

where

.. math::

   M_s = M(M+M^T-A)^{-1}M^T.

From the above exercise, note the relationship
between `M_s` and `\hat{A}` for symmetric methods.

.. proof:definition:: Convergence of stationary methods

   An iterative method based on the splitting `A=M+N` with iteration
   matrix `C=-M^{-1}N` is called {convergent} if

   .. math::

      {y}^k = C^k{y}^0 \to {0}

   for any initial vector `{y}^0`.

.. proof:exercise::

   Show that this implies that `{e}^k={x}^*-{x}^k\to{0}` i.e.
   `{x}^k\to 0` as `k\to\infty`.

.. proof:theorem:: A first convergence criterion

   If `\|C\|_p<1` for `p=1`, 2 or `\infty`, then the iterative method
     converges.

.. proof:proof::

   .. math::

      \|{y}^k\|_p  = & \|C^k{y}^0\|_p \\
      \leq & \|C^k\|_p\|{y}^0\|_p \\
       \leq & \ \left(\|C\|_p\right)^k\|{y}^0\|_p
      \to  0\quad\mbox{as}\,k\to\infty.

This is only a sufficient condition. There may be matrices `C` with
`\|C\|_p>1` for some `p=1,2,\infty`, but the method is still
convergent.

To obtain a necessary condition, we need to use the spectral radius.

.. proof:definition::

   The spectral radius `\rho(C)` of a matrix `C` is
   the maximum of the absolute values of all the eigenvalues `\lambda_i`
   of `C`:

   .. math::

      \rho(C) = \max_{1\leq i\leq n}|\lambda_i|.

.. proof:theorem::

   An iterative method converges `\iff \rho(C)<1`.

.. proof:proof::

   [Proof that `\rho(C)\geq 1\implies` non-convergence]

   If `\rho(C)\geq 1`, then `C` has an eigenvector `{v}` with
   `\|{v}\|_2=1` and eigenvalue `\lambda` with `|\lambda|>1`. Then

   .. math::
      
      \|C^k{v}\|_2 = \|\lambda^k{v}\|_2 = |\lambda|^k\|{v}\|_2\geq 1,
  
   which does not converge to zero.

   [Proof that `\rho(C)< 1\implies` convergence]

   Assume a linearly independent eigenvalue expansion (not necessary
   for the proof but it simplifies things a lot)
   `{z} = \sum_{i=1}^n\alpha_i{v}_i`. Then,
   
   .. math::

      C^k{z} = \sum_{i=1}^n\alpha_iC^k{v}_i
      = \sum_{i=1}^n\alpha_i\lambda^k{v}_i\to 0.

For symmetric matrices `B`, `\rho(B)=\|B\|_2`, so
the two convergence theorems are related.
If `\|C\|_p=c<1`, then

.. math::

   \|{e}^{k+1}\|_p = \|C{e}^k\|_p \leq \|C\|_p\|{e}^k\|_p
   = c\|{e}^k\|_p.

This guarantees that the error will be reduced by a factor of at least
`c` in each iteration. If we only have `\rho(C)<1`, not `\|C\|_p<1`
then the error may not converge monotonically.

.. proof:example:: Range of SOR parameter

   We can use this to analyse the SOR parameter `\omega`.

   .. math::
      
      \left(\frac{1}{\omega}D+L\right){x}^{k+1} =
      {b}+\left(\left(\frac{1}{\omega}-1\right)D-U\right){x}^k
  
   What values of `\omega`? For SOR,iteration matrix `C` is

   .. math::
      
      C = \left(\frac{1}{\omega}D+L\right)^{-1}
      \left(\frac{1-\omega}{\omega}D-U\right) = (D+\omega L)^{-1}
      ((1-\omega)D-\omega U).

   so

  .. math::

     \det(C) & = 
     \det\left((D+\omega L)^{-1}
     ((1-\omega)D-\omega U)\right) \\
     & =  \det\left((D+\omega L)^{-1}\right)\det\left(
     (1-\omega)D - \omega U\right) \\
     & =  \det\left(D^{-1}\right)\det(D)\det\left((I-\omega I) -
     \omega D^{-1}U\right)\\
     & =  \det\left((1-\omega)I\right) = (1-\omega)^n.

  The determinant is the product of the eigenvalues, hence `\rho(C)<1`
  requires `|1-\omega|<1`.

Splitting methods as preconditioners
------------------------------------

Recall that preconditioned GMRES converges well if the eigenvalues
of `\hat{A}^{-1}A` are clustered together.

.. proof:theorem::

   Let `A` be a matrix with splitting `M+N`, such that `\rho(C) < c <
   1`.  Then, the eigenvalues of the left preconditioned matrix
   `\hat{A}^{-1}A` with `\hat{A}=M` are located in a disk of radius
   `c` around `1` in the complex plane.

.. proof:proof::

   .. math::

      C=-M^{-1}N = M^{-1}(M-A) = I-M^{-1}A.

   Then,

   .. math::
      
      1>c>\rho(C)=\rho(I-M^{-1}A),

   and the result follows since `I` and `M^{-1}A` have a simultaneous
   eigendecomposition.

We deduce that good convergence of the GMRES algorithm occurs when `c`
is small.

For symmetric splittings, we have already observed that the iteration
matrix is

.. math::

   C = I-\left(M_s\right)^{-1}A, 

where

.. math::

   M_s = M(M+M^T-A)^{-1}M^T.

For symmetric splittings we can say a little more about the preconditioner.

.. proof:theorem::

   Let `A` be a matrix with splitting
   `M+N`, such that the symmetric splitting has iteration matrix

   .. math::
      \rho(C) = c < 1,

   and assume further that `M_s` is positive definite.

   Then, the eigenvalues of the symmetric preconditioned matrix
   `\hat{A}^{-1}A` are contained in the interval `[1-c,1+c]`.

.. proof:proof::

   We have

   .. math::

      C & = I-\left(M_s\right)^{-1}A, \\
      & = I - \hat{A}^{-1}A,

   so `\rho(I- \hat{A}^{-1}A) = \rho(C) = c`. Further, `M_s` is
   symmetric and positive definite, so there exists a unique symmetric
   positive definite matrix square root `S` such that `SS =
   M_s`. Then,

   .. math::

      M^{S}A = SSA = S(SAS)S^{-1}.

   Thus, `M_sA` is similar to (and therefore has the same eigenvalues as)
   `SAS`, which is symmetric, and therefore has real eigenvalues,
   and the result follows.

Convergence analysis for Richardson
-----------------------------------

First we examine Richardson iteration. In the unscaled case,

.. math::
   {x}^{k+1} = {x}^k - \left(A{x}^k - {b}\right), \quad
   M = I, \, N = A-I, \, \implies C = I-A.

Let `{e}` be an eigenvector of `A` with eigenvalue `\lambda`, so
`A{e}=\lambda{e}`.  Then `(I-A){e}={e}-\lambda{e}=(1-\lambda){e}`.
So, `{e}` is an eigenvector of `I-A` with eigenvalue `1-\lambda`.
Richardson's method will converge if `\rho(C)<1` \emph{i.e.}
`|1-\lambda|<1` for all eigenvalues `\lambda` of `A`.

This is restrictive, which motivates the scaled Richardson iteration,

.. math::
   {x}^{k+1} = {x}^k - \omega\left(A{x}^k - {b}\right), \quad
   M = \frac{I}{\omega}, \, N = A-\frac{I}{\omega}, \, \implies C =
   I-\omega A.

If `A` has eigenvalues `\lambda_1,\lambda_2,\ldots,\lambda_n` then the
iterative matrix `C` has eigenvalues
`1-\omega\lambda_1,1-\omega\lambda_2,\ldots,1-\omega\lambda_n`.  This
requires `|1-\omega\lambda_i|<1`, `i=1,\ldots,n`, for convergence.

If, further, `A` is symmetric positive definite, then all eigenvalues
are real and positive. Then, all of the eigenvalues of `C` lie between
`1-\omega\lambda_{\min}` and `1-\omega\lambda_{\max}`.  We can
minimise `\rho(C)` by choosing
`\omega=2/(\lambda_{\min}+\lambda_{\max})`. The resulting iteration
matrix has spectral radius

.. math::
   
   \rho(C) = 1-2\frac{\lambda_{\min}}{\lambda_{\min}+\lambda_{\max}}
   = \frac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\min}+\lambda_{\max}}.

.. \frame{
   \frametitle{Scaled Richardson iteration}
   The eigenvectors for the model Poisson problem are
   \[
   E_{i,j}^{k,l} = \sin ik\pi h\sin jl\pi h, \quad i,j,k,l=1,\ldots,m-1.
   \]
   with eigenvalues 
   \[
   \lambda_{k,l} = 
   4\left(\sin^2\frac{k\pi h}{2} + \sin^2\frac{l\pi
      h}{2}\right), \quad k,l=1,\ldots,m-1.
   \]
   Min and max eigenvalues are `\lambda_{\min}=8\sin^2(\pi h/2)`
   and `\lambda_{\max}=8\cos^2(\pi h/2)`, so optimal value of `\omega`
   is 
   \[
   \omega = \frac{2}{\lambda_{\min}+\lambda_{\max}}
   = \frac{1}{4}.
   \]
   This is the same as Jacobi's method for our model problem.
   }

   \frame{
   \frametitle{Scaled Richardson iteration}
   The optimal value of `\omega` for the model problem is `\omega=1/4`.\\
   The 2-norm\footnote{Recall that the 2-norm is the ratio of the largest
     and smallest eigenvalues} of `C` is 
   \[
   \|C\|_2 = \frac{\cot^2\left(\frac{\pi h}{2}\right)-1}
   {\cot^2\left(\frac{\pi h}{2}\right)+1} = \cos\pi h 
   = 1-\frac{\pi^2 h^2}{2} + \mathcal{O}(h^4),
   \]
   as `h\to\infty`. \\
   Convergence gets slower and slower as the mesh is refined. \\
   This is a general problem for iterative methods applied to
   discretisations of differential equations.
   }

   %% \frame{
   %% \frametitle{Model problem}
   %% For our model problem:
   %% \begin{eqnarray*}
   %% \rho_J &=& \cos\pi h = 1 - \frac{\pi^2h^2}{2} + \mathcal{O}(h^4).\\
   %% \rho_{GS} &=& \cos^2\pi h = 1 - \pi^2h^2 + \mathcal{O}(h^4).\\
   %% \omega_{opt} & = & \frac{2}{1+\sin\pi h} = 2-2\pi h + \mathcal{O}(h^2)
   %% \\
   %% \rho_{SOR} & = & \frac{1-\sin\pi h}{1+\sin\pi h} = 1-2\pi h + \mathcal{O}(h^2).
   %% \end{eqnarray*}
   %% }

   \frame{
   \frametitle{Knowing when to stop}
   Instead of looking at the error for convergence we can look at:
   \begin{enumerate}
   \item The residual `{r}^k=A{x}^k-{b}`, or
   \item The pseudo-residual `{s}^k = {x}^{k+1}-{x}^k`,
   \end{enumerate}
   which both tend to zero as `{x}^k\to{x}^*`. \\
   How do their sizes relate to the size of `{e}^k={x}^k-{x}^*`?
   }

   \frame{
   \frametitle{From error to residual}
   \begin{eqnarray*}
   {e}^k & = & {x}^* - {x}^k \\
   & = & A^{-1}(A{x}^*-A{x}^k) \\
   & = & A^{-1}({b}-A{x}^k) \\
   & = & A^{-1}{r}^k.
   \end{eqnarray*}
   We need to estimate the size of `A^{-1}`.
   }

   \frame{
   \frametitle{Residual: model problem}
   For our model Poisson problem:
   \[
   \|A^{-1}\|_2 = \frac{1}{8\sin^2(\pi h/2)},
   \]
   so
   \begin{eqnarray*}
   \|{e}^k\|_2 & \leq & \|A^{-1}\|_2\|{r}^k\|_2 \\
   & = & \frac{\|{r}^k\|_2}{8\sin^2(\pi h/2)} \\
   & \approx & \frac{\|{r}^k\|_2}{2\pi^2 h^2},
   \end{eqnarray*}
   so `\|{e}^k\|_2` can be much bigger than `\|{r}^k\|_2` for small `h`.
   }

   \frame{
   \frametitle{Pseudo-residual}
   \begin{eqnarray*}
   {s}^k & = & {e}^{k+1} - {e}^k \\
    & = & (-M^{-1}N-I){e}^k \\
   \end{eqnarray*}
   so `{e}^k=(-M^{-1}N-I)^{-1}{s}^k`. \\
   For convergent methods, `\rho(-M^{-1}N)<1` so `(-M^{-1}N-I)` is non-singular.
   }

   \frame{
   \frametitle{Pseudo-residual: Richardson's method}
   Richardsons method with `A` symmetric positive definite and
   optimal `\omega`:
   \[
   M=\frac{1}{\omega}I, \, N=A-\frac{1}{\omega}I, \quad \omega=\frac{2}{\lambda_{\min}+\lambda_{\max}}.
   \]
   %-M^{-1}N-I = -\omega(A-I/\omega)-I = -\omega A
   In this case `(-M^{-1}N-I)^{-1}=-A^{-1}/\omega`, so
   \begin{eqnarray*}
   \|{e}^k\|_2 & \leq & \|(-M^{-1}A-I)^{-1}\|_2\|{s}^k\|_2 \\
   & = &
   \frac{\lambda_{\min}+\lambda_{\max}}{2\lambda_{\min}}\|{s}^k\|_2 \\
   & = & \frac{1}{2}\left(1+\cond_2(A)\right)\|{s}^k\|_2.
   \end{eqnarray*}
   For the model problem this becomes 
   \[
   \|{e}^k\|+2 \leq \frac{\|{s}\|_2}{2\sin^2(\pi h/2)} 
   \approx 2\frac{\|{s}\|_2}{\pi^2 h^2},\quad\mbox{as}\,h\to 0.
   \]
   }

Convergence analysis for symmetric matrices
-------------------------------------------

For a symmetric positive definite matrix `A`, recall the
Rayleigh Quotient formula,

.. math::
   \lambda_{\max}=\max_{{x}\ne
    0}\frac{{x}^TA{x}}{{x}^T{x}}\equiv \|A\|_2^2, \quad
   \lambda_{\min}=\min_{{x}\ne
    0}\frac{{x}^TA{x}}{{x}^T{x}},

implying that

.. math::
   \lambda_{\min}\|{y}\|_2^2\leq {y}^TA{y} 
   \leq\lambda_{\max}\|{y}\|_2^2

for any non-zero vector `{y}`.

.. proof:definition:: `A`-weighted norm
		      
   For symmetric positive definite `A`, we can define the weighted
   vector norm

   .. math::

      \|{x}\|_A = \sqrt{{x}^TA{x}},

   and the corresponding matrix (operator) norm

   .. math::
      
      \|B\|_A = \|A^{1/2}BA^{-1/2}\|_2.

These norms are useful for studying convergence of iterative methods
for `A{x}={b}` in the symmetric positive definite case.

.. proof:theorem::

   For a splitting `A=M+N`, if the (symmetric) matrix `M+M^T-A` is
   positive definite then

   .. math::
	
      \|I-M^{-1}A\|_A<1.

.. proof:proof::
      
   If `{y}=(I-M^{-1}A){x}`, `{w}=M^{-1}A{x}`, then

   .. math::

      \|{y}\|_A^2 & =  ({x}-{w})^TA({x}-{w}) 
      = {x}^TA{x}-2{w}^TM{w} + {w}^TA{w} \\
      &= {x}^TA{x}-{w}^T(M+M^T){w} + {w}^TA{w} \\
      &= {x}^TA{x}-{w}^T(M+M^T-A){w} \\
      &\leq \|{x}\|_A^2 - \mu_{\min}\|{w}\|^2_2,

   where `\mu_{\min}` is the (positive) minimum eigenvalue of `M^T+M-A`.

   Further,

   .. math::
      \|{w}\|_2^2 & =  {x}^TA\left(M^{-1}\right)^TM^{-1}A{x} \\
      &= \left(A^{1/2}{x}\right)^TA^{1/2}
      \left(M^{-1}\right)^TM^{-1}A^{1/2}\left(A^{1/2}{x}\right) \\
      &\geq \hat{\mu}_{\min}\|A^{{1/2}}{x}\|_2^2 = 
      \hat{\mu}_{\min}\|{x}\|^2_A,  

   where `\hat{\mu}_{\min}` is the minimum eigenvalue of
   `A^{1/2}\left(M^{-1}\right)^TM^{-1}A^{1/2}` i.e.  the square
   of the minimum eigenvalue of `M^{-1}A^{1/2}`, which is invertible so
   `\hat{\mu}_{\min}>0`. If `{y}=(I-M^{-1}A){x}`,
   `{w}=M^{-1}A{x}`, then

   .. math::

      \|{y}\|^2_A 
      \leq \left(1-\mu_{\min}\hat{\mu}_{\min}\right)\|{x}\|_A^2<\|{x}\|_A^2.


This enables us to show the following useful result for symmetric
positive definite matrices.
      
.. proof:theorem::
     Let `A` be a symmetric positive definite matrix with
     splitting `A=M+N`, if `M` is positive definite, then

   .. math::
      \rho(I-M^{-1}A) = \|I-M^{-1}A\|_A=\|I-M^{-1}A\|_M.

.. proof:proof::

   .. math::
      I-A^{1/2}M^{-1}A^{1/2} &= A^{1/2}(I-M^{-1}A)A^{-1/2}, \\
      I-M^{-1/2}AM^{-1/2} &= M^{1/2}(I-M^{-1}A)M^{-1/2}, \\

   so `I-M^{-1}A`, `I-A^{1/2}M^{-1}A^{1/2}`, and `I-M^{-1/2}AM^{-1/2}`
   all have the same eigenvalues, since they are similar matrices.
   Hence,

   .. math::
      \rho(I-M^{-1}A) & =  \rho(I-A^{1/2}M^{-1}A^{1/2}) \\
      & =  \|I-A^{1/2}M^{-1}A^{1/2}\|_2 \\
      & =  \|I-M^{-1}A\|_A,

   and similarly for `I-M^{-1/2}AM^{-1/2}`.

The consequence of this is that if `M+M^T-A` is symmetric
positive definite then there is a guaranteed reduction in the `A`-norm
of the error in each iteration. If `M` is also symmetric
positive definite the there is guaranteed reduction in the `M`-norm of
the error in each iteration.

Now we apply this to the convergence of Jacobi iteration.  In this
case `M=D`, so `M^T+M-A=2D-A` which may not be positive definite.
We generalise to scaled Jacobi iteration with `M=D/\omega`.

.. proof:proposition::

   Let `A` be a symmetric positive definite matrix. Let `\lambda` be
   the (real) maximum eigenvalue of `D^{-1/2}AD^{-1/2}`.  If `\omega <
   2/\lambda` then scaled Jacobi iteration converges.

.. proof:proof::

   For scaled Jacobi iteration with `M=D/\omega`, we have
   `M^T+M-A=2D/\omega-A`. We have

   .. math::
      2D/\omega-A=D^{1/2}(2I/\omega-D^{-1}A)D^{1/2},

   so `2D/\omega-A` and `2I/\omega-D^{-1}A` have the same eigenvalues.
   If `\lambda` is the maximum eigenvalue of `D^{-1}A` (which has real
   eigenvalues because it is similar to `D^{-1/2}AD^{1/2}`, a symmetric
   matrix), then `2/\omega-\lambda` is the minimum eigenvalue of
   `2I/\omega-D^{-1}A` and hence of `2D/\omega-A`. Hence, `2D/\omega-A`
   is positive definite (so scaled Jacobi converges) if
   `2/\omega-\lambda>0` i.e. `\omega<2/\lambda`.

.. proof:proposition::

   Let `A` be a symmetric positive definite matrix. Then Gauss-Seidel
   iteration always converges.

.. proof:proof::

   For Gauss-Seidel,

   .. math::
      M^T+M-A = (D+L)^T + D+L - A = D+U+D+L-A = D,

   which is symmetric-positive definite, so Gauss-Seidel always converges.

.. proof:proposition::

   Let `A` be a symmetric positive definite matrix. Then SOR converges
   provided that `0<\omega 2`.

.. proof:proof::

   For SOR,

   .. math::
      M^T+M-A =& \left(\frac{1}{\omega}D+L\right)^T + 
      \frac{1}{\omega}D+L - A \\
      &= \frac{2}{\omega}D
      +U+L-(L+D+U)=\left(\frac{2}{\omega}-1\right)D,

   which is symmetric positive definite provided that
   `0<\omega<2`.

An example matrix
-----------------

We consider stationary methods for an example arising from the
finite difference discretisation of the two point boundary value
problem

.. math::

   -\frac{d^2u}{dx^2} = f, \quad u(0) = u(1) = 0.

Here, `f` is assumed known and we have to find `u`. We approximate
this problem by writing `u_k = u(k/(n+1))` for `k=0,1,2,\ldots,n+1`.
From the boundary conditions we have `u_0=u_{n+1}=0`, meaning we
just have to find `u_k` with `1\leq k \leq n`, that solve the
finite difference approximation

.. math::

   -u_{k-1} + 2u_k - u_{k+1} = f_k, \quad 1\leq k \leq n,

where `f_k=f(k/n)/n^2`, `1\leq k\leq n`. Taking into account the
boundary conditions `u_0=u_{n+1}=0`, we can write this as a matrix
system `Ax=b` with

.. math::

   A = \begin{pmatrix}
   2 & -1 & \cdots & \cdots & 0 \\
   -1 & 2 & -1 & \cdots & 0 \\
   0 & -1 & 2 & \cdots &
   \vdots \\
   \vdots & & & & \vdots \\
   \vdots & 0 & -1 & 2 & -1 \\
   0 & 0 & \cdots & -1 & 2
   \end{pmatrix},
   \quad
   x = \begin{pmatrix}
   u_1 \\
   u_2 \\
   \vdots \\
   u_n 
   \end{pmatrix},
   \quad
   b =
   \begin{pmatrix}
   f_1 \\
   f_2 \\
   \vdots \\
   f_n 
   \end{pmatrix}.

However, it is possible to evaluate `Ax` and to implement our classic
stationary iterative methods without ever forming `A`. This is critically
important for efficient implementations (especially when extending
to 2D and 3D problems).

We introduce this example matrix because it is possible to compute
spectral radii for all of the matrices arising in the analysis
of classic stationary methods. In the next example we consider
Jacobi.

.. proof:example:: Jacobi iteration for the example matrix

   In this case, `D=2I`. Thus in fact, scaled Jacobi and scaled Richardson
   are equivalent. We have to find the maximum eigenvalue of `K=D^{-1}A`.
   We can compute this by knowing that the eigenvectors `v` of `K`
   are all of the form

   .. math::

      v =
      \left(
      \begin{pmatrix}
      \sin(l\pi/(n+1)) \\
      \sin(2l\pi/(n+1)) \\
      \sin(nl\pi/(n+1)) \\
      \end{pmatrix}
      \right),

   with one eigenvector for each value of `0< l <n+1`. This can be proved
   by considering symmetries of the matrix, but here we just assume this
   form and establish that we have eigenvectors after substituting into the
   definition of an eigenvector `Av=\lambda v`. This is a general approach
   that can be tried for any matrices arising in the analysis of convergence
   of classic stationary methods for this example matrix.

   .. math::

      (D^{-1}A u)_k = 
      -u_{k-1}/2 + u_k - u_{k+1}/2 =
      \lambda u_k

   which becomes

   .. math::

      \lambda\sin(kl\pi/(n+1)) = -\sin((l-1)k\pi/(n+1))/2 + \sin(kl\pi/(n+1)) - \sin((l+1)k\pi/(n+1))/2,

   and you can use trigonometric formulae, or write

   .. math::

      \lambda\sin(kl\pi/(n+1) & = 
      -\sin((l-1)k\pi/(n+1))/2 + \sin(lk\pi/(n+1)) - \sin((l+1)k\pi/(n+1))/2\\
      & = \sin(kl\pi/(n+1)) - \Im\left(\exp(ik(l-1)\pi/(n+1)) + \exp(ik(l+1)\pi/(n+1))\right)/2 \\
      & = \sin(kl\pi/(n+1)) - \Im\left(\left(
      \exp(-ik\pi/(n+1)) + \exp(ik\pi/(n+1))\right)\exp(ikl\pi/(n+1))\right)/2 \\
      & = \sin(kl\pi/(n+1)) - \Im\left(\sin(k\pi/(n+1))\exp(ikl\pi/(n+1))\right) \\
      & = \sin(kl\pi/(n+1))(1 - \sin(k\pi/(n+1)))

and we conclude that `\lambda=1-\sin(k\pi/(n+1))` are the eigenvalues
with `0<k<n+1`. The maximum eigenvalue corresponds to `k=1` and `k=n`,
with `\lambda=1-\sin(\pi/(n+1))`. 

The condition `\omega<2/\lambda` thus requires that

.. math::

   \omega < \frac{2}{1-\sin(\pi/(n+1))}.

.. proof:exercise::

   Find the value of `\omega`
   for scaled Jacobi such that the convergence rate is maximised,
   i.e. so that `\rho(C)` is minimised. What happens to this rate
   as `n\to \infty`?
   
.. \frame{ \frametitle{Optimal scaling for SOR: model problem} For our
     model problem it can be shown that the optimal value of `\omega` is
   \[
   \omega = \frac{2}{1+\sqrt{1-\rho_J^2}},
   \]
   where `\rho_J` is the spectral radius of the Jacobi iteration:
   \[
   \rho_J = \rho\left(I-D^{-1}A\right),
   \]
   for which the convergence rate is
   \[
   \rho_{SOR} = \frac{1-\sqrt{1-\rho_{J}^2}}
   {1+\sqrt{1-\rho_{J}^2}}
   \]
   }

   \frame{\frametitle{Optimal scaling for SOR: model problem}
   \begin{eqnarray*}
   \rho_J &=& \cos\pi h = 1 - \frac{\pi^2h^2}{2} + \mathcal{O}(h^4).\\
   \rho_{GS} &=& \cos^2\pi h = 1 - \pi^2h^2 + \mathcal{O}(h^4).\\
   \omega_{opt} & = & \frac{2}{1+\sin\pi h} = 2-2\pi h + \mathcal{O}(h^2)
   \\
   \rho_{SOR} & = & \frac{1-\sin\pi h}{1+\sin\pi h} = 1-2\pi h + \mathcal{O}(h^2).
   \end{eqnarray*}
   \centerline{\includegraphics[width=6cm]{sorconvergence}}
   }

   \frame{
   \frametitle{Symmetric iterative methods}
   \begin{itemize}
   \item If `M+M^T-A` is positive definite, then we know that the splitting
   `A=M+N` gives a convergent method. \\
   \item ...but we also know that the splitting `A=M^T+N^T` gives a
     convergent method too.
   \end{itemize}
   Let's combine them:
   \[ 
   M{x}^{k+\frac{1}{2}}=-N{x}^k + {b}, \quad 
   M^T{x}^{k+1}=-N^T{x}^{k+\frac{1}{2}} + {b}.
   \] 
   For example, SOR (`L=U^T`):
   \[
   (L+D){x}^{k+\frac{1}{2}}=-U{x}^k + {b}, \quad
   (U+D){x}^{k+1}=-L{x}^{k+\frac{1}{2}} + {b}.
   \]
   This is called {Symmetric Successive Over-Relaxation} (SSOR).
   }

   \frame{
   \frametitle{Symmetric iterative methods}
   \[
   M{x}^{k+\frac{1}{2}}=-N{x}^k + {b}, \quad
   M^T{x}^{k+1}=-N^T{x}^{k+\frac{1}{2}} + {b}.
   \]
   The iteration matrix is
   \begin{eqnarray*}
   C &=& \left(I - \left(M^T\right)^{-1}A\right)
   \left(I - M^{-1}A\right) \\
   &=& I-\left(M_s\right)^{-1}A, 
   \end{eqnarray*}
   where
   \[
   M_s = M(M+M^T-A)^{-1}M^T.
   \]
   }

   \frame{
   \frametitle{Symmetric iterative methods}
   \[
   M{x}^{k+\frac{1}{2}}=-N{x}^k + {b}, \quad
   M^T{x}^{k+1}=-N^T{x}^{k+\frac{1}{2}} + {b}.
   \]
   \[
   C =
   I-\left(M_s\right)^{-1}A, 
   \quad
   M_s = M(M+M^T-A)^{-1}M^T.
   \]
   For convergent methods, `\|I-M^{-1}A\|_A<1` and `\|I-(M^T)^{-1}A\|_A<1`, so
   \[
   \|I-(M_s)^{-1}A\|_A \leq \|I-M^{-1}A\|_A\|I-(M^T)^{-1}A\|_A<1,
   \]
   and the method also converges.
   }

   \frame{
   \frametitle{Symmetric iterative methods}
   \[
   M{x}^{k+\frac{1}{2}}=-N{x}^k + {b}, \quad
   M^T{x}^{k+1}=-N^T{x}^{k+\frac{1}{2}} + {b}.
   \]
   \[
   C =
   I-\left(M_s\right)^{-1}A, 
   \quad
   M_s = M(M+M^T-A)^{-1}M^T.
   \]
   We have
   \begin{eqnarray*}
   & & A^{\frac{1}{2}}\left(I - \left(M^T\right)^{-1}A\right)
   \left(I - M^{-1}A\right)A^{-\frac{1}{2}} \\
   & = & \left(I-A^{\frac{1}{2}}M^{-1}A^{\frac{1}{2}}\right)^T
   \left(I-A^{\frac{1}{2}}M^{-1}A^{\frac{1}{2}}\right),
   \end{eqnarray*}
   so
   \[
   \rho\left(I-(M_s)^{-1}A\right)
   =\|I-A^{\frac{1}{2}}M^{-1}A^{\frac{1}{2}}\|_2^2 = \|I-M^{-1}A\|_A^2,
   \]
   which is the square `\rho\left(I-M^{-1}A\right)^2` of the spectral
   radius of the original scheme.
   }

   \frame{
   \frametitle{Symmetric iterative methods}
   \begin{eqnarray*}
   \rho\left(I-(M_s)^{-1}A\right)
   &=&
   \|I-M^{-1}A\|_A^2,\\
   &=& \frac{\sqrt{4\beta}-\sqrt{1/\alpha}}{\sqrt{4\beta}+\sqrt{1/\alpha}},
   \end{eqnarray*}
   instead of 
   \[
   \sqrt{\frac{\sqrt{4\beta}-\sqrt{1/\alpha}}{\sqrt{4\beta}+\sqrt{1/\alpha}}}.
   \]
   }

Chebyshev acceleration
----------------------

Say we have computed iterates `{x}^0,{x}^1,\ldots,{x}^k` using

.. math::
   M{x}^{k+1} = -N{x}^k + {b}.
   
If the method is convergent, then these iterates are homing in on the
solution. Can we use extrapolation through these iterates to obtain a
better guess for the solution?

.. math::
   \mbox{Find}\,c_{jk},\,j=1,\ldots,k,\,\mbox{with}\,
   {y}^k = \sum_{j=0}^kc_{jk}{x}^j,

with `{y}^k` the best possible approximation to `{x}^*`.

The usual iterative method has `c_{kk}=1`, and `c_{jk}=0` for `j<k`.
If `{x}^i={x}^*`, `i=0,1,\ldots,k` then

.. math::
   {y}^k =
   \sum_{j=0}^kc_{jk}{x}^*={x}^*\sum_{j=0}^kc_{jk},

so we need `\sum_{j=0}^kc_{jk}=1`. Subject to this constraint, we
     seek to minimise `{y}^k-{x}^* =
     \sum_{j=0}^kc_{jk}({x}^j-{x}^*)`.

We can interpret this in terms of matrix polynomials
by writing

.. math::
   {x}^*-{y}^k &= \sum_{j=0}^kc_{jk}({x}^*-{x}^j), \\
    & =  \sum_{j=0}^kc_{jk}\left(-M^{-1}N\right)^j{e}^0, \\
    & =  p_k\left(-M^{-1}N\right){e}^0,

where

.. math::
   p_k(X) = c_{0k} + c_{1k}X + c_{2k}X^2 + \ldots + c_{kk}X^k,

with `p_k(1)=1` (from our condition `\sum_{j=0}^kc_{jk}=1)`.

We want to try to minimise `{y}^k-{x}^*` by choosing `c_{0k}`,
`c_{1k}`, `\ldots`, `c_{kk}` so that the eigenvalues of `p_k` are as
small as possible.  If `\lambda` is an eigenvalue of `C=-M^{-1}N`,
then `p_k(\lambda)` is an eigenvalue of `p_k(C)`.  It is not practical
to know all the eigenvalues of a large matrix, so we will develop
methods that work if we know that all eigenvalues of `C` are real, and
satisfy `-1<\alpha<\lambda<\beta<1`, for some constants `\alpha`
and `\beta` (we know that `|\lambda|<1` otherwise the basic method is
not convergent.

If all eigenvalues of `C` are real, and satisfy
`-1<\alpha<\lambda<\beta<1`,
then we try to make `\rho_{\max} = \max_{\alpha\leq t\leq\beta}|p_k(t)|`
as small as possible.
Then, if `\lambda` is an eigenvalue of `C`, then the corresponding
eigenvalue of `p_k(C)` will satisfy
`|\lambda_{p_k}|  =  |p_k(\lambda)| \leq \rho_{\max}`. 
We have reduced the problem to trying to find polynomials `p(t)` that have the
smallest absolute value in a given range, subject to `p(1)=1`.
The solution to this problem is known: Chebyshev polynomials.

.. proof:definition::
   The Chebyshev polynomial of degree `k`, `T_k(t)` is defined by
   the recurrence

   .. math::
      T_0(t) = 1, \, T_1(t)=t, \, T_k(t)=2tT_{k-1}(t)-T_{k-2}(t).

For example: `T_2(t) = 2tT_1(t)-T_0(t) = 2t^2-1`.

If we search for the `k`-th degree polynomial `p_k(t)` that
minimises

.. math::
   \max_{-1\leq t\leq 1}|p_k(t)|

subject to the constraint that the coefficient of `t^k` is `2^{k-1}`
then we get the `k`-th order Chebyshev polynomial `T_k(t)`. The
maximum value is `1`.

This is not quite what we want, so we change variables, to get

.. math::
   T_k\left(\frac{2t-\beta-\alpha}{\beta-\alpha}\right)\quad\mbox{minimises}
   \quad \max_{\alpha\leq t\leq \beta}|p_k(t)|

subject to the constraint that the coefficient of `t^k` is
`2^{2k-1}/(\beta-\alpha)`.
The maximum value is `1`.

Then we scale the polynomial to reach the condition `p_k(0)=1`.

.. math::
   p_k=\frac{T_k\left(\frac{2t-\beta-\alpha}{\beta-\alpha}\right)}
   {T_k\left(\frac{2-\beta-\alpha}{\beta-\alpha}\right)}
   \quad\mbox{minimises}
   \quad \max_{\alpha\leq t\leq \beta}|p_k(t)|
   
subject to the constraint that `p_k(0)=1`.
The maximum value is 

.. math::
   \frac{1}{T_k\left(\frac{2-\beta-\alpha}{\beta-\alpha}\right)}.

Say we have computed iterates `{x}^0,{x}^1,\ldots,{x}^k` using

.. math::
   M{x}^{k+1} = -N{x}^k + {b}.

Write

.. math::
   p_k=\frac{T_k\left(\frac{2t-\beta-\alpha}{\beta-\alpha}\right)}
   {T_k\left(\frac{2-\beta-\alpha}{\beta-\alpha}\right)}

in the form

.. math::
   p_k(t) = c_{0k} + c_{1k}t + c_{2k}t^2 + \ldots + c_{kk}t^k,

then

.. math::
   {y}^k = \sum_{j=0}^kc_{jk}{x}^k.

.. Since `k`-th order Chebyshev acceleration requires `k` iterations of
   the original method, we compute the {average rate of
     convergence}
   \[
   \left(
   \rho\left(p_k(-M^{-1}N)
   \right)\right)^{\frac{1}{k}} = 
   \left(
   \frac{1}{T_k\left(\frac{2-\beta-\alpha}{\beta-\alpha}\right)}
   \right)^{\frac{1}{k}}.
   \]
   This should hopefully be much smaller than `\rho(-M^{-1}N)=\max\left(|\alpha|,|\beta|\right)`.
   }

   \frame{
   \frametitle{Example: Jacobi's method for model problem}
   For our model problem
   \[
   \alpha = -\cos(\pi h), \, \beta =\cos(\pi h), \quad
   \mbox{so}\,\frac{2-\beta-\alpha}{\beta-\alpha}=\frac{1}
   {\cos(\pi h)}.
   \]
   \begin{center}
   \includegraphics[width=8cm]{jacobitop}\\
   \includegraphics[width=8cm]{jacobibottom}
   \end{center}
   }

There appears to be a practical problem: we need to store `{x}^0`,
`{x}^1`, `\ldots`, `{x}^k` in order to calculate `{y}^k`. However,
we can get a formula for `{y}^k` in terms of `{y}^{k-1}` and
`{y}^{k-2}` by using

.. math::
   T_k(t) = 2tT_{k-1}(t)-T_{k-2}(t).

We get

.. math::
   p_k(t) = 2\frac{2t-\beta-\alpha}{\beta-\alpha}
   \frac{T_{k-1}(s)}{T_k(s)}p_{k-1}(t) -
   \frac{T_{k-2}(s)}{T_k(s)}p_{k-2}(t), 

where `s=\frac{2-\beta-\alpha}{\beta-\alpha}`.

After some manipulations we obtain

.. math::
   {y}^k = \omega_k\left({y}^{k-1}-{y}^{k-2}+\gamma{z}^{k-1}
   \right)+{y}^{k-2},

where

.. math::
   \gamma=\frac{2}{2-\beta-\alpha}, \quad M{z}^{k-1}={b}-A{y}^{k-1}.

with starting formulas

.. math::
   {y}^0 & =  {x}^0 \\
   {y}^1 & =  {x}^0 + \gamma M^{-1}({b}-A{x}^0).

Also,

.. math::

   \omega_k = \frac{1}{1-\omega_{k-1}/(4s^2)}, \, \omega_1=2.

(See Golub and Van Loan for details).

Chebyshev can dramatically accelerate preconditioners provided that
the preconditioned operator is positive definite and upper
and lower bounds on the eigenvalues are known.

.. \frame{
   \frametitle{Accelerated Richardson's method}
   If A is symmetric positive definite with eigenvalues `0\leq
   \lambda_{\min}\leq \lambda \leq \lambda_{\max}`, we use the optimal
   scaling `M=2I/(\lambda_{\min}+\lambda_{\max})`.\\
   \[
   C= I-\frac{2A}{\lambda_{\min}+\lambda_{\max}}, \quad
   \mbox{so}\,
   \alpha=\frac{\lambda_{\min}-\lambda_{\max}}{\lambda_{\min}+\lambda_{\max}},
   \,\beta=\frac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\min}+\lambda_{\max}}.
   \]
   Accelerated Richardson method is
   \[
   {y}^k =
   \omega_k\left({y}^{k-1}-{y}^{k-2}+\frac{2}{\lambda_n+\lambda_1}
   ({b}-A{y}^{k-1})\right) + {y}^{k-2}.
   \]
   }

   \frame{
   \frametitle{Accelerated Jacobi's method}
   We need to check that `M^{-1}N=I-D^{-1}A` has real eigenvalues so we
   can use Chebyshev acceleration.
   \begin{enumerate}
   \item  `D^{-\frac{1}{2}}AD^{-\frac{1}{2}}` is symmetric, so it has
     real eigenvalues.
   \item
     `D^{-1}A=D^{-\frac{1}{2}}D^{-\frac{1}{2}}AD^{-\frac{1}{2}}D^{\frac{1}{2}}`so
       `D^{-1}A` has real eigenvalues.
   \item If `\lambda` is an eigenvalue of `I-D^{-1}A` then `1-\lambda` is
     an eigenvalue of `D^{-1}`, so `I-D^{-1}A` has real eigenvalues. 
   \end{enumerate}
   Need some way of estimating or bounding the maximum and minimum
   eigenvalues of `I-D^{-1}A` for acceleration to work well.
   }

   \frame{
   \frametitle{Accelerated SSOR}
   \begin{itemize}
   \item The eigenvalues of `C` for SOR are not guaranteed to be real.
   \item Even if `A` is symmetric positive definite, `C=-(L+D)^{-1}U` is
     not symmetric.
   \item Instead, use SSOR
   \[
   (\omega L+D){x}^{k+\frac{1}{2}}=\omega(-U{x}^k + {b}), \quad
   (\omega U+D){x}^{k+1}=\omega(-L{x}^{k+\frac{1}{2}} +{b}).
   \]
   \item A good (nearly optimal) choice of `\omega` is 
   \[
   \omega=\frac{2}{1+\sqrt{2(1-\rho_J)}}, \quad
   \implies \rho_{SSOR}\leq \frac{1-\sqrt{(1-\rho_J)/2}}{1+\sqrt{(1-\rho_J)/2}}.
   \]
   \end{itemize}
   }

   \frame{
   \frametitle{Accelerated SSOR}
   For the model problem with accelerated SSOR we can choose (possibly
   not totally optimal)
   \[
   \alpha=0, \quad \beta = \frac{1-\sin(\pi h/2)}{1+\sin(\pi h/2)}
   \]
   and we get the following average convergence rates:\newline
   \centerline{
   \includegraphics[width=8cm]{acceleratedSOR}.
   }\\
   Need some way of estimating or bounding the maximum and minimum
   eigenvalues of `A` for acceleration to work well.
   }
