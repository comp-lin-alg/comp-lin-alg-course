.. default-role:: math

Analysing algorithms
====================

In the previous section we saw three algorithms to compute the QR
factorisation of a matrix. They have a beautiful mathematical
structure based on orthogonal projectors. But are they useful? To
answer this we need to know:

#. Is one faster than others?
#. Is one more sensitive than others to small perturbations due to
   early truncation of the algorithm or due to round-off errors?

In this course we will characterise answers to the first question by
operation count (acknowledging that this is an incomplete evaluation
of speed), and answers to the second question by analysing stability.

In this section we will discuss both of these questions by introducing
some general concepts but also looking at the examples of the QR
algorithms that we have seen so far.

Operation count
---------------

.. details:: Supplementary video

   .. vimeo:: 450203625

Operation count is one aspect of evaluating how long algorithms take.
Here we just note that this is not the only aspect, since transferring
data between different levels of memory on chips can be a serious (and
often dominant) consideration, even more so when we consider
algorithms that make use of large numbers of processors running in
parallel. However, operation count is what we shall focus on here.

In this course, a floating point operation (FLOP) will be any
arithmetic unary or binary operation acting on single numbers (such as
`+`, `-`, `\times`, `\div`, `\sqrt{}`). Of course, in reality, these
different operations have different relative costs, and codes can be
made more efficient by blending multiplications and additions (fused
multiply-adds) for example.  Here we shall simply apologise to
computer scientists in the class, and proceed with this
interpretation, since we are just making relative comparisons between
schemes. We shall also concentrate on asymptotic results in the limit
of large `n` and/or `m`.

Operation count for modified Gram-Schmidt
-----------------------------------------

We shall discuss operation counts through the example of the modified
Gram-Schmidt algorithm. We shall find that the operation count
is `\sim 2mn^2` to compute the QR factorisation, where the `\sim` symbol
means

   .. math::

      \lim_{m,n\to \infty}\frac{N_{\mbox{FLOPS}}}{2mn^2} = 1.

To get this result, we return to the pseudocode for the modified Gram-Schmidt
algorithm, and concentrate on the operations that are happening
inside the inner `j` loop. Inside that loop there are two operations,

#. `r_{ij} \gets q^*_iv_i`. This is the inner product of two vectors
   in `\mathbb{R}^m`, which requires `m` multiplications and `m-1` additions,
   so we count `2m-1` FLOPS per inner iteration.
#. `v_j \gets v_j - r_{ij}q_i`. This requires `m` multiplications and `m`
   subtractions, so we count `2m` FLOPS per inner iteration.

At each iteration we require a combined operation count of `\sim 4m` FLOPS.
There are `n` outer iterations over `i`, and `n-i-1` inner iterations
over `j`, which we can estimate by approximating the sum as an integral,

   .. math::

      N_{\mbox{FLOPS}} \sim \sum_{i=1}^n \sum_{j=i+1}^n 4m
      \sim 4m \int_0^n \int_{x}^n x'\,d x' \,d x
      = 4m\frac{n^2}{2} = 2mn^2,

as suggested above.

Operation count for Householder
-------------------------------

In the Householder algorithm, the computation is dominated by the
transformation

   .. math::

      A_{k:m,k:n} \gets A_{k:m,k:n} -
      \underbrace{2v_k\underbrace{(v_k^*A_{k:m,k:n})}_{1}}_{2},

which must be done for each `k` iteration. To evaluate the part marked
1 requires `n-k` inner products of vectors in `\mathbb{C}^{m-k}`, at a
total cost of `\sim 2(n-k)(m-k)` (we already examined inner products
in the previous example). To evaluate the part marked 2 then requires
the outer product of two vectors in `\mathbb{C}^{m-k}` and
`\mathbb{C}^{n-k}` respectively, at a total cost of `(m-k)(n-k)` FLOPs.
Finally two `(k-m)\times(n-k)` matrices are substracted, at cost
`(k-m)(n-k)`. Putting all this together gives `\sim 4(n-k)(m-k)` FLOPs
per `k` iteration.

Now we have to sum this over `k`, so the total operation count is

   .. math::

      4\sum_{k=1}^n(n-k)(m-k) = 4\sum_{k=1}^n(nm - k(n+m) + k^2)

      \sim 4n^2m - 4(n+m)\frac{n^2}{2} + 4\frac{n^3}{3}
      = 2mn^2 - \frac{2n^3}{3}.

.. proof:exercise::

   Compute FLOP counts for the following operations.

   #. `\alpha = x^*y` for `x,y\in \mathbb{C}^m`.
   #. `y = y + ax` for `x,y \in \mathbb{C}^m`, `a\in\mathbb{C}`.
   #. `y = y + Ax` for `x\in \mathbb{C}^n`, `y\in\mathbb{C}^m`, `A\in\mathbb{C}^{m\times n}`.
   #. `C = C + AB` for `A\in \mathbb{C}^{m\times r}`,
      `B\in \mathbb{C}^{r\times n}`, `C\in \mathbb{C}^{m\times n}`.

.. proof:exercise::

   Suppose `D=ABC` where `A\in \mathbb{C}^{m\times n}`,
   `B\in\mathbb{C}^{n\times p}`, `C\in\mathbb{C}^{p\times q}`. This
   can either be computed as `D=(AB)C` (multiply `A` and `B` first,
   then `C`), or `D=A(BC)` (multiply `B` and `C` first, then `A`).
   Compute the FLOP count for both approaches. For which values
   of `m,n,p,q` would the first approach be more efficient?

.. proof:exercise::

   Suppose `W\in \mathbb{C}^{n\times n}` is defined by

   .. math::

      w_{ij} = \sum_{q=1}^n\sum_{p=1}^n x_{ip}y_{pq}z_{qj},

   where `X,Y,Z\in\mathbb{C}^{n\times n}`. What is the FLOP count for
   computing the entries of `W`?

   The equivalent formula

   .. math::

      w_{ij} = \sum_{p=1}^nx_{ip}\left(\sum_{q=1}^ny_{pq}z_{qj}\right),

   computes the bracket contents first for all `p,j`, before doing
   the sum over `p`. What is the FLOP count for this alternative
   method of computing the entries of `W`?

   Using what you have learned, propose an `\mathcal{O}(n^3)` procedure
   for computing `A\in\mathbb{C}^{n\times n}` with entries

   .. math::

      a_{ij} = \sum_{k=1}^n\sum_{l=1}^n\sum_{m=1}^n
      E_{ki}F_{ki}G_{lk}H_{lm}F_{lm}G_{mj}.

.. proof:exercise::

   Let `L_1,L_2\in\mathbb{C}^{m\times m}` be lower triangular
   matrices.  If we apply the usual formula for multiplying matrices,
   we will waste computation time by multiplying numbers by zero and
   then adding the result to other numbers. Describe a more efficient
   algorithm as pseudo-code and compute the FLOP count, comparing with
   the FLOP count for the standard algorithm.
      
Matrix norms for discussing stability
-------------------------------------

.. details:: Supplementary video

   .. vimeo:: 450204495

In the rest of this section we will discuss another important aspect
of analysing computational linear algebra algorithms, stability. To do
this we need to introduce some norms for matrices in addition to the
norms for vectors that we discussed in Section 1.

If we ignore their multiplication properties, matrices in
`\mathbb{C}^{m\times n}` can be added and scalar multiplied, hence we
can view them as a vector space, in which we can define norms, just
as we did for vectors.

One type of norm arises from simply treating the matrix entries as
entries of a vector and evaluating the 2-norm.

.. proof:definition:: Frobenius norm

   The Frobenius norm is the matrix version of the 2-norm, defined as

      .. math::

	 \|A\|_F = \sqrt{\sum_{i=1}^m\sum_{j=1}^nA_{ij}^2}.

(Exercise: show that `\|AB\|_F \leq \|A\|_F\|B\|_F`.)

Another type of norm measures the maximum amount of stretching the matrix
can cause when multiplying a vector.

.. proof:definition:: Induced matrix norm

   Given an `m\times n` matrix `A` and any chosen vector norms
   `\|\cdot\|_{(n)}` and `\|\cdot\|_{(m)}` on `\mathbb{C}^n` and
   `\mathbb{C}^m`, respectively, the induced norm on `A` is

      .. math::

	 \|A\|_{(m,n)} = \sup_{x\in\mathbb{C}^n, x\neq 0}\frac{\|Ax\|_{(m)}}
	 {\|x\|_{(n)}}.

Directly from the definition we can show

   .. math::

      \frac{\|Ax\|_{(m)}}{\|x\|_{(n)}} \leq \sup_{x\in\mathbb{C}^n, x\neq 0}
      \frac{\|Ax\|_{(m)}}
      {\|x\|_{(n)}} = \|A\|_{(m,n)},

and hence `\|Ax\|\leq \|A\|\|x\|` whenever we use an induced matrix norm.

.. _o2norm:

.. proof:exercise::

   We can reformulate the induced definition as a constrained optimisation
   problem

     .. math::
      
	\|A\|_{(m,n)} = \sqrt{\sup_{x\in\mathbb{C}^n, \|x\|^2=1}\|Ax\|_{(m)}^2}.

   Introduce a Lagrange multiplier `\lambda\in \mathbb{C}` to enforce
   the constraint `\|x\|^2=1`.  Consider the case above where the norms
   on `\mathbb{C}^m` and `\mathbb{C}^n` are both 2-norms. Show that
   `\lambda` must be an eigenvalue of some matrix (which you should
   compute). Hence, given those eigenvalues, provide an expression
   for the operator norm of `A`.

   The :func:`cla_utils.exercises4.operator_2_norm` function has been
   left unimplemented. It takes in an `m\times n` matrix `A` and
   returns the operator norm using the procedure in this exercise.
   You may use the built in function :func:`numpy.linalg.eig` to
   compute the eigenvalues of any matrices that you need. (We will
   discuss algorithms to compute eigenvalues later in the course.) The
   test script ``test_exercises4.py`` in the ``test`` directory will
   test this function.

.. proof:exercise::

   Add a function to :mod:`cla_utils.exercises4` to verify the
   inequality `\|Ax\|\leq \|A\|\|x\|` using
   :func:`cla_utils.exercises4.operator_2_norm`, considering various
   `m` and `n`.

Norm inequalities
-----------------

Often it is difficult to find exact values for norms, so we compute upper
bounds using inequalities instead. Here are a few useful inequalities.

.. proof:definition:: Hölder inequality

   Let `x,y\in \mathbb{C}^m`, and `p,q \in \mathbb{R}+` such that
   `\frac{1}{p}+\frac{1}{q} = 1`. Then

      .. math::

	 |x^*y| \leq \|x\|_p\|y\|_q.

In the case `p=q=2` this becomes the Cauchy-Schwartz inequality.

.. proof:definition:: Cauchy-Schwartz inequality

   Let `x,y\in \mathbb{C}^m`. Then

      .. math::

	 |x^*y| \leq \|x\|_2\|y\|_2.

For example, we can use this to bound the operator norm of the outer
product `A=uv^*` of two vectors.

   .. math::

      \|Ax\|_2 = \|uv^*x\|_2 = \|u(v^*x)\|_2 = |v^*x|\|u\|_2
      \leq \|u\|_2\|v\|_2\|x\|_2,

so `\|A\|_2 \leq \|u\|_2\|v\|_2`.

We can also compute bounds for `\|AB\|_2`.

.. proof:theorem::

   Let `A\in \mathbb{C}^{l\times m}`, `B\in \mathbb{C}^{m\times n}`. Then

      .. math::

	 \|AB|_{(l,n)} \leq \|A\|_{(l,m)}\|B\|_{(m,n)}.

.. proof:proof::

      .. math::

	 \|ABx\|_{(l)} \leq \|A\|_{(l,m)}\|Bx\|_{(m)}
	 \leq \|A\|_{(l,m)}\|B\|_{(m,n)}\|x\|_{(n)},

   so

      .. math::

	 \|AB\|_{(l,n)} = \sup_{x\neq 0}\frac{\|ABx\|_{(l)}}{\|x\|_{(n)}}
	 \leq \|A\|_{(l,m)}\|B\|_{(m,n)},

   as required.

.. proof:exercise::

   Add a function to :mod:`cla_utils.exercises4` to verify this
   theorem for various `l`, `m` and `n`.
   
Condition number
----------------

.. details:: Supplementary video

   .. vimeo:: 450205296

The key tool to understanding numerical stability of computational
linear algebra algorithms is the condition number.  The condition
number is a very general concept that measures the behaviour of a
mathematical problem under perturbations. Here we think of a
mathematical problem as a function `f:X\to Y`, where `X` and `Y` are
normed vector spaces (further generalisations are possible). It is
often the case that `f` has different properties under perturbation
for different values of `x\in X`.

.. proof:definition:: Well conditioned and ill conditioned.

   We say that a problem is well conditioned (at `x`) if small changes
   in `x` lead to small changes in `f(x)`. We say that a problem is
   ill conditioned if small changes in `x` lead to large changes in
   `f(x)`.

These changes are measured by the condition number.

.. proof:definition:: Absolute condition number.

   Let `\delta x` be a perturbation so that `x\mapsto x + \delta x`.
   The corresponding change in `f(x)` is `\delta f(x)`,

   .. math::

      \delta f(x) = f(x + \delta x) - f(x).

   The absolute condition number of `f` at `x` is

   .. math::

      \hat{\kappa} = \sup_{\delta x \neq 0}\frac{\|\delta f\|}{\|\delta x\|},

   i.e. the maximum that `f` can change relative to the size of the
   perturbation `\delta x`.

   It is easier to consider linearised perturbations, defining
   a Jacobian matrix `J(x)` such that
   
   .. math::

      J(x)\delta x = \lim_{\epsilon \to 0}
      \frac{f(x+\epsilon\delta x)-f(x)}{\epsilon}

   and then the linear absolute condition number is

   .. math::

      \hat{\kappa} = \sup_{\delta x \neq 0}\frac{\|J(x)\delta x\|}
      {\|\delta x\|} = \|J(x)\|,

   which is the operator norm of `J(x)`.
      
This definition could be improved by measuring this change relative to the
size of `f` itself.
   
.. proof:definition:: Relative condition number.

   The relative condition number of a problem `f` measures the changes
   `\delta x` and `\delta f` relative to the sizes of `x` and `f`.

   .. math::

      \kappa = \sup_{\delta \neq 0}\frac{\|\delta f\|/\|f\|}
      {\|\delta x\|/\|x\|}.

   The linear relative condition number is

   .. math::

      \kappa = \frac{\|J\|/\|f\|}{1/\|x\|} = \frac{\|J\|\|x\|}{\|f\|}.

Since we use floating point numbers on computers, it makes more sense
to consider relative condition numbers in computational linear
algebra, and from here on we will always use them whenever we mention
condition numbers. If `\kappa` is small (`1-100`, say) then we say that
a problem is well conditioned. If `\kappa` is large (`>10^6`, say),
then we say that a problem is ill conditioned.

.. details:: Supplementary video

   .. vimeo:: 450211558

As a first example, consider the problem of finding the square root,
`f:x\mapsto \sqrt{x}`, a one dimensional problem. In this case,
`J=x^{1/2}/2`. The (linear) condition number is

   .. math::

      \kappa = \frac{|x^{-1/2}/2||x|}{|x^{1/2}|}=1/2.

Hence, the problem is well-conditioned.

As a second example, consider the problem of finding the roots of a
polynomial, given its coefficients. Specifically, we consider the
polynomial `x^2 - 2x +1 = (x-1)^2`, which has two roots equal
to 1. Here we consider the change in roots relative to the coefficient
of `x^0` (which is 1). Making a small perturbation to the polynomial,
`x^2 - 2x + 0.9999 = (x-0.99)(x-1.01)`, so a relative change of `10^{-4}`
gives a relative change of `10^{-2}` in the roots. Using the general formula

   .. math::

      r = 1 \pm\sqrt{1-c} = 1 \pm \sqrt{\delta c} \implies
      \delta r = \pm \sqrt{\delta c},

where `r` returns the two roots with perturbations `\delta r` and `c`
is the coefficient of `x^0` with perturbatino `\delta c`.
is the perturbation to the coefficient of `x^0` (so 1 becomes
`1+\delta c`). The (nonlinear) condition number is then
the sup over `\delta c\neq 0` of 

    .. math::

       \frac{|{\delta r}|/|r|}{|\delta c|/|c|}
       = \frac{|{\delta r}|}{|\delta c|} = \frac{|\delta c|^{1/2}}{|\delta c|}
       = |\delta c|^{-1/2} \to \infty \mbox{ as } \delta c \to 0,

so the condition number is unbounded and the problem is
catastrophically ill conditioned. For an even more vivid example, see
the conditioning of the roots of the Wilkinson polynomial.

Conditioning of linear algebra computations
-------------------------------------------

.. details:: Supplementary video

   .. vimeo:: 450211706

We now look at the condition number of problems from linear algebra.
The first problem we examine is the problem of matrix-vector
multiplication, i.e. for a fixed matrix `A\in \mathbb{C}^{m\times n}`,
the problem is to find `Ax` given `x`. The problem is linear,
with `J=A`, so the condition number is

   .. math::

      \kappa = \frac{\|A\|\|x\|}{\|Ax\|}.

When `A` is non singular, we can write `x = A^{-1}Ax`, and

   .. math::

      \|x\| = \|A^{-1}Ax\| \leq \|A^{-1}\|\|Ax\|,

so

   .. math::

      \kappa \leq \frac{\|A\|\|A^{-1}\|\|Ax\|}{\|Ax\|}
      = \|A\|\|A^{-1}\|.

We call this upper bound the condition number `\kappa(A)` of the matrix `A`.

.. details:: Supplementary video

   .. vimeo:: 450212408

The next problem we consider is the condition number of solving
`Ax=b`, with `b` fixed but considering perturbations to `A`. So, we
have `f:A\mapsto x`. The condition number of this problem measures how
small changes `\delta A` to `A` translate to changes `\delta x` to
`x`. The perturbed problem is

   .. math::

      (A + \delta A)(x + \delta x) = b,

which simplifies (using `Ax=b`) to

   .. math::

      \delta A(x + \delta x) + A\delta x = 0,

which is independent of `b`. If we are considering the linear
condition number, we can drop the nonlinear term, and we get

   .. math::

      \delta A x + A \delta x = 0, \implies \delta x = -A^{-1}\delta Ax,

 from which we may compute the bound

   .. math::

      \|\delta x\| \leq \|A^{-1}\|\|\delta A\|\|x\|.    

Then, we can compute the condition number

   .. math::

      \kappa = \sup_{\|\delta A\|\neq 0}
      \frac{\|\delta x\|/\|x\|}{\|\delta A\|/\|A\|}.
      \leq \sup_{\|\delta A\|\neq 0}
      \frac{\|A^{-1}\|\|\delta A\|\|x\|/\|x\|}{\|\delta A\|/\|A\|}.
      = \|A^{-1}\|\|A\| = \kappa(A),

having used the bound for `\delta x`. Hence the bound on the condition
number for this problem is the condition number of `A`.

.. proof:exercise::

   The :func:`cla_utils.exercises4.cond` function has been left
   unimplemented. It takes in an `m\times m` matrix `A` and returns
   the condition number. You should use a method similar to that in
   :numref:`Exercise {number}<o2norm>`, using the
   :func:`numpy.linalg.eig` to compute the eigenvalues of any matrices
   that you need. Try to think about minimising the number of
   eigenvalue calculations you need to do. The test script
   ``test_exercises4.py`` in the ``test`` directory will test this
   function.
   
Floating point numbers and arithmetic
-------------------------------------

.. details:: Supplementary video

   .. vimeo:: 450212648

Floating point number systems on computers use a discrete and finite
representation of the real numbers. One of the first things we can
deduce from this fact is that there exists a largest and a smallest
positive number.  In "double precision", the standard floating point
number format for scientific computing these days, the largest number
is `N_{\max}\approx 1.79\times 10^{308}`, and the smallest number is
`N_{\min}\approx 2.23 \times 10^{-308}`. The second thing that we can
deduce is that there must be gaps between adjacent numbers in the
number system. In the double precision format, the interval `[1,2]` is
subdivided as `(1,1+2^{-52},1+2\times 2^{-52},1+3\times 2^{-52},
\ldots, 2)`. The next interval `[2,4]` is subdivided as `(2, 2 +
2^{-51}, 2 + 2\times 2^{-51}, \ldots, 4)`.  In general, the interval
`[2^j, 2^{j+1}]` is subdivided by multiplying the set subdividing
`[1,2]` by `2^j`. In this representation, the gaps between numbers
scale with the number size. We call this set of numbers the (double
precision) floating point numbers `\mathbb{F}\subset \mathbb{R}`.

A key aspect of a floating point number system is "machine epsilon"
(`\varepsilon`), which measures the largest relative distance between
two numbers. Considering the description above, we see that
`\varepsilon` is the the distance between 1 and the adjacent number, i.e.

   .. math::

      \varepsilon = 2^{-53} \approx 1.11 \times 10^{-16}.

`\varepsilon` defines the accuracy with which arbitrary real numbers
(within the range of the maximum magnitude above) can be approximated
in `\mathbb{F}`.

   .. math::

      \forall x \in \mathbb{R}, \, \exists x'\in \mathbb{F}
      \mbox{ such that } |x-x'| \leq \varepsilon |x|.

.. details:: Supplementary video

   .. vimeo:: 450213018

.. proof:definition:: Floating point rounding function

   We define `f_L:\mathbb{R}\to \mathbb{F}` as the function that rounds
   `x\in \mathbb{R}` to the nearest floating point number.

The following axiom is just a formal presentation of the properties
of floating point numbers that we discussed below.
   
.. proof:definition:: Floating point axiom I

      .. math::

	 \forall x \in \mathbb{R}, \, \exists \epsilon' \mbox{ with }
	 |\epsilon'| \leq \varepsilon, 

	 \mbox{ such that } f_L(x) = x(1+\epsilon').

The arithmetic operations `+,-,\times,\div` on `\mathbb{R}` have
analogous operations `\oplus,\ominus,\otimes`, etc. In general, binary
operators `\odot` (as a general symbol representing the floating point
version of a real arithmetic operator `\cdot` which could be any of the
above) are constructed such that

   .. math::

      x\odot y = f_L(x\cdot y),

for `x,y\in \mathbb{F}`, with `\cdot` being one of `+,-,\times,\div`.

.. proof:definition:: Floating point axiom II

   .. math::
		      
      \forall x,y \in \mathbb{F}, \exists \epsilon' \mbox{ with }
      |\epsilon'|\leq \varepsilon,\mbox{ such that }

      x\odot y = (x\cdot y)(1 + \epsilon').

.. proof:exercise::

   The formula for the roots of a quadratic equation `x^2 - 2px - q=0`
   is well-known,

      .. math::

	 x = p \pm\sqrt{p^2 + q}.

   Show that the smallest root (with the minus sign above) also
   satisfies

      .. math::

	 x = -\frac{q}{p + \sqrt{p^2 + q}}.

   In the case `p=12345678` and `q=1`, compare the result of these two
   methods for computing the smallest root when using double floating
   point arithmetic (the default floating point numbers in
   Python/NumPy). Which is more accurate? Why is this?
      
Stability
---------

.. details:: Supplementary video

   .. vimeo:: 450213263

Stability describes the perturbation behaviour of a numerical algorithm
when used to solve a problem on a computer. Now we have two problems
`f:X\to Y` (the original problem implemented in the real numbers), and
`\tilde{f}:X\to Y` (the modified problem where floating point numbers
are used at each step).

Given a problem `f` (such as computing the QR factorisation), we are given:

#. A floating point system `\mathbb{F}`,
#. An algorithm for computing `f`,
#. A floating point implementation `\tilde{f}` for `f`.

Then the chosen `x\in X` is rounded to `x'=f_L(x)`, and supplied to
the floating point implementation of the algorithm to obtain
`\tilde{f}(x)\in Y`.

Now we want to compare `f(x)` with `\tilde{f}(x)`. We can measure the
absolute error

   .. math::

      \|\tilde{f}(x)-f(x)\|,

 or the relative error (taking into account the size of `f`),

   .. math::

      \frac{\|\tilde{f}(x)-f(x)\|}{\|f(x)\|}.

An aspiration (but an unrealistic one) would be to aim for an algorithm
to accurate to machine precision, i.e. 

   .. math::

      \frac{\|\tilde{f}(x)-f(x)\|}{\|f(x)\|} = \mathcal{O}(\varepsilon),

by which we mean that `\exists C>0` such that

   .. math::

      \frac{\|\tilde{f}(x)-f(x)\|}{\|f(x)\|} \leq C\varepsilon,

for sufficiently small `\varepsilon`. We shall see below that
we have to lower our aspirations depending on the condition number of `A`.

.. proof:definition:: Stability

   An algorithm `\tilde{f}` for `f` is stable if for each `x\in X`,

      .. math::

	 \frac{\|\tilde{f}(x)-f(\tilde{x})\|}{\|f(\tilde{x})\|} = \mathcal{O}(\varepsilon),

   there exists `\tilde{x}` with

      .. math::

	 \frac{\|\tilde{x}-x\|}{\|x\|} = \mathcal{O}(\varepsilon).	 

We say that a stable algorithm gives nearly the right answer to nearly the
right question.

.. details:: Supplementary video

   .. vimeo:: 450213664

.. details:: Supplementary video

   .. vimeo:: 454094432

.. proof:definition:: Backward stability

   An algorithm `\tilde{f}` for `f` is backward stable if for each `x\in X`,
   `\exists\tilde{x}` such that
   
      .. math::

	 \tilde{f}(x) = f(\tilde{x}),
	 \mbox{ with }
	 \frac{\|\tilde{x}-x\|}{\|x\|} = \mathcal{O}(\varepsilon).

A backward stable algorithm gives exactly the right answer to nearly 
the right answer. The following result shows what accuracy we can expect
from a backward stable algorithm, which involves the condition number
of `f`.

.. _accuracy_backward:

.. proof:theorem:: Accuracy of a backward stable algorithm

   Suppose that a backward stable algorithm is applied to solve problem
   `f:X\to Y` with condition number `\kappa` using a floating point
   number system satisfying the floating point axioms I and II. Then
   the relative error satisfies

      .. math::

	 \frac{\|\tilde{f}(x) - f(x)\|}{\|f(x)\|}
	 = \mathcal{O}(\kappa(x)\epsilon).

.. proof:proof::

   Since `\tilde{f}` is backward stable, we have `\tilde{x}` with
   `\tilde{f}(x)=f(\tilde{x})` and `\|\tilde{x}-x\|/\|x\| =
   \mathcal{O}(\varepsilon)` as above.
   Then,

      .. math::

	 \frac{\|\tilde{f}(x)-f(x)\|}{\|f(x)\|} =
	 \frac{\|f(\tilde{x})-f(x)\|}{\|f(x)\|},

	 = 	 \underbrace{\frac{\|f(\tilde{x})-f(x)\|}{\|f(x)\|}
	 \frac{\|x\|}{\|\tilde{x}-x\|}}_{=\kappa}
	 \underbrace{\frac{\|\tilde{x}-x\|}{\|x\|}}_{=\mathcal{O}(\epsilon)},

   as required.

This type of calculation is known as backward error analysis,
originally introduced by Jim Wilkinson to analyse the accuracy of
eigenvalue calculations using the PILOT ACE, one of the early
computers build at the National Physical Laboratory in the late 1940s
and early 1950s. In backward error analysis we investigate the
accuracy via conditioning and stability. This is usually much easier
than forward analysis, where one would simply try to keep a running
tally of errors committed during each step of the algorithm.

Backward stability of the Householder algorithm
-----------------------------------------------

.. details:: Supplementary video

   .. vimeo:: 450214127

We now consider the example of the problem of finding the QR
factorisation of a matrix `A`, implemented in floating point
arithmetic using the Householder method. The input is `A`, and the
exact output is `Q,R`, whilst the floating point algorithm output is
`\tilde{Q},\tilde{R}`. Here, we consider `\tilde{Q}` as the exact
unitary matrix produced by composing Householder rotations made by
the floating point vectors `\tilde{v}_k` that approximate the `v_k`
vectors in the exact arithmetic Householder algorithm.

For this problem, backwards stability means
that there exists a perturbed input `A+\delta A`, with `\|\delta
A\|/\|A\| =\mathcal{O}(\varepsilon)`, such that `\tilde{Q},\tilde{R}`
are exact solutions to the problem, i.e. `\tilde{Q}\tilde{R}=A+\delta
A`. This means that there is very small backward error,

   .. math::

      \frac{\|A-\tilde{Q}\tilde{R}\|}{\|A\|} = \mathcal{O}(\varepsilon).
      
It turns out that the Householder method is backwards stable.

.. proof:theorem::

   Let the QR factorisation be computed for `A` using a floating point
   implementation of the Householder algorithm. This factorisation is
   backwards stable, i.e. the result `\tilde{Q}\tilde{R}` satisfy

      .. math::

	 \tilde{Q}\tilde{R} = A + \delta A, \quad
	 \frac{\|\delta A\|}{\|A\|} = \mathcal{O}(\varepsilon).

.. proof:proof::

   See the textbook by Trefethen and Bau, Lecture 16.

.. proof:exercise::

   The :func:`cla_utils.exercises5.backward_stability_householder`
   function has been left unimplemented. It generates random `Q_1` and
   `R_1` matrices of dimension `m` provided, and forms `A=QR`. It is
   very important that the two matrices `Q_1` and `R_1` are
   uncorrelated (in particular, computing them as the QR factorisation
   of the same matrix would spoil the experiment). To complete the
   function, pass `A` to the built-in QR factorisation function
   :func:`numpy.linalg.qr` (which uses Householder transformations) to
   get `Q_2` and `R_2`. Print out the value of `\|Q_2-Q_1\|`,
   `\|R_2-R_1\|`, `\|A-Q_2R_2\|`. Explain what you see using what you
   know about the stability of the Householder algorithm.

   
Backward stability for solving a linear system using QR
-------------------------------------------------------

.. details:: Supplementary video

   .. vimeo:: 450214601

The QR factorisation provides a method for solving systems of
equations `Ax=b` for `x` given `b`, where `A` is an invertible
matrix. Substituting `A=QR` and then left-multiplying by `Q^*`
gives

   .. math::
   
      Rx = Q^*b = y.

The solution of this equation is `x=R^{-1}y`, but if there is one
message to take home from this course, it is that you should *never*
form the inverse of a matrix. It is especially disasterous to use
Kramer's rule, which the `m` dimensional extension of the formula for
the inverse of `2\times 2` matrices that you learned at
school. Kramer's rule has an operation count scaling like
`\mathcal{O}(m!)` and is numerically unstable. Hence it is so
disasterous that we won't even show the formula for Kramer's rule
here.

There are some better
algorithms for finding the inverse of a matrix if you really need it,
but in almost every situation it is better to *solve* a matrix system
rather than forming the inverse of the matrix and multiplying it.  It
is particularly easy to solve an equation formed from an upper
triangular matrix.  Written in components, this equation is

  .. math::

     R_{11}x_1 + R_{12}x_2 + \ldots + R_{1(m-1)}x_{m-1} + R_{1m}x_m = y_1,

     0x_1 + R_{22}x_2 + \ldots + R_{2(m-1)}x_{m-1} + R_{2m}x_m = y_2,
     
     \vdots

     0x_1 + 0x_2 + \ldots + R_{(m-1)(m-1)}x_{m-1} + R_{(m-1)m}x_m = y_{m-1},
     
      0x_1 + 0x_2 + \ldots + 0x_{m-1} + R_{mm}x_m = y_{m}.    

The last equation yields `x_m` directly by dividing by `R_{mm}`, then
we can use this value to directly compute `x_{m-1}`. This is repeated
for all of the entries of `x` from `m` down to 1. This procedure is
called back substitution, which we summarise in the following
pseudo-code.

* `x_m  \gets y_m/R_{mm}`
* FOR `i= m-1` TO 1 (BACKWARDS)
  
  * `x_i \gets (y_i - \sum_{k=i+1}^mR_{ik}x_k)/R_{ii}`

In each iteration, there are `m-i-1` multiplications and subtractions
plus a division, so the total operation count is `\sim m^2` FLOPs.

In comparison, the least bad way to form the inverse `Z` of `R` is to
write `RZ = I`. Then, the `k`-th column of this equation is

   .. math::

      Rz_k = e_k,

where `z_k` is the kth column of `Z`. Solving for each column
independently using back substitution leads to an operation count of
`\sim m^3` FLOPs, much slower than applying back substitution directly
to `b`. Hopefully this should convince you to always seek an
alternative to forming the inverse of a matrix.

.. proof:exercise::

   The :func:`cla_utils.exercises5.solve_R` function has been left
   unimplemented. It should implement the `\mathcal{O}(m^2)`
   back-substitution algorithm to solve `Rx=b`, with a single loop
   over the columns.
   The test script ``test_exercises5.py`` in the ``test`` directory
   will test this function.


There are then three steps to solving `Ax=b` using QR factorisation.

#. Find the QR factorisation of `A` (here we shall use the Householder
   algorithm).
#. Set `y=Q^*b` (using the implicit multiplication algorithm).
#. Solve `Rx=y` (using back substitution).

So our `f` here is the solution of `Ax=b` given `b` and `A`, and our
`\tilde{f}` is the composition of the three algorithms above. Now we
ask: "Is this composition of algorithms stable?"

We already know that the Householder algorithm is stable, and a
floating point implementation produces `\tilde{Q},\tilde{R}` such that
`\tilde{Q}\tilde{R}=A+\delta A` with `\|\delta
A\|/\|A\|=\mathcal{O}(\varepsilon)`. It turns out that the implicit
multiplication algorithm is also backwards stable, for similar reasons
(as it is applying the same Householder reflections). This means that
given `\tilde{Q}` (we have already perturbed `Q` when forming it using
Householder) and `b`, the floating point implementation gives
`\tilde{y}` which is not exactly equal to `\tilde{Q}^*b`, but instead
satisfies

   .. math::

      \tilde{y}= (\tilde{Q}+\delta{Q})^*b \implies
      (\tilde{Q} + \delta{Q})\tilde{y} = b,

for some perturbation `\delta Q` with `\|\delta
Q\|=\mathcal{O}(\varepsilon)` (note that `\|Q\|=1` because it is
unitary). Note that here, we are treating `b` as fixed and considering
the backwards stability under perturbations to `\tilde{Q}`.

Finally, it can be shown (see Lecture 17 of Trefethen and Bau for a
proof) that the backward substitution algorithm is backward
stable. This means that given `\tilde{y}` and `\tilde{R}`, the
floating point implementation of backward substitution produces
`\tilde{x}` such that

   .. math::

      (\tilde{R} + \delta \tilde{R})\tilde{x} = \tilde{y},

for some upper triangular perturbation such that `\|\delta
\tilde{R}\|/\|\tilde{R}\|=\mathcal{O}(\varepsilon)`.

.. proof:exercise::

   Complete the function :func:`cla_utils.exercises5.back_stab_solve_R`
   so that it verifies backward stability for back substitution, using
   :func:`cla_utils.exercises5.solve_R`.

Using the individual backward stability of these three algorithms,
we show the following result.

.. proof:theorem::

   The QR algorithm to solve `Ax=b` is backward stable, producing
   a solution `\tilde{x}` such that

      .. math::

	 (A+\Delta A)\tilde{x} = b,

   for some `\|\Delta A\|/\|A\|=\mathcal{O}(\varepsilon)`.

.. proof:proof::

   From backward stability for the calculation of `Q^*b`, we have

      .. math::

	 b = (\tilde{Q}+\delta Q)\tilde{y},

	 = (\tilde{Q} + \delta Q)(\tilde{R} + \delta R)\tilde{x},

   having substituted the backward stability formula for back
   substitution in the second line. Multiplying out the brackets
   and using backward stability for the Householder method gives

      .. math::

	 b = (\tilde{Q}\tilde{R} + (\delta Q)\tilde{R} + \tilde{Q}\delta R
	 + (\delta Q)\delta R)\tilde{x},

	 = (A + \underbrace{\delta A + (\delta Q)\tilde{R} +
	 \tilde{Q}\delta R
	   + (\delta Q)\delta R}_{=\Delta A}\tilde{x}).

   This defines `\Delta A` and it remains to estimate each of these
   terms. We immediately have `\|\delta A\|=\mathcal{O}(\varepsilon)`
   from backward stability of the Householder method.

   Next we estimate the second term. Using `A + \delta A =
   \tilde{Q}\tilde{R}`, we have

      .. math::

	 \tilde{R} = \tilde{Q}^*(A + \delta A),

   we have

      .. math::

	 \frac{\|\tilde{R}\|}{\|A\|} \leq \|\tilde{Q}^*\|
	 \frac{\|A+\delta A\|}{\|A\|} = \mathcal{O}(1), \mbox{ as }
	 \varepsilon \to 0.

   Then we have

      .. math::

	 \frac{\|(\delta Q)\tilde{R}\|}{\|A\|}
	 \leq \|\delta Q\|\frac{\|\tilde{R}\|}{\|A\|}
	 = \mathcal{O}(\varepsilon).

   To estimate the third term, we have

      .. math::

	 \frac{\|\tilde{Q}\delta R\|}{\|A\|} \leq \frac{\|\delta
	 R\|}{\|A\|}\underbrace{\|\tilde{Q}\|}_{=1} =
	 \underbrace{\frac{\|\delta
	 R\|}{\|\tilde{R}\|}}_{\mathcal{O}(\varepsilon)}
	 \underbrace{\frac{\|\tilde{R}\|}{\|A\|}}_{\mathcal{O}(1)}
	 = \mathcal{O}(\varepsilon).

   Finally, the fourth term has size

   .. math::

      \frac{\|\delta Q\delta R\|}{\|A\|} \leq
      \underbrace{\|\delta Q\|}_{\mathcal{O}(\varepsilon)}
      \underbrace{\frac{\|\delta R\|}{\|\tilde{R}\|}}_{\mathcal{O}(\varepsilon)}
      \underbrace{\frac{\|\tilde{R}\|}
      {\|A\|}}_{\mathcal{O}(1)} = \mathcal{O}(\epsilon^2),

   hence `\|\Delta A\|/\|A\|=\mathcal{O}(\varepsilon)`.

.. details:: Supplementary video

   .. vimeo:: 450215261

   
.. proof:Corollary::

   When solving `Ax=b` using the QR factorisation procedure above, the
   floating point implementation produces an approximate solution
   `\tilde{x}` with
   
      .. math::

	 \frac{\|\tilde{x}-x\|}{\|{x}\|} = \mathcal{O}(\kappa(A)\varepsilon).
   
.. proof:proof::
   
   From :numref:`Theorem {number}<accuracy_backward>`, using the
   backward stability that we just derived, we know that
   
      .. math::

	 \frac{\|\tilde{x}-x\|}{\|{x}\|} = \mathcal{O}(\kappa\varepsilon),

   where `\kappa` is the condition number of the problem of solving
   `Ax=b`, which we have shown is bounded from above by `\kappa(A)`.

.. proof:exercise::

   Complete the function :func:`cla_utils.exercises5.back_stab_householder_solve`
   so that it verifies backward stability for solving `m\times m` dimensional
   square systems `Ax=b` using :func:`cla_utils.exercises3.householder_solve`.
