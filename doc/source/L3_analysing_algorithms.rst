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
of large $n$ and/or $m$.

Operation count for modified Gram-Schmidt
-----------------------------------------

We shall discuss operation counts through the example of the modified
Gram-Schmidt algorithm. We shall find that the operation count
is `\sim mn^2` to compute the QR factorisation, where the `\sim` symbol
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
      \sim 4m \sum_{i=1}^n i\int_1^n x\,d x
      \sim 4m\frac{n^2}{2} = 2mn^2,

as suggested above.

Operation count for Householder
-------------------------------

In the Householder algorithm, the computation is dominated by the
transformation

   .. math::

      A_{k:m,k:n} \gets A_{k:m,k:n} -
      \underbrace{2v_k\underbrace{(v_k^*A_{k:m,k:n})}_{1}}_{2},

which must be done for each 'k' iteration. To evaluate the part marked
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

Matrix norms for discussing stability
-------------------------------------

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
      {\|x\|_{(n)} = \|A\|_{(m,n)}},

and hence `\|Ax\|\leq \|A\|\|x\|` whenever we use an induced matrix norm.

Norm inequalities
-----------------

Often it is difficult to find exact values for norms, so we compute upper
bounds using inequalities instead. Here are a few useful inequalities.

.. proof:definition:: H\"older inequality

   Let `x,y\in \mathbb{C}^m`, and `p,q \in \mathbb{R}+` such that
   `\frac{1}{p}+{1}{q} = 1`. Then

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

Condition number
----------------

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

      \kappa = \frac{\|J\|/\|f\|}{\|x\|} = \frac{\|J\|\|x\|}{\|f\|}.

Since we use floating point numbers on computers, it makes more sense
to consider relative condition numbers in computational linear
algebra, and from here on we will always use them whenever we mention
condition numbers. If `\kappa` is small (`1-100`, say) then we say that
a problem is well conditioned. If `\kappa` is large (`>10^6`, say),
then we say that a problem is ill conditioned.

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

Floating point numbers and arithmetic
-------------------------------------

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
