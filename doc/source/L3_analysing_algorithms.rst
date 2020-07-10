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
linear algebra algorithms is the condition number of matrices.


   
