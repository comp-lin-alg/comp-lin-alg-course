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

In this course, a floating point operation (FLOP) will be any arithmetic
unary or binary operation acting on single numbers (such as `+`, `-`,
`\times`, `\div`, `\sqrt`). Of course, in reality, these different operations
have different relative costs, and codes can be made more efficient by
blending multiplications and additions (fused multiply-adds) for example.
Here we shall simply apologise to computer scientists in the class, and


