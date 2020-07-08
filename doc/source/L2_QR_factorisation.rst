.. default-role:: math

QR Factorisation
================

A common theme in computational linear algebra is transformations
of matrices and algorithms to implement them. A transformation is
only useful if it can be computed efficiently and sufficiently
free of pollution from truncation errors (either due to finishing
an iterative algorithm early, or due to round-off errors). A particularly
powerful and insightful transformation is the QR factorisation.
In this section we will introduce the QR factorisation and some
good and bad algorithms to compute it.

What is the QR Factorisation?
-----------------------------

We start with another definition.

.. proof:definition:: Upper triangular matrix

   An `m\times n` upper triangular matrix `R` has coefficients satisfying
   `r_{ij}=0` when `i\geq j`.

   It is called upper triangular because the nonzero rows form a triangle
   on and above the main diagonal of `R`.
   
Now we can describe the QR factorisation.

.. proof:definition:: QR factorisation

   A QR factorisation of an `m\times n` matrix `A` consists of an `m\times m`
   unitary matrix `Q` and an `m\times n` upper triangular matrix `R` such
   that `A=QR`.

The QR factorisation is a key tool in analysis of datasets, and
polynomial fitting. It is also at the core of one of the most widely
used algorithms for finding eigenvalues of matrices. We shall discuss
all of this later during this course.

When `m > n`, `R` must have all zero rows after the `n`th row. Hence,
it makes sense to only work with the top `n\times n` block of `R`
consisting of the first `n` rows, which we call `\hat{R}`. Similarly,
in the matrix vector product `QR`, all columns of `Q` beyond the `n`th
column get multiplied by those zero rows in `R`, so it makes sense to
only work with the first `n` columns of `Q`, which we call `\hat{Q}`.
We then have the reduced QR factorisation, `\hat{Q}\hat{R}`.
