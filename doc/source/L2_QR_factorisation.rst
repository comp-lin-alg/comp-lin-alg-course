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

In the rest of this section we will examine some algorithms for computing
the QR factorisation, before discussing the application to least squares
problems. We will start with a bad algorithm, before moving on to some
better ones.

QR factorisation by classical Gram-Schmidt algorithm
----------------------------------------------------

The classical Gram-Schmidt algorithm for QR factorisation is motivated
by the column space interpretation of the matrix-matrix multiplication
`A = QR`, namely that the jth column `a_j` of `A` is a linear
combination of the orthonormal columns of `Q`, with the coefficients
given by the jth column `r_j` of `R`. 

The first column of `R` only has a non-zero entry in the first row, so
the first column of `Q` must be proportional to `A`, but normalised
(i.e. rescaled to have length 1). The scaling factor is this first row
of the first column of `R`. The second column of `R` has only non-zero
entries in the first two rows, so the second column of `A` must be
writeable as a linear combination of the first two columns of
`Q`. Hence, the second column of `Q` must by the second column of `A`
with the first column of `Q` projected out, and then normalised. The
first row of the second column of `R` is then the coefficient for this
projection, and the second row is the normalisation scaling
factor. The third row of `Q` is then the third row of `A` with the
first two columns of `Q` projected out, and so on.

Hence, finding a QR factorisation is equivalent to finding an
orthonormal spanning set for the columns of `A`, where the span of the
first 'j' elements of the spanning set and of the first `j' columns of
'A' is the same, for 'j=1,\ldots, n'.

Hence we have to find `R` coefficients such that

.. math::

   q_1 = \frac{a_1}{r_11},

   q_2 = \frac{a_2-r_{12}q_1}{r_{22}}

   \vdots

   q_n = \frac{q_n - \sum_{i=1}^{n-1}r_{in}q_i}{r_{nn}},
   
with `(q_1,q_2,\ldots,q_n)` an orthonormal set. The non-diagonal
entries of `R` are found by inner products, i.e.,

.. math::

   r_{ij} = q_i^*a_j, \, i > j,

and the diagonal entries are chosen so that `\|q_i\|=1`, for
`i=1,2,\ldots,n`, i.e.

.. math::

   |r_{jj}| = \left\| a_j - \sum_{i=1}^{j-1} r_{ij} q_i \right\|.

Note that this absolute value does leave a degree of nonuniqueness
in the definition of `R`. It is standard to choose the diagonal entries
to be real and non-negative.

We now present the classical Gram-Schmidt algorithm as pseudo-code.

* FOR `j = 1` TO `n`
  
  * `v_j \gets a_j`
  * FOR `i = 1` TO `j-1`
    
    * `r_{ij} \gets q_i^*a_j`
    * `v_j \gets v_j - r_{ij}q_i`
  * END FOR
  * `r_{jj} \gets \|v_j\|_2`
  * `q_j \gets v_j/r_{jj}`
* END FOR

(Remember that Python doesn't have END FOR statements, but instead
uses indentation to terminate code blocks. We'll write END statements
for code blocks in pseudo-code in these notes.)
