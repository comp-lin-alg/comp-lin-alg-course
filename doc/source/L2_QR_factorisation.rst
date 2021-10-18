

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

What is the QR factorisation?
-----------------------------

.. details:: Supplementary video

   .. vimeo:: 450191857

We start with another definition.

.. proof:definition:: Upper triangular matrix

   An `m\times n` upper triangular matrix `R` has coefficients satisfying
   `r_{ij}=0` when `i > j`.

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

.. proof:exercise::

   The :func:`cla_utils.exercises2.orthog_space` function has been
   left unimplemented. Given a set of vectors `v_1,v_2,\ldots,v_n`
   that span the subspace `U \subset \mathbb{C}^m`, the function
   should find an orthonormal basis for the orthogonal complement
   `U^{\perp}` given by

      .. math::

	 U^{\perp} = \{x \in \mathbb{C}^m: x^*v = 0, \, \forall v \in U\}. 

   It is expected that it will only compute this up to a tolerance.
   You should make use of the built in QR factorisation routine
   :func:`numpy.linalg.qr`. The test script ``test_exercises2.py`` in
   the ``test`` directory will test this function.

In the rest of this section we will examine some algorithms for computing
the QR factorisation, before discussing the application to least squares
problems. We will start with a bad algorithm, before moving on to some
better ones.

QR factorisation by classical Gram-Schmidt algorithm
----------------------------------------------------


   

.. details:: Supplementary video

   .. vimeo:: 450192200

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
first `j` elements of the spanning set and of the first `j` columns of
`A` is the same, for `j=1,\ldots, n`.

Hence we have to find `R` coefficients such that

.. math::

   q_1 = \frac{a_1}{r_{11}},

   q_2 = \frac{a_2-r_{12}q_1}{r_{22}}

   \vdots

   q_n = \frac{a_n - \sum_{i=1}^{n-1}r_{in}q_i}{r_{nn}},
   
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

.. proof:exercise::

   The :func:`cla_utils.exercises2.GS_classical` function has been
   left unimplemented. It should implement the classical Gram-Schmidt
   algorithm above, using Numpy slice notation so that only one Python
   for loop is used. The function should work "in place" by changing
   the values in `A`, without introducing additional intermediate
   arrays (you will need to create a new array to store `R`). The test
   script ``test_exercises2.py`` in the ``test`` directory will test
   this function.

Projector interpretation of Gram-Schmidt
----------------------------------------


   

.. details:: Supplementary video

   .. vimeo:: 450192723

At each step of the Gram-Schmidt algorithm, a projector is applied to
a column of `A`. We have

.. math::

   q_1 = \frac{P_1a_1}{\|P_1a_1\|},

   q_2 = \frac{P_2a_2}{\|P_2a_2\|},

   \vdots

   q_n = \frac{P_na_n}{\|P_na_n\|},

where `P_j` are orthogonal projectors that project out the first `j-1`
columns `(q_1,\ldots,q_{j-1})` (`P_1` is the identity as this set is
empty when `j=1`). The orthogonal projector onto the first `j-1` columns
is `\hat{Q}_{j-1}\hat{Q}_{j-1}^*`, where

.. math::

   \hat{Q}_{j-1} =
   \begin{pmatrix} q_1 & q_2 & \ldots & q_{j-1} \end{pmatrix}.

Hence, `P_j` is the complementary projector, `P_j=I -
\hat{Q}_{j-1}\hat{Q}_{j-1}^*`.

Modified Gram-Schmidt
---------------------


   

.. details:: Supplementary video

   .. vimeo:: 450193303

There is a big problem with the classical Gram-Schmidt algorithm. It
is unstable, which means that when it is implemented in inexact
arithmetic on a computer, round-off error unacceptably pollutes the
entries of `Q` and `R`, and the algorithm is not useable in
practice. What happens is that the columns of `Q` are not quite
orthogonal, and this loss of orthogonality spoils everything. We will
discuss stability later in the course, but right now we will just
discuss the fix for the classical Gram-Schmidt algorithm, which is
based upon the projector interpretation which we just discussed.

To reorganise Gram-Schmidt to avoid instability, we decompose `P_j`
into a sequence of `j-1` projectors of rank `m-1`, that each project
out one column of `Q`, i.e.

.. math::

   P_j = P_{\perp q_{j-1}}\ldots P_{\perp q_2} P_{\perp q_1},

where

.. math::

   P_{\perp q_j} = I - q_jq_j^*.

Then, 

.. math::

   v_j = P_ja_j = P_{\perp q_{j-1}}\ldots P_{\perp q_2}P_{\perp q_1}a_j.

Here we notice that we must apply `P_{\perp q_1}` to all but one
columns of `A`, and `P_{\perp q_2}` to all but two columns of `A`,
`P_{\perp q_3}` to all but three columns of `A`, and so on.

By doing this, we gradually transform `A` to a unitary matrix, as follows.

   .. math::

      A = 
      \begin{pmatrix}
      a_1 & a_2 & a_3 & \ldots & a_n \\
      \end{pmatrix}
      
      \begin{pmatrix}
      q_1 & v_2^1 & v_3^1 & \ldots & v_n^1 \\
      \end{pmatrix}

      \to
      \begin{pmatrix}
      q_1 & q_2 & v_3^2 & \ldots & v_n^2 \\
      \end{pmatrix}

      \ldots
      \to 
      \begin{pmatrix}
      q_1 & q_2 & q_3 & \ldots & q_n \\
      \end{pmatrix}.

Then it is just a matter of keeping a record of the coefficients
of the projections and normalisation scaling factors and storing
them in `R`.

This process is mathematically equivalent to the classical Gram-Schmidt
algorithm, but the arithmetic operations happen in a different order,
in a way that turns out to reduce accumulation of round-off errors.

We now present this modified Gram-Schmidt algorithm as pseudo-code.

* FOR `i = 1` TO `n`

  * `v_i \gets a_i`
* END FOR
* FOR `i = 1` TO `n`
  
  * `r_{ii} \gets \|v_i\|_2`
  * `q_i = v_i/r_{ii}`
    
  * FOR `j = i+1` TO `n`

    * `r_{ij} \gets q_i^*v_j`
    * `v_j \gets v_j - r_{ij}q_i`
  * END FOR
* END FOR

This algorithm can be applied "in place", overwriting the entries
in `A` with the `v` s and eventually the `q` s.

.. proof:exercise::

   The :func:`cla_utils.exercises2.GS_modified` function has been
   left unimplemented. It should implement the modified Gram-Schmidt
   algorithm above, using Numpy slice notation where possible.
   What is the minimal number of Python
   for loops possible?

   The function should work "in place" by changing the values of `A`,
   without introducing additional intermediate arrays. The test script
   ``test_exercises2.py`` in the ``test`` directory will test this
   function.

.. proof:exercise::

   Investigate the mutual orthogonality of the `Q` matrices that are
   produced by your classical and modified Gram-Schmidt
   implementations. Is there a way to test mutual orthogonality
   without writing a loop? Round-off typically causes problems for
   matrices with large condition numbers and large off-diagonal
   values. You could also try the opposite of what was done in
   ``test_GS_classical``: instead of ensuring that all of the entries
   in the diagonal matrix `D` are `\mathcal{O}(1)`, try making some of
   the values small and some large. See if you can find a matrix that
   illustrates the differences in orthogonality between the two
   algorithms.

Modified Gram-Schmidt as triangular orthogonalisation
-----------------------------------------------------


   

.. details:: Supplementary video

   .. vimeo:: 450193575

This iterative transformation process can be written as
right-multiplication by an upper triangular matrix. For
example, at the first iteration,

   .. math::

      \underbrace{
      \begin{pmatrix}
      v_1^0 & v_2^0 & \ldots & v_n^0
      \end{pmatrix}}_{A}
      \underbrace{
      \begin{pmatrix}
      \frac{1}{r_{11}} & -\frac{r_{12}}{r_{11}} & \ldots &
      \ldots & -\frac{r_{11}}{r_{11}} \\
      0 & 1 & 0 & \ldots & 0 \\    
      0 & 0 & 1 & \ldots & 0 \\
      \vdots & \ddots & \ddots & \ldots & \vdots \\
      0 & 0 & 0 & \ldots & 1 \\
      \end{pmatrix}}_{R_1}
      =
      \underbrace{
      \begin{pmatrix}
      q_1 & v_2^1 & \ldots & v_n^1
      \end{pmatrix}}_{A_1}.

To understand this equation, we can use the column space
interpretation of matrix-matrix multiplication. The columns of `A_1`
are linear combinations of the columns of `A` with coefficients
given by the columns of `R_1`.  Hence, `q_1` only depends on `v_1^0`,
scaled to have length 1, and `v_i^1` is a linear combination of
`(v_1^0,v_i^0)` such that `v_i^1` is orthogonal to `q_1`, for `1<i\leq
n`. 

Similarly, the second iteration may be written as

   .. math::

      \underbrace{
      \begin{pmatrix}
      v_1^1 & v_2^1 & \ldots & v_n^1
      \end{pmatrix}}_{A_1}
      \underbrace{
      \begin{pmatrix}
      1 & 0 & 0 &
      \ldots & 0 \\
      0 & \frac{1}{r_{22}} & -\frac{r_{23}}{r_{22}} & \ldots & -\frac{r_{2n}}{r_{nn}} \\      0 & 0 & 1 & \ldots & 0 \\
      \vdots & \ddots & \ddots & \ldots & \vdots \\
      0 & 0 & 0 & \ldots & 1 \\
      \end{pmatrix}}_{R_2}
      =
      \underbrace{
      \begin{pmatrix}
      q_1 & q_2 & v_3^2 \ldots & v_n^2
      \end{pmatrix}}_{A_2}.

It should become clear that each transformation from `A_i` to `A_{i+1}`
takes place by right multiplication by an upper triangular matrix `R_{i+1}`,
which is an identity matrix plus entries in row i. By combining these
transformations together, we obtain

   .. math::

      A\underbrace{R_1R_2\ldots R_n}_{\hat{R}^{-1}} = \hat{Q}.

Since upper triangular matrices form a group, the product of the `R_i`
matrices is upper triangular. Further, all the `R_i` matrices have
non-zero determinant, so the product is invertible, and we can write
this as `\hat{R}^{-1}`. Right multiplication by `\hat{R}` produces the
usual reduced QR factorisation. We say that modified Gram-Schmidt
implements triangular orthogonalisation: the transformation of `A` to
an orthogonal matrix by right multiplication of upper triangular
matrices.

This is a powerful way to view the modified Gram-Schmidt process from
the point of view of understanding and analysis, but of course we do not
form the matrices `R_i` explicitly (we just follow the pseudo-code given
above).

.. proof:exercise::

   In a break from the format so far, the
   :func:`cla_utils.exercises2.GS_modified_R` function has been
   implemented. It implements the modified Gram-Schmidt algorithm in
   the form describe above using upper triangular matrices. This is
   not a good way to implement the algorithm, because of the inversion
   of `R` at the end, and the repeated multiplication by zeros in
   multiplying entries of the `R_k` matrices, which is a
   waste. However it is important as a conceptual tool for
   understanding the modified Gram-Schmidt algorithm as a triangular
   orthogonalisation process, and so it is good to see this in a code
   implementation. Study this function to check that you understand
   what is happening.

   However, the :func:`cla_utils.exercises2.GS_modified_get_R`
   function has not been implemented. This function computes the `R_k`
   matrices at each step of the process. Complete this code. The test
   script ``test_exercises2.py`` in the ``test`` directory will also
   test this function.


Householder triangulation
-------------------------


   

.. details:: Supplementary video

   .. vimeo:: 450199222

This view of the modified Gram-Schmidt process as triangular
orthogonalisation gives an idea to build an alternative algorithm.
Instead of right multiplying by upper triangular matrices to transform
`A` to `\hat{Q}`, we can consider left multiplying by unitary
matrices to transform `A` to `R`,

   .. math::

      \underbrace{Q_n\ldots Q_2Q_1}_{=Q^*}A = R.

Multiplying unitary matrices produces unitary matrices, so we obtain
`A=QR` as a full factorisation of `A`.


   

.. details:: Supplementary video

   .. vimeo:: 450199366

To do this, we need to work on the columns of `A`, from left to right,
transforming them so that each column has zeros below the
diagonal. These unitary transformations need to be designed so that they
don't spoil the structure created in previous columns. The easiest
way to ensure this is construct a unitary matrix `Q_k` with an identity
matrix as the `(k-1)\times (k-1)` submatrix,

   .. math::

      Q_k =
      \begin{pmatrix}
      I_{k-1} & 0 \\
      0 & F \\
      \end{pmatrix}.

This means that multiplication by `Q_k` won't change the first `k-1`
rows, leaving the previous work to remove zeros below the diagonal
undisturbed. For `Q_k` to be unitary and to transform all below
diagonal entries in column `k` to zero, we need the
`(n-k+1)\times(n-k+1)` submatrix `F` to also be unitary, since

   .. math::

      Q_k^* = 
      \begin{pmatrix}
      I_{k-1} & 0 \\
      0 & F^* \\
      \end{pmatrix}, \,
      Q_k^{-1} = 
      \begin{pmatrix}
      I_{k-1} & 0 \\
      0 & F^{-1} \\
      \end{pmatrix}.

We write the `k` th column `v_k^k` of `A_k` as

   .. math::

      v_k^k =
      \begin{pmatrix}
      \hat{v}_k^k \\
      x
      \end{pmatrix},

where `\hat{v}_k^k` contains the first `k-1` entries of `v_k^k`. The column
gets transformed according to

   .. math::

      Q_kv_k^k = \begin{pmatrix}
      \hat{v}_k^k \\
      Fx
      \end{pmatrix}.

and our goal is that `Fx` is zero, except for the first entry (which
becomes the diagonal entry of `Q_kv_k^k`). Since `F` is unitary, we must
have `\|Fx\|=\|x\|`. For now we shall specialise to
real matrices, so we choose to have

   .. math::

      Fx = \pm\|x\|e_1,

where we shall consider the sign later. Complex matrices have a more
general formula for Householder transformations which we shall not
discuss here.

We can achieve this by using a Householder reflector for `F`, which is
a unitary transformation that does precisely what we
need. Geometrically, the idea is that we consider a line joining `x`
and `Fx=\pm\|x\|e_1`, which points in the direction `v=\pm\|x\|e_1-x`. We can
transform `x` to `Fx` by a reflection in the hyperplane `H` that is
orthogonal to `v`. Since reflections are norm preserving, `F` must be
unitary. Applying the projector `P` given by

   .. math::

      Px = \left(I - \frac{vv^*}{v^*v}\right)x,

does half the job, producing a vector in `H`. To do a reflection we
need to go twice as far,

   .. math::

      Fx = \left(I - 2\frac{vv^*}{v^*v}\right)x.

We can check that this does what we want,

   .. math::

      Fx = \left(I - 2\frac{vv^*}{v^*v}\right)x,

         = x - 2\frac{(\pm\|x\|e_1 - x)}{\|\pm\|x\|e_1 - x\|^2}
	 (\pm\|x\|e_1 - x)^*x,

	 = x - 2\frac{(\pm\|x\|e_1 - x)}{\|\pm\|x\|e_1 - x\|^2}
	 \|x\|(\pm x_1 - \|x\|),

	 = x + (\pm\|x\|e_1 - x) = \pm\|x\|e_1,

as required, having checked that (assuming `x` is real)

   .. math::

      \|\pm \|x\|e_1 - x\|^2 = \|x\|^2 \mp 2\|x\|x_1 + \|x\|^2
      = -2\|x\|(\pm x_1 - \|x\|).

We can also check that `F` is unitary. First we check that `F`
is Hermitian,

   .. math::

      \left(I - 2\frac{vv^*}{v^*v}\right)^*
      = I - 2\frac{(vv^*)^*}{v^*v},

      = I - 2\frac{(v^*)^*v^*}{v^*v},

      = I - 2\frac{vv^*}{v^*v} = F.

Now we use this to show that `F` is unitary,
      
   .. math::

      F^*F = \left(I - 2\frac{vv^*}{v^*v}\right)
      \left(I - 2\frac{vv^*}{v^*v}\right)

      = I - 4\frac{vv^*}{v^*v}\frac{vv^*}{v^*v} +
      4 \frac{vv^*}{v^*v}\frac{vv^*}{v^*v} = I,

so `F^*=F^{-1}`. In summary, we have constructed a unitary
matrix `Q_k` that transforms the entries below the diagonal
of the kth column of `A_k` to zero, and leaves the previous
`k-1` columns alone.


   

.. details:: Supplementary video

   .. vimeo:: 450200163

Earlier, we mentioned that there is a choice of sign in `v`.  This
choice gives us the opportunity to improve the numerical stability of
the algorithm. In the case of real matrices, to avoid unnecessary
numerical round off, we choose the sign that makes `v` furthest from
`x`, i.e.

   .. math::

      v = \mbox{sign}(x_1)\|x\|e_1 + x.

(Exercise, show that this choice of sign achieves this.) It is critical
that we use a definition of `\mbox{sign}` that always returns a number
that has magnitude 1, so we conventionally choose `\mbox{sign}(0)=1`.

.. hint::

   Note that the ``numpy.sign`` function has `\mbox{sign}(0)=0`, so
   you need to take care of this case separately in your Python
   implementation.

For complex valued matrices, the Householder reflection uses `x_1/|x_1|`
(except for `x_1=0` where we use 1 as above).
   
We are now in a position to describe the algorithm in
pseudo-code. Here it is described an "in-place" algorithm, where the
successive transformations to the columns of `A` are implemented as
replacements of the values in `A`. This means that we can allocate
memory on the computer for `A` which is eventually replaced with the
values for `R`. To present the algorithm, we will use the "slice"
notation to describe submatrices of `A`, with `A_{k:l,r:s}` being
the submatrix of `A` consisting of the rows from `k` to `l` and
columns from `r` to `s`.

* FOR `k = 1` TO `n`

  * `x = A_{k:m,k}`
  * `v_k \gets \mbox{sign}(x_1)\|x\|_2e_1 + x`
  * `v_k \gets v_k/\|v_k\|`
  * `A_{k:m,k:n} \gets A_{k:m,k:n} - 2v_k(v_k^*A_{k:m,k:n})`.
* END FOR

.. proof:exercise::

   The :func:`cla_utils.exercises3.householder` function has been left
   unimplemented. It should implement the algorithm above, using only
   one loop over `k`. It should return the resulting `R` matrix. The
   test script ``test_exercises3.py`` in the ``test`` directory will
   test this function.

.. hint::

   Don't forget that Python numbers from zero, which will be important
   when implementing the submatrices using Numpy slice notation. 


   

.. details:: Supplementary video

   .. vimeo:: 450201578

Note that we have not explicitly formed the matrix `Q` or the product
matrices `Q_i`. In some applications, such as solving least squares
problems, we don't explicitly need `Q`, just the matrix-vector product
`Q^*b` with some vector `b`. To compute this product, we can just
apply the same operations to `b` that are applied to the columns of
`A`. This can be expressed in the following pseudo-code, working
"in place" in the storage of `b`.

* FOR `k = 1` TO `n`

  * `b_{k:m} \gets b_{k:m} - 2v_k(v_k^*b_{k:m})`
* END FOR

We call this procedure "implicit multiplication".


   

.. details:: Supplementary video

   .. vimeo:: 450202242

.. proof:exercise::

   Show that the implicit multiplication procedure is equivalent to computing
   an extended array

      .. math::

	 \hat{A} = \begin{pmatrix}
	 a_1 & a_2 & \ldots & a_n & b
	 \end{pmatrix}

   and performing Householder on the first `n` rows. Transform the
   equation `Ax=b` into `Rx=\hat{b}` where `QR=A`, and find the form
   of `\hat{b}`, explaining how to get `\hat{b}` from Householder
   applied to `\hat{A}` above. Solving systems with upper triangular
   matrices is much cheaper than solving general matrix systems as
   we shall discuss later.

   Now, say that we want to solve multiple equations

   .. math::

      Ax_i =b_i, i=1,2,\ldots,k,

   which have the same matrix `A` but different right hand sides
   `b=b_i`, `i=1,2,\ldots,k`. Extend this idea above to the case
   `k>1`, by describing an extended `\hat{A}` containing all the `b_i`
   vectors.

   The :func:`cla_utils.exercises3.householder_solve` function has
   been left unimplemented. It takes in a set of right hand side
   vectors `b_1,b_2,\ldots,b_k` and returns a set of solutions
   `x_1,x_2,\ldots,x_k`.  It should construct an extended array
   `\hat{A}`, and then pass it to
   :func:`cla_utils.exercises3.householder`.  If you have not already
   done so, you will need to modified
   :func:`cla_utils.exercises3.householder` to use the ``kmax``
   argument. You may make use of the built-in triangular solve
   algorithm :func:`scipy.linalg.solve_triangular` (we shall consider
   triangular matrix algorithms briefly later). The test script
   ``test_exercises3.py`` in the ``test`` directory will also test this
   function.

If we really need `Q`, we can get it by matrix-vector products with
each element of the canonical basis `(e_1,e_2,\ldots,e_n)`.  This
means that first we need to compute a matrix-vector product `Qx` with
a vector `x`. One way to do this is to apply the Householder
reflections in reverse, since

   .. math::

      Q = (Q_n\ldots Q_2Q_1)^* = Q_1Q_2\ldots Q_n,

having made use of the fact that the Householder reflections are
Hermitian. This can be expressed in the following pseudo-code.

* FOR `k = n` TO `1` (DOWNWARDS)

  * `x_{k:m} \gets x_{k:m} - 2v_k(v_k^*x_{k:m})`
* END FOR

Note that this requires to record all of the history of the `v` vectors,
whilst the `Q^*` application algorithm above can be interlaced with the
steps of the Householder algorithm, using the `v` values as they are
needed and throwing them away. Then we can compute `Q` via

    .. math::

       Q = \begin{pmatrix}
       Qe_1 & Qe_2 & \ldots & Qe_n
       \end{pmatrix},

with each column using the `Q` application algorithm described above.

.. proof:exercise::

   Show that the implicit multiplication procedure applied to the
   columns of `I` produces `Q^*`, from which we can easily obtain `Q`,
   explaining how. Show how to implement this by applying Householder
   to an augmented matrix `\hat{A}` of some appropriate form.

   The :func:`cla_utils.exercises3.householder_qr` function has been
   left unimplemented. It takes in the `m\times n` array `A` and
   returns `Q` and `R`. It should use the method of this exercise to
   compute them by forming an appropriate `\hat{A}`, calling
   :func:`cla_utils.exercises3.householder` and then extracting
   appropriate subarrays using slice notation. The test script
   ``test_exercises3.py`` in the ``test`` directory will also test
   this function.
   

Application: Least squares problems
-----------------------------------


   

.. details:: Supplementary video

   .. vimeo:: 450202726

Least square problems are relevant in data fitting problems,
optimisation and control, and are also a crucial ingredient of modern
massively parallel linear system solver algorithms such as GMRES,
which we shall encounter later in the course. They are a way of
solving "long thin" matrix vector problems `Ax=b` where we want to
obtain `x\in \mathbb{C}^n` from `b\in\mathbb{C}^m` with `A` an
`m\times n` matrix.  Often the problem does not have a solution as it
is overdetermined for `m>n`. Instead we just seek `x` that minimises
the 2-norm of the residual `r=b-Ax`, i.e. `x` is the minimiser of

   .. math::

      min_x \|Ax - b\|^2.

This residual will not be zero in general, when `b` is not in the
range of `A`. The nearest point in the range of `A` to `b` is `Pb`,
where `P` is the orthogonal projector onto the range of `A`. From
:numref:`Theorem {number}<orthogonal_projector>`, we know that
`P=\hat{Q}\hat{Q}^*`, where `\hat{Q}` from the reduced `QR`
factorisation has the same column space as `A` (but with orthogonal
columns).

Then, we just have to solve

   .. math::

      Ax = Pb,

which is now solveable since `Pb` is in the column space of `A` (and
hence can be written as a linear combination of the columns of `A` i.e.
as a matrix-vector product `Ax` for some unknown `x`).

Now we have the reduced `QR` factorisation of `A`, and we can write

   .. math::

      \hat{Q}\hat{R}x = \hat{Q}\hat{Q}^*b.

Left multiplication by `\hat{Q}^*` then gives

   .. math::

      \hat{R}x = \hat{Q}^*b.

This is an upper triangular system that can be solved efficiently using
back-substitution (which we shall come to later.)

.. proof:exercise::

   The :func:`cla_utils.exercises3.householder_ls` function has been
   left unimplemented. It takes in the `m\times n` array `A` and a
   right-hand side vector `b` and solves the least squares problem
   minimising `\|Ax-b\|` over `x`. It should do this by forming an
   appropriate augmented matrix `\hat{A}`, calling
   :func:`cla_utils.exercises3.householder` and extracting appropriate
   subarrays using slice notation, before using
   :func:`scipy.linalg.solve_triangular` to solve the resulting upper triangular
   system, before returning the solution `x`. The test script
   ``test_exercises3.py`` in the ``test`` directory will also test this
   function.

.. hint::

   You will need to do extract the appropriate submatrix to obtain the
   square (and invertible) reduced matrix `\hat{R}`.
