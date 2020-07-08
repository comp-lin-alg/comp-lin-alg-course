.. default-role:: math

Preliminaries
=============

In this preliminary section we revise a few key linear algebra
concepts that will be used in the rest of the course, emphasising
the column space of matrices.

Matrices, vectors and matrix-vector multiplication
--------------------------------------------------

We will consider the multiplication of a vector

   .. math::

      x = \begin{pmatrix} x_1 \\
      x_2 \\
      \vdots
      x_n \\
      \end{pmatrix}, \quad x_i \in \mathbb{C}, \, i=1,2,\ldots,n,
      \mbox{ i.e. } x \in \mathbb{C}^n,

by a matrix

   .. math::

      A = \begin{pmatrix}
      a_{11} & a_{12} & \ldots & a_{1n} \\
      a_{21} & a_{22} & \ldots & a_{2n} \\
      \vdots & \vdots & \ddots & \vdots \\
      a_{m1} & a_{m2} & \ldots & a_{mn} \\
      \end{pmatrix},

i.e. `A\in \mathbb{C}^{m\times n}`. `A` has `m` rows and `n` columns
so that the product

   .. math::

      b = Ax

produces `b \in \mathbb{C}^m`. In this course it is important to
consider the general case where `m \neq n`, which has many applications
in data analysis, curve fitting etc. We will usually state generalities
in this course for vectors over the field `\mathbb{C}`, noting where things
specialise to `\mathbb{R}`.

We can quickly check that the map `x \to Ax` given by matrix
multiplication is a linear map from `\mathbb{C}^n \to \mathbb{C}^m`, since

   .. math::

      A(\alpha x + y) = \alpha Ax + Ay,

for all `x,y \in \mathbb{C}^n` and `\alpha\in \mathbb{C}`.
