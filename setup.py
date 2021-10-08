from distutils.core import setup
from glob import glob

setup(name='comp-lin-alg-course',
      version=str(2021.0),
      description="""Computational linear algebra course notes and exercises in support of Computational Linear Algebra 3/4/5 at
Imperial College London""",
      author="Colin Cotter",
      author_email="colin.cotter@imperial.ac.uk",
      url="https://comp-lin-alg.github.io/",
      packages=["cla_utils"],
      scripts=glob('scripts/*'))
