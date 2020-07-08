from distutils.core import setup
from glob import glob

setup(name='finite-element-course',
      version=2018.0,
      description="""Computational linear algebra tasks in support of Computational Linear Algebra 3/4/5 at
Imperial College London""",
      author="Colin Cotter",
      author_email="colin.cotter@imperial.ac.uk",
      url="XXX",
      packages=["cla_utils"],
      scripts=glob('scripts/*'))
