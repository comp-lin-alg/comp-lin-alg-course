=========================
 Computational exercises
=========================

In the course notes you will encounter computational exercises for you
to complete.  The object of these exercise is to build up
understanding about how computational linear algebra algorithms
actually work. Along the way you will have the opportunity to pick up
valuable scientific computing skills in coding, software engineering
and rigorous testing. They involve completing unfinished "skeleton"
code.

In this section, we will explain how to get set up to do the
computational exercises.

You can either use the machines in the maths department labs, or
you can use your own laptop running Windows, Linux, or macOS. 

There is a lot of information below, but here is a summary checklist to
check that you have everything ready to do your work.

   1. Install Python, Git and a text editor using the instructions below.
   2. Check out your course repository from Github Classroom.
   3. Create a virtual environment (venv) in your course repository.
   4. Activate the venv.
   5. Install the course module to the venv.
   6. Install `numpy` to the venv.
   7. Remember to activate the venv every time you work on the course module.
      
To follow these steps read the information below.
   
Editing code
============

In order to write the code required for the implementation exercise,
you'll need to use a Python-aware text editor. There are many such
editors available and you can use any you like, but you haven't used
Python and/or Git before, it is a good idea to use Visual Studio Code
(VSCode) which is a Python-aware text editor, since VSCode also
provides a command line and an interface to Git.

Getting the software that you need
==================================

The
core requirements are Python 3, git, and a Python-aware text editor.

Using your own machine (recommended)
------------------------------------

Up to date information on how to install Python, Git and VSCode on a
Windows, Linux or Mac machine is available at the `FONS Computing
<https://imperial-fons-computing.github.io/>`_ page.


Using the Windows machines in the computer labs
-----------------------------------------------

Imperial's lab machines have the software you need installed, in some
cases via the `Software Hub`. To get started, double click the
`Software Hub` icon on the desktop or visit the `Software Hub
<https://softwarehub.imperial.ac.uk/>`_ page.

Alternatively, from Software Hub you can run `git`, which will give
you a terminal to run git and Python commands, and `atom`, which will
provide you with an editor. This is the link to the `VSCode
installation instructions
<https://imperial-fons-computing.github.io/vscode.html>`_.

The command line
================

A lot of the routine activity involved in this module revolves around
executing commands on the Bash terminal (sometimes referred to as the
"command line"). For example you use the terminal to work with the
revision control system. If you're not familiar with the Linux
terminal, then you can read this `brief guide to the terminal
<http://www.tuxarena.com/static/intro_linux_cli.php>`_. That guide
focusses on the Bash shell, which is the one we will use.

.. hint::

   In VSCode you can get a terminal by selecting New Terminal from the
   Terminal menu. This should open a Terminal window on your VS Code
   screen.  To the top right of this window is a pulldown menu to
   select the interpreter, which needs to be Bash. The image to the
   right shows this pulldown menu.

   .. image:: _static/bash.png
      :align: right
      :width: 20ex


Python
======

Your implementation will be written in Python based on a code skeleton
provided. This means that you'll need a certain familiarity with the
Python language. But don't panic! Python is a very easy language to
work with. This module will use Python 3. 

If you haven't done any Python before, then go through `the official
Python tutorial <https://docs.python.org/3/tutorial/index.html>`_
(Sections 1 to 5 are a good start).

The Matlab-like array features of Python are provided by `Numpy
<http://www.numpy.org/>`_ for which there is a `helpful tutorial
<http://wiki.scipy.org/Tentative_NumPy_Tutorial>`_. There is also a
handy `guide for Matlab users
<http://wiki.scipy.org/NumPy_for_Matlab_Users>`_. In that context, the
code provided in this course will always use Numpy arrays, and never
Numpy matrices.

.. hint::

   In VSCode, to ensure you are using the correct Python interpreter,

   1. Go to the View menu and select the Command Palette.
   2. Start typing `Python: Select Interpreter`, and click on it when it
      comes up.
   3. Select the correct Python interpreter from the pulldown menu (on
      Mac, the one you installed from Homebrew or Anaconda, on
      Windows, the one that you installed from Anaconda).
   
.. _bitbucket-git:

GitHub and git
==============

Revision control is a mechanism for recording and managing different
versions of changing software. This enables changes to be tracked and
helps in the process of debugging code, and in managing conflicts when
more than one person is working on the same project. Revision control
can be combined with online hosting to provide secure backups and to
enable you to work on code from different locations.

In this module, you'll use revision control to access the skeleton
files, and to update those files if and when they change. You'll also
use the same revision control system to record the edits you make over
time and to submit your work for feedback and, eventually, marking.

We will be using the revision control system `git
<http://git-scm.com/>`_, which is the current state of the art and is
widely adopted. We'll be combining git with the online hosting service GitHub.

Getting started with git and GitHub
-----------------------------------

The very first thing you'll need is a GitHub account. Navigate to
`GitHub <https://github.com/>`_ and sign up.

.. note::

   Make sure you use your Imperial College email address on
   GitHub. This enables you to request unlimited free private GitHub
   repositories and other goodies by `applying here
   <https://education.github.com/pack>`_. You don't strictly need this
   for this module, but there are some nice things in there that you
   might want anyway.

Next you need to do just a little Git setup. At the Terminal, type the
following::
  
  git config --global user.name "Jane Bloggs"

Obviously you put in your own name rather than "Jane Bloggs". Similarly, you need to set your email::

  git config --global user.email "Jane.Bloggs12@imperial.ac.uk"

Once again, you obviously use your own email address. Now there is a
small setting which makes the output of git colourful and therefore a
lot easier to read::
  
  git config --global color.ui "auto"

.. hint::

   If you are a more confident computer user, you could go ahead and
   set up git to work with ssh, the secure shell. This will save a lot
   of password typing but it's not essential so if you are not so
   confident with computers, you can skip this bit. You can follow these `ssh key generating instructions
   <https://help.github.com/articles/generating-an-ssh-key/>`_.
   
If you haven't used Git before, it might be a good idea to look at the
excellent `git tutorial <https://swcarpentry.github.io/git-novice/>`_
over at Software Carpentry.

Setting up your repository
==========================

We're using a tool called `GitHub classroom
<https://classroom.github.com>`_ to automate the creation of your
copies of the repository. This classroom will be updated for the
2022/23 academic year.

Cloning a local copy
--------------------

At the Terminal on your working machine type::

  git clone <url> comp-lin-alg-course

Substituting your git repository url for <url>. Your git repository
url can be found by clicking on `clone or download` at the top right
of your repository page on GitHub. You have to select the `ssh` version
of the repository, and it may be necessary to set up "ssh keys" for this.


.. hint::

   If you are using VSCode, you can do this by:

   1. Opening the Command Palette using the View menu.
   2. Type `git clone` into the Command Palette prompt and paste in
      the repository URL.



.. hint::

   If you get stuck cloning your repository, try reading the `FONS help on git
   <https://imperial-fons-computing.github.io/git.html>`_.

Setting up your venv
--------------------

We're going to use a Python Virtual Environment (venv). This is a
private Python environment in which we'll install the packages we
need, including our own implementation exercise. This minimises
interference between this project and anything else which might be
using Python on the system.  You need to get this right or we won't be
able to mark your code correctly.

In your Terminal, change folder to the repository that you just
checked out (this should contain folders called `doc`, `cla_utils`,
`test`, etc.). Then, create the venv by typing::

  python3 -m venv clavenv

This creates a venv called "clavenv" (you can choose another name).

In VSCode, you will be asked if you want to make this venv the default
for your project. Select "yes" as this will help to ensure that it is
activated.

.. hint::

   To change folder in the terminal, type `cd <path>` where `<path>`
   is the path to the folder you want to change to. Paths can be
   "absolute" e.g. `/home/users/jbloggs/comp-lin-alg/` or "relative"
   e.g. if you are currently in `/home/users/jbloggs` then you can use
   `comp-lin-alg`.  Typing `pwd` shows the current path, and typing
   `ls` shows the contents of the current folder.  Typing `cd ..`
   changes to the enclosing folder, and typing `cd -` changes back to
   the previous folder. For more information see the "brief guide to
   the terminal" linked above.

.. hint::

   If you get stuck with your venv, try reading the `FONS help on venvs
   <https://imperial-fons-computing.github.io/python.html#python-virtual-environments>`_.

Activating your venv
--------------------

**Every time** you want to work on the implementation exercise, you need
to activate the venv. On Linux or Mac do this in the Terminal with::

  source venv/bin/activate

This assumes that you have already changed folder to the repository
that you just checked out (this should contain folders called `doc`,
`cla_utils`, `test`, etc.). Otherwise, you need to provide the full
path to `venv/bin/activate`.
  
On Windows the command is::

  source venv/Scripts/activate

Obviously if you are typing this in a folder other than the one
containing the venv, you need to modify the path accordingly.

Installing the course package to the venv
-----------------------------------------

In this course we will be working on skeleton code stored as a Python
package in the repository. This means that we will be able to import
everything as a module using `from cla_utils import *` without needing
to be in a particular directory. This is what makes the tests work,
for example.

To do this:
   1. Activate the venv as above.
   2. Change folder to the repository that you just checked out (this
should contain folders called `doc`, `cla_utils`, `test`, etc.).
   3. Type `python -m pip install -e .`

The package has a dependency on Numpy (Numeric Python). To install
Numpy, follow steps 1-2 above and then type `python -m pip install
numpy`.

.. hint::

   Here is some useful information on `Modules
   <https://object-oriented-python.github.io/2_programs_in_files.html#modules>`_
   and `Packages
   <https://object-oriented-python.github.io/2_programs_in_files.html#packages>`_
   that might be useful later.

How to do the computational exercises
=====================================

For the computational exercises, quite a lot of the coding
infrastructure you will need is provided already. Your task is to
write the crucial mathematical operations at key points, as described
on this website.

The code on which you will build is in the ``cla_utils`` folder of
your repository. The code has embedded documentation which is used to
build the :doc:`cla_utils` web documentation.

As you do the exercises, **commit your code** to your repository. This
will build up your computational exercise solution sets. You should
commit code early and often - small commits are easier to understand
and debug than large ones. 

.. hint::

   In Git, we use the Terminal to commit changes and push them to the
   remote repository on Github Classroom. A repository is a record of
   the history of the code as you are working.  To add a file to the
   list of files whose changes will be committed to the repository,
   type `git add <filename> -m <log message>`, where `<log message>`
   is a short description of the changes you made.  To commit those
   changes, type `git commit`. They will now be saved locally.  To
   push these changes to the "remote" repository on Github Classroom,
   type `git push` (you may be asked to set the name of the remote,
   just paste the suggested command into the Terminal). To pull
   changes from the remote repository on Github Classroom, type `git
   pull`. For further features and better explanation, please take
   a look at the Github Tutorial linked above.

.. warning::

   Never use `git add *`, since this will add unwanted files to the
   repository which shouldn't be there. You should never add machine
   specific files such as your venv, or `.pyc` files which are
   temporary machine specific files generated by the Python
   interpreter.  This really slows down the marking process and makes
   the markers grumpy. You should only add the `.py` files that you
   are working on.

.. warning::
   
   Do not commit to the feedback branch.  This branch is just there so
   that we can provide feedback on your changes to the main branch,
   and if you commit there, it will mess up our marking system.

Running your work
=================

If you want to execute your code written in `cla_utils`, this can be
imported into IPython (in the terminal, or using a Jupyter notebook),
or in a script.

To use IPython, type `ipython` in the Terminal (when the venv is
activated). You may need to install it first using `python -m install
ipython`. Then you can import `cla_utils` interactively using `from
cla_utils import *`. To exit IPython type Ctrl-D.

.. hint::

   Here is some `Information about IPython <https://object-oriented-python.github.io/2_programs_in_files.html#ipython>`_.

If you also import `numpy` then you can create
example `numpy` arrays and pass them to `cla_utils` functions to try
them out. You can also do this in a script, e.g.::

  from cla_utils import *
  from numpy import *
  A = numpy.array([[1.0,2.0,0.,0.,1.0+1.0j],
                  [0.0,1.0,3.,0.,0.],
		  [0.0,0.0,1.,0.,0.],
		  [0.0,0.0,0.,1.,0.],
		  [0.0,0.0,0.,0.,1.]])
  xr = numpy.array([1.,2.,1.,0.5,0.3])
  xi = numpy.array([1.1,0.2,0.,1.5,-0.7])
  ABiC(A, xr, xi)

After saving your text to a script with a filename ending in `.py`,
e.g.  `run_ABiC.py`, you can execute the script in the Terminal by
typing `python run_ABiC.py` (remember to change to the folder where
the file is located). Scripts are better because you can run the whole
thing again more easily if you make a mistake, and you can save them.

.. hint::

   Here is some `information about Python scripts <https://object-oriented-python.github.io/2_programs_in_files.html#python-scripts-and-text-editors>`_.

.. warning::

   Don't clutter up your repository by adding these experimental scripts
   with `git add`. If you want to store them it is best to use another
   separate git repository for that.

Testing your work
=================

As you complete the exercises, there will often be test scripts which
check the code you have just written. These are located in the
``test`` folder and employ the `pytest <http://pytest.org/>`_
testing framework. You run the tests with:: 

   py.test test_script.py

from the bash Terminal, replacing ``test_script.py`` with the appropriate
test file name. The ``-x`` option to ``py.test`` will cause the test
to stop at the first failure it finds, which is often the best place
to start fixing a problem. For those familiar with debuggers, the
``--pdb`` option will drop you into the Python debugger at the first
error.

You can also run all the tests by running ``py.test`` on the tests
folder. This works particularly well with the -x option, resulting
in the tests being run in course order and stopping at the first
failing test::

  py.test -x tests/

You should make sure that your code passes tests before moving on
to the next exercise.

Coding style and commenting
===========================

Computer code is not just functional, it also conveys information to
the reader. It is important to write clear, intelligible code. **The
readability and clarity of your code will count for marks**.

The Python community has agreed standards for coding, which are
documented in `PEP8
<https://www.python.org/dev/peps/pep-0008/>`_. There are programs and
editor modes which can help you with this. The skeleton implementation
follows PEP8 quite closely. You are encouraged, especially if you are
a more experienced programmer, to follow PEP8 in your
implementation. However nobody is going to lose marks for PEP8
failures.

Skeleton code documentation
===========================

There is web documentation for the complete :doc:`cla_utils`. There is
also an :ref:`alphabetical index <genindex>` and a :ref:`search page<search>`.
