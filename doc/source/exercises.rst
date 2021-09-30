Computational exercises
=======================

In the course notes you will encounter computational exercises for you
to complete.  The object of these exercise is to build up
understanding about how computational linear algebra algorithms
actually work. Along the way you will have the opportunity to pick up
valuable scientific computing skills in coding, software engineering
and rigorous testing. They involve completing unfinished "skeleton"
code.


Obtaining the skeleton code
---------------------------

This section assumes you've already done the :ref:`Git tutorial <bitbucket-git>`.

Setting up your repository
~~~~~~~~~~~~~~~~~~~~~~~~~~

We're using a tool called `GitHub classroom <https://classroom.github.com>`_ to automate the creation of your
copies of the repository.
To create your repository, click `here <https://classroom.github.com/a/OY_NHMeU>`_.


Cloning a local copy
~~~~~~~~~~~~~~~~~~~~

At the command line on your working machine type::

  git clone <url> comp-lin-alg-course

Substituting your git repository url for <url>. Your git repository
url can be found by clicking on `clone or download` at the top right of your repository page on GitHub. 

Setting up your venv
~~~~~~~~~~~~~~~~~~~~

We're going to use a Python venv. This is a private Python environment
in which we'll install the packages we need, including our own
implementation exercise. This minimises interference between this
project and anything else which might be using Python on the
system. We can run a script from the git repository to make the venv::

  ./comp-lin-alg-course/scripts/cla_install_venv venv

This has to install several packages in the venv, so it might take a
few minutes to run.

On Windows, the set of commands is somewhat different. In this case
you would run::

  ./comp-lin-alg-course/scripts/cla_install_venv_win venv

Activating your venv
~~~~~~~~~~~~~~~~~~~~

**Every time** you want to work on the implementation exercise, you need
to activate the venv. On Linux or Mac do this with::

  source venv/bin/activate

while on Windows the command is::

  source venv/Scripts/activate

Obviously if you are typing this in a directory other than the one
containing the venv, you need to modify the path accordingly.
   
Skeleton code documentation
---------------------------

There is web documentation for the complete :doc:`cla_utils`. There is
also an :ref:`alphabetical index <genindex>` and a :ref:`search page<search>`.

How to do the computational exercises
-------------------------------------

For the computational exercises, quite a lot of the coding
infrastructure you will need is provided already. Your task is to
write the crucial mathematical operations at key points, as described
on this website.

The code on which you will build is in the ``cla_utils`` directory of
your repository. The code has embedded documentation which is used to
build the :doc:`cla_utils` web documentation.

As you do the exercises, **commit your code** to your repository. This
will build up your computational exercise solution sets. You should
commit code early and often - small commits are easier to understand
and debug than large ones. **Never** commit back to the ``master``
branch of your fork, that should always remain a clean copy of the
main repository.

Pull requests for feedback
--------------------------

#. Click on the ``New pull request`` button at the top of your
   repository page on GitHub.
#. Make sure **left** dropdown box ("base") is set to ``master``.
#. Make sure **right** dropdown box ("compare") is set to ``implementation``.
#. Type a suitable title in the title box. For example 
   ``Request for feedback 30/1/19``.
#. If you have any comments you would like to pass on to the lecturer
   (for example questions about how you should have done a particular
   exercise) then type these in the ``Description`` box.
#. Click ``Create pull request``.


Testing your work
-----------------

As you complete the exercises, there will often be test scripts which
exercise the code you have just written. These are located in the
``test`` directory and employ the `pytest <http://pytest.org/>`_
testing framework. You run the tests with:: 

   py.test test_script.py

from the bash command line, replacing ``test_script.py`` with the appropriate
test file name. The ``-x`` option to ``py.test`` will cause the test
to stop at the first failure it finds, which is often the best place
to start fixing a problem. For those familiar with debuggers, the
``--pdb`` option will drop you into the Python debugger at the first
error.

You can also run all the tests by running ``py.test`` on the tests
directory. This works particularly well with the -x option, resulting
in the tests being run in course order and stopping at the first
failing test::

  py.test -x tests/


Coding style and commenting
---------------------------

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
