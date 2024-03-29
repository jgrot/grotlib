KSP Calculator
==============

Setup
-----

1. Source the setup_env.sh file.

2. If you have already created a Python virtual environment with Python 3.8 or newer, then skip to step 3.  Otherwise create a virtual environment, activate it, then:

   a. Update pip

        ``$ pip install -U pip``

   b. Install requirements

        ``$ pip install -r requirements.txt``

3. Activate the python virtual environment


Command Line Interface (CLI)
----------------------------
Run the CLI with
  ``$ ksp.py <command> <args>``

Use the argparse help system to get more details.
  ``$ ksp.py -h``
  
  ``$ ksp.py <command> -h``

The current list of commands is:

**craft**
  Analyze a craft's DV.

**dvMap**
  Compute a DV map based on a maneuvers file.
  
**dvOrbit**
  Compute DV to get into a circular orbit at an altitude.

**g**
  Compute gravitational acceleration of a body at an altitude.

**jump**
  Solve analytical (experimental) expression for the thrust or amount of fuel needed to launch a craft to a certain altitude.
  
**orbitV**
  Compute speed at a circular orbit around a body.

**plotfuncs**
  Plot the functors currently used in the various models of kspcalc.

**uconv**
  Perform a unit conversion.

API
---

.. automodule:: ksp
   :members:

Data Directory
--------------

The ``data`` directory holds various types of input files and
databases.  Each subdirectory name indicates the contents.  The
``ksp.py`` module features ``pth*()`` functions in the API that build
absolute paths to locations in the data directory.

Physical Quantities and Units Databases
---------------------------------------

Unless units are explicitly specified for a quantity input, the format
will (in general) be a tuple or list of the form ``( <value>, "<name
of unit>" )``.  Please consult the various units databases in the ``#
Databases for uconv()`` section of ``ksp.py``.
