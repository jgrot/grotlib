KSP - Main Library Module
=========================

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
  ``$ python -m ksp <command> <args>``

Use the argparse help system to get more details.
  ``$ python -m ksp -h``
  
  ``$ python -m ksp <command> -h``

The current list of commands is:

**craft**
  Analyze a craft's DV

**dvMap**
  Compute a DV map based on a maneuvers file
  
**dvOrbit**
  Compute DV to get into a circular orbit at an altitude

**g**
  Compute gravitational acceleration of a body at an altitude
  
**orbitV**
  Compute speed at a circular orbit around a body

API
---

.. automodule:: ksp
   :members:
