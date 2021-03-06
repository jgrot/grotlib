Drag Finder
===========

Command Line Interface (CLI)
----------------------------
Run the CLI with
  ``$ drag_finder.py <command> <args> <options>``

Use the argparse help system to get more details.
  ``$ drag_finder.py -h``
  
  ``$ drag_finder.py <command> -h``

The current list of commands is:

**inst**
  Prints brief instructions on how to use the drag finder.

**test**
  Runs the drag finder test.

Determining Drag Term
---------------------

Currently, KSPCalc lumps the 1/2*dragco*A into one term.

**If simultaneously determining drag divergence:**

If simultaneously determining drag divergence (see :doc:`optim_dd`),
try to stay well out of the drag divergence regime by running a "slow
test," and supply ``--nodd`` to ``drag_finder.py.`` This transition
point occurs somewhere around Mach 0.7.  Keep in mind that the speed
of sound increases with decreasing air density, so the higher the
altitude, the lower the density, the faster the speed of sound.

Run an experiment and record the max altitude and time to crash. Tips:

* Construct a craft stage file and put thrust limiting in this file.  Run with full throttle for liquid engines.  Setting thrust limit is exact while trying to set throttle is not.
* Fly straight up (not tracking velocity) with SAS enabled.
* When vertical speed is zero, hit escape and record altitude.
* While paused
  * check to make sure you are in surface mode and not orbit mode
  * make sure you are tracking prograde direction
* Unpause and let the craft crash keeping an eye on the time (instead of watching the crash, which is fun) because sometimes the clock keeps going depending on what survives the crash.
* Record the max altitude and crash time.  An example is shown below.

Experiment File (ancAslow_test.json)::

  { "nosecone name" : "Aerodynamic Nose Cone - Type A",
    "stage file"    : "ancAslow_stage.json",
    "htarget"       : 2241,
    "crash time"    : 46
  }

Stage File - note: the ``"dragco"`` term is ignored in this test - (ancAslow_stage.json)::
  
  {"m0": [1.564, "t"], "elist": [["RT-10", 1, 37, [37.5,"su"]]], "dragco": 0.0 }

Once the experiment is recorded, run the drag finder.  Note, this is a slow test to avoid drag divergence, so we use the ``--nodd`` option::

  $ drag_finder.py test ancAslow_test.json --nodd
  iter:    1, drag term: 10.00000000, herr: -1.6864e+03, terr:     -16.90
  iter:    2, drag term: 10.50000000, herr: -1.6985e+03, terr:     -17.00
  iter:    3, drag term: 9.50000000, herr: -1.6733e+03, terr:     -16.80
  iter:    4, drag term: 9.00000000, herr: -1.6592e+03, terr:     -16.60
  iter:    5, drag term: 8.00000000, herr: -1.6270e+03, terr:     -16.30
  ...
  iter:   47, drag term: 0.51994324, herr: 1.5159e-01, terr:       0.70
  iter:   48, drag term: 0.51994324, herr: 1.5159e-01, terr:       0.70
  iter:   49, drag term: 0.51993179, herr: 3.9577e-03, terr:       0.70
  iter:   50, drag term: 0.51992798, herr: 3.1249e-02, terr:       0.70
  iter:   51, drag term: 0.51993370, herr: -5.8204e-03, terr:       0.70
  [0.51993179]
