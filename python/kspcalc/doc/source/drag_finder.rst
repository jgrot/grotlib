Drag Finder
===========

Command Line Interface (CLI)
----------------------------
Run the CLI with
  ``$ drag_finder.py <command> <args>``

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

The idea here is to stay out of the drag divergence regime.  This
transition point occurs somewhere around Mach 0.7.  Keep in mind that
the speed of sound increases with decreasing air density, so the
higher the altitude, the lower the density, the faster the speed of
sound.

Run an experiment with decreased fuel and throttle and record the max
altitude and time to crash. Tips:

* Fly straight up (not tracking velocity) with SAS enabled.
* When vertical speed is zero, hit escape and record altitude.
* While paused, check to make sure you are tracking ground velocity.
* Unpause and let the craft crash keeping an eye on the time (instead of watching the crash, which is fun) because sometimes the clock keeps going depending on what survives the crash.
* Record the max altitude and crash time, along with the throttle setting in an experiment file.  An example is shown below.
* If you modified the SRM fuel amount (haven't worked liquid engines into the code yet), add the fuel amount in the engine record in "elist".

Experiment File (ancAslow_test.json)::

  { "nosecone name" : "Aerodynamic Nose Cone - Type A",
    "stage file"    : "ancAslow_stage.json",
    "htarget"       : 2262,
    "crash time"    : 45,
    "throttle"      : 0.5
  }

Stage File - note: the ``"dragco"`` term is ignored in this test - (ancAslow_stage.json)::
  
  {"m0": [1.564, "t"], "elist": [[1, "RT-10", [37.5,"su"]]], "dragco": 0.45102539 }

Once the experiment is recorded, run the drag finder::

  $ drag_finder.py test ancAslow_test.json 
  iter:    1, drag term: 0.10000000, herr: 1.3108e+03, terr:      11.60
  iter:    2, drag term: 0.10500000, herr: 1.2717e+03, terr:      11.30
  iter:    3, drag term: 0.11000000, herr: 1.2340e+03, terr:      11.00
  iter:    4, drag term: 0.11500000, herr: 1.1974e+03, terr:      10.70
  iter:    5, drag term: 0.12500000, herr: 1.1279e+03, terr:      10.20
  iter:    6, drag term: 0.13500000, herr: 1.0629e+03, terr:       9.70
  iter:    7, drag term: 0.15500000, herr: 9.4419e+02, terr:       8.70
  iter:    8, drag term: 0.17500000, herr: 8.3820e+02, terr:       7.90
  iter:    9, drag term: 0.21500000, herr: 6.5582e+02, terr:       6.40
  iter:   10, drag term: 0.25500000, herr: 5.0341e+02, terr:       5.10
  iter:   11, drag term: 0.33500000, herr: 2.5993e+02, terr:       3.00
  iter:   12, drag term: 0.41500000, herr: 7.1249e+01, terr:       1.30
  iter:   13, drag term: 0.57500000, herr: -2.0726e+02, terr:      -1.20
  iter:   14, drag term: 0.49500000, herr: -8.1100e+01, terr:      -0.10
  iter:   15, drag term: 0.33500000, herr: 2.5993e+02, terr:       3.00
  iter:   16, drag term: 0.45500000, herr: -8.6597e+00, terr:       0.60
  iter:   17, drag term: 0.49500000, herr: -8.1100e+01, terr:      -0.10
  iter:   18, drag term: 0.43500000, herr: 3.0260e+01, terr:       1.00
  iter:   19, drag term: 0.47500000, herr: -4.5723e+01, terr:       0.30
  iter:   20, drag term: 0.44500000, herr: 1.0586e+01, terr:       0.80
  iter:   21, drag term: 0.46500000, herr: -2.7399e+01, terr:       0.40
  iter:   22, drag term: 0.45000000, herr: 8.9411e-01, terr:       0.70
  iter:   23, drag term: 0.44500000, herr: 1.0586e+01, terr:       0.80
  iter:   24, drag term: 0.45250000, herr: -3.9009e+00, terr:       0.70
  iter:   25, drag term: 0.44750000, herr: 5.7027e+00, terr:       0.70
  iter:   26, drag term: 0.45125000, herr: -1.5225e+00, terr:       0.70
  iter:   27, drag term: 0.44875000, herr: 3.3697e+00, terr:       0.70
  iter:   28, drag term: 0.45062500, herr: -3.2844e-01, terr:       0.70
  iter:   29, drag term: 0.45125000, herr: -1.5225e+00, terr:       0.70
  iter:   30, drag term: 0.45031250, herr: 2.7200e-01, terr:       0.70
  iter:   31, drag term: 0.45000000, herr: 8.9411e-01, terr:       0.70
  iter:   32, drag term: 0.45046875, herr: -1.0001e-02, terr:       0.70
  iter:   33, drag term: 0.45062500, herr: -3.2844e-01, terr:       0.70
  iter:   34, drag term: 0.45039063, herr: 1.7104e-01, terr:       0.70
  iter:   35, drag term: 0.45054688, herr: -1.6926e-01, terr:       0.70
  iter:   36, drag term: 0.45050781, herr: -8.9620e-02, terr:       0.70
  iter:   37, drag term: 0.45042969, herr: 9.1402e-02, terr:       0.70
  iter:   38, drag term: 0.45048828, herr: -4.9790e-02, terr:       0.70
  iter:   39, drag term: 0.45044922, herr: 2.9810e-02, terr:       0.70
  iter:   40, drag term: 0.45045898, herr: 9.8966e-03, terr:       0.70
  [0.45045898]

