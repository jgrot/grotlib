Drag Divergence Optimizer
=========================

This optimizer is a start on trying to tease out the exact drag
divergence model that KSP is using.  I suppose I could just Google it,
but doing it experimentally is more fun.

Command Line Interface (CLI)
----------------------------
Run the CLI with
  ``$ optim_dd.py <command> <args>``

Use the argparse help system to get more details.
  ``$ optim_dd.py -h``
  
  ``$ optim_dd.py <command> -h``

The current list of commands is:

**inst**
  Prints brief instructions on how to use the drag divergence optimizer.

**test**
  Runs the drag optimizer test.


