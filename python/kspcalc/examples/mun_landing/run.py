#! /usr/bin/env python

# Copyright © 2021 Jonathan Grot

'''Mun landing
'''

import math

import ksp
import cli_tools as ct
import mks_polar_motion as mpm
import mpl_tools as mpt

if __name__ == "__main__" :
    import matplotlib
    import matplotlib.pyplot as plt
    
    ksp.augmentBodyDbs()
    ksp.processBodyDbs()

    rw = ct.RowWriter([40])

    plots = []
    plot_opts = []

    # Stage we will be flying
    s2 = ksp.Stage()
    s2.loadJSON("stage_2.json")
    s2.dumpInfo()

    rw.print("Analysis...")

    # A FlyingStage is used when there are non-central forces acting on the stage.
    
    

    # Plot everything
    fig, ax = plt.subplots()
    bbox = mpt.square_plot_bounds(plots)
    mpt.square_plots(ax, plots, bbox, plot_opts)
    
    plt.show()
