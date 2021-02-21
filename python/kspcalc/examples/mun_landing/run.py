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

    # Mun
    mun_r_m = ksp.dInterp("R('Mun')")
    # Plot Mun
    mun_XY = mpt.sample_circle(mun_r_m, 100)
    plots.append(mun_XY)
    plot_opts.append(None)
    
    # A FlyingStage is used when there are non-central forces acting on the stage.
    def fthrottle(t, y) :
        m, r, th, vr, om = y

        v, n = mpm.v_and_dir(y)

        if v > 10.0 :
            return min(1.0, 0.0003*(v-10.0))
        else :
            return 0.0
    
    def falpha(t, y, fs) :
        v, n = mpm.v_and_dir(y)
        retro = (-n[0], -n[1])
        return retro
    
    fs2 = ksp.FlyingStage(s2, "fs2", "Mun", fthrottle, falpha)
    fs2_h0 = 8000.0 # meters
    fs2_r0 = mun_r_m + fs2_h0
    fs2_v0 = ksp.orbitV("Mun", (fs2_h0, "m"))
    y_init = [s2.m0_kg, fs2_r0, 0.0, 0.0, fs2_v0/fs2_r0]
    fs2.launch(y_init)

    # Plot trajectory
    fs2_XY = fs2.sample(0.0, 100*60.0, 5.0)
    plots.append(fs2_XY)
    plot_opts.append({"marker":"o"})
    # Bbox surrounds fs2 traj
    bbox = mpt.square_plot_bounds([fs2_XY])

    # Plot everything
    if len(plots) > 0 :
        fig, ax = plt.subplots()
        # bbox = mpt.square_plot_bounds(plots)
        mpt.square_plots(ax, plots, bbox, plot_opts)
    
        plt.show()
