#! /usr/bin/env python
#
# Copyright © 2021 Jonathan Grot
# This file is part of grotlib. The latest version of grotlib is
# available at https://github.com/jgrot/grotlib
#
# Grotlib is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Grotlib is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Grotlib.  If not, see <https://www.gnu.org/licenses/>.
#

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
    mun_XY = mpt.sample_circle(mun_r_m, 10000)
    plots.append(mun_XY)
    plot_opts.append(None)
    
    # A FlyingStage is used when there are non-central forces acting on the stage.
    def fthrottle(t, y) :
        m, r, th, vr, om = y

        v, n = mpm.v_and_dir(y)

        h = r - 200E3

        vtarget = h/10.0

        if t < 5.0 :
            # Initial kick out of circ orbit
            return 0.1
        else :
            return min(1.0, 0.01*(max(v-vtarget,0.0)))
    
    def fthustdir(t, y, fs) :
        v, n = mpm.v_and_dir(y)
        retro = (-n[0], -n[1])
        return retro
    
    fs2 = ksp.FlyingStage(s2, "fs2", "Mun", fthrottle, fthustdir)
    fs2_h0 = 30000.0 # meters
    fs2_r0 = mun_r_m + fs2_h0
    fs2_v0 = ksp.orbitV("Mun", (fs2_h0, "m"))
    y_init = [s2.m0_kg, fs2_r0, 0.0, 0.0, fs2_v0/fs2_r0]
    fs2.launch(y_init)

    tend = 100*60.0
    dt = 5.0
    
    # Plot trajectory
    fs2_XY = fs2.sample(0.0, tend, dt)
    plots.append(fs2_XY)
    plot_opts.append({"marker":"o"})
    # Bbox surrounds fs2 traj
    bbox = mpt.square_plot_bounds([fs2_XY])

    # Dump trajectory
    fs2.dumpTraj(0.0, tend, dt)
    
    # Plot everything
    if len(plots) > 0 :
        fig, ax = plt.subplots()
        # bbox = mpt.square_plot_bounds(plots)
        mpt.square_plots(ax, plots, bbox, plot_opts)
    
        plt.show()
