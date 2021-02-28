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
    s1 = ksp.Stage()
    s1.loadJSON("TestLFE1.json")
    s1.dumpInfo()

    def fthrottle( t, y ) :
        return 1.0
    def falpha( t, y, flyer ) :
        m, r, th, vr, om = y
        vth = om/r
        #if r < (flyer.R+7.5E3) :
        n = (1.0, 0.0)
        #else :
        #    v,n = mpm.v_and_dir(y)
        return n
    
    flyer = ksp.FlyingStage(s1, "Stage 1", "Kerbin", fthrottle, falpha)
    flyer.launch()
    # Fly until crash
    flyer.flyTo(30000)

    import matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    
    plots = [flyer.sample()]
    plot_opts = [{"marker":"o"}]
    
    bbox = mpt.square_plot_bounds(plots)
    
    plots.append(mpt.sample_circle(flyer.R, 100))
    plot_opts.append(None)
    plots.append(mpt.sample_circle(flyer.R+70, 100))
    plot_opts.append(None)
    
    mpt.square_plots(ax, plots, bbox, plot_opts)
    
    plt.show()

    flyer.dumpTraj()
