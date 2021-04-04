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

    s1 = ksp.Stage()
    s1.loadJSON("minmus_express_2_s1.json")
    s1.dumpInfo()

    s2 = ksp.Stage()
    s2.loadJSON("minmus_express_2_s2.json")
    s2.dumpInfo()

    s3 = ksp.Stage()
    s3.loadJSON("minmus_express_2_s3.json")
    s3.dumpInfo()

    t1 = 300.0
    t2 = 400.0
    t3 = 10000.0
    
    def fthrottle_s1( t, y ) :
        return 1.0
    def falpha_s1( t, y, flyer ) :
        m, r, th, vr, om = y
        # vth = om/r
        #if r < (flyer.R+7.5E3) :
        n = (1.0, 0.0)
        #else :
        #    v,n = mpm.v_and_dir(y)
        return n
    
    flyer_s1 = ksp.FlyingStage(s1, "Stage 1", "Kerbin", fthrottle_s1, falpha_s1)
    flyer_s1.launch()
    flyer_s1.flyTo(t1)

    def fthrottle_s2( t, y ) :
        return 1.0
    def falpha_s2( t, y, flyer ) :
        m, r, th, vr, om = y
        # vth = om/r
        v, n = mpm.v_and_dir(y)        

        return n
    
    flyer_s2 = ksp.FlyingStage(s2, "Stage 2", "Kerbin", fthrottle_s2, falpha_s2)
    flyer_s2.launch(sm1=flyer_s1, t0=t1)
    flyer_s2.flyTo(t2)

    def fthrottle_s3( t, y ) :
        if t < t2 + 60.0 :
            return 1.0
        else :
            return 0.0
    def falpha_s3( t, y, flyer ) :
        m, r, th, vr, om = y
        # vth = om/r
        v, n = mpm.v_and_dir(y)        

        return n
    
    flyer_s3 = ksp.FlyingStage(s3, "Stage 3", "Kerbin", fthrottle_s3, falpha_s3)
    flyer_s3.launch(sm1=flyer_s2, t0=t2)
    flyer_s3.flyTo(t3)
    
    import matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    
    plots = [flyer_s3.sample(t0 = 0.0)]
    plot_opts = [{"marker":"o"}]
    
    bbox = mpt.square_plot_bounds(plots)
    
    plots.append(mpt.sample_circle(flyer_s3.R, 100))
    plot_opts.append(None)
    plots.append(mpt.sample_circle(flyer_s3.R+70, 100))
    plot_opts.append(None)
    
    mpt.square_plots(ax, plots, bbox, plot_opts)
    
    plt.show()

    flyer_s3.dumpTraj(t0 = 0.0)
