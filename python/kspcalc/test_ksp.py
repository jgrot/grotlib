#!/usr/bin/env python
#
# Copyright © 2021 Jonathan Grot
#
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
import sys

import compare as cmp
import ksp
import mks_polar_motion as mpm
import mpl_tools as mpt

import json
import math
import os
import pickle

def T2DS1ME() :
    
    print("Testing traj2d with one stage with mixed engines (T2DS1ME)")
    
    tname = "T2DS1ME"
    tfile = ksp.pthut(tname + ".pkl")
    do_generate = not os.access(tfile, os.F_OK)

    stage_1 = ksp.Stage()
    stage_1.loadJSON(ksp.pthut(tname + "_stage1.json"))
    # stage_1.dumpInfo( )

    def fthrottle( t, y ) :
        return 1.0
    def falpha( t, y, flyer ) :
        m, r, th, vr, om = y
        vth = om/r
        if r < (flyer.R+7.5E3) :
            n = (1.0, 0.0)
        elif vth == 0.0 :
            n = (1.0, 0.0)
        else :
            v,n = mpm.v_and_dir(y)
        return n
    
    flyer = ksp.FlyingStage(stage_1, "Stage 1", "Kerbin", fthrottle, falpha)
    flyer.launch( )
    # Fly until crash
    flyer.flyTo( 30000 )

    if do_generate :
        with open( tfile, 'wb' ) as f :
            pickle.dump( flyer.soln, f )
    else :
        with open( tfile, 'rb' ) as f :
            soln_compare = pickle.load( f )
            cmp.compare_datasets( flyer.soln, soln_compare )      
        print("SUCCESS")

    if True :
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

        flyer.dumpTraj( )    


def T2DS2ME() :
    # KSP Flight
    #
    # 50s stage
    # 30km track orbit
    #
    # Max alt:
    # Crash Time:
    #
    print("Testing traj2d with two stages with mixed engines (T2DS2ME)")
    
    tname = "T2DS2ME"
    tfile = ksp.pthut(tname + ".pkl")
    do_generate = not os.access(tfile, os.F_OK)

    stage_1 = ksp.Stage()
    stage_1.loadJSON(ksp.pthut(tname + "_stage1.json"))
    # stage_1.dumpInfo( )

    stage_2 = ksp.Stage()
    stage_2.loadJSON(ksp.pthut(tname + "_stage2.json"))
    # stage_2.dumpInfo( )
    
    def fthrottle( t, y ) :
        return 1.0
    def falpha( t, y, flyer ) :
        r = y[1]
        if r < (flyer.R+30.0E3) :
            n = (1.0, 0.0)
        else :
            v, n = mpm.v_and_dir(y)
        return n
    
    fly_s1 = ksp.FlyingStage( stage_1, "Stage 1", "Kerbin", fthrottle, falpha )
    fly_s1.launch( )
    fly_s1.flyTo( 84.00 )

    fly_s2 = ksp.FlyingStage( stage_2, "Stage 2",  "Kerbin", fthrottle, falpha )
    fly_s2.launch( sm1 = fly_s1, t0 = 84.00 )
    fly_s2.flyTo( 30000 )

    t = 0.0
    DV = []
    while t <= fly_s2.solnt[-1] :
        Y, crashed, flyer = fly_s2.flyTo( t )
        m, r, th, vr, om = Y
        craft_asl_dvremain = flyer.dv_at_m( m, ksp.C_p0 )
        craft_dvremain = flyer.dv_at_m( m, flyer.fpress.call(r-flyer.R)[0] )
        stage_asl_dvremain = flyer.stage.dv_at_m( m, ksp.C_p0 )
        DV.append( (craft_asl_dvremain, craft_dvremain, stage_asl_dvremain) )
        t += 10
    
    if do_generate :
        with open( tfile, 'wb' ) as f :
            data = {
                "stage1" : fly_s1.soln,
                "stage2" : fly_s2.soln,
                "DV"     : DV
            }
            pickle.dump( data, f )
    else :
        with open( tfile, 'rb' ) as f :
            soln_compare = pickle.load( f )
            cmp.compare_datasets( fly_s1.soln, soln_compare["stage1"] )
            cmp.compare_datasets( fly_s2.soln, soln_compare["stage2"] )
            cmp.compare_datasets( DV, soln_compare["DV"] )
        print("SUCCESS")

    if True :
        import matplotlib
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        
        plots = [fly_s2.sample(t0=0.0, dt=10.0)]
        plot_opts = [{"marker":"o"}]
        
        bbox = mpt.square_plot_bounds(plots)
        
        plots.append(mpt.sample_circle(fly_s2.R, 100))
        plot_opts.append(None)
        plots.append(mpt.sample_circle(fly_s2.R+70, 100))
        plot_opts.append(None)

        mpt.square_plots(ax, plots, bbox, plot_opts)
        
        plt.show()

        flyer.dumpTraj()


def TORBIT() :
    print("Testing Orbit object (TORBIT)")
    
    tname = "TORBIT"
    tfile = ksp.pthut(tname + ".pkl")
    do_generate = not os.access(tfile, os.F_OK)

    GM = ksp.bodies_db["Kerbin"]["GM"][0]
    r = 670E3
    
    # Test circular orbit
    print("Checking forced circular orbit (v_r ignored)")
    o = mpm.Orbit( [ 2000, r, 0, 100.0, 0.01 ], GM, force="circle" )
    # Eccentricity must be identically zero, so tolfrac is zero
    if not cmp.fsame(o.e, 0.0, tolfrac=0.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.k, 706E13, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.h, 1537888162.383728, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.E, -5268656716.417911, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.v0, 2295.3554662443707, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.r1, r, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    tau_circ = 2.0*math.pi/math.sqrt(GM/math.pow(r,3.0))
    if not cmp.fsame(o.tau, tau_circ, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    t1 = 0.25*o.tau
    t2 = 0.5*o.tau
    t3 = 0.75*o.tau
    t4 = -1.25*o.tau

    phi1 = o.phi_at_t(t1)
    phi2 = o.phi_at_t(t2)
    phi3 = o.phi_at_t(t3)
    phi4 = o.phi_at_t(t4)
    
    print("  ...checking phi at 1/4 tau")
    if not cmp.fsame(phi1, 0.5*math.pi, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 1/2 tau")
    if not cmp.fsame(phi2, math.pi, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 3/4 tau")
    if not cmp.fsame(phi3, 1.5*math.pi, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at -1.25 tau")
    if not cmp.fsame(phi4, -1.25*2.0*math.pi, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    phi1 = -1.25*2.0*math.pi
    t1_truth = -1.25*o.tau
    t1 = o.t_at_phi(phi1)
    print("  ...checking t at -1.25*2pi phi")
    if not cmp.fsame(t1, t1_truth, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    # Test ellipse
    print("Checking elliptical orbit")
    o = mpm.Orbit( [ 2000, r, 0, 0.0, 0.004 ], GM )
    if not cmp.fsame(o.e, 0.3632317280453257, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.k, 706E13, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.h, 1795600000.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.E, -3354913432.835821, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.v0, 2680.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.r1, 1434376.2056275664, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.tau, 3609.379904084419, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    t1 = 0.25*o.tau
    t2 = 0.5*o.tau
    t3 = 0.75*o.tau
    t4 = 4.0*o.tau
    t5 = -4.0*o.tau

    phi1 = o.phi_at_t(t1)
    phi2 = o.phi_at_t(t2)
    phi3 = o.phi_at_t(t3)
    phi4 = o.phi_at_t(t4)
    phi5 = o.phi_at_t(t5)

    phi1_truth = 2.2431742381218593
    phi3_truth = 2.0*math.pi - phi1_truth
    phi4_truth = 8.0*math.pi
    phi5_truth = -8.0*math.pi
    print("  ...checking phi at 1/4 tau")
    if not cmp.fsame(phi1, phi1_truth, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 1/2 tau")
    if not cmp.fsame(phi2, math.pi, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 3/4 tau")
    if not cmp.fsame(phi3, phi3_truth, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 4 tau")
    if not cmp.fsame(phi4, phi4_truth, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at -4 tau")
    if not cmp.fsame(phi5, phi5_truth, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    

    # Test parabola
    print("Checking forced parabolic orbit (v_r ignored)")
    o = mpm.Orbit( [ 2000, r, 0, 100.0, 0.01 ], GM, force="parabola" )
    # e must be identically 1.0
    if not cmp.fsame(o.e, 1.0, tolfrac=0.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.k, 706E13, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.h, 2174902296.656105, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    # E must be identically 0.0
    if not cmp.fsame(o.E, 0.0, tolfrac=0.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.v0, 3246.1228308300074, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.phi_i, 0.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    t1 = 10.0
    t2 = 100.0
    t3 = 1000.0

    phi1 = o.phi_at_t(t1)
    phi2 = o.phi_at_t(t2)
    phi3 = o.phi_at_t(t3)

    print("  ...checking phi at 10 seconds")
    if not cmp.fsame(phi1, 0.0484306562739271, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 100 seconds")
    if not cmp.fsame(phi2, 0.46686597787401796, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 1000 seconds")
    if not cmp.fsame(phi3, 1.9248853641334023, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    # Test hyperbola (Note, positive vr entry point)
    print("Checking hyperbolic orbit with insertion not at phi=0 (v_r > 0)")
    o = mpm.Orbit( [ 2000, r, 0.0, 100.0, 0.01 ], GM)
    if not cmp.fsame(o.e, 7.521273426540002, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.k, 706E13, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.h, 4489000000.0, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.E, 34362686567.16418, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.v0, 6700.845443458724, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    if not cmp.fsame(o.phi_i, 0.016908466324052317, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    t1 = 10.0
    t2 = 100.0
    t3 = 1000.0

    phi1 = o.phi_at_t(t1)
    phi2 = o.phi_at_t(t2)
    phi3 = o.phi_at_t(t3)

    print("  ...checking phi at 10 seconds")
    if not cmp.fsame(phi1, 0.09973240289205285, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 100 seconds")
    if not cmp.fsame(phi2, 0.8035063096303883, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)
    print("  ...checking phi at 1000 seconds")
    if not cmp.fsame(phi3, 1.5791723428116158, tolfrac=1E-5, report_to=sys.stderr, report=cmp.CMP_ONLY_NOT_SAME) : exit(1)

    print("SUCCESS")

def TSOIISECT() :
    print("Testing SOI intersections")

    ins = [
        # m, r, th, vr, om, c1, c2, c3
        [2000.0, 670E3, 0.0, 0.0, 0.0048, 1.0, 40.0, 6.162],
        [2000.0, 670E3, 0.0, 0.0, 0.0048, 1.0, 40.0, 6.17],
        [2000.0, 670E3, 0.0, 0.0, 0.0048, 0.95, 40.0, 6.17],
        [2000.0, 670E3, 0.0, 0.0, 0.005, 0.6, 5.0, 2.0],
    ]

    outs = [
        [2.9962290251689883, 3.2869562820106277],
        [2.992304492616973, 2.9978778327411484, 3.285307474421102, 3.290880814566911],
        [2.934401028222881, 3.0310486370486043],
        [1.832548507281095, 2.1857173338290483],
    ]

    for i, w in enumerate(ins) :
        x = outs[i]

        m, r, th, vr, om, c1, c2, c3 = w

        y = [m, r, th, vr, om]
        o = mpm.Orbit(y, ksp.bodies_db["Kerbin"]["GM"][0] )
        
        phi_soi = c1*math.pi
        d_soi = c2*o.r0
        r_soi = c3*o.r0
        
        phi_x = o.intersect_soi(phi_soi, d_soi, r_soi)

        cmp.compare_datasets([phi_x], [x])

    print("SUCCESS")

    
if __name__ == "__main__" :

    ksp.augmentBodyDbs( )
    ksp.processBodyDbs( )
    
    print("Testing g at Kerbin at 0 m")
    truth = 9.805555555555555
    test = ksp.g("Kerbin", (0.0,"m"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing g at Kerbin at 100 km")
    truth = 7.204081632653061
    test = ksp.g("Kerbin", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing g at Kerbin at 100E3 m (unit conversion test)")
    truth = 7.204081632653061
    test = ksp.g("Kerbin", (100E3,"m"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing g at Mun at 0 m")
    truth = 1.6285
    test = ksp.g("Mun", (0.0,"m"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing g at Mun at 100 km")
    truth = 0.7237777777777777
    test = ksp.g("Mun", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
    
    print("Testing g at Minmus at 0 m")
    truth = 0.4905555555555556
    test = ksp.g("Minmus", (0.0,"m"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing g at Minmus at 100 km")
    truth = 0.068984375
    test = ksp.g("Minmus", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
    
    print("Testing dvOrbit at at Kerbin at 100 km")
    truth = 916.7748853869145
    test = ksp.dvOrbit("Kerbin", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing orbitV at Kerbin at 100 km")
    truth = 2245.6306781964713
    test = ksp.orbitV("Kerbin", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dvHohmannPeri at Kerbin r_peri=100km, r_apo=30Mm")
    truth = 894.457725062739
    test = ksp.dvHohmannPeri("Kerbin", (100,"km"), (30,"Mm"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dvHohmannApo at Kerbin r_peri=100km, r_apo=30Mm")
    truth = 267.8140180540594
    test = ksp.dvHohmannApo("Kerbin", (100,"km"), (30,"Mm"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Mun and Kerbin")
    truth = 12000000.0
    test = ksp.dInterpDist("D('Mun','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Kerbin and Mun (flip order)")
    truth = 12000000.0
    test = ksp.dInterpDist("D('Mun','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Minmus and Kerbin")
    truth = 47000000.0
    test = ksp.dInterpDist("D('Minmus','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Kerbin and Minmus (flip order)")
    truth = 47000000.0
    test = ksp.dInterpDist("D('Minmus','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpRadius Kerbin")
    truth = 600000.0
    test = ksp.dInterpRadius("R('Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpRadius Mun")
    truth = 200000.0
    test = ksp.dInterpRadius("R('Mun')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpRadius Minmus")
    truth = 60000.0
    test = ksp.dInterpRadius("R('Minmus')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterp on [100,\"pc\"]")
    truth = 3.0857e+18
    test = ksp.dInterp([100,"pc"])
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
    
    print("Testing dInterp on (100,\"pc\")")
    truth = 3.0857e+18
    test = ksp.dInterp((100,"pc"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    s = "(0.5,'pc')"
    print("Testing dInterp on \"%s\"" % s)
    truth = 1.54285e+16
    test = ksp.dInterp(s)
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    s = "[0.5,'pc']"
    print("Testing dInterp on \"%s\"" % s)
    truth = 1.54285e+16
    test = ksp.dInterp(s)
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    s = "D('Kerbin','Minmus') + (100,'km') + R('Minmus') - (50000,'m')"
    print("Testing dInterp on \"%s\"" % s)
    truth = 47110000.0
    test = ksp.dInterp(s)
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dvInterp with dvmap_kerbin_to_mun.json")
    with open(ksp.pthex("dvmap_kerbin_to_mun.json"), "rt") as f :
        maneuvers = json.load(f)
    truth = 3674.084081820914
    test = ksp.dvInterp( maneuvers )
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
    
    print("Testing jump")
    truth = 2.4543073334711605
    test = ksp.jump((44,"km"), "Kerbin", 175, (3,"t"), (200,"kN"), "mfuel")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    T2DS1ME()
    T2DS2ME()
    TORBIT()
    TSOIISECT()
