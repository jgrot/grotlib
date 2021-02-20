#! /usr/bin/env python

# Copyright © 2021 Jonathan Grot

'''DV to go from kerbin circ orbit to mun landing
'''

import math

import ksp
import cli_tools as ct
import mks_polar_motion as mpm
import mpl_tools as mpt

if __name__ == "__main__" :

    ksp.augmentBodyDbs()
    ksp.processBodyDbs()

    plots = []
    plot_opts = []
    
    import matplotlib
    import matplotlib.pyplot as plt
    
    rw = ct.RowWriter([40])
    
    s2 = ksp.Stage()
    s2.loadJSON("stage_2.json")
    s2.dumpInfo()

    rw.print("Analysis...")

    # Available stage DV
    dv = s2.dv_at_m(s2.m0_kg, 0.0)
    rw.write("DV remaining in vacuum",dv)

    GM, GMu = ksp.bodies_db["Kerbin"]["GM"]
    rw.write("GM Kerbin", GM)

    # Plot Kerbin
    plots.append(mpt.sample_circle(6E5, 100))
    plot_opts.append(None)

    # Mun
    mun_th = math.pi
    d_mun = ksp.dInterp("D('Kerbin','Mun')")
    r_mun = ksp.dInterp("R('Mun')")
    XYmun = mpt.sample_circle(r_mun, 100)
    x_mun, y_mun = mpm.xy_at_rth(d_mun, mun_th)
    XYmun = mpt.offset_xy(XYmun, x_mun, y_mun)
    # Plot Mun
    plots.append(XYmun)
    plot_opts.append(None)
    # Mun SOI
    rsoi_mun, rsoi_munu = ksp.bodies_db["Mun"]["SOI"]
    rsoi_mun *= ksp.uconv(ksp.dist_db, rsoi_munu, "m")
    rw.write("Mun SOI", rsoi_mun)
    XYmun_soi = mpt.sample_circle(rsoi_mun, 100)
    XYmun_soi = mpt.offset_xy(XYmun_soi, x_mun, y_mun)
    plots.append(XYmun_soi)
    plot_opts.append(None)
    
    rw.write("Distance from Kerbin to Mun center (Mm)", 1E-6*d_mun)

    # Initial (circular) orbit
    r0 = ksp.dInterp("R('Kerbin') + (800,'km')")
    v0 = mpm.v_circular_orbit(GM, r0)
    rw.write("Circular orbit radius", r0)
    rw.write("Circular orbit speed", v0)

    # Hohmann transfer out to Mun
    dv1 = mpm.dv_r0_hohmann(GM, r0, 0.85*d_mun)
    rw.write("DV to reach Mun", dv1)

    # New stage mass.
    m = s2.m_at_dv(dv - dv1, 0.0)
    rw.write("M after DV", m)

    # Stage orbit after Hohmann transfer
    th_i = 10.0*math.pi/180.0
    y = [m, r0, th_i, 0.0, (v0+dv1)/r0]
    o = mpm.OrientedOrbit(y, GM)
    TH=[i*math.pi*2/100 for i in range(101)]
    R=o.sample_th(TH)
    # Plot stage orbit
    plots.append(mpt.polar_to_xy(R,TH))
    plot_opts.append(None)

    # Intersect elliptical orbit with Mun SOI
    xth = o.intersect_soi(mun_th, d_mun, rsoi_mun)
    if len(xth) > 0 :
        Yx = o.y_th(xth[0])
        mx, rx, thx, vrx, omx = Yx
        vthx = rx*omx
        # Intersection point
        x_x, y_x = mpm.xy_at_rth(Yx[1], xth[0])
        # Speed of stage in grand sys at SOI intersection
        s2v_x, s2v_y = mpm.vxy_at_th(vrx, vthx, thx)
        # Plot intersection point
        plots.append([[x_x],[y_x]])
        plot_opts.append({"marker":"o"})

        # Set up Mun reference frame
        mun_v, mun_vu = ksp.bodies_db["Mun"]["VORB"]
        rw.write("Mun orbit speed (m/s)", mun_v)
        munfrm = mpm.Frame2D()
        munfrm.set_rth(d_mun, mun_th)
        munfrm.set_vrth(0.0, mun_v)
        rw.write("Mun orbit speed vector", repr((munfrm.vx, munfrm.vy)))
        rw.write("Intersection speed vector", repr((s2v_x, s2v_y)))
        # Intersection velocity in the Mun reference frame
        vx_mun, vy_mun = munfrm.xform_vel_to(s2v_x, s2v_y)
        rw.write("Intersection velocity in Mun frame", repr((vx_mun, vy_mun)))
        # Intersection point in Mun frame
        r_mun, th_mun = munfrm.xform_pos_to('c', x_x, y_x, 'p')
        rw.write("Intersection point in Mun frame", repr((r_mun, th_mun)))
        # Intersection velocity in polar coordinates
        vr_mun, vth_mun = mpm.vrth_at_th(vx_mun, vy_mun, th_mun)
        rw.write("Intersection Vr, Vth in Mun frame", repr((vr_mun, vth_mun)))

        # Stage orbit in Mun reference frame
        GM_mun, GM_munu = ksp.bodies_db["Mun"]["GM"]
        y_mun = [m, r_mun, th_mun, vr_mun, vth_mun/r_mun]
        o_mun = mpm.OrientedOrbit(y_mun, GM_mun)
        rw.write("New orbit in Mun frame is", o_mun.ORB_ORBITS[o_mun.classify()])
        rw.write("New orbit in Mun frame r0 is", o_mun.r0)
        
        # Plot stage orbit in Mun reference frame
        phi_interval = 0.8*(o_mun.phi_max - o_mun.phi_i)
        rw.write("phi_interval", phi_interval)
        dphi = phi_interval/100.0
        TH = [th_mun + dphi*i for i in range(101)]
        R = o_mun.sample_th(TH)
        o_mun_plt = mpt.polar_to_xy(R,TH)
        o_mun_plt2 = mpt.offset_xy(o_mun_plt, munfrm.x, munfrm.y)
        plots.append(o_mun_plt2)
        plot_opts.append({"marker":"x"})

        if o_mun.phi_i < 0.0 :
            # Before r0
            dv_circ = o_mun.r0_dv_to_e(0.0)
            rw.write("DV to circ orbit about mun is", dv_circ)
        else :
            rw.write("DV to circ orbit about mun is", "NONE: PAST R0!")
        
    # Plot everything
    fig, ax = plt.subplots()
    bbox = mpt.square_plot_bounds(plots)
    mpt.square_plots(ax, plots, bbox, plot_opts)
    
    plt.show()
