#! /usr/bin/env python

# Copyright © 2021 Jonathan Grot

'''DV to go from kerbin circ orbit to mun circ orbit
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
    dv = []
    dv_desc = []

    # Stage we will be flying
    s2 = ksp.Stage()
    s2.loadJSON("stage_2.json")
    s2.dumpInfo()

    rw.print("Analysis...")

    # Available stage DV
    dv_stage = s2.dv_at_m(s2.m0_kg, 0.0)
    dv.append(dv_stage)
    dv_desc.append("Initial")
    rw.write("Initial stage DV", dv_stage)

    # Kerbin
    GMkerbin, GMkerbinu = ksp.bodies_db["Kerbin"]["GM"]
    rw.write("GM Kerbin", GMkerbin)
    # Plot Kerbin
    plots.append(mpt.sample_circle(6E5, 100))
    plot_opts.append(None)

    # Mun
    mun_th = math.pi
    mun_d = ksp.dInterp("D('Kerbin','Mun')")
    mun_r = ksp.dInterp("R('Mun')")
    mun_XY = mpt.sample_circle(mun_r, 100)
    mun_x, mun_y = mpm.xy_at_rth(mun_d, mun_th)
    mun_XY = mpt.offset_xy(mun_XY, mun_x, mun_y)
    # Plot Mun
    plots.append(mun_XY)
    plot_opts.append(None)
    # Mun SOI
    mun_soi_r, mun_soi_ru = ksp.bodies_db["Mun"]["SOI"]
    mun_soi_r *= ksp.uconv(ksp.dist_db, mun_soi_ru, "m")
    rw.write("Mun SOI", mun_soi_r)
    mun_soi_XY = mpt.sample_circle(mun_soi_r, 100)
    mun_soi_XY = mpt.offset_xy(mun_soi_XY, mun_x, mun_y)
    plots.append(mun_soi_XY)
    plot_opts.append(None)
    
    rw.write("Distance from Kerbin to Mun center (Mm)", 1E-6*mun_d)

    # Initial (circular) orbit
    o1_r0 = ksp.dInterp("R('Kerbin') + (800,'km')")
    o1_v0 = mpm.v_circular_orbit(GMkerbin, o1_r0)
    rw.write("Circular orbit radius", o1_r0)
    rw.write("Circular orbit speed", o1_v0)

    # Hohmann transfer out to Mun
    dv_o1_o2 = mpm.dv_r0_hohmann(GMkerbin, o1_r0, 0.85*mun_d)
    dv.append(dv[-1] - dv_o1_o2)
    dv_desc.append("After Hohmann to Mun")
    rw.write("DV to reach Mun", dv_o1_o2)

    # New stage mass.
    s2_m_o2 = s2.m_at_dv(dv_stage - dv_o1_o2, 0.0)
    rw.write("M after DV", s2_m_o2)

    # Stage orbit after Hohmann transfer
    o2_th0 = 10.0*math.pi/180.0 # Inserting orbiter at r0
    o2_y = [s2_m_o2, o1_r0, o2_th0, 0.0, (o1_v0 + dv_o1_o2)/o1_r0]
    o2 = mpm.OrientedOrbit(o2_y, GMkerbin)
    o2_TH=[i*math.pi*2/100 for i in range(101)]
    o2_R=o2.sample_th(o2_TH)
    # Plot stage orbit
    plots.append(mpt.polar_to_xy(o2_R, o2_TH))
    plot_opts.append(None)

    # Intersect elliptical orbit with Mun SOI
    TH_o2_x_mun_soi = o2.intersect_soi(mun_th, mun_d, mun_soi_r)
    if len(TH_o2_x_mun_soi) > 0 :
        Y_o2_x_mun_soi = o2.y_th(TH_o2_x_mun_soi[0])
        m_o2_x_mun_soi, r_o2_x_mun_soi, th_o2_x_mun_soi, vr_o2_x_mun_soi, om_o2_x_mun_soi = Y_o2_x_mun_soi
        vth_o2_x_mun_soi = r_o2_x_mun_soi*om_o2_x_mun_soi
        # Intersection point
        x_o2_x_mun_soi, y_o2_x_mun_soi = mpm.xy_at_rth(r_o2_x_mun_soi, th_o2_x_mun_soi)
        # Speed of stage in grand sys at SOI intersection
        s2_vx_o2_x_mun_soi, s2_vy_o2_x_mun_soi = mpm.vxy_at_th(vr_o2_x_mun_soi, vth_o2_x_mun_soi, th_o2_x_mun_soi)
        # Plot intersection point
        plots.append([[x_o2_x_mun_soi],[y_o2_x_mun_soi]])
        plot_opts.append({"marker":"o"})

        # Set up Mun reference frame
        mun_v, mun_vu = ksp.bodies_db["Mun"]["VORB"]
        rw.write("Mun orbit speed (m/s)", mun_v)
        mun_frm = mpm.Frame2D()
        mun_frm.set_rth(mun_d, mun_th)
        mun_frm.set_vrth(0.0, mun_v)
        rw.write("Mun orbit speed vector", repr((mun_frm.vx, mun_frm.vy)))
        rw.write("Intersection speed vector", repr((s2_vx_o2_x_mun_soi, s2_vy_o2_x_mun_soi)))
        # Intersection velocity in the Mun reference frame
        s2_vx_mun, s2_vy_mun = mun_frm.xform_vel_to(s2_vx_o2_x_mun_soi, s2_vy_o2_x_mun_soi)
        rw.write("Intersection velocity in Mun frame", repr((s2_vx_mun, s2_vy_mun)))
        # Intersection point in Mun frame
        s2_r_mun, s2_th_mun = mun_frm.xform_pos_to('c', x_o2_x_mun_soi, y_o2_x_mun_soi, 'p')
        rw.write("Intersection point in Mun frame", repr((s2_r_mun, s2_th_mun)))
        # Intersection velocity in polar coordinates
        s2_vr_mun, s2_vth_mun = mpm.vrth_at_th(s2_vx_mun, s2_vy_mun, s2_th_mun)
        rw.write("Intersection Vr, Vth in Mun frame", repr((s2_vr_mun, s2_vth_mun)))

        # Stage orbit in Mun reference frame
        GM_mun, GM_munu = ksp.bodies_db["Mun"]["GM"]
        Y_mun = [s2_m_o2, s2_r_mun, s2_th_mun, s2_vr_mun, s2_vth_mun/s2_r_mun]
        o3 = mpm.OrientedOrbit(Y_mun, GM_mun)
        rw.write("New orbit in Mun frame is", o3.ORB_ORBITS[o3.classify()])
        rw.write("New orbit in Mun frame r0 is", o3.r0)
        
        # Plot stage orbit in Mun reference frame
        o3_phi_interval = 0.8*(o3.phi_max - o3.phi_i)
        rw.write("o3_phi_interval", o3_phi_interval)
        o3_dphi = o3_phi_interval/100.0
        o3_TH = [s2_th_mun + o3_dphi*i for i in range(101)]
        o3_R = o3.sample_th(o3_TH)
        o3_plt = mpt.polar_to_xy(o3_R,o3_TH)
        o3_plt2 = mpt.offset_xy(o3_plt, mun_frm.x, mun_frm.y)
        plots.append(o3_plt2)
        plot_opts.append({"marker":"x"})

        if o3.phi_i < 0.0 :
            # Before r0
            dv_o3_o4 = o3.r0_dv_to_e(0.0)
            dv.append(dv[-1]-abs(dv_o3_o4))
            dv_desc.append("After Hohmann to circ around mun")
            rw.write("DV to circ orbit about mun is", dv_o3_o4)

            s2_v0_o4 = o3.v0 + dv_o3_o4
            rw.write("V for circ orbit", s2_v0_o4)

            s2_m_o4 = s2.m_at_dv(dv_stage - dv_o1_o2 + dv_o3_o4, 0.0)
            
            rw.write("M after DV circ", s2_m_o4)
            o4_y = [s2_m_o4, o3.r0, o3.th_at_phi(0.0), 0.0, s2_v0_o4/o3.r0]
            o4 = mpm.OrientedOrbit(o4_y, GM_mun)
            rw.write("Circular orbit eccentricity", o4.e)

            o4_TH=[i*2.0*math.pi/100.0 for i in range(101)]
            o4_R=o4.sample_th(o4_TH)
            o4_pltc = mpt.polar_to_xy(o4_R, o4_TH)
            o4_plt2c = mpt.offset_xy(o4_pltc, mun_frm.x, mun_frm.y)
            plots.append(o4_plt2c)
            plot_opts.append({"marker":"s"})
            
        else :
            rw.write("DV to circ orbit about mun is", "NONE: PAST R0!")

        dvdata = zip(dv_desc, dv)
        rw.tabulate(["Desc","DV"], dvdata)

    # Plot everything
    fig, ax = plt.subplots()
    bbox = mpt.square_plot_bounds(plots)
    mpt.square_plots(ax, plots, bbox, plot_opts)
    
    plt.show()
