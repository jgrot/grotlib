import sys

import compare as cmp
import json
import ksp
import math
import os
import pickle


def compareSolutions( soln, soln_ref ) :
    
    for irow, row in enumerate( soln ) :
        try :
            row_ref = soln_ref[ irow ]
        except :
            row_ref = None
                
        if row_ref is None :
            print("Ran out of reference data")
            exit(1)
                
        for jelem, elem in enumerate( row ) :
            elem_ref = row_ref[ jelem ]
            if not cmp.fsame(elem, elem_ref, report_to = sys.stderr, report = cmp.CMP_ONLY_NOT_SAME) : exit(1)


def T2DS1ME( ) :
    # KSP Flight:
    #
    # Max alt: 138247
    # Crash time: 6:52 (412 s)
    #
    # Model:
    #
    # Drag   Max Alt    Crash
    # 0.130  112690     365
    # 0.100  125243     390
    # 0.090  129713     395
    # 0.080  134325     405
    # 0.070  139099     415
    # 0.072  138121     410
    
    print("Testing traj2d with one stage with mixed engines (T2DS1ME)")
    
    tname = "T2DS1ME"
    tfile = tname + ".pkl"    

    do_generate = not os.access(tfile, os.F_OK)

    stage_1 = ksp.Stage()
    stage_1.loadJSON(tname + "_stage1.json")
    stage_1.dumpInfo( )

    def fthrottle( t, y ) :
        return 1.0
    def falpha( t, y, flyer ) :
        m, r, th, vr, om = y
        vth = r*om
        if r < (flyer.R+7.5E3) :
            return 0.5 * math.pi
        if vth == 0.0 :
            a = 0.5 * math.pi
        else :
            a = math.atan( vr / vth )
        return a
    
    flyer = ksp.FlyingStage( stage_1, "Stage 1", "Kerbin", fthrottle, falpha )
    flyer.launch( )
    # Fly until crash
    flyer.flyTo( 30000 )

    if do_generate :
        with open( tfile, 'wb' ) as f :
            pickle.dump( flyer.soln, f )
    else :
        with open( tfile, 'rb' ) as f :
            soln_compare = pickle.load( f )
            compareSolutions( flyer.soln, soln_compare )      
        print("SUCCESS")

    flyer.plot( )
    flyer.dumpTraj( )    


def T2DS2ME( ) :
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
    tfile = tname + ".pkl"    

    do_generate = not os.access(tfile, os.F_OK)

    stage_1 = ksp.Stage()
    stage_1.loadJSON(tname + "_stage1.json")
    stage_1.dumpInfo( )

    stage_2 = ksp.Stage()
    stage_2.loadJSON(tname + "_stage2.json")
    stage_2.dumpInfo( )
    
    def fthrottle( t, y ) :
        return 1.0
    def falpha( t, y, flyer ) :
        m, r, th, vr, om = y
        vth = r*om
        if r < (flyer.R+30.0E3) :
            return 0.5 * math.pi
        else :
            a = math.atan( vr / vth )
        return a
    
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
        craft_asl_dvremain = flyer.dvRemain( m, ksp.C_p0 )
        craft_dvremain = flyer.dvRemain( m, flyer.fpress.call(r-flyer.R)[0] )
        stage_asl_dvremain = flyer.stage.dvRemain( m, ksp.C_p0 )
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
            compareSolutions( fly_s1.soln, soln_compare["stage1"] )
            compareSolutions( fly_s2.soln, soln_compare["stage2"] )
            compareSolutions( DV, soln_compare["DV"] )
        print("SUCCESS")
            
    fly_s2.plot( t0 = 0.0, dt = 10.0 )
    fly_s2.dumpTraj( t0 = 0, dt = 10.0 )

    
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
    truth = 11400000.0
    test = ksp.dInterpDist("D('Mun','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Kerbin and Mun (flip order)")
    truth = 11400000.0
    test = ksp.dInterpDist("D('Mun','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Minmus and Kerbin")
    truth = 46400000.0
    test = ksp.dInterpDist("D('Minmus','Kerbin')")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dInterpDist Kerbin and Minmus (flip order)")
    truth = 46400000.0
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
    truth = 46510000.0
    test = ksp.dInterp(s)
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dvInterp with dvmap_kerbin_to_mun.json")
    with open("dvmap_kerbin_to_mun.json", "rt") as f :
        maneuvers = json.load(f)
    truth = 3671.3933492693336
    test = ksp.dvInterp( maneuvers )
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
    
    print("Testing jump")
    truth = 2.4543073334711605
    test = ksp.jump((44,"km"), "Kerbin", 175, (3,"t"), (200,"kN"), "mfuel")
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    T2DS1ME( )
    T2DS2ME( )
