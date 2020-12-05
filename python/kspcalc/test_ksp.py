import sys

import compare as cmp
import json
import ksp

if __name__ == "__main__" :

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
    truth = 2646.501134322117 
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
    truth = 6215.674952948572
    test = ksp.dvInterp( maneuvers )
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
    
    print("Testing mFuelToReachAlt")
    truth = 726.6922810772608
    # Note: 0.0075 t per unit of solid fuel
    test = ksp.mFuelToReachAlt((44,"km"), "Kerbin", 175, (250E3,"N"), (6.3,"t"))/0.0075
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
