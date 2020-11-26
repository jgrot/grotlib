import sys

import compare as cmp
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

    print("Testing dv_orbit at at Kerbin at 100 km")
    truth = 2646.501134322117
    test = ksp.dv_orbit("Kerbin", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing orbit_v at Kerbin at 100 km")
    truth = 2245.6306781964713
    test = ksp.orbit_v("Kerbin", (100,"km"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dv_hohmann_peri at Kerbin r_peri=100km, r_apo=30Mm")
    truth = 894.457725062739
    test = ksp.dv_hohmann_peri("Kerbin", (100,"km"), (30,"Mm"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)

    print("Testing dv_hohmann_apo at Kerbin r_peri=100km, r_apo=30Mm")
    truth = 267.8140180540594
    test = ksp.dv_hohmann_apo("Kerbin", (100,"km"), (30,"Mm"))
    if not cmp.fsame(test, truth, report_to = sys.stderr) : exit(1)
