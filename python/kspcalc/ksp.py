#! /usr/bin/env python

# UNITS
#
# For now all functions should output MKS units
#

import argparse
import json
import math
import os
import sys

# Verify Python version
if ( sys.version_info.major < 3 and sys.version_info.minor < 8 ) :
    print("PYTHON VERSION ", sys.version)
    sys.stderr.write("*** Wrong python version.  Make sure virtual environment is activated.\n")
    exit(1)

bodies_db = {
    "Kerbin" : {
        "GM" : (3.53E12, "m3/s2"),
        "R"  : (600E3, "m"),
    },
    "Mun" : {
        "GM" : (6.514E10, "m3/s2"),
        "R"  : (200E3, "m")
    }
}


force_db = {
    "N" : 1.0,
    "kN" : 1E3
}


isp_db = {
    # Number of s in ISP unit
    "s" : 1.0,
    "m/s" : 1.0/9.81
}


length_db = {
    # number of meters in unit
    "m" : 1.0,
    "km" : 1E3,
    "Mm" : 1E6
}


mass_db = {
    # Number of kgs in unit
    "kg" : 1.0,
    "t"  : 907.185
}


template_craft = [
    ( "rocket", {
        "name" : "Templ 1",
        "stages" : [
            {
	        "m0" : (0,"t"),
	        "m"  : (0,"t"),
	        "T"  : (0,"kN"),
	        "Isp" : (0, "s")
            },
        ]
    }
    )
]


def analyzeRocket(rocket_def) :

    craft_type, craft_def = rocket_def
    
    if craft_type != "rocket" :
        raise Exception("Not a rocket record")

    print("\nRocket \"%s\"" % craft_def["name"])

    last_stage_m = None
    
    for istage, stage in enumerate(craft_def['stages']) :
        print("\nStage %i" % (istage+1))
        
        m0, m0u = stage["m0"]
        m0 *= uconv(mass_db, m0u, "kg")

        m, mu = stage["m"]
        m *= uconv(mass_db, mu, "kg")

        if last_stage_m is not None :
            if m0 != last_stage_m :
                print("\nWARNING: insconsistent inter-stage masses\n")
        last_stage_m = m

        isp, ispu = stage["Isp"]
        isp *= uconv(isp_db, ispu, "m/s")

        thr, thru = stage["T"]
        thr *= uconv(force_db, thru, "N")
        
        print("m0: %s kg" % m0)
        print("m: %s kg" % m)
        print("Isp: %s m/s" % isp)
        print("T: %s N" % thr)

        dv = isp * math.log(m0/m)

        print("DV: %s m/s" % dv)

        if istage == 0 :
            kerbin_ground_g = g("Kerbin", (0, "km"))
            print("Kerbin ground g: %s m/s2" % kerbin_ground_g)
            # Gravity and thrust correction for first stage
            dv_corr = isp * kerbin_ground_g / thr * (m - m0)
            print("First stage DV correction on Kerbin: %s m/s" % dv_corr)
            dv_final = dv + dv_corr
            print("First stage DV on Kerbin: %s m/s" % dv_final)

            mfuel = m0-m
            x1 = -0.5*kerbin_ground_g*((isp*mfuel/thr)**2)
            x2 = (isp**2)/thr*(m*math.log(m/m0)+mfuel)
            h_at_stage = x1+x2
            h_at_stage *= uconv(length_db, "m", "km")
            print("Alt at end of stage: %s km" % h_at_stage)


def dv_hohmann_apo(body, r_peri, r_apo) :
    '''DV computed when thrusting at apoapsis to change periapsis.
    '''
    
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    (R, Ru) = body_rec["R"]
    
    rperi, rperiu = r_peri
    rperi *= uconv(length_db, rperiu, "m")
    
    rapo, rapou = r_apo
    rapo *= uconv(length_db, rapou, "m")

    A = math.sqrt(GM/(rapo + R))
    B = 1.0 - math.sqrt(2.0*(rperi + R)/(rperi + rapo + 2.0*R))

    return A*B
    

def dv_hohmann_peri(body, r_peri, r_apo) :
    '''DV computed when thrusting at periapsis to change apoapsis.
    '''
    
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    (R, Ru) = body_rec["R"]
    
    rperi, rperiu = r_peri
    rperi *= uconv(length_db, rperiu, "m")
    
    rapo, rapou = r_apo
    rapo *= uconv(length_db, rapou, "m")

    A = math.sqrt(GM/(rperi + R))
    B = math.sqrt(2.0*(rapo + R)/(rperi + rapo + 2.0*R)) - 1.0

    return A*B
    

def dv_orbit(body, alt) :
    
    try :
        alt, altu = alt
    except :
        sys.stderr.write("*** Bad form for altitude.  Required: (value, \"unit\")\n")
        exit(1)

    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    (R, Ru) = body_rec["R"]

    alt *= uconv(length_db, altu,"m")
    
    r = R + alt

    g_ground = g(body, (0,"m"))

    dv_orbit = math.sqrt((GM/r) + 2.0*g_ground*alt)

    return dv_orbit
    

def g(body, alt) :
    '''Computes acceleration due to gravity near body "body" at altitude "alt"

    body: name of body
    alt: tuple, (value, "units")

    TODO: full blown unit conversion
    '''

    try :
        alt, altu = alt
    except :
        sys.stderr.write("*** Bad form for altitude.  Required: (value, \"unit\")\n")
        exit(1)

    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    (R, Ru) = body_rec["R"]

    alt *= uconv(length_db, altu,"m")
    
    r = R + alt

    g = GM / (r*r)

    return g



def listOfBodyNames(sep=" | ") :
    body_names = bodies_db.keys()
    return sep.join(body_names)


def orbit_v(body, alt) :
    
    try :
        alt, altu = alt
    except :
        sys.stderr.write("*** Bad form for altitude.  Required: (value, \"unit\")\n")
        exit(1)

    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    (R, Ru) = body_rec["R"]

    alt *= uconv(length_db, altu, "m")
    
    r = R + alt

    return math.sqrt(GM/r)


def uconv(db, ufrom, uto) :
    '''Computes multiplication factor to convert value of units
    "ufrom" to value of units "uto"

    '''
    return db[ufrom] / db[uto]


if __name__ == "__main__" :
    
    parser = argparse.ArgumentParser(description="KSP Calculator")

    subparsers = parser.add_subparsers(help="Commands", dest="command")

    # Command "g"
    cmd_g = subparsers.add_parser("g", help="Compute acceleration due to gravity for a specified body")
    cmd_g.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_g.add_argument("alt", help="\"(alt, 'unit')\"")

    # Command "craft"
    cmd_craft = subparsers.add_parser("craft", help="Analyzes a set of craft")
    cmd_craft.add_argument("fpath", help="Name of a craft JSON file.  If file can't be found, then generates a template to file to the name")

    # Command "dv_orbit"
    cmd_dv_orbit = subparsers.add_parser("dv_orbit", help="Computes circular orbital velocity")
    cmd_dv_orbit.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_dv_orbit.add_argument("alt", help="\"(alt, 'unit')\"")

    # Command "orbit_v"
    cmd_orbit_v = subparsers.add_parser("orbit_v", help="Computes circular orbital velocity")
    cmd_orbit_v.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_orbit_v.add_argument("alt", help="\"(alt, 'unit')\"")
    
    args = parser.parse_args()

    if args.command == "craft" :

        if not os.access( args.fpath, os.F_OK ) :
            print("Generating template craft file \"%s\"" % args.fpath)
            with open(args.fpath, "wt") as craft_file :
                json.dump(template_craft, craft_file, indent=4)
            exit(0)

        with open(args.fpath, "rt") as craft_file :
            crafts = json.load(craft_file)

        for craft in crafts :
            if craft[0] == "rocket" :
                analyzeRocket(craft)

                
    if args.command == "dv_orbit" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        print("DV to reach circular orbit: %s" % dv_orbit(args.body, altitude))

                
    if args.command == "g" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        print("Accel of gravity: %s" % g(args.body, altitude))

            
    if args.command == "orbit_v" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        print("Circular orbit speed: %s" % orbit_v(args.body, altitude))

