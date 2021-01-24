#! /usr/bin/env python

# UNITS
#
# For now all functions should output MKS units
#

import argparse
import json
import math
import os
import scipy.optimize as spop
import sys
import tabulate

# Verify Python version
if ( sys.version_info.major < 3 and sys.version_info.minor < 8 ) :
    print("PYTHON VERSION ", sys.version)
    sys.stderr.write("*** Wrong python version.  Make sure virtual environment is activated.\n")
    exit(1)

angle_db = {
    # Number of radians in angle unit
    "radians" : 1.0,
    "degrees" : math.pi / 180.0
}
    
bodies_db = {
    "Kerbin" : {
        "GM" : (3.53E12, "m3/s2"),
        "R"  : (600E3, "m"),
        "satellites" : {
            "Mun"    : (11.4,"Mm"),
            "Minmus" : (46.4,"Mm")
        }
    },
    "Mun" : {
        "GM" : (6.514E10, "m3/s2"),
        "R"  : (200E3, "m"),
    },
    "Minmus" : {
        "GM" : (1.766E9, "m3/s2"),
        "R"  : (60, "km"),
    }
}


force_db = {
    "N" : 1.0,
    "kN" : 1E3
}


isp_db = {
    # Number of s in ISP unit
    "s" : 1.0, # KSP Isp units are seconds
    "m/s" : 1.0/9.81
}


dist_db = {
    # number of meters in unit
    "m" : 1.0,
    "km" : 1E3,
    "Mm" : 1E6,
    "pc" : 3.0857E16 # Parsec
}


mass_db = {
    # Number of kgs in unit
    "kg" : 1.0,
    "t"  : 907.185,
    "su" : 6.8038875
}


# Unit DB's
udbs = [ dist_db, isp_db, force_db, mass_db ]


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


template_maneuvers = [
    {
        "type"  : "launch",
        "desc"  : "Launch from Kerbin to circular orbit",
	"body"  : "Kerbin",
	"orbit" : [100,"km"]
    },
        {
	"type" : "hperi",
        "desc" : "Extend orbit to the Mun",
	"body" : "Kerbin",
	"peri" : [100,"km"],
	"apo"  : "D('Kerbin','Mun') + R('Mun') + (100,'km')"
    },
    {
        "type"  : "launch",
        "desc"  : "Land on Mun",
	"body"  : "Mun",
	"orbit" : [100,"km"]
    },
    {
        "type"  : "launch",
        "desc"  : "Launch from Mun",
	"body"  : "Mun",
	"orbit" : [100,"km"]
    },
    {
	"type" : "hapo",
        "desc" : "Pull orbit into Kerbin",
	"body" : "Kerbin",
	"peri" : [45,"km"],
	"apo"  : "D('Kerbin','Mun') + R('Mun') + (100,'km')"
    },
    {
	"type" : "hperi",
        "desc" : "Reduce orbital energy",
	"body" : "Kerbin",
	"peri" : [45,"km"],
	"apo"  : "D('Kerbin','Mun') + R('Mun') + (100,'km')"
    },
]


def analyzeRocket(rocket_def) :
    '''Analyze a rocket DV specified in a JSON file

    Takes into account gravity when computing the DV for the first stage.

    :param string rocket_def: Path to JSON file containing the rocket definition.
    '''

    craft_type, craft_def = rocket_def
    
    if craft_type != "rocket" :
        raise Exception("Not a rocket record")

    print("\nRocket \"%s\"" % craft_def["name"])

    last_stage_m = None

    dvtot = 0.0
    
    for istage, stage in enumerate(craft_def['stages']) :
        print("\nStage %i" % (istage+1))
        
        m0, m0u = stage["m0"]
        m0 *= uconv(mass_db, m0u, "kg")

        m, mu = stage["m"]
        m *= uconv(mass_db, mu, "kg")

        if last_stage_m is not None :
            if m0 != last_stage_m :
                print("\nWARNING: inconsistent inter-stage masses\n")
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
            h_at_stage *= uconv(dist_db, "m", "km")
            print("Alt at end of stage: %s km" % h_at_stage)
        else :
            dv_final = dv

        dvtot += dv_final

        print("DV Total:", dvtot)

            
def dInterp(term, dunit="m") :
    '''General distance interpreter.  Allowed terms are:

    1. Regular value, units list or tuple pair
    2. String: D('body1', 'body2')
    3. String: R('body')
    4. String containing all valid terms above separated by operators

    '''
    if isinstance(term, str) :
        terms = term.split()
        # Currently only allowing an odd number of terms where even terms must be + or -
        nterms = len(terms)
        if nterms % 2 != 1 :
            raise Exception("Bad number of terms")
        # Convert odd terms (0,2,...) to meters
        terms2 = []
        for iterm, term in enumerate(terms) :
            if iterm % 2 == 1 :
                terms2.append(term)
            else :
                if term[0] == "D" :
                    # This is a distance calculation
                    d = dInterpDist(term, dunit)
                    terms2.append(d)
                elif term[0] == "R" :
                    r = dInterpRadius(term, dunit)
                    terms2.append(r)
                elif term[0] in ("(","[") :
                    # Expecting regular distance term
                    d, du = eval(term)
                    d *= uconv(dist_db, du, dunit)
                    terms2.append(d)
                    
        # Compine the distance terms

        dtot = 0.0
        op = None
        for i, x in enumerate(terms2) :
            if i % 2 == 0 :
                if op == '+' or op is None :
                    dtot += x
                elif op == '-' :
                    dtot -= x
            else :
                op = x

        return dtot
        
    elif isinstance(term, list) or isinstance(term, tuple) :
        d, du = term
        d *= uconv(dist_db, du, dunit)
        return d


def dInterpDist(term, dunit="m") :
    '''Operates on a string "D('body1', 'body2')"
    '''
    
    if term[0] != "D" :
        raise Exception("Not a distance term")

    b1, b2 = eval(term[1:])

    b1_rec = bodies_db[b1]
    b2_rec = bodies_db[b2]

    try :
        if "satellites" in b1_rec :
            d, du = b1_rec["satellites"][b2]
        else :
            d, du = b2_rec["satellites"][b1]
    except KeyError :
        raise Exception("Could not look up distance between %s and %s" % (b1, b2))

    d *= uconv(dist_db, du, dunit)

    return d


def dInterpRadius(term, dunit="m") :

    if term[0] != "R" :
        raise Exception("Not a radius term")

    body = eval(term[1:])

    body_rec = bodies_db[body]

    R, Ru = body_rec["R"]

    R *= uconv(dist_db, Ru, dunit)

    return R


def dvHohmannApo(body, r_peri, r_apo) :
    '''DV computed when thrusting a) prograde at apo to attain circular orbit, or b) retrograde from apo to attain elliptical orbit.
    '''
    
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    
    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    rperi, rperiu = r_peri
    rperi *= uconv(dist_db, rperiu, "m")
    
    rapo, rapou = r_apo
    rapo *= uconv(dist_db, rapou, "m")

    A = math.sqrt(GM/(rapo + R))
    B = 1.0 - math.sqrt(2.0*(rperi + R)/(rperi + rapo + 2.0*R))

    return A*B
    

def dvHohmannPeri(body, r_peri, r_apo) :
    '''DV computed when thrusting a) prograde at peri to obtain elliptical orbit, or b) retrograde at peri to obtain circular orbit.
    '''
    
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]

    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    rperi, rperiu = r_peri
    rperi *= uconv(dist_db, rperiu, "m")
    
    rapo, rapou = r_apo
    rapo *= uconv(dist_db, rapou, "m")

    A = math.sqrt(GM/(rperi + R))
    B = math.sqrt(2.0*(rapo + R)/(rperi + rapo + 2.0*R)) - 1.0

    return A*B
    

def dvInterp(maneuvers) :
    '''Interpret a set of maneuvers to compute DV map'''

    dv_tot = 0.0

    tabrows = []
    
    for m in maneuvers :
        tabrow = []

        try :
            tabrow.append(m['desc'])
        except :
            tabrow.append("Maneuver")
        
        mtype = m["type"]

        tabrow.append(mtype)
        
        if mtype == "launch" :
            dv = dvOrbit(m["body"], m["orbit"])

        elif mtype == "hperi" :
            peri = m["peri"]
            apo = m["apo"]
            p = dInterp(peri, "m")
            a = dInterp(apo, "m")
            dv = dvHohmannPeri(m["body"], (p,"m"), (a,"m"))

        elif mtype == "hapo" :
            peri = m["peri"]
            apo = m["apo"]
            p = dInterp(peri, "m")
            a = dInterp(apo, "m")
            dv = dvHohmannApo(m["body"], (p,"m"), (a,"m"))

        elif mtype == "turn" :
            speed = m["speed"]
            ang = m["angle"]
            dv = dvTurn(speed, ang)

        dv_tot += dv

        tabrow.append(dv)
        tabrow.append(dv_tot)

        tabrows.append(tabrow)
        
    print(tabulate.tabulate(tabrows, headers=["Description","Type","DV","DV Sum"]))
    
    return dv_tot


def dvOrbit(body, alt) :
    
    try :
        alt, altu = alt
    except :
        sys.stderr.write("*** Bad form for altitude.  Required: (value, \"unit\")\n")
        exit(1)

    alt *= uconv(dist_db, altu,"m")
            
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]

    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    r = R + alt

    dvOrbit = math.sqrt((GM/R) - (GM/r))

    return dvOrbit


def dvTurn(speed_xpr, theta) :
    '''Computes an on DV needed for a particular turn.

    :param speed_xpr float_or_python_expr: orbital speed of craft
    :param theta (ang, 'u_ang'): turn angle
    '''

    if isinstance(speed_xpr, float) or isinstance(speed_xpr, int) :
        speed = float(speed_xpr)
    else :
        speed = eval(speed_xpr)
    
    theta, u_theta = theta
    theta *= uconv(angle_db, u_theta, "radians")

    return speed * theta


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

    alt *= uconv(dist_db, altu,"m")

    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]

    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    r = R + alt

    g = GM / (r*r)

    return g


def listOfBodyNames(sep = " | ") :
    body_names = bodies_db.keys()
    return sep.join(body_names)


def jump(alt, body, Isp, Me, x, solve_for="thrust", vf=0.0, Edrag=0.0, ) :
    '''Compute fuel or thrust needed for a first stage to reach a certain
       altitude

    Solves for fuel mass or thrust in the equation below, such that
    the kinetic energy after stage 1 exhausts is eaten up by the
    remaining potential energy to get to the final altitude with a
    final velocity (v_f)

    .. math::
       1/2 (\Delta v(m_F))^2 - g(h-h_{\Delta v}(m_F)) - 1/2 {v_f}^2 = 0

    :param (alt,"altu") alt: Altitude
    :param str body: Body launching from
    :param float Isp:  Units of seconds
    :param (m,"mu") Me: Mass (empty) after fuel is spent
    :param (x,"Xu") x: Mass Fuel (if solve_for=="thrust"), Thrust (if solve_for=="mfuel")
    :param string solve_for: "thrust" (default) or "mass"
    :param float vf: speed at target altitude
    :param float Edrag: (J/kg) (g*Dh) Energy per mass lost due to drag.  Have to do an experiment to get it.

    '''

    alt, altu = alt
    alt *= uconv(dist_db, altu, "m")
    Isp *= 9.81
    x, Xu = x
    g_ground = g(body, (0,"m"))
    Me, Meu = Me
    Me *= uconv(mass_db, Meu, "kg")
    

    if solve_for == "thrust":
        # Easier to solve for inv thrust
        
        Mf = x * uconv(mass_db, Xu, "kg")

        def dv(R) :
            return Isp*(math.log((Me + Mf)/Me) - g_ground*Mf*R)
        def h(R) :
            return -0.5*g_ground*((Isp*R)**2) + (Isp**2)*R*(Me*math.log(Me/(Me+Mf)) + Mf)
        def F(R) :
            return 0.5*(dv(R)**2) - g_ground*(alt - h(R)) - 0.5*(vf**2) - Edrag

        # The guess is 100 kN.  fsolve returns a list of zeros, so take first element.
        inv_thrust = spop.fsolve( F, 1E-5 )[0]
    
        return 1.0/inv_thrust
        
    elif solve_for == "mfuel" :
        T = x * uconv(force_db, Xu, "N")

        def dv(Mf) :
            return Isp*(math.log((Me + Mf)/Me) - g_ground*Mf/T)
        def h(Mf) :
            return -0.5*g_ground*((Isp/T)**2) + (Isp**2)/T*(Me*math.log(Me/(Me+Mf)) + Mf)
        def F(Mf) :
            return 0.5*(dv(Mf)**2) - g_ground*(alt - h(Mf)) - 0.5*(vf**2) - Edrag

        # The guess is zero fuel.  fsolve returns a list of zeros, so take first element.
        mfuel_kg = spop.fsolve( F, 0.0 )[0]
    
        mfuel_t = mfuel_kg*uconv(mass_db, "kg", "t")

        return mfuel_t
        
    else :
        raise Exception("Bad value for solve_for")


def orbitV(body, alt) :
    
    try :
        alt, altu = alt
    except :
        sys.stderr.write("*** Bad form for altitude.  Required: (value, \"unit\")\n")
        exit(1)

    alt *= uconv(dist_db, altu, "m")
    
    body_rec = bodies_db[body]

    (GM, GMu) = body_rec["GM"]

    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    r = R + alt

    return math.sqrt(GM/r)


def uconv(db, ufrom, uto) :
    '''Computes multiplication factor to convert value of units
    "ufrom" to value of units "uto"

    '''
    return db[ufrom] / db[uto]


def main() :
    '''Command Line Interface

    Commands: g, craft, dvOrbit, orbitV
    '''
    parser = argparse.ArgumentParser(description="KSP Calculator")

    subparsers = parser.add_subparsers(help="Commands", dest="command")

    # Command "craft"
    cmd_craft = subparsers.add_parser("craft", help="Analyzes a set of craft")
    cmd_craft.add_argument("fpath", help="Name of a craft JSON file.  If file can't be found, then generates a template to file to the name")

    # Command "dvMap"
    cmd_dvMap = subparsers.add_parser("dvMap", help="Computes DV map for a set of maneuvers")
    cmd_dvMap.add_argument("fpath", help="Name of a maneuvers JSON file.  If file can't be found, then generates a template to file to the name")
    
    # Command "dvOrbit"
    cmd_dvOrbit = subparsers.add_parser("dvOrbit", help="Computes DV to get into circular orbit")
    cmd_dvOrbit.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_dvOrbit.add_argument("alt", help="\"(alt, 'unit')\"")

    # Command "jump"
    cmd_jump = subparsers.add_parser("jump", help="Computes first stage fuel or thrust to reach a target altitude and speed.  Tip: MINIMIZE DRAG ON VEHICLE.")
    cmd_jump.add_argument("alt", help="\"(alt, 'unit')\"")
    cmd_jump.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_jump.add_argument("Isp", type=float, help="Isp (s)")
    cmd_jump.add_argument("Me", help="Mass when empty \"(Mass, 'unit')\"")
    cmd_jump.add_argument("X", help="\"(X, 'unit of X')\" where X is thrust if solving for mass, and vice-versa")
    cmd_jump.add_argument("solve_for", help="\"thrust\" or \"mass\"")
    cmd_jump.add_argument("speed", type=float, help="Speed at target altitude (m/s).")
    cmd_jump.add_argument("edrag", type=float, help="Measured drag energy lost (J/kg) = g*(H_no_drag - H_drag) done with target speed = 0.")

    # Command "g"
    cmd_g = subparsers.add_parser("g", help="Compute acceleration due to gravity for a specified body")
    cmd_g.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_g.add_argument("alt", help="\"(alt, 'unit')\"")
    cmd_g.add_argument("--mass", help="Optionally compute force on a mass at the altitude\"(mass, 'unit')\"")

    # Command "orbitV"
    cmd_orbitV = subparsers.add_parser("orbitV", help="Computes circular orbital velocity")
    cmd_orbitV.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_orbitV.add_argument("alt", help="\"(alt, 'unit')\"")

    # Command "uconv"
    cmd_uconv = subparsers.add_parser("uconv", help="Call the units converter function")
    cmd_uconv.add_argument("vin", help="\"(value, 'unit')\"")
    cmd_uconv.add_argument("uout", help="Unit out")
    
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

    if args.command == "dvMap" :
        if not os.access( args.fpath, os.F_OK ) :
            print("Generating template maneivers file \"%s\"" % args.fpath)
            with open(args.fpath, "wt") as man_file :
                json.dump(template_maneuvers, man_file, indent=4)
            exit(0)
            
        with open(args.fpath, "rt") as man_file :
            maneuvers = json.load(man_file)

        dvInterp(maneuvers)
                
    if args.command == "dvOrbit" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        print("DV to reach circular orbit: %s" % dvOrbit(args.body, altitude))

    if args.command == "jump" :
        alt = eval(args.alt)
        body = args.body
        Isp = args.Isp
        Me = eval(args.Me)
        X = eval(args.X)
        solve_for = args.solve_for
        speed = args.speed
        edrag = args.edrag

        if solve_for == "mfuel" :
            m_t = jump(alt, body, Isp, Me, X, solve_for, speed, edrag)
            solid_units = m_t/.0075
            liq_units = m_t/.005
            tabrows=[
                ["Fuel mass (t)", m_t],
                ["Solid units", solid_units],
                ["Liquid units", liq_units]
            ]
            print(tabulate.tabulate(tabrows))
        elif solve_for == "thrust" :
            thrust = jump(alt, body, Isp, Me, X, solve_for, speed, edrag)
            tabrows=[
                ["Thrust", thrust],
            ]
            print(tabulate.tabulate(tabrows))
        else :
            raise Exception("Bad value for solve_for")
        
    if args.command == "g" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        accel_g = g(args.body, altitude)
        print("Accel of gravity: %s" % accel_g)

        if args.mass is not None :
            mass, umass = eval( args.mass )
            mass2 = mass * uconv( mass_db, umass, "kg" )
            force_kN = accel_g * mass2 / 1E3
            print("Force on mass of %f %s is %f kN" % ( mass, umass, force_kN ))

            
    if args.command == "orbitV" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        print("Circular orbit speed: %s" % orbitV(args.body, altitude))


    if args.command == "uconv" :
        value, unit = eval( args.vin )

        # Find database with input unit.  If the unit is not found,
        # then udb will be the last one in the list and we'll get an
        # error downstream.
        for udb in udbs :
            if unit in udb :
                break

        try :
            value2 = value * uconv( udb, unit, args.uout )
        except :
            value2 = None

        if value2 is not None :
            print("%f %s -> %f %s" % ( value, unit, value2, args.uout ))
        else :
            print("Could not convert %f %s" % (value, unit) )


if __name__ == "__main__" :
    main()
