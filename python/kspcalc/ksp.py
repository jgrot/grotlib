#! /usr/bin/env python
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

import argparse
import copy
import json
import math
import numpy as np
import os
import scipy.optimize as spop
from scipy.integrate import ode
import sys
import tabulate

# Grotlib modules
import compare as cmp
import mks_polar_motion as mpm
import moremath as mm
import mpl_tools as mpt

# Verify Python version
if ( sys.version_info.major < 3 and sys.version_info.minor < 8 ) :
    print("PYTHON VERSION ", sys.version)
    sys.stderr.write("*** Wrong python version.  Make sure virtual environment is activated.\n")
    exit(1)


##
## PROJECT PATHS
##

KSPCALCDIR = os.path.dirname(os.path.realpath(__file__))
KSPDBDIR = os.path.join( KSPCALCDIR, "data/databases" )
KSPEXAMPLESDIR = os.path.join( KSPCALCDIR, "data/examples" )
KSPUTDIR = os.path.join( KSPCALCDIR, "data/unit_tests" )

def pthdat( filename ) :
    '''Convenience: generate a full path to a relative location in the data directory'''
    return os.path.join( KSPCALCDIR, "data", filename )

def pthdb( filename ) :
    '''Convenience: generate a full path to a database file.
    '''
    return os.path.join( KSPDBDIR, filename )

def pthex( filename ) :
    '''Convenience: generate a full path to an example file.
    '''
    return os.path.join( KSPEXAMPLESDIR, filename )

def pthut( filename ) :
    '''Convenience: generate a full path to a unit test file.
    '''
    return os.path.join( KSPUTDIR, filename )


##
## PHYSICAL CONSTANTS
##


C_Rgas = 8.3145           # J/mol/K
C_air_mol_mass = 28.84E-3 # kg/mol
C_p0 = 100174.2           # Pa (ground pressure on Kerbin)
kerbin_day_s = 21549.0
mun_rot_period_s = 6.0*kerbin_day_s + 2.0*3600.0 + 36*60.0
minmus_rot_period_s = 1.0*kerbin_day_s + 5.0*3600.0 + 13*60.0
# Fuel to ox burn rate ratio
f_to_o_eng = 0.8181842495887883

##
## DATABASES
##

# Solar System Data
bodies_db = {
    "Kerbin" : {
        "GM" : (3.53E12, "m3/s2"),
        "R"  : (600E3, "m"),
        "satellites" : {
            "Mun"    : (11.4+.6, "Mm"),
            "Minmus" : (46.4+.6, "Mm")
        },
        "trot" : (kerbin_day_s, "s"), # Rotational period
    },
    "Mun" : {
        "GM" : (6.514E10, "m3/s2"),
        "R"  : (200E3, "m"),
        "VORB"  : (542.5, "m/s"),
        "SOI" : (2430.0, "km"),
        "trot" : (mun_rot_period_s, "s")
    },
    "Minmus" : {
        "GM" : (1.766E9, "m3/s2"),
        "R"  : (60, "km"),
        "VORB" : (274.1, "m/s"),
        "SOI" : (2247.0, "km"),
        "trot" : (minmus_rot_period_s, "s")
    }
}

# Databases for uconv()
angle_db = {
    # Number of radians in angle unit
    "radians" : 1.0,
    "degrees" : math.pi / 180.0
}
    
dist_db = {
    # Number of meters in unit
    "m" : 1.0,
    "km" : 1E3,
    "Mm" : 1E6,
    "pc" : 3.0857E16 # Parsec
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

mass_db = {
    # Number of kgs in unit
    "kg" : 1.0,
    "t" : 1000.0,
    "su" : 7.5,
    "lu" : 5.0
}

time_db = {
    # Number of seconds in unit
    "s" : 1.0,
    "m" : 60.0,
    "h" : 3600.0
}

udbs = [ angle_db, dist_db, force_db, isp_db, mass_db, time_db ]

# KSP Parts Databases 
engine_db = {
    "BACC" : {
        "model"  : "BACC",
        "name"   : "Thumper",
        "type"   : "s",
        "rate"   : (19.423, "su"), # Per second (found under Propellant)
        "amount" : (820.0, "su"),
        "ispsl"  : (175.0, "s"), # Isp at sea level
        "ispvac" : (210.0, "s")  # Isp at vacuum
    },
    "F3S0" : {
        "model"  : "F3S0",
        "name"   : "Shrimp",
        "type"   : "s",
        "rate"   : (1.897, "su"), # Per second (found under Propellant)
        "amount" : (90.0, "su"),
        "ispsl"  : (190.0, "s"), # Isp at sea level
        "ispvac" : (215.0, "s")  # Isp at vacuum
    },
    "S2-17" : {
        "model"  : "S2-17",
        "name"   : "Thoroghbred",
        "type"   : "s",
        "rate"   : (100.494, "su"),
        "amount" : (8000.0, "su"),
        "ispsl"  : (205.0, "s"),
        "ispvac" : (230.0, "s")
    },
    "RE-M3" : {
        "model"  : "RE-M3",
        "name"   : "Mainsail",
        "type"   : "l",
        "rate"   : ((44.407, "lu"),(54.275,"lu")),
        "amount" : ((0.0, "lu"),(0.0, "lu")),
        "ispsl"  : (285.0, "s"),
        "ispvac" : (310.0, "s")
    },
    "RT-5" : {
        "model"  : "RT-5",
        "name"   : "Flea",
        "type"   : "s",
        "rate"   : (15.821, "su"),
        "amount" : (140.0, "su"),
        "ispsl"  : (140.0, "s"),
        "ispvac" : (165.0, "s")
    },
    "RT-10" : {
        "model"  : "RT-10",
        "name"   : "Hammer",
        "type"   : "s",
        "rate"   : (15.827, "su"),
        "amount" : (375.0, "su"),
        "ispsl"  : (170.0, "s"),
        "ispvac" : (195.0, "s")        
    },
    "LFB KR-1x2" : {
        "model"  : "LFB KR-1x2",
        "name"   : "Twin-Boar",
        "type"   : "l",
        "rate"   : ((61.183, "lu"),(74.779, "lu")), # (fuel, ox)
        "amount" : ((2880.0, "lu"),(3520.0, "lu")),
        "ispsl"  : (280.0, "s"),
        "ispvac" : (300.0, "s")
    },
    "LV-T45" : {
        "model"  : "LV-T45",
        "name"   : "Swivel",
        "type"   : "l",
        "rate"   : ((6.166, "lu"),(7.536, "lu")), # (fuel, ox)
        "amount" : ((0.0, "lu"),(0.0, "lu")),
        "ispsl"  : (250.0, "s"),
        "ispvac" : (320.0, "s")
    },
}

tanks_db = {
    "FL-T400" : {
        "model" : "FL-T400",
        "name"  : "FL-T400 Fuel Tank",
        "amount": ((180.0, "lu"),(220.0, "lu")) # (fuel, ox)
    },
    "RJ-64" : {
        "model" : "RJ-64",
        "name"  : "Rockomax Jumbo-64 Fuel Tank",
        "amount": ((2880.0, "lu"),(3520.0, "lu"))
    },
    "FL-A215" : {
        "model" : "FL-A215",
        "name"  : "FL-A215 Fuel Tank Adapter",
        "amount": ((540.0, "lu"),(660.0, "lu"))
    },
    "FL-A151L": {
        "model" : "FL-A151L",
        "name"  : "FL-A151L Fuel Tank Adapter",
        "amount": ((270.0, "lu"),(330.0, "lu"))
    },
}

drag_db = {
    "ANC" : {
        "name" : "Aerodynamic Nose Cone",
        "drag" : 0.51047516,
        "dd_terms" : (2.84071867,  9.55284246, 0.77831561)
    },
    "ANC-A" : {
        "name" : "Aerodynamic Nose Cone Type A",
        "drag" : 0.51993179,
        "dd_terms" : (2.82551746,  9.54428153, 0.77963552)
    },
    "Mk5A" : {
        "name" : "",
        "drag" : 0.80194092,
        "dd_terms" : (2.71342946, 10.38393587, 0.83268328)
    },
    "KV-3" : {
        "name" : "'Pomegranate' Reentry Module",
        "drag" : 1.27636719,
        "dd_terms" : (3.16721787,  9.75074905, 0.67653654)
    },
    "Mk7" : {
        "name" : "Protective Rocket Nose Cone Mk7",
        "drag" : 1.10253906,
        "dd_terms" : (2.66280059, 10.467076,   0.83946844)
    }
}

def augmentBodyDbs() :
    '''Looks for dictionary (JSON) files named <Body Key>.json and adds data to the bodies_db dictionary'''
    
    for body_key in bodies_db :
        fname = pthdb( "%s.json" % body_key )
        exists = os.access( fname, os.F_OK )
        if exists :
            with open( fname, 'rt' ) as f :
                extra_dat = json.load(f)
                bodies_db[body_key].update( extra_dat )

def listOfBodyNames(sep = " | ") :
    body_names = bodies_db.keys()
    return sep.join(body_names)

def processBodyDbs() :
    '''Pre-processes bodies_db data.
    '''
    
    for body_key in bodies_db :
        bodyrec = bodies_db[body_key]
        
        if "apt" in bodyrec :
            # apt - Atmospleric pressure and temperature
            alts, ps, ts = zip( *bodyrec["apt"] )

            dens = [ C_air_mol_mass*p/ts[i]/C_Rgas for i,p in enumerate(ps) ]
            snd  = [ math.sqrt(142E3 / di) for di in dens ]

            bodyrec["fdens"] = mm.Interp1DFunctor(alts, dens, kind="quadratic")
            bodyrec["fpress"] = mm.Interp1DFunctor(alts, ps, kind="quadratic")
            bodyrec["fsnd"] = mm.Interp1DFunctor(alts, snd, kind="quadratic", low_fill = snd[0], high_fill=1E40)
        else :
            bodyrec["fdens"] = mm.Const1DFunctor(0.0)
            bodyrec["fpress"] = mm.Const1DFunctor(0.0)
            bodyrec["fsnd"] = mm.Const1DFunctor(math.inf)

        if "satellites" in bodyrec :

            for sat in bodyrec["satellites"] :
                d, du = bodyrec["satellites"][sat]
                d_m = d * uconv(dist_db, du, "m")
                
                satrec = bodies_db[sat]

                try :
                    vo, vou = satrec["VORB"]
                    if vou != "m/s" :
                        raise Exception("TODO: speed_db")
                    om = vo/d_m
                    satrec["rev_om"] = om
                except :
                    pass

def uconv(db, ufrom, uto) :
    '''Computes multiplication factor to convert value of units "ufrom" to value of units "uto"

    :param string db: name of units db
    :param string ufrom: units converting from
    :param string uto: units converting to

    :returns: float conversion factor
    '''
    return db[ufrom] / db[uto]


##
## PROCEDURES
##


def dInterp(term, dunit="m") :
    '''General distance interpreter.

    Terms are separated by spaces and therefore must not contain
    spaces.  Valid terms are:

    1. Regular value, units list or tuple pair
    2. String: D('body1','body2')
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
    '''Interprets a distance expression "D('body1', 'body2')"

    :param string dunit: destination unit

    :returns: float distance in units of dunit
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
    '''Interprets a body radius expression "R('body')"

    :returns: float radius in units of dunit
    '''

    if term[0] != "R" :
        raise Exception("Not a radius term")

    body = eval(term[1:])

    body_rec = bodies_db[body]

    R, Ru = body_rec["R"]

    R *= uconv(dist_db, Ru, dunit)

    return R

def dvHohmannApo(body, h_peri, h_apo) :
    '''DV computed when firing a) prograde at apo to attain circular orbit, or b) retrograde from apo to attain elliptical orbit.

    :param (d,"du") h_peri: periapsis *altitude* (not radius)
    :param (d,"du") h_apo:  apoapsis  *altitude* (not radius)

    :returns: float Delta-V for maneuver
    '''
    
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]
    
    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    hperi, hperiu = h_peri
    hperi *= uconv(dist_db, hperiu, "m")
    
    hapo, hapou = h_apo
    hapo *= uconv(dist_db, hapou, "m")

    r0 = R+hperi
    r1 = R+hapo

    return mpm.dv_r1_hohmann(GM, r0, r1)
    
def dvHohmannPeri(body, h_peri, h_apo) :
    '''DV computed when firing a) prograde at peri to attain elliptical orbit, or b) retrograde at peri to attain circular orbit.

    :param (d,"du") h_peri: periapsis *altitude* (not radius)
    :param (d,"du") h_apo:  apoapsis  *altitude* (not radius)

    :returns: float Delta-V for maneuver
    '''
    
    body_rec = bodies_db[body]
    (GM, GMu) = body_rec["GM"]

    (R, Ru) = body_rec["R"]
    R *= uconv(dist_db, Ru, "m")
    
    hperi, hperiu = h_peri
    hperi *= uconv(dist_db, hperiu, "m")
    
    hapo, hapou = h_apo
    hapo *= uconv(dist_db, hapou, "m")

    r0 = R+hperi
    r1 = R+hapo

    return mpm.dv_r0_hohmann(GM, r0, r1)

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

def dvTurn(speed_xpr, theta) :
    '''Computes DV needed for orbit inclination change

    :param float_or_python_expr speed_xpr: orbital speed of craft
    :param (ang,'u_ang') theta: turn angle
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

    :param string body: name of body
    :param (d,"du"() alt: altitude

    :returns: acceleration of gravity (m/s^2)
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

def orbitV(body, alt) :
    '''Computes speed of circular orbit.

    :param string body: name of body about which to orbit
    :param (d,"du") alt: altitude

    :returns: float orbit speed (m/s)
    '''
    
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

##
## CLASSES AND MODULE OBJECTS
##

class DragDivergence(mm.Functor) :
    '''Drag divergence is a phenomenon where the drag coefficient
    increases as vehicle speed approaches Mach 1.

    This model is in early stages of dev and testing, but a
    parameterized tanh form seems to best approximate the effect.

    :See also: optim_dd.py

    :Todo: cite references on tanh form.
    '''

    def __init__(self, c0=None, c1=None, c2=None) :
        rangemin = [0.0]
        rangemax = [20.0]
        
        super().__init__(rangemin, rangemax)

        if c0 is not None :
        
            self.c0 = c0
            self.c1 = c1
            self.c2 = c2
            
        else :

            # Collect all the drag divergence data
            experiments = []
            for part in drag_db :
                record = drag_db[part]
                if "dd_terms" in record :
                    experiments.append(record["dd_terms"])
            
            nexp = float(len(experiments))

            C = [0.0, 0.0, 0.0]
            for Cex in experiments :
                C = [C[i] + Cex[i] for i in range(3)]
            C = [c/nexp for c in C]

            # Formula parameters
            self.c0 = C[0] # Currently, like amplitude
            self.c1 = C[1] # Currently, like sharpness
            self.c2 = C[2] # Currently, like divergence Mach number
            
    def call( self, *X ) :
        self.checkRange( X )

        M = X[0]

        y = 1.0 + self.c0*0.5*(1.0 + math.tanh( self.c1*(M - self.c2) ))

        # print("C0 = %s C1 = %s C2 = %s M = %s DD = %s" % (self.c0, self.c1, self.c2, M, y))
        return [y]
        
    def nDep( self ) :
        return 1
    
dd = DragDivergence()


class PropulsionPhase :
    '''Given a collection of engines and tanks (with the same burn time)
    represents an equivalent single engine with a specific burn time,
    burn rate, Isp SL and Isp vac.

    The term "phase" refers to a DT section of the Mfuel v Time graph.

    CURRENT ASSUMPTIONS:

    - SRM phases contain only one type of SRM.  This is enforced in Stage._processEngines.

    - Fuel and ox are always in the correct proportions and fuel_level
      for liquid engines is the sum of the two.  Burn rate is the sum
      of the two.

    '''
    
    def __init__(self) :
        # Currently only supporting either solid or liquid since it
        # would be a little nuts to mix the two such that the burn
        # time is equal.
        self.type = None        # "s" or "l"
        self.burn_time = None   # seconds
        self.burn_rate = None   # kg/s
        self.ispsl = None       # m/s
        self.ispvac = None      # m/s
        self.engine_specs = []  # [(name, n, thrust_limit, fuel_level)]
        self.tank_specs = []    # [(name, n, fuel_level)]

    def _compute_parameters(self) :
        
        # Currently only needed for liquid systems.
        if self.type != "l" :
            raise Exception("This method does not yet support SRMs")

        # Add up all the fuel and oxidizer
        mfuel_sum = 0.0 # kg
        mox_sum   = 0.0 # kg
        
        for tank_spec in self.tank_specs :
            
            tank_model, ntanks, fuel_level = tank_spec

            if not tank_model in tanks_db :
                raise Exception("Tank %s is not in tanks_db" % tank_model)

            tank_rec = tanks_db[tank_model]

            max_fuel, u = tank_rec["amount"][0]
            max_fuel *= uconv(mass_db, u, "kg") # kg

            max_ox, u = tank_rec["amount"][1]
            max_ox *= uconv(mass_db, u, "kg") # kg

            if fuel_level is not None :
                
                liq_level, u = fuel_level
                liq_level *= uconv(mass_db, u, "kg")

                if liq_level <= 0.0 or liq_level > max_fuel + max_ox :
                    raise Exception("Modified fuel level out of range")

                liq_frac = liq_level / (max_fuel + max_ox)

            else :

                liq_frac = 1.0

            mfuel = ntanks*liq_frac*max_fuel
            mox   = ntanks*liq_frac*max_ox
            
            fo = mfuel / mox
            if not cmp.fsame(fo, f_to_o_eng, tolfrac=1E-4, report_to=None) :
                raise Exception("Unsupported: fuel and ox ratios must be in engine burn rate proportions for tank %s" % tank_model)

            mfuel_sum += mfuel # kg
            mox_sum   += mox # kg

        # Add up all the engine burn rates and also additional fuel amounts

        burn_rate_sum = 0.0 # kg/s
        thrust_vac_sum = 0.0 # N
        thrust_sl_sum = 0.0 # N

        for engine_spec in self.engine_specs :
            
            engine_model, neng, thrust_limit, fuel_level = engine_spec

            if not engine_model in engine_db :
                raise Exception("Engine %s is not in engine_db" % engine_model)

            engine_rec = engine_db[engine_model]

            #
            # Fuel and ox
            #

            max_fuel, u = engine_rec["amount"][0]
            max_fuel *= uconv(mass_db, u, "kg") # kg

            max_ox, u = engine_rec["amount"][1]
            max_ox *= uconv(mass_db, u, "kg") # kg

            if (max_fuel + max_ox) > 0.0 :
            
                if fuel_level is not None :
                    liq_level, u = fuel_level
                    liq_level *= uconv(mass_db, u, "kg")

                    if liq_level <= 0.0 or liq_level > max_fuel + max_ox :
                        raise Exception("Modified fuel level out of range")

                    liq_frac = liq_level / (max_fuel + max_ox)

                else :

                    liq_frac = 1.0

                mfuel = neng*liq_frac*max_fuel
                mox   = neng*liq_frac*max_ox

                if mox > 0.0 :
                    fo = mfuel / mox
                else :
                    if mfuel != 0.0 :
                        raise Exception('Zero mox but non-zero mfuel for engine %s' % engine_model)
                    else :
                        fo = f_to_o_eng
                
                if not cmp.fsame(fo, f_to_o_eng, tolfrac=1E-4, report_to=None) :

                    sys.stderr.write("Error: M oxidizer is zero\n")
                    sys.stderr.write("     : Engine model: %s\n" % engine_model)
                    sys.stderr.write("     : Number of engines specified: %i\n" % neng)
                    sys.stderr.write("     : Thrust limit specified: %s\n" % thrust_limit)
                    sys.stderr.write("     : Fuel level specified: %s\n" % repr(fuel_level))
                    sys.stderr.write("     : Max fuel from DB (kg): %s\n" % max_fuel)
                    sys.stderr.write("     : Max ox from DB (kg): %s\n" % max_ox)
                    sys.stderr.write("     : Computed liquid fraction from fuel level: %s\n" % liq_frac)
                    sys.stderr.write("     : Fuel to ox ratio: %s\n" % fo)
                    sys.stderr.write("     : Expected fuel to ox ratio: %s\n" % f_to_o_eng)
                    raise Exception("Unsupported: fuel and ox ratios must be in engine burn rate proportions for engine %s" % engine_model)

            else :

                mfuel = 0.0
                mox   = 0.0

            mfuel_sum += mfuel
            mox_sum   += mox

            #
            # Burn rate
            #
            
            max_fburn_rate, u = engine_rec["rate"][0]
            max_fburn_rate *= uconv(mass_db, u, "kg") # kg/s

            max_oburn_rate, u = engine_rec["rate"][1]
            max_oburn_rate *= uconv(mass_db, u, "kg") # kg/s
            
            if thrust_limit > 100.0 or thrust_limit <= 0.0 :
                raise Exception("Thrust limiter out of range")

            fburn_rate = neng*0.01*thrust_limit*max_fburn_rate # kg/s
            oburn_rate = neng*0.01*thrust_limit*max_oburn_rate # kg/s

            fo = fburn_rate / oburn_rate
            if not cmp.fsame(fo, f_to_o_eng, tolfrac=1E-4, report_to=None) :
                raise Exception("Logic error: fuel and ox burn rate ratio is incorrect for engine %s.  Expected: %s, actual: %s" % (engine_model, fo, f_to_o_eng))

            burn_rate = fburn_rate + oburn_rate
            burn_rate_sum += burn_rate
            
            #
            # Equivalent Isp
            #
            
            ispsl, u = engine_rec["ispsl"]
            ispsl *= uconv(isp_db, u, "m/s")

            ispvac, u = engine_rec["ispvac"]
            ispvac *= uconv(isp_db, u, "m/s")

            thrust_sl = ispsl*burn_rate
            thrust_vac = ispvac*burn_rate

            thrust_sl_sum += thrust_sl
            thrust_vac_sum += thrust_vac

        mliq = mfuel_sum + mox_sum
        
        self.burn_time = mliq / burn_rate_sum
        self.burn_rate = burn_rate_sum
        self.ispsl = thrust_sl_sum / burn_rate_sum
        self.ispvac = thrust_vac_sum / burn_rate_sum
            
    def _add_le(self, engine_rec, nengine, thrust_limit, fuel_level) :

        if self.type is None :
            # New propulsion phase
            self.type = "l"

        if engine_rec["type"] != self.type :
            raise Exception("Logic error: attempt to add liquid engine to solid motor propulsion phase")

        # Currently, can always add a liquid burner.
        self.engine_specs.append((engine_rec["model"], nengine, thrust_limit, fuel_level))
        self._compute_parameters()
        return True

    def _add_srm(self, engine_rec, nmotor, thrust_limit, fuel_level) :
        '''Add a group of identical SRMs

        :param dict engine_rec: a record from engines_db
        :param integer nmotor: number of SRMs in this group
        :param thrust_limit: >0 - 100
        :param tuple fuel_level: (float <amount>, string <unit>) alternate amount of fuel (None -> max)

        :returns: True if succeeded, False otherwise
        '''

        if self.type is not None :
            # print("Only one model of SRM is allowed at this time")
            return False

        self.type = "s"

        max_burn_rate, u = engine_rec["rate"]
        max_burn_rate *= uconv(mass_db, u, "kg") # kg/s

        if thrust_limit > 100.0 or thrust_limit <= 0.0 :
            raise Exception("Thrust limiter out of range")

        # Burn rate for one motor of this group (they are all the same)
        burn_rate = 0.01 * thrust_limit * max_burn_rate # kg/s

        # Fuel amount for one motor of this group (they are all the same)
        max_fuel, u = engine_rec["amount"]
        max_fuel *= uconv(mass_db, u, "kg") # kg

        if fuel_level is None :
            fuel_level = max_fuel # kg
        else :
            fuel_level, u = fuel_level
            fuel_level *= uconv(mass_db, u, "kg") # kg

            if fuel_level > max_fuel :
                raise Exception("Input fuel level exceeds the design level of %s" % srm_name)
        
        ispsl, u = engine_rec["ispsl"]
        ispsl *= uconv(isp_db, u, "m/s")

        ispvac, u = engine_rec["ispvac"]
        ispvac *= uconv(isp_db, u, "m/s")
            
        self.burn_time = fuel_level / burn_rate
        self.burn_rate = nmotor*burn_rate # kg/s
        self.ispsl = ispsl
        self.ispvac = ispvac
        self.engine_specs.append((engine_rec["model"], nmotor, thrust_limit, fuel_level))
        return True

    def add_engine(self, ename, nengines, thrust_limit=100.0, fuel_level=None) :
        '''Adds a group of identical engines

        :param string ename: name of engine in engines_db
        :param integer nengines: number of engines in this group (with the same specs)
        :param thrust_limit: >0 - 100
        :param tuple fuel_level: (float <amount>, string <unit>) alternate amount of fuel (None -> max)
        '''

        if ename not in engine_db :
            raise Exception("Could not find %s in the engine database" % ename)
            
        engine_rec = engine_db[ename]

        added_engine = False
        
        if engine_rec["type"] == "s" :
            added_engine = self._add_srm(engine_rec, nengines, thrust_limit, fuel_level)
        elif engine_rec["type"] == "l" :
            added_engine = self._add_le(engine_rec, nengines, thrust_limit, fuel_level)
        else :
            # There is also a helium engine we can add later
            raise Exception("Logic error: engine database contains an invalid engine type")

        return added_engine
        
    def add_tank(self, tank_model, ntanks, fuel_amount=None) :
        self.tank_specs.append((tank_model, ntanks, fuel_amount))
        self._compute_parameters()
    
class Stage :
    '''Model of an isolated stage.

    :param tuple m0: (<mass value>,'<mass unit>') mass of stage with engines completely full.  If None, assumes init via loadJSON().
    :param list engines: [("Name of engine 1", <count of engine 1>, thrust_limiter(>0-100), (<alternate fuel level>,"<units>") | None), ...]
    :param list tanks: [("Name of tank 1", <count of tank 1>, (<alternate fuel level>,"<units>")* | None)]
    :param list drags: [("Name of drag-adding part", <count of part>)] | None **
    :param float dragco: currently lumped term of 1/2*Cd*A **


    :Notes:
    
    \* The alternate fuel level for fuel/ox tanks is the sum of fuel
    and ox, and the ratio is assumed to be ``ksp.f_to_o_eng``. See
    also ``PropulsionPhase._compute_parameters()``.

    \*\* If ``drags`` is set, then ``dragco`` is ignored.

    '''
    
    def __init__(self, m0=None, engines=[], tanks=[], drags=None, dragco=0.0) :
        self.m0 = m0
        self.engines = engines
        self.tanks = tanks
        self.drags = drags
        self.dragco = dragco

        # Function cacheing
        self._dv_v_m_p = None

        if m0 is not None :
            self._assimilate()

    def _assimilate( self ) :
        '''(private) Perform all pre-processing based on input parameters.

        Call this function when top level input params change.
        '''
        self._processEngines()
        self._processDrags()

    def _dv_v_m(self, p_Pa) :
        '''Generates a DV node for each node in self.mass_nodes

        TODO: ascii art for the dv v m graph from rocket book note "RB 2021-01-30 04.15.00.pdf"

        
        '''

        if (p_Pa == self._dv_v_m_p) :
            # Already have the table for p (most likely at p=0)
            return
        
        vac = (C_p0 - p_Pa) / C_p0

        # One for each mass node
        self.dv_nodes = [0.0]
        
        for imass in range(1, len(self.mass_nodes)) :
            
            rate_rec = self.rate_edges[imass-1]

            sum_rate_kgps = 0.0
            sum_thrust = 0.0

            for rate_kgps, isl, iva in rate_rec :
                isp = (1.0 - vac)*isl + vac*iva
                sum_rate_kgps += rate_kgps
                sum_thrust += (rate_kgps * isp)

            isp_eqv = sum_thrust / sum_rate_kgps

            dv_phase = isp_eqv * math.log(self.mass_nodes[imass] / self.mass_nodes[imass-1])
            
            self.dv_nodes.append(dv_phase + self.dv_nodes[imass-1])
        
    def _processDrags(self) :

        self.total_drag = 0.0
        
        if self.drags is None :
            self.total_drag = self.dragco
            return

        for drag in self.drags :
            drag_part, ndrags = drag
            drag_rec = drag_db[drag_part]
            drag = drag_rec["drag"]
            self.total_drag += ndrags*drag

    def _processEngines(self) :
        '''(private) Preprocess engine data.'''

        # Here, we sort engine configurations by burn time.

        phases = []

        for eng_rec in self.engines :

            added_engine = False
            
            for phase in phases :
                
                # Try to add engine_rec to the phase.  If everything
                # checks out in terms of engine type and burn time then
                # the add succeeds and we are done.

                added_engine = phase.add_engine(*eng_rec)
                if added_engine : break

            if not added_engine :

                # If the engine could not be added to existing phases, then
                # create a new phase
                
                prop_phase = PropulsionPhase()
                prop_phase.add_engine(*eng_rec)

                phases.append(prop_phase)


        # Add fuel tanks.  New phases are not created for fuel tanks,
        # they must find a home in a liquid propulsion phase,
        # otherwise it is an error.
        if self.tanks is not None :
            for tank_rec in self.tanks :

                added_tank = False

                for phase in phases :

                    # Add tank to the first (and only for now) liquid prop phase.

                    if phase.type == "l" :

                        phase.add_tank(*tank_rec)
                        added_tank = True

                        break

                    if added_tank : break

                if not added_tank :

                    raise Exception("Could not find a propulsion phase for fuel tank")

        # Sort phase by burn time
        phases_by_t = [(p.burn_time, p) for p in phases]
        phases_by_t.sort()
        phases_by_t.reverse()

        #
        #       |o M_full
        #       | .
        #       |  . dM/dt = rate_i : all engines running
        #       |   .
        # Mfuel |    o M_i
        #       |    |   . dM/dt = rate0 : longest burning engines running
        #       |    |       .
        #       +----|-----------o M_0-------------> time
        #       |DTi |    DT0    |
        #
        #  Build a rate vs mass function using the notions in the plot above.
        #
        #  rate record: [ (<engine group max mass rate>, <engine group ASL Isp>, <engine group Vac Isp>), ... ]

        self.fmass_nodes = [0.0]
        self.rate_edges = []
        rate_kgps_sum = 0.0
        rate_rec = []

        # process engines starting from the longest burning engines
        # for irec, eng_rec in enumerate(engines_by_t_burn) :
        for iphase, rec in enumerate(phases_by_t) :
            
            t_burn, phase = rec
            
            # Note: propulsion phases are sorted by *descending* burn
            # time Try to get previous burn time, if we are out of
            # phase records, then end with t=0.
            try :
                t_burn_prev, phase_prev = phases_by_t[iphase+1]
            except :
                t_burn_prev = 0.0

            # Compute burn time interval
            DT = t_burn - t_burn_prev # DTi in the plot

            # This is used to compute DMass to get the mass nodes.
            rate_kgps_sum += phase.burn_rate

            # Add this engine group to the list of currently burning engine groups
            rate_rec.append((phase.burn_rate, phase.ispsl, phase.ispvac))

            # Add the collection of running engine groups to the rate
            # edges list.  Maybe come up with a new name for this?
            # Probably don't need this?  Although it does simplify things.
            self.rate_edges.append( copy.copy(rate_rec) )

            # Build list of masses at time nodes.
            self.fmass_nodes.append( self.fmass_nodes[-1] + DT*rate_kgps_sum )

        # self.engines_by_tburn = engines_by_t_burn
        self.burn_phases = phases_by_t
        
        m0, m0u = self.m0
        m0_kg = m0 * uconv( mass_db, m0u, "kg" )
        self.m0_kg = m0_kg
        # Empty stage mass (kg)
        self.me_kg = m0_kg - self.fmass_nodes[-1]
        # Stage mass nodes
        self.mass_nodes = [ self.me_kg + m for m in self.fmass_nodes ]
        
    def dmdt_or_thrust(self, mstage_kg, throttle, p_Pa=None) :
        '''Compute stage dm/dt OR thrust (if p_Pa is not none) vs firing phase
        (via stage mass) and throttle level

        :param float mstage_kg: mass of stage in kg
        :param float throttle: 0 - 1
        :param float p_Pa: ambient pressure in Pascals, or None to get dmdt
        '''

        # This might happen with roundoff error or solver overshoot.
        if mstage_kg <= self.me_kg :
            return 0.0
        
        mfuel_kg = mstage_kg - self.me_kg

        if p_Pa is not None :
            is_thrust = True
            vac = (C_p0 - p_Pa) / C_p0
        else :
            is_thrust = False
        
        # Find largest left-bracketing record
        for i in reversed(range(len(self.rate_edges))) :
            if mfuel_kg >= self.fmass_nodes[i] :
                rec = self.rate_edges[i]
                break

        out_sum = 0.0
        # rate_kgps_sum = 0.0
        for rate_kgps, isl, iva in rec : # rate, Isp-ASL, Isp-Vac
            if is_thrust :
                out_sum += throttle*rate_kgps*((1.0-vac)*isl + vac*iva)
            else :
                out_sum += throttle*rate_kgps

        if not is_thrust :
            out_sum = -out_sum # dmdt
            
        return out_sum

    def dumpInfo( self ) :

        tabrows = [ ["Initial mass (kg)", self.m0_kg],
                    ["Empty mass (kg)", self.me_kg],
                    ["Total drag", self.total_drag]
                   ]
        print("\n###")
        print("### STAGE INFO")
        print("###\n")
        print(tabulate.tabulate(tabrows, headers=["Parameter", "Value"]))
        
        tabrows = []
        for iburn, rate_rec in enumerate(self.rate_edges) :
            row = []
            m1 = self.fmass_nodes[iburn]
            m0 = self.fmass_nodes[iburn+1]

            ilast_eng = len(rate_rec)-1
            t_burn = self.burn_phases[ilast_eng][0]
            
            row.append(t_burn)
            row.append(m0)
            row.append(m1)

            for ieng, rate in enumerate(rate_rec) :

                rate_tot, isp_sl, isp_vac = rate
                
                t_burn, phase = self.burn_phases[ieng]

                # Build engine list string
                words = []
                for eng_name, n_eng, thrust_limit, fuel_amount in phase.engine_specs :
                    words.append("%i@%s@%i%%" % (n_eng,eng_name,thrust_limit))
                engines_str = " ".join(words)
                row.append(engines_str)
                row.append(rate_tot)
                row.append(isp_sl)
                row.append(isp_vac)
                
            tabrows.append(row)
            
        tabrows.reverse()

        maxlen = 0
        for tabrow in tabrows :
            maxlen = max(maxlen, len(tabrow))

        headers = [""]*maxlen
        headers[0] = "Burn Time(s)"
        headers[1] = "Fuel M0"
        headers[2] = "Fuel M1"

        nengrec = int((maxlen-3)/4)
        for ieng in range(nengrec) :
            i = 3+(ieng*4)
            headers[i]   = "Engines"
            headers[i+1] = "Rate Tot"
            headers[i+2] = "ISP SL"
            headers[i+3] = "ISP Vac"

        print("\n###")
        print("### STAGE BURN PHASES")
        print("###\n")
        print(tabulate.tabulate(tabrows, headers=headers))
        
        print()

    def dumpJSON( self, fname ) :
        '''Save initial parameters to a JSON file.
        '''
        
        data = { "m0" : self.m0,
                 "elist" : self.engines,
                 "drags" : self.drags,
                 "dragco" : self.dragco
                 }
        
        with open( fname, "wt" ) as f :
            json.dump( data, f )

    def dv_at_m( self, mstage_kg, p_Pa ) :
        '''Compute remaining Delta-V given firing phase (via stage mass) and ambient pressure.

        :param float mstage_kg: Stage mass in kg
        :param float p_pA: ambient pressure in Pascals
        '''

        if mstage_kg <= self.me_kg :
            return 0.0

        # Regenerate DV table (if needed)
        self._dv_v_m(p_Pa)

        if mstage_kg >= self.m0_kg :
            return self.dv_nodes[-1]

        return mm.bisect_interp(mstage_kg, self.mass_nodes, self.dv_nodes)

    def loadJSON( self, fname ) :
        '''Load initial parameters from a JSON file and pre-process.
        '''

        with open( fname, "rt" ) as f :
            data = json.load( f )
            self.m0 = data["m0"]
            self.engines = data["elist"]
            try :
                self.tanks = data["tanks"]
            except :
                self.tanks = None
            try :
                self.drags = data["drags"]
            except :
                self.drags = None
            self.dragco = data["dragco"]
            self._assimilate( )        

    def m_at_dv(self, dv_stage, p_Pa) :
        '''Compute stage mass given firing phase (via dv) and ambient pressure.

        :param float dv_stage: Remaining DV in m/s
        :param float p_pA: Ambient pressure in Pascals
        '''

        if dv_stage <= 0.0 :
            return self.me_kg

        # Regenerate table if necessary
        self._dv_v_m(p_Pa)
        
        return mm.bisect_interp(dv_stage, self.dv_nodes, self.mass_nodes)    

class FlyingStage :
    '''A Stage in the context of a body which provides functions needed for plane-polar trajectory.

    The main inputs are the throttle function and the orientation function.

    :param Stage stage: A Stage object
    :param string stage_name: Arbitrary name for flying stage
    :param string body_name: name of solar system body about which stage is flying
    :param function fthrottle: function( t, y ) -> 0-1
    :param function fthrustdir: function( t, y, FlyingStage ) -> (n_r, n_th)
    '''
    def __init__( self, stage, stage_name, body_name, fthrottle, fthrustdir ) :
        self.stage = stage
        self.stage_name = stage_name
        self.body_name  = body_name
        self.fthrottle = fthrottle
        self.fthrustdir = fthrustdir
        self.sm1 = None # Previous stage
        self.sp1 = None # Next stage

        self.body = bodies_db[ body_name ]
        (self.GM, GMu) = self.body["GM"] # units of m3/s2

        (R, Ru) = self.body["R"]
        self.R = R*uconv(dist_db, Ru, "m")

        self.fpress = self.body["fpress"]
        self.fdens  = self.body["fdens"]
        self.fsnd   = self.body["fsnd"]
 
        # Body rotational period
        tday, tday_u = self.body["trot"]
        tday_s = tday * uconv( time_db, tday_u, "s" )
        self.body_omega = 2.0*math.pi / tday_s
        
    def _a_r( self, t, y ) :
        '''(private) Compute radial component of acceleration

        :param float t: trajectory time
        :param list y: solver dependent values array (See makePolarMotionSolver)
        '''
        m, r, th, vr, om = y
        alt = r - self.R

        omb = om - self.body_omega
        vthb = r*omb
        
        v = math.sqrt(vr*vr + vthb*vthb)
        M = v/self.fsnd.call(alt)[0]

        # print("DEBUG M %f THRUST %f DMDT %f" % (m, self._thrust(t,y), self._dmdt(t,y) ))
        
        return ( ( self._thrust(t, y) * self.fthrustdir(t, y, self)[0] / m )
                 - self.GM/(r*r) - self.stage.total_drag * dd.call(M)[0] * abs(vr) * vr * self.fdens.call(alt)[0]/m )

    def _a_th( self, t, y ) :
        '''(private) Compute polar component of acceleration

        :param float t: trajectory time
        :param list y: solver dependent values array (See makePolarMotionSolver)
        '''
        m, r, th, vr, om = y
        alt = r - self.R
        
        # Body relative theta motion for drag.  R motion is the same.
        omb = om - self.body_omega
        vthb = r*omb

        v = math.sqrt(vr*vr + vthb*vthb)
        M = v/self.fsnd.call(alt)[0]
                
        vth = r*om
        alt = r - self.R
        return ( ( self._thrust(t, y) * self.fthrustdir(t, y, self)[1] / m )
                 - self.stage.total_drag * dd.call(M)[0] * abs(vthb) * vthb * self.fdens.call(alt)[0]/m )
    
    def _dmdt( self, t, y ) :
        '''(private) Compute dmdt

        :param float t: trajectory time
        :param list y: solver dependent values array (See makePolarMotionSolver)
        '''
        
        m, r, th, vr, om = y
        # Calling with no pressure term means compute dmdt
        return self.stage.dmdt_or_thrust(m, self.fthrottle(t, y))

    def _thrust( self, t, y ) :
        '''(private) Compute thrust at a trajectory point

        :param float t: trajectory time
        :param list y: solver dependent values array (See makePolarMotionSolver)
        '''
        
        m, r, th, vr, om = y
        alt = r - self.R
        p = self.fpress.call( alt )[0]
        return self.stage.dmdt_or_thrust( m, self.fthrottle(t, y), p )

    def dumpTraj( self, t0 = None, t1 = None, dt = 5.0 ) :
        '''Dump trajectory information.  Calls flyTo() and so takes into account the entire craft.

        :param float t0: time to start trajectory dump.  If None, will
                         start with the stage for which it is called.
        :param float t1: time to end trajectory dump.  If None, will
                         end with the last time stamp of the current
                         stage.
        :param float dt: time interval between dumped trajectory
                         points.  Default is 5.0 seconds.

        '''

        if t0 is None :
            t0 = self.solnt[0]

        if t1 is None :
            t1 = self.solnt[-1]
        
        t = t0

        rowdat = []

        while t <= t1 :
            
            Y, crashed, flyer = self.flyTo(t)

            if crashed :
                break
            
            m, r, th, vr, om = Y

            throttle = self.fthrottle(t, Y)
            n_r, n_th = self.fthrustdir(t, Y, self)
            
            om_gnd = om - flyer.body_omega
            
            # Compute other things
            x = r*math.cos(th)
            y = r*math.sin(th)

            craft_asl_dvremain = flyer.dv_at_m( m, C_p0 )
            craft_dvremain = flyer.dv_at_m( m, flyer.fpress.call(r-flyer.R)[0] )
            stage_dvremain = flyer.stage.dv_at_m( m, C_p0 )

            vom_gnd = r*om_gnd
            spd_gnd = math.sqrt( vr*vr + vom_gnd*vom_gnd )

            vom = r*om
            spd = math.sqrt( vr*vr + vom*vom )
            
            row = [ flyer.stage_name, t, m, throttle, n_r, n_th, r, th, vr, om, om_gnd,
                    stage_dvremain, craft_asl_dvremain, craft_dvremain, spd_gnd, spd ]

            rowdat.append(row)

            t += dt

        headers = [ "stage", "time", "mass", "throttle", "nr", "nth", "r", "theta", "v_r", "omega",
                    "rel omega", "Stage DV (ASL)", "Craft DV (ASL)",
                    "Craft DV", "Ground speed", "Orbit Speed" ]
        
        print( tabulate.tabulate( rowdat, headers = headers) )

    def dv_at_m( self, mcraft_kg, p_Pa ) :
        '''Compute remaining Delta-V given firing phase (via craft mass) and ambient pressure.

        :param float mcraft_kg: Stage mass in kg
        :param float p_pA: ambient pressure in Pascals
        '''
        
        mydv = self.stage.dv_at_m( mcraft_kg, p_Pa )

        if self.sp1 is not None :
            return mydv + self.sp1.dv_at_m( mcraft_kg, p_Pa )
        else :
            return mydv        

    def flyTo(self, t) :
        '''Advance and / or sample a trajectory.  Will sample from earlier stages if they exist.

        :param float t: trajectory time.

        :returns: ( y(t), <Bool: crashed at t>, FlyingStage at t )
        '''

        # Currently a fixed ODE time interval.  Needs more research as
        # to what the ODE solvers are doing with time stepping, but a
        # larger value of dt will give less accurate results.
        dt = 0.1
        
        if t < self.solnt[0] :
            if self.sm1 is not None :
                return self.sm1.flyTo( t )
            else :
                raise Exception("input time is before launch time")
        
        if t >= self.solnt[-1] :
            
            if not self.crashed :
                while self.solv.successful() and t >= self.solv.t:
                    y = self.solv.integrate( self.solv.t + dt )
                    if y[1] <= self.R : ## y[1] := r
                        self.crashed = True
                        break
                    else :
                        self.maxr = max( y[1], self.maxr )
                        self.solnt.append( self.solv.t )
                        self.soln.append( list(y) )

            if self.crashed :
                return ( copy.copy(self.soln[-1]), True, self )

        y = mm.bisect_interp(t, self.solnt, self.soln)
        
        return (y, False, self)

    def launch( self, y0 = None, sm1 = None, t0 = None ) :
        '''Launch the stage.  The stage can be anywhere in space.  Can add a previous FlyingStage at this point.

        :param list y0: (optional) intitial makePolarMotionSolver variables, defaults to ground launch, *** OR ***
        :param FlyingStage sm1: (optional) previous FlyngStage object
        :param float t0: (required if sm1 is set) initial solver time, defaults to zero.
        '''
        
        self.crashed = False

        # Link stages
        self.sm1 = sm1
        if sm1 is not None :
            sm1.sp1 = self

        if y0 is not None and sm1 is not None :
            raise Exception("Cannot set both y0 and sm1")
        
        if y0 is None and sm1 is None :
            # When nothing is specified, stage is launched from Kerbin base.
            y0 = [ self.stage.m0_kg, self.R+71.0, 0.0, 0.0, self.body_omega ]

        elif sm1 is not None :
            if t0 is None :
                raise Exception("Must set t0 for previous stage")
            y0, crashed, flyer = sm1.flyTo( t0 )
            if crashed :
                raise Exception("Cannot stage after a crash")
            # Replace m with M0 of this stage
            y0[0] = self.stage.m0_kg

        if t0 is None :
            t0 = 0.0

        dmdt  = lambda t, y : self._dmdt(t, y)
        a_r   = lambda t, y : self._a_r(t, y)
        a_th  = lambda t, y : self._a_th(t, y)
        
        self.solv = mpm.makePolarMotionSolver( dmdt, a_r, a_th, y0, t0 )

        self.solnt = [ t0 ]
        self.soln = [ y0 ]
        self.maxr = 0.0
        
    def sample(self, t0 = None, t1 = None, dt = 5.0) :
        '''Samples a trajectory for plotting.

        :param float t0: time in seconds to start trajectory dump.  If None, will
                         start with the stage for which it is called.
        :param float t1: time in seconds to end trajectory dump.  If None, will
                         end with the last time stamp of the current
                         stage.
        :param float dt: time interval between dumped trajectory
                         points.  Default is 5.0 seconds.

        :returns: tuple ([<x values>], [<y values>])
        '''

        if t0 is None :
            t0 = self.solnt[0]

        if t1 is None :
            t1 = self.solnt[-1]
        
        x = []
        y = []
        
        t = t0

        while t <= t1 :
            Y, crashed, flyer = self.flyTo( t )
            m, r, th, vr, om = Y
            x.append( r * math.cos( th ) )
            y.append( r * math.sin( th ) )
            t += dt
            
        return(x,y)

######################################################
##                                                  ##
##                                                  ##
##           EARLY EXPERIMENTAL CODE                ##
##                                                  ##
##        stuff written before Stage, FlyingStage   ##
##        Orbit, OrientedOrbit, etc.                ##
##                                                  ##
######################################################

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
    '''Analyze a rocket DV specified in a JSON file.

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
            
def dvOrbit(body, alt) :
    '''(Experimental) Compute DV to reach circular orbit about a body from launch.
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

    dvOrbit = math.sqrt((GM/R) - (GM/r))

    return dvOrbit

def jump(alt, body, Isp, Me, x, solve_for="thrust", vf=0.0, Edrag=0.0, ) :
    '''(Experimental) Compute fuel or thrust needed for a first stage to reach a certain
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
    :param string solve_for: "thrust" (default) or "mfuel"
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

def main() :
    '''Command Line Interface
    '''
    
    augmentBodyDbs( )
    processBodyDbs( )

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

    # Command "experiment"
    cmd_experiment = subparsers.add_parser("experiment", help="Run experimental code")

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

    # Command "plotfuncs"
    cmd_plotfuncs = subparsers.add_parser("plotfuncs", help="Plot up the model functions")

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
                
    if args.command == "orbitV" :
        print("Body: %s" % args.body)
        altitude = eval(args.alt)
        alt, unit = altitude
        print("Altitude: %s %s" % (alt, unit))

        print("Circular orbit speed: %s" % orbitV(args.body, altitude))

    if args.command == "plotfuncs" :

        import matplotlib
        import matplotlib.pyplot as plt

        fig1, ax1 = plt.subplots()
        dd.plot( ax1, [(0.0, 7.5)] )
        ax1.set_title("Drag Divergence")
        ax1.set_xlabel("Mach number")
        ax1.set_ylabel("Multiplier")

        fig2, ax2 = plt.subplots()
        bodies_db["Kerbin"]["fdens"].plot( ax2, [True] )
        ax2.set_title("Kerbin Density")
        ax2.set_xlabel("Altitude (m)")
        ax2.set_ylabel("Density (kg/m^3)")
        
        fig3, ax3 = plt.subplots()
        bodies_db["Kerbin"]["fpress"].plot( ax3, [True] )
        ax3.set_title("Kerbin Pressure")
        ax3.set_xlabel("Altitude (m)")
        ax3.set_ylabel("Pressure (Pa)")

        fig4, ax4 = plt.subplots()
        bodies_db["Kerbin"]["fsnd"].plot( ax4, [True] )
        ax4.set_title("Kerbin Sound Speed")
        ax4.set_xlabel("Altitude (m)")
        ax4.set_ylabel("Sound Speed (m/s)")
        
        plt.show()
                
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

    if args.command == "experiment" :

        print("No new experiments at this time.  See examples directory.")
            

if __name__ == "__main__" :
    main()
