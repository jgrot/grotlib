#! /usr/bin/env python

# UNITS
#
# For now all functions should output MKS units
#

import argparse
import bisect
import copy
import json
import math
import numpy as np
import os
import scipy.optimize as spop
from scipy.integrate import ode
from scipy.interpolate import interp1d
import sys
import tabulate

# Verify Python version
if ( sys.version_info.major < 3 and sys.version_info.minor < 8 ) :
    print("PYTHON VERSION ", sys.version)
    sys.stderr.write("*** Wrong python version.  Make sure virtual environment is activated.\n")
    exit(1)

KSPCALCDIR = os.path.dirname(os.path.realpath(__file__))

Rgas = 8.3145 # J/mol/K
AirMmass = 28.84E-3 # kg/mol
p0_Pa = 100174.2 # For Isp scaling

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
        },
        "tday" : ( 6, "h" ),
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
    "t" : 1000.0,
    "su" : 7.5
}


time_db = {
    # number of seconds in unit
    "s" : 1.0,
    "m" : 60.0,
    "h" : 3600.0
}

# Unit DB's
udbs = [ dist_db, isp_db, force_db, mass_db, time_db ]

engine_db = {
    "BACC" : {
        "model" : "BACC",
        "name"  : "Thumper",
        "rate" : (19.423, "su"), # Per second (found under Propellant)
        "amount" : (820, "su"),
        "ispsl" : (175, "s"), # Isp at sea level
        "ispvac" : (210, "s")
        },
    "S2-17" : {
        "model" : "S2-17",
        "name" : "Thoroghbred",
        "rate" : (100.494, "su"), # Per second (found under Propellant)
        "amount" : (8000, "su"),
        "ispsl" : (205, "s"), # Isp at sea level
        "ispvac" : (230, "s")
    },
    "RT-5" : {
        "model" : "RT-5",
        "name" : "Flea",
        "rate" : (15.821, "su"),
        "amount" : (140, "su"),
        "ispsl" : (140, "s"),
        "ispvac" : (165, "s")
    },
    "RT-10" : {
        "model" : "RT-10",
        "name" : "Hammer",
        "rate" : (15.827, "su"),
        "amount" : (375, "su"),
        "ispsl" : (170, "s"),
        "ispvac" : (195, "s")        
    }
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


class Stage :
    ''' 

    @param engine_list: [ (N1, "Name of engine 1"), etc. ]
    '''
    
    def __init__( self, m0=None, engine_list=[], dragco=0.0 ) :
        self.m0 = m0
        self.engine_list = engine_list
        self.dragco = dragco

        if m0 is not None :
            self._assimilate( )

    def _assimilate( self ) :
        self._processEngines()
        
    def _processEngines( self ) :
        '''Compute engine info'''

        self.engines_by_t = []

        for erec in self.engine_list :
            n, ename = erec
            engine = engine_db[ename]
            
            fuel, fuelu = engine["amount"]
            fuel_kg = fuel * uconv( mass_db, fuelu, "kg" )

            rate, rateu = engine["rate"]
            rate_kgps = rate * uconv( mass_db, rateu, "kg" )

            burn_time = fuel_kg / rate_kgps

            self.engines_by_t.append( (burn_time, n, engine) )

        self.engines_by_t.sort()
        
        self.engines_by_t.reverse()

        self.mass_nodes = [0.0]
        self.rate_edges = []
        rate_total = 0.0
        rate_rec = []

        for irec, erec in enumerate(self.engines_by_t) :
            t, n, e = erec
            try :
                tm1, x1, x2 = self.engines_by_t[irec+1]
            except :
                tm1 = 0.0
            DT = t - tm1

            rate, rateu = e["rate"]
            rate_kg = rate * uconv( mass_db, rateu, "kg" )
            rate = n*rate_kg
            rate_total += rate

            ispsl, ispslu = e["ispsl"]
            ispsl_mps = ispsl * uconv( isp_db, ispslu, "m/s" )

            ispvac, ispvacu = e["ispvac"]
            ispvac_mps = ispvac * uconv( isp_db, ispvacu, "m/s" )

            rate_rec.append( (rate, ispsl_mps, ispvac_mps) )

            self.rate_edges.append( copy.copy(rate_rec) )
            self.mass_nodes.append( self.mass_nodes[-1] + DT*rate_total )

        self.engines_by_t.reverse()

        m0, m0u = self.m0
        m0_kg = m0 * uconv( mass_db, m0u, "kg" )
        
        self.m0_kg = m0_kg
        self.me_kg = m0_kg - self.mass_nodes[-1]
        
        
    def dmdt( self, m_kg, throttle ) :
        '''
        @param m_kg: current mass of stage in kg
        @param throttle: 0 - 1
        '''

        mfuel_kg = m_kg - self.me_kg

        if mfuel_kg < 1.0E-7 :
            return 0.0

        if mfuel_kg >= self.mass_nodes[-1] :
            rec = self.rate_edges[-1]
        else :
            rec = None
            for i in range(len(self.rate_edges)) :
                if mfuel_kg >= self.mass_nodes[i] and mfuel_kg < self.mass_nodes[i+1] :
                    rec = self.rate_edges[i]
                    break

        if rec is None :
            return 0.0
        else :
            
            r_tot = 0.0
            for r, isl, iva in rec : # rate, Isp-SL Isp-Vac
                r_tot += throttle*r

            return -r_tot

        
    def dvRemain( self, m_kg, p_Pa ) :
        ### See rocket book note "RB 2021-01-30 04.15.00.pdf"
        vac = (p0_Pa - p_Pa) / p0_Pa

        if m_kg <= self.me_kg :
            return 0.0
        
        # Mass of vehicle (at fuel mass nodes)
        M = [ self.me_kg + m for m in self.mass_nodes ]
        # One for each mass node
        DV = [0.0]
        
        for i in range(1,len(M)) :
            
            rate_rec = self.rate_edges[i-1]

            sr = 0.0
            st = 0.0

            for r, isl, iva in rate_rec :
                isp = (1.0 - vac)*isl + vac*iva
                sr += r
                st += (r*isp)

            isp_eqv = st / sr

            dv = isp_eqv * math.log( M[i] / M[i-1] )
            DV.append( dv+DV[i-1] )

            if m_kg >= M[i-1] and m_kg <= M[i] :
                a = ( m_kg - M[i-1] ) / ( M[i] - M[i-1] )
                return DV[i-1]*(1.0-a) + DV[i]*a

        # If we get here then it probably means this is a later
        # (ligther) stage being queried for its DV
        return DV[-1]

        
    def thrust( self, m_kg, throttle, p_Pa ) :
        '''
        @param m_kg: current mass of stage
        @param throttle: 0 - 1
        @param p_Pa: ambient pressure in Pa
        '''
        
        vac = (p0_Pa - p_Pa) / p0_Pa
        
        mfuel_kg = m_kg - self.me_kg
        
        if mfuel_kg < 1.0E-7 :
            return 0.0

        if mfuel_kg >= self.mass_nodes[-1] :
            rec = self.rate_edges[-1]
        else :
            rec = None
            for i in range(len(self.rate_edges)) :
                if mfuel_kg >= self.mass_nodes[i] and mfuel_kg < self.mass_nodes[i+1] :
                    rec = self.rate_edges[i]
                    break

        if rec is None :
            return 0.0
        else :
            
            th = 0.0
            for r, isl, iva in rec : # rate, Isp-SL Isp-Vac
                th += throttle*r*((1.0-vac)*isl + vac*iva)
                
            return th

    
    def dumpInfo( self ) :

        print()
        
        tabrows = [ [ "Initial mass (kg)", self.m0_kg ],
                    [ "Empty mass (kg)", self.me_kg] ]
        print(tabulate.tabulate(tabrows, headers=["Parameter", "Value"]))
        print()
        
        tabrows = []
        for t, n, engine in self.engines_by_t :
            tabrows.append([t, n, engine["model"], engine["name"]])
        print(tabulate.tabulate(tabrows, headers=["Burn Time (s)", "Count", "Model", "Name"]))
        print()
        
        tabrows = []
        for i in range(len(self.rate_edges)) :
            tabrows.append( [ self.mass_nodes[i], self.mass_nodes[i+1], self.rate_edges[i] ] )
        print(tabulate.tabulate(tabrows, headers=["M0", "M1", "Rate Record"]))
        print()

    def loadJSON( self, fname ) :
        data = None
        with open( fname, "rt" ) as f :
            data = json.load( f )
        self.m0 = data["m0"]
        self.engine_list = data["elist"]
        self.dragco = data["dragco"]
        self._assimilate( )        
        
    def dumpJSON( self, fname ) :
        data = { "m0" : self.m0,
                 "elist" : self.engine_list,
                 "dragco" : self.dragco
                 }
        with open( fname, "wt" ) as f :
            json.dump( data, f )



class FlyingStage :
    '''A Stage in the context of a body which provides functions needed for 2D trajectory.

    The main inputs are the throttle function and the orientation function

    @param fthrottle = function( t, y ) -> 0-1
    @param falpha = function( t, y, FlyingStage ) -> angle relative to body normal, RHR applies
    '''
    def __init__( self, stage, stage_name,  body_name, fthrottle, falpha ) :
        self.stage = stage
        self.stage_name = stage_name
        self.body_name  = body_name
        self.fthrottle = fthrottle
        self.falpha = falpha
        self.sm1 = None # Previous stage
        self.sp1 = None # Next stage

        self.body = bodies_db[ body_name ]
        (self.GM, GMu) = self.body["GM"] # units of m3/s2

        (R, Ru) = self.body["R"]
        self.R = R*uconv(dist_db, Ru, "m")

        self.fpress = self.body["fpress"]
        self.fdens = self.body["fdens"]

        tday, tday_u = self.body["tday"]
        tday_s = tday * uconv( time_db, tday_u, "s" )
        self.body_omega = 2.0*math.pi / tday_s
        
    def _thrust( self, t, y ) :
        m, r, th, vr, om = y
        alt = r - self.R
        p = self.fpress( alt )
        return self.stage.thrust( m, self.fthrottle(t, y), p )

    def a_r( self, t, y ) :
        m, r, th, vr, om = y
        alt = r - self.R
        
        return ( ( self._thrust(t, y) * math.sin( self.falpha(t, y, self) ) / m )
                 - self.GM/(r*r) - self.stage.dragco*abs(vr)*vr*self.fdens(alt)/m )

    def a_th( self, t, y ) :
        m, r, th, vr, om = y

        # Body relative theta motion for drag.  R motion is the same.
        omb = om - self.body_omega
        vthb = r*omb
        
        vth = r*om
        alt = r - self.R
        return ( ( self._thrust(t, y) * math.cos( self.falpha(t, y, self) ) / m )
                 - self.stage.dragco*abs(vthb)*vthb*self.fdens(alt)/m )
    
    def dmdt( self, t, y ) :
        m, r, th, vr, om = y
        return self.stage.dmdt( m, self.fthrottle(t, y) )


    def dumpTraj( self, t0 = None, t1 = None, dt = 5.0 ) :

        if t0 is None :
            t0 = self.solnt[0]

        if t1 is None :
            t1 = self.solnt[-1]
        
        t = t0

        rowdat = []

        flyer_last = None
        DV_last = 0.0
        DV = 0.0
        
        while t <= t1 :
            Y, crashed, flyer = self.flyTo( t )

            if flyer != flyer_last :
                flyer_last = flyer
                DV_last += DV
                print("DVLAST IS ", DV_last)
            
            m, r, th, vr, om = Y

            om_gnd = om - flyer.body_omega
            
            # Compute other things
            x = ( r * math.cos( th ) )
            y = ( r * math.sin( th ) )

            dvremain = flyer.dvRemain( m, p0_Pa )
            print("DEBUG FLYER PRESSURE ", flyer.fpress(r - flyer.R))
            act_dvremain = flyer.dvRemain( m, flyer.fpress(r-flyer.R) )
            sdvremain = flyer.stage.dvRemain( m, p0_Pa )

            vom_gnd = r*om_gnd
            spd_gnd = math.sqrt( vr*vr + vom_gnd*vom_gnd )

            vom = r*om
            spd = math.sqrt( vr*vr + vom*vom )
            
            row = [ flyer.stage_name, t, m, r, th, vr, om, om_gnd, sdvremain, dvremain, act_dvremain, spd_gnd, spd ]

            rowdat.append(row)

            t += dt

        headers = [ "stage", "time", "mass", "r", "theta", "v_r", "omega", "rel omega", "Stage DV (SL)", "Craft DV (SL)", "Craft DV", "Ground speed", "Orbit Speed" ]
        print( tabulate.tabulate( rowdat, headers = headers) )

        
    def dvRemain( self, m_kg, p_Pa ) :
        
        mydv = self.stage.dvRemain( m_kg, p_Pa )

        if self.sp1 is not None :
            return mydv + self.sp1.dvRemain( m_kg, p_Pa )
        else :
            return mydv
        

    def launch( self, y0 = None, sm1 = None, t0 = None ) :
        '''Launch the solver.  The stage can be anywhere in space.

        @param y0: intitial solver variables  ***OR***
        @param sm1: previous FlyngStage object
        @param t0: initial solver time.  Also, extracts the motion of the prvious stage if sm1 is set.
        '''
        
        self.crashed = False

        # Link the stages
        self.sm1 = sm1
        if sm1 is not None :
            sm1.sp1 = self

        if y0 is not None and sm1 is not None :
            raise Exception("Cannot set both y0 and sm1")
        
        if y0 is None and sm1 is None :
            y0 = [ self.stage.m0_kg, self.R, 0.0, 0.0, self.body_omega ]
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

        dmdt  = lambda t, y : self.dmdt(t, y)
        a_r   = lambda t, y : self.a_r(t, y)
        a_th  = lambda t, y : self.a_th(t, y)
        
        self.solv = trajectory2D_solver( dmdt, a_r, a_th, y0, t0 )

        self.solnt = [ t0 ]
        self.soln = [ y0 ]
        self.maxr = 0.0
        

    def flyTo( self, t ) :
        '''Advance and / or sample a trajectory.  Will sample from earlier
        stages if they exist.

        Returns: ( y(t), <Bool: crashed at t>, FlyingStage at t ) 

        '''

        if t < self.solnt[0] :
            if self.sm1 is not None :
                return self.sm1.flyTo( t )
            else :
                raise Exception("input time is before launch time")
        
        if t >= self.solnt[-1] :
            if not self.crashed :
                while self.solv.successful() and t >= self.solv.t:
                    y = self.solv.integrate( self.solv.t + 0.1 )
                    if y[1] <= self.R : ## y[1] := r
                        self.crashed = True
                        break
                    else :
                        self.maxr = max( y[1], self.maxr )
                        self.solnt.append( self.solv.t )
                        self.soln.append( list(y) )
            if self.crashed :
                print( "WARNING: FlyingStage is crashed at time %f" % t )
                return ( copy.copy(self.soln[-1]), True, self )

        y = interp2Dtraj( t, self.solnt, self.soln )
        
        return ( y, False, self )

    
    def plot( self, t0 = None, t1 = None, dt = 5.0 ) :

        import matplotlib
        import matplotlib.pyplot as plt

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

        xb = [ min(x), max(x) ]
        yb = [ min(y), max(y) ]

        xrng = xb[1] - xb[0]
        yrng = yb[1] - yb[0]

        maxrng = max( xrng, yrng )

        xpad = 0.5*(maxrng - xrng)
        ypad = 0.5*(maxrng - yrng)

        xmin = xb[0] - xpad
        xmax = xb[1] + xpad
        ymin = yb[0] - ypad
        ymax = yb[1] + ypad

        fig, ax = plt.subplots()
        ax.set_xlim( xmin, xmax )
        ax.set_ylim( ymin, ymax )
        ax.set_aspect(1.0)
        ax.plot(x, y, marker = "o")

        Nth = 100
        dth = 2.0*math.pi / float(Nth - 1)

        x = []
        y = []
        for ith in range(Nth) :
            th = dth * ith
            x.append(self.R*math.cos(th))
            y.append(self.R*math.sin(th))
        ax.plot(x, y)

        x = []
        y = []
        for ith in range(Nth) :
            th = dth * ith
            x.append((self.R+70E3)*math.cos(th))
            y.append((self.R+70E3)*math.sin(th))
        ax.plot(x, y)

        plt.show()

        
def augmentBodyDbs( ) :
    '''Looks for dictionary (JSON) file named <Body Key>.json and adds data to the dictionary'''
    for body_key in bodies_db :
        fname = os.path.join( KSPCALCDIR, "%s.json" % body_key )
        exists = os.access( fname, os.F_OK )
        if exists :
            with open( fname, 'rt' ) as f :
                extra_dat = json.load(f)
                bodies_db[body_key].update( extra_dat )


def processBodyDbs( ) :

    for body_key in bodies_db :
        bodyrec = bodies_db[body_key]
        if "apt" in bodyrec :
            alts, ps, ts = zip( *bodyrec["apt"] )

            dens = [ AirMmass*p/ts[i]/Rgas for i,p in enumerate(ps) ]

            bodyrec["fdens"] = interp1d( alts, dens, kind="quadratic", bounds_error=False, fill_value = 0.0 )
            bodyrec["fpress"] = interp1d( alts, ps, kind="quadratic", bounds_error=False, fill_value = 0.0 )


def interp2Dtraj( t, solnt, soln ) :
    '''
    @param soln: array of trajectory2D_solver y 
    '''
    
    # Rationale: I reserve the right to use a variable time
    # step interval in the future.
    i = (bisect.bisect( solnt, t ) - 1)
    t0 = solnt[i]
    t1 = solnt[i+1]
    a = (t - t0)/(t1 - t0)
    y0 = soln[i]
    y1 = soln[i+1]
    y = [ (1.0 - a)*y0[j] + a*y1[j] for j in range(len(y0)) ]
    
    return y

def trajectory2D_solver( dmdt, a_r, a_th, y0, t0 ) :
    '''Numerically integrates the planar orbit equations as a function of time.
    '''
    # y := [ m(kg), r(m), th(radians), vr(m/s), om(rad/s) ]
    #
    # dmdt:   function( t, y ) -> ks/s
    # a_r:    function( t, y ) -> acceleration along r axis
    # a_th:   function( t, y ) -> acceleration along th axis
    # y0:     initial conditions
    # t0:     starting time (s)

    def f(t, y) :
        m, r, th, vr, om = y
        return [ dmdt(t, y),                   # dmdt 
                 vr,                           # drdt
                 om,                           # dThdt
                 a_r(t, y) + r*om*om,          # dvrdt
                 ( a_th(t, y) - 2*vr*om ) / r  # domdt 
        ]

    solv = ode( f ).set_integrator( "vode", method="adams" )

    solv.set_initial_value( y0, t0 )

    return solv


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
    '''Computes acceleration due to gravity (in m/s^2) near body "body" at altitude "alt"

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
    
    augmentBodyDbs( )
    processBodyDbs( )

    if False :
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.use('TkAgg')

        kerbin_dens = bodies_db["Kerbin"]["fdens"]
        
        y = []
        for i in range(60000) :
            y.append( kerbin_dens(i) )
            
        fig, ax = plt.subplots()
        ax.plot(range(60000), y)
        plt.show()

    
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

    # Command "stage"
    cmd_stage = subparsers.add_parser("stage", help="Test stage")

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

        
    if args.command == "stage" :
        s2 = Stage((7.573,'t'), [(2,"RT-5"),(1,"RT-10")], .130)
        s2.dumpJSON( "T2DS2ME_stage2.json" )

        s1 = Stage((15.263,'t'), [(1,"BACC")], .130)
        s1.dumpJSON( "T2DS2ME_stage1.json" )
        
        
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
