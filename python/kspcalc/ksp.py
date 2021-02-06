#! /usr/bin/env python
# Copyright © 2021 Jonathan Grot

import argparse
import bisect
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
import functor

# Verify Python version
if ( sys.version_info.major < 3 and sys.version_info.minor < 8 ) :
    print("PYTHON VERSION ", sys.version)
    sys.stderr.write("*** Wrong python version.  Make sure virtual environment is activated.\n")
    exit(1)

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



#
# Constants
#

C_Rgas = 8.3145           # J/mol/K
C_air_mol_mass = 28.84E-3 # kg/mol
C_p0 = 100174.2           # Pa (ground pressure on Kerbin)

#
# Solar system data
#

bodies_db = {
    "Kerbin" : {
        "GM" : (3.53E12, "m3/s2"),
        "R"  : (600E3, "m"),
        "satellites" : {
            "Mun"    : (11.4, "Mm"),
            "Minmus" : (46.4, "Mm")
        },
        "tday" : (6.0, "h"),
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

#
# Databases for uconv()
#

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
    "su" : 7.5
}

time_db = {
    # Number of seconds in unit
    "s" : 1.0,
    "m" : 60.0,
    "h" : 3600.0
}

udbs = [ angle_db, dist_db, force_db, isp_db, mass_db, time_db ]

#
# KSP Parts databases 
#

engine_db = {
    "BACC" : {
        "model"  : "BACC",
        "name"   : "Thumper",
        "rate"   : (19.423, "su"), # Per second (found under Propellant)
        "amount" : (820, "su"),
        "ispsl"  : (175, "s"), # Isp at sea level
        "ispvac" : (210, "s")  # Isp at vacuum
    },
    "F3S0" : {
        "model"  : "F3S0",
        "name"   : "Shrimp",
        "rate"   : (1.897, "su"), # Per second (found under Propellant)
        "amount" : (90, "su"),
        "ispsl"  : (190, "s"), # Isp at sea level
        "ispvac" : (215, "s")  # Isp at vacuum
    },
    "S2-17" : {
        "model"  : "S2-17",
        "name"   : "Thoroghbred",
        "rate"   : (100.494, "su"),
        "amount" : (8000, "su"),
        "ispsl"  : (205, "s"),
        "ispvac" : (230, "s")
    },
    "RT-5" : {
        "model"  : "RT-5",
        "name"   : "Flea",
        "rate"   : (15.821, "su"),
        "amount" : (140, "su"),
        "ispsl"  : (140, "s"),
        "ispvac" : (165, "s")
    },
    "RT-10" : {
        "model"  : "RT-10",
        "name"   : "Hammer",
        "rate"   : (15.827, "su"),
        "amount" : (375, "su"),
        "ispsl"  : (170, "s"),
        "ispvac" : (195, "s")        
    }
}

nosecone_db = {
    
}

#
# Prototype code
#

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


class DragDivergence( functor.Functor ) :
    # Original defaults by eye - looking at graphs on the internet: 1.0, 80.0, 0.15
    # 2.98903399 14.56279572 0.12133894 anc  test
    # 3.29937834 11.29413372 0.09386697 ancA test
    # 3.11236594  5.35477702 0.12315374 chj3 test
    def __init__( self, cpeak=2.98903399, cdiv=14.56279572, ctail=0.12133894 ) :
        rangemin = [0.0]
        rangemax = [20.0]
        
        super().__init__( rangemin, rangemax )

        # Divergence peak
        self.c0 = cpeak
        
        # Divergence coeff
        self.c1 = cdiv
        
        # Tail coeff
        self.c2 = ctail

    def call( self, *X ) :
        tr = 1.0
        
        self.checkRange( X )
        
        M = X[0]
        if M <= tr :
            y = 1.0 + self.c0*math.exp( -self.c1*math.pow(M - tr, 2.0) )
        else :
            y= 1.0 + self.c0*math.exp( -self.c2*math.pow(M - tr, 2.0) )
            # y= 1.0 + self.c0*math.exp( - self.c2*(M - tr) )
            # y = 1.0 + self.c0*math.exp( -self.c1*math.pow(M - tr, 2.0) )
        # print("C0 = %s C1 = %s C2 = %s M = %s DD = %s" % (self.c0, self.c1, selM, y))
        return [y]
        
    def nDep( self ) :
        return 1
    
dd = DragDivergence()

#
# Stage and FlyingStage classes used to model craft for
# makePolarMotionSolver.
#

class Stage :
    '''Model of an isolated stage.

    :param float m0: mass of ready stage.  If None, assumes init via loadJSON().
    :param list engines: [ (N1, "Name of engine 1"), etc. ]
    :param float dragco: currently lumped term of 1/2*Cd*A
    '''
    
    def __init__( self, m0 = None, engines = [], dragco = 0.0 ) :
        self.m0 = m0
        self.engines = engines
        self.dragco = dragco

        if m0 is not None :
            self._assimilate( )

    def _assimilate( self ) :
        '''(private) Perform all pre-processing based on input parameters.

        Call this function when top level input params change.
        '''
        self._processEngines()
        
    def _processEngines( self ) :
        '''(private) Preprocess engine data.'''

        engines_by_t_burn = []

        for eng_rec in self.engines :

            if len(eng_rec) == 3 :
                n_eng, ename, fuel_mass = eng_rec
            elif len(eng_rec) == 2 :
                n_eng, ename = eng_rec
                fuel_mass = None

            engine = engine_db[ename]

            if fuel_mass is None:
                # Get from database
                mfuel, mfuelu = engine["amount"]
            else :
                mfuel, mfuelu = fuel_mass
                
            mfuel_kg = mfuel * uconv( mass_db, mfuelu, "kg" )

            rate, rateu = engine["rate"]
            rate_kgps = rate * uconv( mass_db, rateu, "kg" )

            t_burn = mfuel_kg / rate_kgps

            engines_by_t_burn.append( (t_burn, n_eng, engine) )

        engines_by_t_burn.sort()
        engines_by_t_burn.reverse()

        self.mass_nodes = [0.0]
        self.rate_edges = []
        rate_kgps_sum = 0.0
        rate_rec = []

        #
        #       |o M_full
        #       | .
        #       |  . dM/dt = rate_i : all engines running
        #       |   .
        # Mfuel |    o Mi
        #       |    |   . dM/dt = rate0 : longest burning engines running
        #       |    |       .
        #       +----|-----------o M0-------------> time
        #       |DTi |    DT0    |
        #
        #  Build a rate vs mass function using the notions in the plot above.
        #
        #  rate record: [ (<engine group max mass rate>, <engine group ASL Isp>, <engine group Vac Isp>), ... ]
        for irec, eng_rec in enumerate(engines_by_t_burn) :
            
            t_burn, n_eng, engine = eng_rec
            
            # Note: engines list is sorted by *descending* burn time
            try :
                t_burn_prev, x, x = engines_by_t_burn[irec+1]
            except :
                t_burn_prev = 0.0
                
            DT = t_burn - t_burn_prev # DTi in the plot

            rate, rateu = engine["rate"]
            rate_kgps = n_eng * rate * uconv( mass_db, rateu, "kg" )
            rate_kgps_sum += rate_kgps

            ispsl, ispslu = engine["ispsl"]
            ispsl_mps = ispsl * uconv( isp_db, ispslu, "m/s" )

            ispvac, ispvacu = engine["ispvac"]
            ispvac_mps = ispvac * uconv( isp_db, ispvacu, "m/s" )

            rate_rec.append( (rate_kgps, ispsl_mps, ispvac_mps) )

            self.rate_edges.append( copy.copy(rate_rec) )
            self.mass_nodes.append( self.mass_nodes[-1] + DT*rate_kgps_sum )

        m0, m0u = self.m0
        m0_kg = m0 * uconv( mass_db, m0u, "kg" )
        
        self.m0_kg = m0_kg
        # Empty stage mass (kg)
        self.me_kg = m0_kg - self.mass_nodes[-1]
        
    def dmdt( self, mstage_kg, throttle ) :
        '''Compute stage dm/dt vs firing phase (via stage mass) and throttle level

        :param float mstage_kg: mass of stage in kg
        :param float throttle: 0 - 1
        '''

        # This might happen with roundoff error or solver overshoot.
        if mstage_kg <= self.me_kg :
            return 0.0
        
        mfuel_kg = mstage_kg - self.me_kg

        # Find largest left-bracketing record
        for i in reversed(range(len(self.rate_edges))) :
            if mfuel_kg >= self.mass_nodes[i] :
                rec = self.rate_edges[i]
                break
            
        rate_kgps_sum = 0.0
        for rate_kgps, isl, iva in rec : # rate, Isp-ASL, Isp-Vac
            rate_kgps_sum += throttle * rate_kgps

        return -rate_kgps_sum

    def dumpInfo( self ) :

        print()
        
        tabrows = [ [ "Initial mass (kg)", self.m0_kg ],
                    [ "Empty mass (kg)", self.me_kg] ]
        print(tabulate.tabulate(tabrows, headers=["Parameter", "Value"]))
        print()
        
        tabrows = [[self.mass_nodes[0], None]]
        for i in range(len(self.rate_edges)) :
            tabrows.append( [ None, self.rate_edges[i] ] )
            tabrows.append( [ self.mass_nodes[i+1], None ] )
        print(tabulate.tabulate(tabrows, headers=["Stage Mass", "Rate Record"]))
        print()

    def dvRemain( self, mstage_kg, p_Pa ) :
        '''Compute remaining Delta-V given firing phase (via stage mass) and ambient pressure.

        :param float mstage_kg: Stage mass in kg
        :param float p_pA: ambient pressure in Pascals
        '''
        ### See rocket book note "RB 2021-01-30 04.15.00.pdf"

        if mstage_kg <= self.me_kg :
            return 0.0
        
        vac = (C_p0 - p_Pa) / C_p0

        # Mass of vehicle (at fuel mass nodes)
        MSTAGE = [ self.me_kg + m for m in self.mass_nodes ]
        # One for each mass node
        DV = [0.0]
        
        for i in range(1, len(MSTAGE)) :
            
            rate_rec = self.rate_edges[i-1]

            sum_rate_kgps = 0.0
            sum_thrust = 0.0

            for rate_kgps, isl, iva in rate_rec :
                isp = (1.0 - vac)*isl + vac*iva
                sum_rate_kgps += rate_kgps
                sum_thrust += (rate_kgps * isp)

            isp_eqv = sum_thrust / sum_rate_kgps

            dv_phase = isp_eqv * math.log( MSTAGE[i] / MSTAGE[i-1] )
            
            DV.append( dv_phase + DV[i-1] )

            # If input stage mass falls within the current firing
            # phase, then interpolate the answer and we're done.
            if mstage_kg >= MSTAGE[i-1] and mstage_kg <= MSTAGE[i] :
                a = ( mstage_kg - MSTAGE[i-1] ) / ( MSTAGE[i] - MSTAGE[i-1] )
                return DV[i-1]*(1.0-a) + DV[i]*a

        # If we get here then it probably means this is a later
        # (ligther) stage being queried for its available DV.
        return DV[-1]

    def thrust( self, mstage_kg, throttle, p_Pa ) :
        '''Thrust of stage given firing phase (via stage mass), throttle, and ambient pressure.

        :param float mstage_kg: Stage mass in kg
        :param float throttle: 0 - 1
        :param float p_Pa: ambient pressure in Pascals
        '''
        
        if mstage_kg < self.me_kg :
            return 0.0

        vac = (C_p0 - p_Pa) / C_p0
        
        mfuel_kg = mstage_kg - self.me_kg

        # Find largest left-bracketing record
        for i in reversed(range(len(self.rate_edges))) :
            if mfuel_kg >= self.mass_nodes[i] :
                rec = self.rate_edges[i]
                break
        
        thrust = 0.0
        for rate_kgps, isl, iva in rec : # rate, Isp-ASL Isp-Vac
            thrust += throttle * rate_kgps * ((1.0-vac)*isl + vac*iva)
                
        return thrust
    
    def loadJSON( self, fname ) :
        '''Load initial parameters from a JSON file and pre-process.
        '''

        with open( fname, "rt" ) as f :
            data = json.load( f )
            self.m0 = data["m0"]
            self.engines = data["elist"]
            self.dragco = data["dragco"]
            self._assimilate( )        
        
    def dumpJSON( self, fname ) :
        '''Save initial parameters to a JSON file.
        '''
        
        data = { "m0" : self.m0,
                 "elist" : self.engines,
                 "dragco" : self.dragco
                 }
        
        with open( fname, "wt" ) as f :
            json.dump( data, f )

class FlyingStage :
    '''A Stage in the context of a body which provides functions needed for plane-polar trajectory.

    The main inputs are the throttle function and the orientation function.

    :param Stage stage: A Stage object
    :param string stage_name: Arbitrary name for flying stage
    :param string body_name: name of solar system body about which stage is flying
    :param function fthrottle: function( t, y ) -> 0-1
    :param function falpha: function( t, y, FlyingStage ) -> angle relative to body normal, RHR applies
    '''
    def __init__( self, stage, stage_name, body_name, fthrottle, falpha ) :
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
        self.fdens  = self.body["fdens"]
        self.fsnd   = self.body["fsnd"]
 
        # Body rotational period
        tday, tday_u = self.body["tday"]
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
        
        return ( ( self._thrust(t, y) * math.sin( self.falpha(t, y, self) ) / m )
                 - self.GM/(r*r) - self.stage.dragco * dd.call(M)[0] * abs(vr) * vr * self.fdens.call(alt)[0]/m )

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
        return ( ( self._thrust(t, y) * math.cos( self.falpha(t, y, self) ) / m )
                 - self.stage.dragco * dd.call(M)[0] * abs(vthb) * vthb * self.fdens.call(alt)[0]/m )
    
    def _dmdt( self, t, y ) :
        '''(private) Compute dmdt

        :param float t: trajectory time
        :param list y: solver dependent values array (See makePolarMotionSolver)
        '''
        
        m, r, th, vr, om = y
        return self.stage.dmdt( m, self.fthrottle(t, y) )

    def _thrust( self, t, y ) :
        '''(private) Compute thrust at a trajectory point

        :param float t: trajectory time
        :param list y: solver dependent values array (See makePolarMotionSolver)
        '''
        
        m, r, th, vr, om = y
        alt = r - self.R
        p = self.fpress.call( alt )[0]
        return self.stage.thrust( m, self.fthrottle(t, y), p )

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
            
            Y, crashed, flyer = self.flyTo( t )

            if crashed :
                break
            
            m, r, th, vr, om = Y

            om_gnd = om - flyer.body_omega
            
            # Compute other things
            x = ( r * math.cos( th ) )
            y = ( r * math.sin( th ) )

            craft_asl_dvremain = flyer.dvRemain( m, C_p0 )
            craft_dvremain = flyer.dvRemain( m, flyer.fpress.call(r-flyer.R)[0] )
            stage_dvremain = flyer.stage.dvRemain( m, C_p0 )

            vom_gnd = r*om_gnd
            spd_gnd = math.sqrt( vr*vr + vom_gnd*vom_gnd )

            vom = r*om
            spd = math.sqrt( vr*vr + vom*vom )
            
            row = [ flyer.stage_name, t, m, r, th, vr, om, om_gnd,
                    stage_dvremain, craft_asl_dvremain, craft_dvremain, spd_gnd, spd ]

            rowdat.append(row)

            t += dt

        headers = [ "stage", "time", "mass", "r", "theta", "v_r", "omega",
                    "rel omega", "Stage DV (ASL)", "Craft DV (ASL)",
                    "Craft DV", "Ground speed", "Orbit Speed" ]
        
        print( tabulate.tabulate( rowdat, headers = headers) )

    def dvRemain( self, mcraft_kg, p_Pa ) :
        '''Compute remaining Delta-V given firing phase (via craft mass) and ambient pressure.

        :param float mcraft_kg: Stage mass in kg
        :param float p_pA: ambient pressure in Pascals
        '''
        
        mydv = self.stage.dvRemain( mcraft_kg, p_Pa )

        if self.sp1 is not None :
            return mydv + self.sp1.dvRemain( mcraft_kg, p_Pa )
        else :
            return mydv        

    def flyTo( self, t ) :
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
                # print( "WARNING: FlyingStage is crashed at time %f" % t )
                return ( copy.copy(self.soln[-1]), True, self )

        y = interp2Dtraj( t, self.solnt, self.soln )
        
        return ( y, False, self )

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
        
        self.solv = makePolarMotionSolver( dmdt, a_r, a_th, y0, t0 )

        self.solnt = [ t0 ]
        self.soln = [ y0 ]
        self.maxr = 0.0
        
    def plot( self, t0 = None, t1 = None, dt = 5.0 ) :
        ''' Plots a trajectory.

        :param float t0: time to start trajectory dump.  If None, will
                         start with the stage for which it is called.
        :param float t1: time to end trajectory dump.  If None, will
                         end with the last time stamp of the current
                         stage.
        :param float dt: time interval between dumped trajectory
                         points.  Default is 5.0 seconds.
        '''

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
    '''Looks for dictionary (JSON) files named <Body Key>.json and adds data to the bodies_db dictionary'''
    
    for body_key in bodies_db :
        fname = pthdb( "%s.json" % body_key )
        exists = os.access( fname, os.F_OK )
        if exists :
            with open( fname, 'rt' ) as f :
                extra_dat = json.load(f)
                bodies_db[body_key].update( extra_dat )

def processBodyDbs( ) :
    '''Pre-processes bodies_db data.
    '''
    
    for body_key in bodies_db :
        bodyrec = bodies_db[body_key]
        if "apt" in bodyrec :
            alts, ps, ts = zip( *bodyrec["apt"] )

            dens = [ C_air_mol_mass*p/ts[i]/C_Rgas for i,p in enumerate(ps) ]
            snd  = [ math.sqrt(142E3 / di) for di in dens ]

            bodyrec["fdens"] = functor.Interp1DFunctor( alts, dens, kind="quadratic" )
            bodyrec["fpress"] = functor.Interp1DFunctor( alts, ps, kind="quadratic" )
            bodyrec["fsnd"] = functor.Interp1DFunctor( alts, snd, kind="quadratic", low_fill = snd[0], high_fill=1E40 )

def interp2Dtraj( t, solnt, soln ) :
    '''Utility function for makePolarMotionSolver.
    
    :param float t: query time of point along trajectory
    :param list solnt: array of times corresponding to each element of soln.
    :param list soln: array of makePolarMotionSolver y.

    :Rationale: bisect is used in case a non-uniform time interval is used.
    :TODO: consider using the scipy interpolator.
    '''
    
    i = (bisect.bisect( solnt, t ) - 1)
    t0 = solnt[i]
    t1 = solnt[i+1]
    a = (t - t0)/(t1 - t0)
    y0 = soln[i]
    y1 = soln[i+1]
    y = [ (1.0 - a)*y0[j] + a*y1[j] for j in range(len(y0)) ]
    
    return y

def makePolarMotionSolver( dmdt, a_r, a_th, y0, t0 ) :
    '''Creates an ODE solver that numerically integrates the planar orbit equations as a function of time.
    
    Dependent variable list:

        y := [ m(kg), r(m), th(radians), vr(m/s), om(rad/s) ]

    :param function dmdt: function( t, y ) -> ks/s
    :param function a_r:  function( t, y ) -> acceleration along r axis
    :param function a_th: function( t, y ) -> acceleration along th axis
    :param list y0:       initial conditions
    :param float t0:      starting time (s)

    :returns: scipy.integrate.ode object.
    '''

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
            
def dInterp(term, dunit="m") :
    '''General distance interpreter.

    Allowed terms are:

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

def dvHohmannApo(body, r_peri, r_apo) :
    '''DV computed when thrusting a) prograde at apo to attain circular orbit, or b) retrograde from apo to attain elliptical orbit.

    :param (d,"du") r_peri: periapsis
    :param (d,"du") r_apo:  apoapsis

    :returns: float Delta-V for maneuver
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

    :param (d,"du") r_peri: periapsis
    :param (d,"du") r_apo:  apoapsis

    :returns: float Delta-V for maneuver
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

def listOfBodyNames(sep = " | ") :
    body_names = bodies_db.keys()
    return sep.join(body_names)

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

def uconv(db, ufrom, uto) :
    '''Computes multiplication factor to convert value of units "ufrom" to value of units "uto"

    :param string db: name of units db
    :param string ufrom: units converting from
    :param string uto: units converting to

    :returns: float conversion factor
    '''
    return db[ufrom] / db[uto]

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

    # Command "fly"
    cmd_fly = subparsers.add_parser("fly", help="Fly a craft")
    cmd_fly.add_argument("body", help="Name of body. Available: [%s]" % listOfBodyNames())
    cmd_fly.add_argument("craft", help="Name of craft")

    # Command "plotfuncs"
    cmd_plotfuncs = subparsers.add_parser("plotfuncs", help="Plot up the model functions")

    # Command "uconv"
    cmd_uconv = subparsers.add_parser("uconv", help="Call the units converter function")
    cmd_uconv.add_argument("vin", help="\"(value, 'unit')\"")
    cmd_uconv.add_argument("uout", help="Unit out")

    args = parser.parse_args()

    print ( dd.call(7.259003981175488e-05) )
    
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

    if args.command == "fly" :
        stage = Stage()
        stage.loadJSON( args.craft )
        
        def fthrottle( t, y ) :
            return 1.0
        def falpha( t, y, flyer ) :
            return 0.5 * math.pi

        flyer = FlyingStage( stage, "Stage 1", args.body, fthrottle, falpha )
        flyer.launch( )

        flyer.dumpTraj(t1 = 30000, dt = 1.0)
        
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


if __name__ == "__main__" :
    main()
