# Copyright © 2021 Jonathan Grot

'''Low level polar motion solvers where all units are assumed to be in MKS units.
'''

import math
from scipy.integrate import ode
from scipy.optimize import minimize

import compare as cmp
import moremath as mm
import mpl_tools as mpt

def dv_r1_hohmann(GM, r0, r1) :
    '''DV computed when firing a) prograde at apo (r1) to attain circular orbit, or b) retrograde from (r1) to attain elliptical orbit.

    :param float GM: G*M for celestial body (m^3/s^2)
    :param float r0: periapsis distance from celestial body center (m)
    :param float r1: apoapsis distance from celestial body center (m)

    :returns: float Delta-V (m/s) for maneuver
    '''
    
    A = math.sqrt(GM/r1)
    B = 1.0 - math.sqrt(2.0*r0/(r0 + r1))
    return A*B
    
def dv_r0_hohmann(GM, r0, r1) :
    '''DV computed when firing a) prograde at peri (r0) to attain elliptical orbit, or b) retrograde at peri (r0) to attain circular orbit.

    :param float GM: G*M for celestial body (m^3/s^2)
    :param float r0: periapsis distance from celestial body center (m)
    :param float r1: apoapsis distance from celestial body center (m)

    :returns: float Delta-V (m/s) for maneuver
    '''
    
    A = math.sqrt(GM/r0)
    B = math.sqrt(2.0*r1/(r0 + r1)) - 1.0
    return A*B

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

class Orbit :
    '''Parameterized orbit under a gravitational force initialized by a
       trajectory point anywhere along the orbit.

    Parameterizes by "phi" such that phi=0 is at the periapsis.

    self.phi0 represents phi at the initialization point.

    :param list y: makePolarMotionSolver variables [ m(kg), r(m), th(radians), vr(m/s), om(rad/s) ]
    :param float GM: body GM

    Some orbit equations:

    .. math::

       k = G M m

       h = r^2 \omega = r v_{\phi} = r_0 v_0 = r_1 v_1

       E = { m v^2 \over 2 } - { k \over r } = \mathit{const}

       e = \sqrt{1 + 2Emh^2k^{-2}}

       r_0 = { mh^2 \over { k(1+e) } }

       r = r_0 { {1+e} \over { 1 + e \cos \phi } }

       \phi_0 = \mathrm{arccos} \\left ( { 1 \over e} \\left ( {r_0 \over r} (1+e) - 1 \\right ) \\right ); \: \mathrm{if} \: \dot{r} < 0, \; \phi_0 \leftarrow 2\pi - \phi_0

       \\tau = 2 \pi \\left( {m \over k} \\right)^{1/2} a^{3/2}

       a = { {m h^2} \over { k(1-e^2) } }

       \\tau = 2 \pi \\left( 1 \over {GM} \\right)^2 { h^3 \over {(1-e^2)^{3/2}}}
 
       {d\phi \over dt } = h \\left( { 1 + e\cos \phi } \over { r_0 (1+e) } \\right) ^2 = A(1+e\cos\phi)^2; \: A= { h \over { ( r_0(1+e) )^2 } }

    '''

    ORB_CIRCLE = 0
    ORB_ELLIPSE = 1
    ORB_PARABOLA = 2
    ORB_HYPERBOLA = 3
    ORB_ORBITS = ( 'circle', 'ellipse', 'parabola', 'hyperbola' )
    
    def __init__(self, y, GM, force=None) :
        
        m, r, th, vr, om = y

        if force is not None :
            if force == "circle" :
                om = math.sqrt(GM/r)/r
                vr = 0.0
            elif force == "parabola" :
                om = math.sqrt(2.0*GM/r)/r
                vr = 0.0

        k = GM*m
        h = r*r*om
        vsq = vr*vr + h*om
        E = 0.5*m*vsq - (k/r)
        e = math.sqrt(1.0 + 2.0*E*m*h*h/(k*k))
        r0 = m*h*h/(k*(1.0 + e))
        v0 = h/r0
        if e < 1.0 :
            r1 = r0*(1.0 + e)/(1.0 - e)
        else :
            r1 = math.inf

        self.GM = GM
        self.m = m
        self.k = k
        self.h = h
        self.E = E
        self.e = e
        self.r0 = r0
        self.v0 = v0
        self.r1 = r1

        # PHI0
        if e > 0.0 :
            # Not a circle
            x1 = 1.0 / e
            x2 = (r0/r)*(1.0 + e)
            x3 = x1*(x2 - 1.0)
            if x3 >= 1.0 :
                phi0 = 0.0
            elif x3 <= -1.0 :
                phi0 = -math.pi
            else :
                phi0 = math.acos(x3)
                if vr < 0.0 :
                    phi0 = 2.0*math.pi - phi0
            
            self.phi0 = phi0
        else :
            # Circle
            self.phi0 = 0.0

        # PERIOD
        if e < 1.0 :
            x1 = 1.0/GM/GM
            x2 = math.pow(h,3.0)/math.pow((1.0 - e*e),1.5)
            self.tau = 2.0*math.pi*x1*x2


        # Set up time integrator for phi.  Works for all orbits.
            
        A = h / math.pow(r0*(1.0 + e), 2.0)

        def f(t, y) :
            x1 = 1.0 + e*math.cos(y[0])
            ans = A*x1*x1
            return [ ans ]

        solv = ode(f).set_integrator("vode", method="adams")

        orbtype = self.classify()
        
        if orbtype in [self.ORB_CIRCLE, self.ORB_ELLIPSE] :
            
            solv.set_initial_value( 0.0, 0.0 )

            self.phivt = [ [0.0] ]
            self.solnt = [ 0.0 ]

            dt = self.tau / 1000.0
            phi = 0.0
            
            while solv.successful() and phi < 2.0*math.pi :
                phi = solv.integrate( solv.t + dt )

                self.phivt.append( list(phi) )
                self.solnt.append( solv.t )
                
        elif orbtype in [self.ORB_HYPERBOLA, self.ORB_PARABOLA] :
            
            self.phi_max = math.acos(-1.0/self.e)
            self.phi_min = -self.phi_max

            solv.set_initial_value( self.phi0, 0.0 )

            self.phivt = [ [self.phi0] ]
            self.solnt = [ 0.0 ]

            dt = 0.1
            phi = 0.0
            
            while solv.successful() and solv.t < 1000.0 :
                phi = solv.integrate( solv.t + dt )

                self.phivt.append( list(phi) )
                self.solnt.append( solv.t )

    def classify( self ) :
        if self.e < 0 :
            raise Exception("Logic error: e < 0")
        if self.e == 0 :
            return self.ORB_CIRCLE
        if self.e < 1 :
            return self.ORB_ELLIPSE
        if self.e == 1 :
            return self.ORB_PARABOLA
        return self.ORB_HYPERBOLA
        
    def intersect_soi(self, phi_soi, d_soi, r_soi) :
        '''Intersect this orbit with a sphere of influence.  Finds all
        intersections, so it's up to the caller to sort out what is
        relevant.

        Procedure: starting at the leading edge (angle) of the SOI,
        looks for alternating minima and maxima, putting the starting
        simplex of the Nelder-Mead optimization algorithm, just a hair
        past the last extrema.  Once the optimizer falls off the SOI,
        we're done.  Testing has shown that a glancing intersection
        will generally count as one intersection, so it is possible to
        have anywhere from one to four intersections.

        :param float phi_soi: angle of line to SOI in the orbit frame.
        :param float d_soi: distance of new SOI from center of current SOI body.
        :param float r_soi: radius of new SOI.

        '''

        # Angle from center of SOI to edge of SOI
        a = math.asin(r_soi/d_soi)

        phistop_1 = phi_soi - 1.5*a
        phistop_2 = phi_soi + 1.5*a
        
        def objf(X, sign) :
            phi = X[0]
            
            # Searching for maxima will get stuck at these extremes.
            if phi < phistop_1 :
                return math.inf
            if phi > phistop_2 :
                return math.inf
            
            r = self.y_phi(phi)[1]
            # Bad value of phi for the orbit
            if r is None :
                return math.inf

            # Distance of point from SOI center
            rs = mm.r_soi(r, phi, d_soi, phi_soi)

            # Objective function looking for minima if sign = +1, maxima if sign = -1
            fr = sign*math.pow((rs - r_soi),2.0)
            
            return fr

        # Only care about the mimina
        minima = []

        # Leading boundary of SOI
        phia = phi_soi - a
        phib = phia + .01*a

        sign = 1.0
        while True :
            # Start the search for the next extremum
            init_simplex = [[phia], [phib]]
            result = minimize(objf, None, args=sign, method="Nelder-Mead", options={"initial_simplex":init_simplex, "fatol": 1E-8})

            fobj = result.fun
            phi = result.x[0]
            
            # If no close call for a minimum, then there is no
            # intersection.  result.fun is the value of the objective
            # function.
            if sign > 0.0 and math.sqrt(fobj) > 1E-8*r_soi :
                return minima

            if ( phi < (phi_soi - a) or phi > (phi_soi + a) ) :
                return minima
            
            if sign > 0.0 :
                minima.append(phi)

            # Advance the starting simplex to over the next hump
            phia = phi
            phib = phi + 0.01*a

            # Look for the opposite type of extremum
            sign *= -1.0

    def phi_t(self, t) :
        return mm.bisect_interp(t, self.solnt, self.phivt)[0]
    
    def plot(self, ax, use_t=False, Rbody=None) :
        orbtype = self.classify()
        
        if orbtype in [self.ORB_HYPERBOLA, self.ORB_PARABOLA] :
            if use_t :
                phi_of_t = lambda t: self.phi_t(t)
                r_of_t = lambda t : self.y_phi( self.phi_t(t) )[1]
                dt = self.solnt[-1]/100.0
                maxt = 99.0*dt
                mpt.plot_polar(ax, r_of_t, phi_of_t, 0.0, dt, maxt, centered=True, marker="o")
            else :
                phi_of_phi = None
                r_of_phi = lambda phi: self.y_phi(phi)[1]
                maxphi = 0.9*self.phi_max
                dphi = (maxphi - self.phi0)/100.0
                mpt.plot_polar(ax, r_of_phi, phi_of_phi, self.phi0, dphi, maxphi, centered=True, marker="o")
                
        elif orbtype in [self.ORB_CIRCLE, self.ORB_ELLIPSE] :
            if use_t :
                phi_of_t = lambda t: self.phi_t(t)
                r_of_t = lambda t : self.y_phi( self.phi_t(t) )[1]
                dt = self.solnt[-1]/100.0
                maxt = 99.0*dt
                mpt.plot_polar(ax, r_of_t, phi_of_t, 0.0, dt, maxt, centered=True, marker="o")
            else :
                phi_of_phi = None
                r_of_phi = lambda phi: self.y_phi(phi)[1]
                dphi = 2.0*math.pi/100.0
                maxphi = dphi*99.0
                mpt.plot_polar(ax, r_of_phi, phi_of_phi, 0.0, dphi, maxphi, centered=True, marker="o")

        if Rbody is not None :
            mpt.plot_circle(ax, Rbody, 100)
            

    def r0_dv_to_e(self, e) :
        '''Computes *signed* delta V to change eccentricity at r0

        * This works for any e.
        * Returns signed value, so take the absolute value when making DV maps.

        Use Hohmann transfer for elliptical orbits.
        '''
        v0_new = math.sqrt(self.GM*(1.0 + e)/self.r0)
        return (v0_new - self.v0)

    def y_phi( self, phi ) :
        '''Returns a polar motion state vector.
        '''
        
        num = self.r0*(1.0 + self.e)
        denom = 1.0 + self.e*math.cos(phi)

        if denom <= 0.0 :
            return [None, None, None, None, None]
        
        r = num / denom

        vphi = self.h / r
        om = vphi / r
        
        vsq = (2.0/self.m) * (self.E + self.k/r)

        vrsq = vsq - vphi*vphi

        if vrsq < 0.0 :
            print("WARNING VR^2 IS < 0.  HOPEFULLY, JUST A ROUNDING ERROR: ", vrsq, " PHI ", phi)
            vrsq = 0.0
        
        vr = math.sqrt( vrsq )

        if phi > math.pi :
            vr = -vr
        
        return [ self.m, r, phi, vr, om ]
