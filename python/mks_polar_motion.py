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

def rth_at_xy(x, y) :
    '''Converts from Cartesian to polar coordinates x,y -> r,th

    See also: xy_at_rth()
    '''

    r = math.sqrt(x*x + y*y)
    th = math.atan2(y,x)
    
    return (r, th)
    
def v_circular_orbit(GM, r) :
    '''Returns v for a circular orbit of radius r.
    '''
    return math.sqrt(GM/r)

def v_and_dir(y) :
    '''Given polar trajectory point y, compute velocity direction vector
    in polar coordinates'''
    
    m, r, th, vr, om = y

    vth = r*om

    v = math.sqrt(vr*vr + vth*vth)

    return (v, (vr/v, vth/v))
    
def vrth_at_th(v_x, v_y, th) :
    '''Converts a vector in Cartesian coordinates to polar coordinates.

    The Cartesian coordinate system is defined such that th is wrt to the x axis.

    :param float v_x: a vector x-component
    :param float v_th: a vector y-component
    :param float th: (radians) the angle of the location of the point at which the vector originates
    '''
    v_r  =  v_x*math.cos(th) + v_y*math.sin(th)
    v_th = -v_x*math.sin(th) + v_y*math.cos(th)

    return (v_r, v_th)
    
def vxy_at_th(v_r, v_th, th) :
    '''Converts a vector in polar coordinates to Cartesian coordinates.

    The Cartesian coordinate system is defined such that th is wrt to the x axis.

    :param float v_r: a vector component in the r direction
    :param float v_th: a vector component in the th direction
    :param float th: (radians) the angle of the location of the point at which the vector originates
    '''
    v_x = v_r*math.cos(th) - v_th*math.sin(th)
    v_y = v_r*math.sin(th) + v_th*math.cos(th)

    return (v_x, v_y)

def xy_at_rth(r, th) :
    '''Polar coordinates r,th -> x,y

    See also rth_at_xy()
    '''
    x = r*math.cos(th)
    y = r*math.sin(th)
    
    return (x, y)

class Frame2D :
    '''A 2D reference frame

    :param Frame2D parent: A parent Frame2D object, or None

    A parent frame of None implies the "grand" coordinate system.

    Currently assumes all reference frames are axis aligned:

    - Velocities in x,y don't need transforming
    - To transform a velocity (vx,vy) to polar in the ref frame:
      1. Call r,th = frm.xform_pos_to('c'|'p', x|r, y|th, 'p') where (vx,vy) is the velocity at the pos
      2. Call vfx, vfy = frm.xform_vel_to(vx, vy)
      2. Call vr, vth = vrth_at_th(vfx, vfy, th)

    '''

    def __init__(self, parent=None):
        self.parent = parent
        # Coordinates and velocity in parent frame
        self.x = None
        self.y = None
        self.r = None
        self.th = None
        self.vx = None
        self.vy = None

    def set_rth(self, r, th) :
        self.r = r
        self.th = th
        self.x, self.y = xy_at_rth(r, th)
        
    def set_xy(self, x, y) :
        self.x = x
        self.y = y
        self.r, self.th = rth_at_xy(x, y)

    def set_vrth(self, vr, vth) :
        '''Given vr and vth in the parent frame, computes vx and vy in the parent frame'''
        if self.th is None :
            raise Exception("Must first set position of reference frame.")
        
        self.vx, self.vy = vxy_at_th(vr, vth, self.th)

    def set_vxy(self, vx, vy) :
        self.vx = vx
        self.vy = vy

    def xform_pos_to(self, c_in, z1, z2, c_out) :
        '''Transforms position from parent to this frame

        :param string c_in: Input coordinate sys: 'p' - polar, 'c' - Cartesian
        :param float z1: First component of input coordinate
        :param float z2: Second component of input coordinate
        :param string c_out: Output coordinate sys: 'p' - polar, 'c' - Cartesian
        '''

        if self.parent is not None :
            raise Exception("TODO: incorporate parent frame")

        # Convert to parent x,y

        if c_in == 'p' :
            x_in, y_in = xy_at_rth(z1, z2)
        else :
            x_in = z1
            y_in = z2

        # Convert to this ref frame
        
        x = x_in - self.x
        y = y_in - self.y

        # Return requested coordinates

        if c_out == 'p' :
            z1_out, z2_out = rth_at_xy(x, y)
        else :
            z1_out = x
            z2_out = y

        return (z1_out, z2_out)

    def xform_vel_to(self, vx, vy) :
        if self.vx is None :
            raise Exception("Must set frame velocity")
        
        vfx = vx - self.vx
        vfy = vy - self.vy

        return (vfx, vfy)
                     
class Orbit :
    '''Parameterized orbit under a gravitational force initialized by a
       trajectory point anywhere along the orbit.

    :param list y: makePolarMotionSolver variables [ m(kg), r(m), th(radians), vr(m/s), om(rad/s) ]
    :param float GM: body GM

    :Conventions:

    phi is the polar angle in the coordinate system of the orbit with r0 at phi=0

    phi_i is the orbit angle of the "intertion point", i.e. the point
    which is used to initialize the orbit

    :See also: OrientedOrbit

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
    
    def __init__(self, y, GM, phi_resolution=5000.0, force=None) :

        self.phi_resolution = phi_resolution

        # Note: th is ignored here, but used in OrientedOrbit
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
        x1 = 1.0 + 2.0*E*m*h*h/(k*k)
        try :
            e = math.sqrt(x1)
        except :
            if abs(x1) < 1E-14 :
                e = 0
            else :            
                raise Exception("Attempt to take sqrt(%e)" % x1)
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

        # Compute phi_i (-pi <= phi_i <= pi)
        
        if e > 0.0 :
            # Not a circle
            x1 = 1.0 / e
            x2 = (r0/r)*(1.0 + e)
            x3 = x1*(x2 - 1.0)
            if x3 >= 1.0 :
                phi_i = 0.0
            elif x3 <= -1.0 :
                phi_i = -math.pi
            else :
                phi_i = math.acos(x3)
                if vr < 0.0 :
                    phi_i = mm.home_angle(2.0*math.pi - phi_i)
            
            self.phi_i = phi_i

        else :
            # Circle
            self.phi_i = 0.0

        # PERIOD OR PHI_MAX
        
        if e < 1.0 :
            
            x1 = 1.0/GM/GM
            x2 = math.pow(h,3.0)/math.pow((1.0 - e*e),1.5)
            self.tau = 2.0*math.pi*x1*x2
            self.phi_max = math.inf
            self.phi_min = -math.inf
            
        else :

            self.phi_max = math.acos(-1.0/self.e)
            self.phi_min = -self.phi_max

            if self.phi_i > self.phi_max :
                raise Exception("Logic Error:  phi_i > phi_max")
            if self.phi_i < self.phi_min :
                raise Exception("Logic Error:  phi_i < phi_min")
            
        # Set up time integrator for phi.  Works for all orbits.
            
        A = h / math.pow(r0*(1.0 + e), 2.0)
        def f(t, y) :
            x1 = 1.0 + e*math.cos(y[0])
            ans = A*x1*x1
            return [ ans ]
        self.solv = ode(f).set_integrator("vode", method="adams")

        # t=0 -> phi=0
        self.solv.set_initial_value(0.0, 0.0)

        # Compute phi vs time with t=0 at phi=0 (r=r0)

        self.path_phi = [0.0]
        self.path_t = [0.0]
        self.orbtype = self.classify()
        
    def _to_path_param(self, x_in, which) :
        '''Converts t_in or phi_in to periodic path values -p/2 <= x <= p/2 where p is period 

        :param float x_in: input t or phi
        :param int which: 0=t 1=phi
        '''

        if self.e >= 1 :
            raise Exception("Called for non-periodic orbit")

        if which == 0 :
            p = self.tau
        else :
            p = 2.0 * math.pi

        # Only need to work out the math on the + side
        x2 = abs(x_in)
            
        if x2 > p :
            f,i = math.modf(x2/p)
            x2 = f*p
            
        if x2 > 0.5*p :
            x2 = x2 - p

        # Now account for negative x_in
        if x_in < 0 :
            x2 = -x2

        return x2
            
    def _propagate_t_phi(self, to, which) :
        '''Extends the phi vs t arrays out to "to"

        - Does not take into account the type of orbit, only
          propagates.  It is up to other code to manage periodicity.

        - Appends to head and tail of path_t and path_phi

        :param float to: value to extend to
        :param int which: 0=t, 1=phi

        '''

        # Only need to deal with the + side of things
        to = abs(to)

        # Check if there is anything to do
        if which == 0 : # t
            if to < self.path_t[-1] : return
        else :
            if to < self.path_phi[-1] : return
                        
        # Extend

        if which == 0 :
            to_last = self.path_t[-1]
        else :
            to_last = self.path_phi[-1]

        # phi resolution = 2pi / self.phi_resolution
        while self.solv.successful() and to_last < to :
            
            phi_last = self.path_phi[-1]
            m, r, phi, vr, om = self.y_at_phi(phi_last)

            dt = 2.0 * math.pi / float(self.phi_resolution) / om

            next_t   = self.solv.t + dt
            next_phi = self.solv.integrate(next_t)[0]

            # print("DEBUG DT, NEXT_T, NEXT_PHI, OM", dt, next_t, next_phi, om)
            
            self.path_t.append(next_t)
            self.path_phi.append(next_phi)
            self.path_t.insert(0, -next_t)
            self.path_phi.insert(0, -next_phi)

            if which == 0 :
                to_last = next_t
            else :
                to_last = next_phi

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

        :returns: list [<intersection phi's>]
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
            
            r = self.y_at_phi(phi)[1]
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

    def phi_at_t(self, t, phi_range=2) :
        '''Returns orbit angle at specified time

        :param float t: time from phi=0, r=r0
        :param int phi_range: For periodic motion, else ignored.  See notes.

        :Notes:

        phi_range for periodic motion:

          - 0 -> -pi to pi
          - 1 -> 0 to 2pi
          - 2 -> -inf to inf
        '''

        # Keep it simple; if t=0, then phi=0
        if t == 0.0 :
            return 0.0

        if self.e < 1.0 : # closed orbit
            t2 = self._to_path_param(t, which=0) # which=0 -> adjust time
        else :
            t2 = t

        self._propagate_t_phi(t2, which=0) # which=0 -> propagate time
        
        phi_out = mm.bisect_interp(t2, self.path_t, self.path_phi)

        if self.e >= 1.0 :
            # For open orbits, just return phi
            return phi_out

        #
        # Periodic orbit
        #
        
        if phi_range == 0 :
            return phi_out

        # phi_range is 1 or 2

        # Solve things on the plus side.  So if t < 0, solve with phi = -phi and then negate at the end.

        sign = 1.0 if t >= 0.0 else -1.0
        
        phi_out = sign*phi_out
        
        if phi_out < 0 :
            phi_out += 2.0*math.pi

        if phi_range == 1 :
            return sign*phi_out

        # phi_range is 2

        # How many time periods out are we (stored in i)?
        f,i = math.modf(sign*t/self.tau)

        # Add that number of phi periods
        phi_out += 2.0*math.pi*i

        return sign*phi_out
    
    
    def r0_dv_to_e(self, e) :
        '''Computes *signed* delta V to change eccentricity at r0

        * This works for any e.
        * Returns signed value, so take the absolute value when making DV maps.

        Use Hohmann transfer for elliptical orbits.
        '''
        v0_new = math.sqrt(self.GM*(1.0 + e)/self.r0)
        return (v0_new - self.v0)

    def sample_phi(self, PHI) :
        '''Generates a series of r values for an input series of polar angles.

        :param iterable PHI: A list of polar angles

        :returns: A list of r values
        '''
        R = []
        for phi in PHI :
            if phi > self.phi_max or phi < -self.phi_max :
                raise Exception("PHI out of range")
            R.append(self.y_at_phi(phi)[1])

        return R

    def sample_t(self, T) :
        '''(UNDER CONSTRUCTION) Generates a series of r, phi values for an input series of time values.

        :param iterable T: A list of times

        :returns: ([<r values>], [<phi values>])

        :todo: make sure to experiment with t > tau (for e<1) and t > solnt[-1] (for e>=1)
        '''
        R = []
        PHI = []

        for t in T :
            phi = self.phi_at_t(t)
            if phi > self.phi_max or phi < -self.phi_max :
                raise Exception("PHI out of range")
            PHI.append(phi)
            R.append(self.y_at_phi(phi)[1])

        return (R, PHI)

    def t_at_phi(self, phi, t_range=2) :
        '''Returns orbit t at specified angle.

        :param float phi: orbit angle
        :param int t_range: For periodic motion, else ignored.  See notes.

        :Notes:

        t_range for periodic motion:

          - 0 -> -tau/2 to tau/2
          - 1 -> 0 to tau
          - 2 -> -inf to inf
        '''
        
        if self.e < 1.0 : # closed orbit
            phi2 = self._to_path_param(phi, which=1) # which=1 -> adjust phi
        else :
            phi2 = phi

        self._propagate_t_phi(phi2, which=1) # which=1 -> propagate phi
        
        t_out = mm.bisect_interp(phi2, self.path_phi, self.path_t)

        if self.e >= 1.0 :
            return t_out
        
        #
        # Periodic orbit
        #

        if t_range == 0 :
            return t_out

        # t_range is 1 or 2
        
        # Solve things on the plus side.  So if phi < 0, solve with t = -t and then negate at the end.

        sign = 1.0 if phi >= 0.0 else -1.0
        
        t_out = sign*t_out
        
        if t_out < 0 :
            t_out += self.tau
            
        if t_range == 1 :
            return sign*t_out

        # t_range is 2

        # How many phi periods out are we (stored in i)?
        f,i = math.modf(0.5*sign*phi/math.pi)

        # Add that number of phi periods
        t_out += self.tau*i

        return sign*t_out

        
    def y_at_phi(self, phi) :
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
            # print("WARNING VR^2 IS < 0.  HOPEFULLY, JUST A ROUNDING ERROR: ", vrsq, " PHI ", phi)
            vrsq = 0.0
        
        vr = math.sqrt( vrsq )

        if phi > math.pi :
            vr = -vr
        
        return [ self.m, r, phi, vr, om ]


class OrientedOrbit(Orbit) :
    '''An Orbit, oriented in the reference frame of the insertion point.
    Currently still in the same plane as an Orbit, but with r0
    positioned at th0 in the outer frame. th0 is automatically
    computed from the insertion point.

    :Conventions:
        
    phi is the polar angle in the coordinate system of the orbit with r0 at phi=0

    theta (or th) is the polar angle in the "outer" coordinate system
        
    phi_i is the orbit angle of the "intertion point", i.e. the point
    which is used to initialize the orbit
        
    th_i is the "outer" coordinate system angle of the intertion point
        
    th0 is the angle of r0 in the "outer" coordinate system.
        
    In general th = phi + th0
        
    th0 = th_i - phi_i
        
    So, r(th) = Orbit.r(th - th0)

    :See also: Orbit

    '''
    
    def __init__(self, y, GM, phi_resolution=5000.0, force=None) :
        
        super().__init__(y, GM, phi_resolution, force)

        m, r, th_i, vr, om = y

        self.th_i = th_i
        self.th0 = th_i - self.phi_i

    def intersect_soi(self, th_soi, d_soi, r_soi):
        phi_x = super().intersect_soi(self.phi_at_th(th_soi), d_soi, r_soi)
        th_x = [phi + self.th0 for phi in phi_x]
        return th_x

    def sample_t(self, T) :
        R, PHI = super().sample_t(T)
        TH = [self.th_at_phi(phi) for phi in PHI]
        return (R, TH)
    
    def sample_th(self, TH) :
        PH = [self.phi_at_th(th) for th in TH]
        return self.sample_phi(PH)
        
    def phi_at_th(self, th) :
        return (th - self.th0)

    def t_at_th(self, th) :
        return self.t_at_phi(self.phi_at_th(th))
    
    def th_at_phi(self, phi) :
        return (phi + self.th0)

    def y_at_th(self, th):
        Y = self.y_at_phi(self.phi_at_th(th))
        Y[2] = self.th_at_phi(Y[2])
        return Y

