# Copyright © 2021 Jonathan Grot

# Tools for working with Matplotlib

import math

def plot_circle(ax, r, n, **plot_opts) :
    '''Plots a circle to axes

    :param object ax: Matplotlib Axes object
    :param float   r: Radius of circle
    :param int     n: Number of points to plot
    '''
    x = []
    y = []

    # Last point will overlap to close the circle
    dth = 2.0*math.pi / (n+1)
    maxth = 2.0*math.pi
    
    plot_polar(ax, lambda th: r, lambda th: th, 0, dth, maxth, **plot_opts)

def plot_polar(ax, r_of_z, th_of_z, z0, dz, max_z, centered=False, **plot_opts) :
    '''Polar plot

    :param object ax: Matplotlib Axes object
    :param float  r_of_z: Radius as a function of parameter z
    :param float  th_of_z: Theta as a function of parameter z.  If None, then th = z
    :param float  z0: Parameter z initial value
    :param float  dz: Parameter z intervals
    :param float max_z: Maximum value of z
    :param bool centered: Calls square_plot if true
    :param args plot_opts: named arguments to pass to the plot function.
    '''
    x = []
    y = []

    iz = 0.0

    while True:
        
        z = z0 + dz*iz

        if z > max_z :
            break
        else :
            iz += 1.0

        r = r_of_z(z)

        if th_of_z is None :
            th = z
        else :
            th = th_of_z(z)
            
        x.append(r*math.cos(th))
        y.append(r*math.sin(th))

    if centered :
        square_plot(ax, x, y, **plot_opts)
    else :
        ax.plot(x, y, **plot_opts)

def square_plot(ax, xs, ys, **plot_opts) :
    '''Centers the plot in square axes

    :param object ax: Matplotlib Axes object
    :param list   xs: List of x values
    :param list   ys: List of y values
    :param args   plot_opts: other named arguments to pass to the Matplotlib plot command
    '''
    xb = [ min(xs), max(xs) ]
    yb = [ min(ys), max(ys) ]
    
    xrng = xb[1] - xb[0]
    yrng = yb[1] - yb[0]
    
    maxrng = max( xrng, yrng )
    
    xpad = 0.5*(maxrng - xrng)
    ypad = 0.5*(maxrng - yrng)
    
    xmin = xb[0] - xpad
    xmax = xb[1] + xpad
    ymin = yb[0] - ypad
    ymax = yb[1] + ypad
    
    ax.set_xlim( xmin, xmax )
    ax.set_ylim( ymin, ymax )
    ax.set_aspect(1.0)
    ax.plot(xs, ys, **plot_opts)

