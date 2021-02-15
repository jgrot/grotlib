# Copyright © 2021 Jonathan Grot

# Tools for working with Matplotlib

import math

def offset_xy(XY, xoff, yoff) :

    X, Y = XY

    X2 = [x+xoff for x in X]
    Y2 = [y+yoff for y in Y]

    return ( tuple(X2), tuple(Y2) )


def sample_circle(r, n) :
    '''Plots a circle to axes

    :param float   r: Radius of circle
    :param int     n: Number of points to plot
    '''

    # Last point will overlap to close the circle
    dth = 2.0*math.pi / (n+1)
    maxth = 2.0*math.pi
    
    return(sample_polar(lambda th: r, lambda th: th, 0, dth, maxth))

def sample_polar(r_of_z, th_of_z, z0, dz, max_z) :
    '''Sample a polar function.

    :param float r_of_z: Radius as a function of parameter z
    :param float th_of_z: Theta as a function of parameter z.  If None, then th = z
    :param float z0: Parameter z initial value
    :param float dz: Parameter z intervals
    :param float max_z: Maximum value of z

    :returns: tuple( tuple(<x values>),tuple(<y values>) )
    '''

    X=[]
    Y=[]
    
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
            
        X.append(r*math.cos(th))
        Y.append(r*math.sin(th))

    return (X,Y)

def square_plots(ax, plots, xmin, xmax, ymin, ymax, plot_opts=None) :
    ax.set_xlim( xmin, xmax )
    ax.set_ylim( ymin, ymax )
    ax.set_aspect(1.0)

    for iplot, plot in enumerate(plots) :
        xs, ys = plot
        if plot_opts is not None and plot_opts[iplot] is not None :
            ax.plot(xs, ys, **plot_opts[iplot])
        else :
            ax.plot(xs, ys)

def square_plot_bounds(plots) :
    '''Computes params for centering a bunch of plots in a square (presumable aspect ratio 1.0) plot.

    :param [[X,Y]] plots: A collection of X,Y arrays for plotting.

    :returns: (xmin, xmax, ymin, ymax)
    '''

    xb = [math.inf, -math.inf]
    yb = [math.inf, -math.inf]
    
    for X,Y in plots :
        xb[0] = min(xb[0], min(X))
        yb[0] = min(yb[0], min(Y))

        xb[1] = max(xb[1], max(X))
        yb[1] = max(yb[1], max(Y))
    
    xrng = xb[1] - xb[0]
    yrng = yb[1] - yb[0]
    
    maxrng = max( xrng, yrng )
    
    xpad = 0.5*(maxrng - xrng)
    ypad = 0.5*(maxrng - yrng)
    
    xmin = xb[0] - xpad
    xmax = xb[1] + xpad
    ymin = yb[0] - ypad
    ymax = yb[1] + ypad
    
    return (xmin, xmax, ymin, ymax)

