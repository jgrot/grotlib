# Copyright © 2021 Jonathan Grot

'''Tools for working with Matplotlib'''

import math

def offset_xy(XY, xoff, yoff) :
    '''Add an x-y offset to an array of x values and an array of y values.

    :param iterable XY: [ [<x values>], [<y values>] ]
    :param float xoff: x-offset value
    :param float yoff: y-offset value

    :returns: ``tuple( tuple(<x values>), tuple(<y values>) )``
    '''

    X, Y = XY

    X2 = [x+xoff for x in X]
    Y2 = [y+yoff for y in Y]

    return ( tuple(X2), tuple(Y2) )

def sample_circle(r, n) :
    '''Generate x-y points along a circle.

    :Tips:

    - Add an offset with offset_xy()

    :param float   r: Radius of circle
    :param int     n: Number of points to plot

    :returns: ``tuple( tuple(<x values>), tuple(<y values>) )``
    '''

    # Last point will overlap to close the circle
    dth = 2.0*math.pi / (n-1)

    TH = [ i*dth for i in range(n) ]
    R  = [ r for i in range(n) ]

    return polar_to_xy(R, TH)
    
def polar_to_xy(R, TH) :
    X = []
    Y = []

    for i, r in enumerate(R) :
        th = TH[i]
        X.append(r*math.cos(th))
        Y.append(r*math.sin(th))

    return (X, Y)

def square_plot_bounds(plots) :
    '''Computes params for centering a bunch of plots in a square (presumably aspect ratio 1.0) plot.

    :param [[X,Y]] plots: A collection of X,Y arrays for plotting.

    :returns: (xmin, xmax, ymin, ymax)

    :Tips:

    - Use in conjunction with square_plots()
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

def square_plots(ax, plots, bbox, plot_opts=None) :
    '''Writes an aspect ratio 1.0 plot to the axes "ax"

    :param object ax: Matplotlib Axes object
    :param iterable plots: A collection of plot XY lists: [ [[<x values>], [y values]], ...]
    :param iterable bbox: (xmin, xmax, ymin, ymax)
    :param iterable plot_opts: Either None, or a collection of corresponding Matplotlib options for each plot: [{<options>} | None, ...]

    :Tips:

    - Use square_plot_bounds() with a suitable collection of plot XYs to set bbox
    '''

    xmin, xmax, ymin, ymax = bbox
    
    ax.set_xlim( xmin, xmax )
    ax.set_ylim( ymin, ymax )
    ax.set_aspect(1.0)

    for iplot, plot in enumerate(plots) :
        xs, ys = plot
        if plot_opts is not None and plot_opts[iplot] is not None :
            ax.plot(xs, ys, **plot_opts[iplot])
        else :
            ax.plot(xs, ys)
