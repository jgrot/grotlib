# Copyright © 2021 Jonathan Grot

from scipy.interpolate import interp1d

# Grotlib imports
import compare as cmp
import type_tools

class BadDimError( Exception ) :
    def __init__( self, msg ) :
        super().__init__( msg )

class BadValueError( Exception ) :
    def __init__( self, msg ) :
        super().__init__( msg )

class InitError( Exception ) :
    def __init__( self, msg ) :
        super().__init__( msg )
        
class RangeError( Exception ) :
    def __init__( self, msg ) :
        super().__init__( msg )

class Functor :
    '''Abstract base class for a multi-variable function with built-in plotting capability.

    Independent variables are referred to as X: [ x1, x2, ... ]
    Dependent variables are referred to as Y: [ y1, y2, ... ]

    :param iterable rangemin: (x1min, x2min, ...)
    :param iterable rangemax: (x1max, x2max, ...)

    :Usage Notes:
      * Derived class must provide nDep() method to return number of independent variables returned in call()
      * Derived class definition of call()
        * call() MUST have the signature call(self, *X)
        * should usually call checkRange(X) before anything else
        * return Y: [ y1, y2, ... ]

    :Throws:
      * BadDimError
      * RangeError
      * InitError
    '''
    def __init__( self, rangemin, rangemax ) :
        
        self.rangemin = rangemin
        self.n_indep_var = len( rangemin )

        if len( rangemax ) != self.n_indep_var :
            raise BadDimError( "Inconsistent range dimensions" )
        else :
            self.rangemax = rangemax

        for i, xmin in enumerate( rangemin ) :
            if rangemax[i] <= rangemin[i] :
                msg = "Functor: range maximum, %f, is <= range minimum, %f for independent variable %i" % ( rangemax[i], rangemin[i], i )
                raise RangeError( msg )
            
        try :
            self.n_dep_var = self.nDep( )
        except :
            raise InitError( "Attempt to call derived class method nDep() failed" )

    def call( self, *X ) :
        raise Exception("Please override the call() method for your derived class.")

    def checkRange( self, X ) :
        '''Checks dimensionality of X and the elements of X against ranges.

        :Throws:
          * BadDimError
          * RangeError
        '''
        
        if len(X) != self.n_indep_var :
            raise BadDimError( "Bad number of input independent variables" )
        
        for i, x in enumerate( X ) :
            if cmp.flessthan(x, self.rangemin[i]) or cmp.fgreaterthan(x, self.rangemax[i]) :
                msg = "Functor: independent variable %i is out of range with value %f" % ( i, x )
                raise RangeError( msg )

    def plot( self, mpl_axes, range_specs, resolution = None, dep_var = 0 ) :
        '''Makes 1D or 2D plots of slices through the range space of the function.

        :param object: Matplotlib Axes object to plot to.  Tip: for 1D plots, call repeatedly for the same Axes object to put mutiple plots on the same Axes.
        :param list range_specs: [ {<True> | (min,max) | <Fixed value>}, ... ].  See discussion below.
        :param list resolution: [ <Number of points along range axis for indep var>, ... ], one for each indep var, even the ones that are fixed, to avoid confusion.  A value of None yields [100, 100, ...].
        :param int dep_var: Index into result Y to plot.

        The 'range_specs' argument contains instructions for each independent variable:
          * <True> indicates that the function should be plotted over the full range (specified in the constructor) for that independent variable
          * (min,max) overrides the specified range.  Exceeding the specified ranges will trigger a RangeError if call() is using checkRange() and not catching the exception.
          * <Fixed value> indicates that this independent variable should have a fixed value.  Exceeding the specified ranges will trigger a RangeError if call() is using checkRange() and not catching the exception.  
        
        :Matplotlib Tips:
        Annotate the Axes object with:
           * set_title()
           * set_xlable()
           * set_ylabel()
       
        :Throws:
          * BadDimError
          * BadValueError
        '''
        import numpy
        
        if len(range_specs) != self.n_indep_var :
            raise BadDimError( "Bad number of elements in argument range_specs" )

        # Process range_specs
        plot_ndim = 0
        dimrange = []
        plot_axes_vars = [] # Indexes of indep variables to use for plot axes
        
        for ispec, spec in enumerate(range_specs) :
            if spec == True : # i.e Explcitly equals True
                # Use the ranges specified in the constructor
                dimrange.append( (self.rangemin[ispec], self.rangemax[ispec]) )
                plot_axes_vars.append( ispec )
                plot_ndim += 1
            elif type_tools.isIterableNonString( spec ) :
                # For the plot, override the ranges specified in the constructor
                if len(spec) != 2 :
                    raise BadDimError("Range limit specification must be (min, max)")
                if spec[1] <= spec[0] :
                    raise BadValueError("Bad order for limited range spec %s" % repr(spec))
                dimrange.append( spec )
                plot_axes_vars.append( ispec )
                plot_ndim += 1
            elif isinstance( spec, float ) or isinstance( spec, int ) :
                # Fixed value for indep var
                dimrange.append( spec )
            else :
                raise BadValueError( "Bad expression for specification: " + repr(spec) )
            
        if plot_ndim > 2 :
            raise BadValueError( "Too many independent variables are used for plot axes in range specification" )

        if plot_ndim == 0 :
            raise BadValueError( "Could not identify independent variables to use for plot dimensions in specifications" )

        if resolution is None :
            # Default of 100 for each plot dimension
            resolution = [100] * self.n_indep_var
        elif len( resolution ) != self.n_indep_var :
            raise BadDimError( "Bad dimensionality for resolution argument" )

        # Make sure elements of resolution argument are integers.
        # Note: the value of the elements are checked during plot
        # pre-processing as some will be ignored.
        for r in resolution :
            if not isinstance( r, int ) :
                raise BadValueError( "Bad resolution term: %s (must be int)" % r )

        # Make sure dep_var is sane.
        if dep_var < 0 or dep_var >= self.n_dep_var :
            raise BadValueError( "Bad value for dep_var" )

        # Computed indep var value lists (and fixed values) for the points being plotted.
        xvalues = []

        # Evaluate independent variable ranges and fixed values.
        for ispec, spec in enumerate(dimrange) :
            try :
                # Try to parse a range
                xmin, xmax = spec
                xfixed = None
            except :
                # Could not get range, this element is a fixed value.
                xmin = None
                xmax = None
                xfixed = spec

            if xmin is not None :
                # Expand a range into a list of values for plotting
                r = resolution[ispec]
                if r < 5 :
                    raise BadValueError( "Resolution term, %i, is too low (must be 10 or more)" % r )
                dx = (xmax - xmin)/float(r-1)
                xx = [ xmin + dx*i for i in range(r) ]
                xvalues.append( xx )
            else :
                # Record a fixed value
                xvalues.append( xfixed )

        # Generate Y values and plot
        if plot_ndim == 1 :
            # Indep var used for horizontal axis
            ii = plot_axes_vars[0]
            # Get resolution
            r = resolution[ii]

            y = []
            x = []
            for i in range(r) :
                X = []
                for l in range( self.n_indep_var ) :
                    if l==ii :
                        X.append( xvalues[l][i] )
                        x.append( xvalues[l][i] )
                    else :
                        X.append( xvalues[l] )
                y.append( self.call( *X )[dep_var] )
                
            mpl_axes.plot(x, y, marker = "o")

        if plot_ndim == 2 :

            ii1 = plot_axes_vars[0]
            ii2 = plot_axes_vars[1]

            r1 = resolution[ii1]
            r2 = resolution[ii2]

            y = numpy.zeros( (r1, r2) )
            for j in range( r2 ) :
                for i in range( r1 ) :
                    X = []
                    for l in range( self.n_indep_var ) :
                        if l == ii1 :
                            x = xvalues[l][i]
                            X.append( x )
                        elif l == ii2 :
                            x = xvalues[l][j]
                            X.append( x )
                        else :
                            X.append( xvalues[l] )
                    y[i,j] = self.call( *X )[dep_var]

            extent = [ self.rangemin[ii1], self.rangemax[ii1], self.rangemin[ii2], self.rangemax[ii2] ]
            mpl_axes.imshow( y.transpose(), origin='lower', extent=extent )

class Interp1DFunctor( Functor ) :
    '''A Functor based on scipy.interpolate.interp1d

    :param list xs: in values
    :param list ys: out values
    '''
    def __init__( self, xs, ys, kind="linear", bounds_error=False, fill_value=0.0 ) :
        rangemin = [xs[0]]
        rangemax = [xs[-1]]
        
        self.f = interp1d( xs, ys, kind=kind, bounds_error=bounds_error, fill_value=fill_value )
        super().__init__( rangemin, rangemax )

    def call( self, *X ) :
        # Note: checkRange is not called here, since the interp1d
        # range handling is used.
        return [ self.f(X[0]) ]

    def nDep( self ) :
        return 1

#
# Unit test (and examples)
#
if __name__ == "__main__" :
    # Unit test and demo

    import math

    class Norm3D( Functor ) :
        
        def __init__( self, rangemin, rangemax ) : 
            super().__init__( rangemin, rangemax )

        def call( self, *X ) :
            self.checkRange( X )

            r = math.sqrt( X[0]*X[0] + 0.5*X[1]*X[1] + 0.25*X[2]*X[2] )
            ans = math.exp( -(r*r) )
            
            return [ ans ]

        def nDep( self ) :
            return 1

    rangemin = (-1, -1, -1)
    rangemax = ( 1,  1,  1)
    badrangemax = (-1, -1, -1)

    try :
        fbad = Norm3D( rangemin, badrangemax )
    except RangeError :
        print("SUCCESS, trapped bad initialization")
    else :
        print("FAIL, did not trap bad initialization")
        exit(1)
    
    f = Norm3D( rangemin, rangemax )

    try :
        Ybad = f.call( -1.1, 1.1, -3.3 )
    except RangeError :
        print("SUCCESS, trapped out of range inputs")
    except :
        print("FAIL, did not trap out of range inputs")
        exit(1)
    else :
        print("FAIL, did not trap out of range inputs")
        exit(1)

    try :
        Ybad = f.call( 0.0, 0.0 )
    except BadDimError :
        print("SUCCESS, trapped bad dimensions")
    except :
        print("FAIL, did not trap bad dimensions")
        exit(1)
    else :
        print("FAIL, did not trap bad dimensions")
        exit(1)
        
    import matplotlib
    import matplotlib.pyplot as plt

    if True :
        fig1, ax1 = plt.subplots()
        f.plot(ax1, [True, 0.0, 0.0] )
        f.plot(ax1, [0.0, True, 0.0] )
        f.plot(ax1, [0.0, 0.0, True] )
        
        fig2, ax2 = plt.subplots()
        f.plot(ax2, [True, 0.5, 0.5] )
        f.plot(ax2, [0.5, True, 0.5] )
        f.plot(ax2, [0.5, 0.5, True] )
        
        fig3, ax3 = plt.subplots()
        f.plot(ax3, [True, True, 0.0] )
        
    if True :
        class Axis( Functor ) :

            def __init__( self, rangemin, rangemax, whichaxis ) : 
                super().__init__( rangemin, rangemax )
                self.whichaxis = whichaxis

            def call( self, *X ) :
                self.checkRange( X )

                return [ X[self.whichaxis] ]

            def nDep( self ) :
                return 1

        fx = Axis(rangemin, rangemax, 0)
        fy = Axis(rangemin, rangemax, 1)
        fz = Axis(rangemin, rangemax, 2)

        # X increasing
        fig4, ax4 = plt.subplots()
        fx.plot(ax4, [True, True, 0.0] )

        # Y increasing
        fig4, ax4 = plt.subplots()
        fy.plot(ax4, [True, True, 0.0] )

        # Z should increase in the Y direction
        fig4, ax4 = plt.subplots()
        fz.plot(ax4, [0.0, True, True] )

    plt.show()

    
