# Copyright © 2021 Jonathan Grot

'''Toleranced floating point comparisons.
'''

import math
import sys

# Control of reporting behavior
CMP_ONLY_NOT_SAME = 1
CMP_ONLY_SAME = 2
CMP_BOTH = 3

def fsame(a, b, tolfrac=1E-15, report_to=sys.stderr, report=CMP_BOTH) :
    '''Compare floats.  Returns true if within fractional tolerance.

    report_to is a file object
    '''
    
    tol = tolfrac*0.5*(abs(a) + abs(b))

    if tol == math.inf :
        same = False
    else :
        same = ( abs(a-b) <= tol )

    if report_to is not None :
    
        if same and report >= 2 :

            report_to.write("++++: %s and %s are within tolerance of %s\n" %
                            (a, b, tol) )

        elif not same and report >= 1 :

            report_to.write("----: %s and %s are outside tolerance of %s\n" %
                            (a, b, tol) )

    return same

def flessthan(a, b, tolfrac=1E-15) :
    ''' a <? b '''

    tol = tolfrac*0.5*(abs(a) + abs(b))
    return ( a-b < -tol )

def fgreaterthan(a, b, tolfrac=1E-15) :
    ''' a >? b '''

    tol = tolfrac*0.5*(abs(a) + abs(b))
    return ( a-b > tol )
