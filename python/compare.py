'''Recommended usage

import compare as cmp

if cmp.fsame(a,b) :
   print("Same!")

'''

import sys


CMP_ONLY_NOT_SAME = 1
CMP_ONLY_SAME = 2
CMP_BOTH = 3

def fsame(a, b, tolfrac=1E-15, report_to = sys.stderr, report = CMP_BOTH ) :
    '''Compare floats.  Returns true if within fractional tolerance.

    report_to is a file object
    '''
    
    tol = tolfrac*0.5*(abs(a) + abs(b))

    same = ( abs(a-b) <= tol )

    if report_to is not None :
    
        if same and report >= 2 :

            report_to.write("++++: %s and %s are within tolerance of %s\n" %
                            (a, b, tol) )

        elif not same and report >= 1 :

            report_to.write("----: %s and %s are outside tolerance of %s\n" %
                            (a, b, tol) )

    return same

