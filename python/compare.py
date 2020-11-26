'''Recommended usage

import compare as cmp

if cmp.fsame(a,b) :
   print("Same!")

'''

def fsame(a, b, tolfrac=1E-15, report_to = None) :
    '''Compare floats.  Returns true if within fractional tolerance.

    report_to is a file object
    '''
    
    tol = tolfrac*0.5*(abs(a) + abs(b))

    same = ( abs(a-b) <= tol )

    if report_to is not None :
    
        if same :

            report_to.write("++++: %s and %s are within tolerance of %s\n" %
                            (a, b, tol) )

        else :

            report_to.write("----: %s and %s are outside tolerance of %s\n" %
                            (a, b, tol) )

    return same

