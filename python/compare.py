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
'''Toleranced floating point comparisons.
'''

import math
import sys

# Control of reporting behavior
CMP_ONLY_NOT_SAME = 1
CMP_ONLY_SAME = 2
CMP_BOTH = 3

def compare_datasets(data, reference, tolfrac=1E-15) :
    '''Compares two lists of lists.
    '''
    for irow, row in enumerate( data ) :
        try :
            row_ref = reference[ irow ]
        except :
            row_ref = None
                
        if row_ref is None :
            print("Ran out of reference data")
            exit(1)
                
        for jelem, elem in enumerate( row ) :
            elem_ref = row_ref[ jelem ]
            if not fsame(elem, elem_ref, tolfrac=tolfrac, report_to=sys.stderr, report=CMP_ONLY_NOT_SAME) : exit(1)


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
