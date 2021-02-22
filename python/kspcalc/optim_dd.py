#! /usr/bin/env python
#
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
import argparse
import json
import math
from scipy.optimize import minimize

import ksp
from ksp import DragDivergence

def falpha( t, y, flyer ) :
    return 0.5 * math.pi

instructions = [
    "Drag Divergence Optimizer",
    "Create a missile and estimate the (slow) drag term using drag_finder.py",
    "Store the (slow) drag term in the stage file",
    "Fly the same missile high and fast and record results in a test file",
    "Run the test subcommand to find drag divergence parameters"
]

if __name__ == "__main__" :

    import matplotlib
    import matplotlib.pyplot as plt

    ksp.augmentBodyDbs( )
    ksp.processBodyDbs( )

    parser = argparse.ArgumentParser(description="Drag Divergence Optmizer")
    subparsers = parser.add_subparsers(help="Commands", dest="command")
    
    cmd_inst = subparsers.add_parser("inst", help="Instructions")

    cmd_test = subparsers.add_parser("test", help="Run the test")
    cmd_test.add_argument("case_file", help="Name of the case file")

    args = parser.parse_args()

    if args.command == "inst" :
        for line in instructions :
            print(line)
        exit(0)

    if args.command == "test" :

        case_data = None
        with open( args.case_file, "rt" ) as f :
            case_data = json.load( f )

        def fthrottle( t, y ) :
            return case_data["throttle"]

        stage = ksp.Stage( )
        stage.loadJSON( ksp.pthdat("drag_tests/" + case_data["stage file"]) )
        flyer = ksp.FlyingStage( stage, "stage", "Kerbin", fthrottle, falpha )

        htarget = case_data["htarget"]
        crash_time = case_data["crash time"]

        itrial = 1
        def trial( X ) :
            global itrial

            if X[0] < 1.0 or \
               X[0] > 10.0 or \
               X[1] < 5.0 or \
               X[1] > 200.0 or \
               X[2] < 0.01 or \
               X[2] > 100.0 :
               return math.inf
            
            ksp.dd = DragDivergence(X[0], X[1], X[2])
            
            flyer.launch( )
            try :
                flyer.flyTo(30000)
            except :
                return math.inf
                
            hmax = flyer.maxr - flyer.R

            z1 = hmax - htarget
            z2 = flyer.solnt[-1] - crash_time

            print( "iter: %4i, c0: %10.4f, c1: %10.4f, c2: %10.4f, herr: %10.4e, terr: %10.2f" % (itrial, X[0], X[1], X[2], z1, z2) )
            itrial += 1
            
            return ( z1*z1 + z2*z2 )

        # Initial drag divergence coefficients
        x0 = [3.0, 10.0, 0.7]
        
        result = minimize(trial, x0, method="Nelder-Mead", options={"maxiter":100})

        print(result.x)
    
    fig, ax = plt.subplots()
    dd = DragDivergence( *result.x )
    dd.plot( ax, [(0,4)], resolution=[1000] )
    plt.show()
