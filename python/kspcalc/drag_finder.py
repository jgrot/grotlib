#! /usr/bin/env python
#
# Copyright © 2021 Jonathan Grot
#
# This file is part of grotlib.
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

def falpha( t, y, flyer ) :
    return 0.5 * math.pi

instructions = [
    "Drag Term Finder",
    "Create a missile to test.",
    "Fly the missile straight up and let it burn all the fuel.",
    "Record the maximum altitude and the crash time in TALO seconds.",
    "Create an input record (nosecone_test.json)",
    "Run the test subcommand to find the drag term by minmizing the error in target altitude and crash time."
]

if __name__ == "__main__" :

    ksp.augmentBodyDbs( )
    ksp.processBodyDbs( )

    parser = argparse.ArgumentParser(description="Drag Term Finder")
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
        stage.loadJSON( ksp.pthdat("drag_tests/" + case_data["stage file"] ))
        flyer = ksp.FlyingStage( stage, "stage", "Kerbin", fthrottle, falpha )

        htarget = case_data["htarget"]
        crash_time = case_data["crash time"]

        
        itrial = 1
        def trial( X ) :

            
            global itrial
            
            stage.dragco = X[0]
            flyer.launch( )
            flyer.flyTo(30000)
            
            hmax = flyer.maxr - flyer.R

            z1 = hmax - htarget
            z2 = flyer.solnt[-1] - crash_time

            print( "iter: %4i, drag term: %10.8f, herr: %10.4e, terr: %10.2f" % (itrial, X[0], z1, z2) )
            itrial += 1
            
            return ( z1*z1 + z2*z2 )

        # Initial drag term
        x0 = [0.1]
        
        result = minimize(trial, x0, method="Nelder-Mead", options={"maxiter":1000})

        print(result.x)
