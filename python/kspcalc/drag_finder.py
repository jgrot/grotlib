#! /usr/bin/env python
# Copyright © 2021 Jonathan Grot

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

        ksp.dd = ksp.DragDivergence(0.99998138, 79.96128745, 0.14998027)
        
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
