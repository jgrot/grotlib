#! /usr/bin/env python
# Copyright © 2021 Jonathan Grot

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
    "Create a missile and estimate the drag term using nosecone_drag.py",
    "Make sure the most recent drag term is stored with the stage",
    "Specify the corresponding nosecone test file (e.g. nosecone_test.json)",
    "Run the test subcommand optmize the drag divergence function."
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
