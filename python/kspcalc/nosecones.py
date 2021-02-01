import ksp
import math
from scipy.optimize import minimize

def fthrottle( t, y ) :
    return 1.0
def falpha( t, y, flyer ) :
    return 0.5 * math.pi

if __name__ == "__main__" :
    ksp.augmentBodyDbs( )
    ksp.processBodyDbs( )

    ZZ = []
    
    def f(X) :
        dragco = X[0]
        s = ksp.Stage(m0=(5,'t'), engines=[(1,"RT-10")], dragco=dragco)
        flyer = ksp.FlyingStage( s, "S1", "Kerbin", fthrottle, falpha )
        flyer.launch()
        flyer.flyTo(30000)
        hmax = flyer.maxr - flyer.R

        z = hmax - 50000.0

        ZZ.append( (z*z, dragco, hmax) )
        
        return z*z

    x0 = [0.1]
    r = minimize(f, x0, method="Nelder-Mead", options={"maxiter":1000})

    ZZ.sort()

    ZZ.reverse()
    
    for zz in ZZ :
        print( zz )

    print(r.x)
                      
