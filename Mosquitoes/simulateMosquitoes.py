#!/usr/bin/env python

import numpy as np
import environment
import mosquito

def setHosts():
    pass

def setMosqs():
    pass

def saveResults():
    # save AND print results
    pass

xc,yc = setHosts()
environ = environment.environment(x,y)
initPosx = setMosqs()
mosqUpwindPop = mosquito.upwind(environ,initPosx)
mosqDownwindPop = mosquito.downwind(environ,initPosx)
mosqCrosswindPop = mosquito.crosswind(environ,initPosx)
dt = environ.simsParams['dt']

for t in np.arange(environ.simsParams['initialTime'],environ.simsParams['finalTime'],dt):

    environ.updateEnvironment()
    if t%1.0 < dt/2.0:
        mosqUpwindPop.updatePosition(environ,t)
        mosqDownwindPop.updatePosition(environ,t)
        mosqCrosswindPop.updatePosition(environ,t)
        if mosqUpwindPop.stopSimulation and mosqDownwindPop.stopSimulation and mosqCrosswindPop.stopSimulation:
            saveResults()
            stopsim = True
            break
        else:
            stopsim = False

if not stopsim:
    print('Not all mosquitoes are out of the domain.')
    saveResults()



