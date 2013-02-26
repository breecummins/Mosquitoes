import numpy as np

class mosquitoPopulation(object):

    '''
    This class represents a group of mosquito agents with the same parameters. 
    It is mandatorily subclassed to choose a plume tracking strategy (subclasses
    tropotaxis and klinotaxis) and subclassed again to choose a plume finding 
    strategy (upwind, downwind, and crosswind).

    '''

    def __init__(self,initPosx,**kwargs):
        '''
        initPosx is a numpy array of initial x positions of the mosquito population. 
        The length of the list is the number of mosquitoes in the population.
        kwargs are optional keyword arguments that may overwrite any of the 
        default parameter values assigned below.

        '''
        self.initPosx = initPosx
        self.currentPosx = initPosx
        # placeholders for subclass assignment
        self.currentPosy = None
        self.currentCO2 = None 
        # construct parameter dictionary, CO2 associated with speed calculations
        self.mosqParams = {'startTime':350.0,'decisionInterval':1.0,'hostRadius':5,'CO2Thresh':0.01,'CO2Sat':1.0,'CO2Kappa':0.0,'CO2WindowMin':0.4,'CO2WindowMax':1.5,'windThresh':0.0,'windSat':0.5,'windKappa':0.0,'windWindowMin':np.pi/6,'windWindowMax':np.pi/2}
        self.mosqParams.update(kwargs)
        self.mosqParams['windScaledThresh'] = self.mosqParams['windThresh']/self.mosqParams['windSat']
        self.mosqParams['CO2ScaledThresh'] = self.mosqParams['CO2Thresh']/self.mosqParams['CO2Sat']
        if self.mosqParams['windScaledThresh'] != 0 and self.mosqParams['windKappa'] <= -1.0/self.mosqParams['windScaledThresh']:
            raise ValueError('windKappa must be > -1.0 / %0.3f' %self.mosqParams['windScaledThresh'])
        if self.mosqParams['CO2ScaledThresh'] != 0 and self.mosqParams['CO2Kappa'] <= -1.0/self.mosqParams['CO2ScaledThresh']:
            raise ValueError('CO2Kappa must be > -1.0 / %0.3f' %self.mosqParams['CO2ScaledThresh'])

    def _updatePosition(self,environ,simulation):
        '''
        environ is an object containing odor plume information, an instance of
        class environment
        simulation is an object containing numerical simulation information, an
        instance of class numericalSims
        
        '''
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)
        mosqsinplume = self._inPlume()
        dx,dy = self._respondInPlume(mosqsinplume)
        self.currentPosx[mosqsinplume] = self.currentPosx[mosqsinplume] + dx
        self.currentPosy[mosqsinplume] = self.currentPosy[mosqsinplume] + dy
        dxw,dyw = self._respondWindOnly(~mosqsinplume) 
        self.currentPosx[~mosqsinplume] = self.currentPosx[~mosqsinplume] + dxw
        self.currentPosy[~mosqsinplume] = self.currentPosy[~mosqsinplume] + dyw
        self._atHost()

    def _inPlume(self):
        return self.currentCO2 >= self.mosqParams['CO2Thresh']

    def _respondInPlume(self,boolarray):
        '''
        Stub for subclass function

        '''
        pass

    def _respondWindOnly(self,boolarray):
        '''
        Stub for subclass function

        '''
        pass


    def _responseCurve(self,responseStr,currentVal):
        '''
        This method determines how the mosquito responds to a signal.
        responseStr is the string that identifies which thresh, sat, and kappa 
        parameter set to use. Examples: 'CO2' or 'wind' or 'diffCO2'.
        currentVal is an array of the numerical values of the wind or CO2 etc. 
        at the location of each mosquito.

        '''
        unscaledSat = self.mosqParams[responseStr+'Sat']
        val = currentVal/unscaledSat
        kappa = self.mosqParams[responseStr+'Kappa']
        thresh = self.mosqParams[responseStr+'ScaledThresh']
        fMax = self.mosqParams[responseStr+'WindowMax']
        fMin = self.mosqParams[responseStr+'WindowMin']
        def response(v):
            if v <= thresh:
                return 0.0
            elif v >= 1.0: # saturation scaled to 1.0
                return 1.0
            else:
                return (1.0+kappa*thresh)*(v - thresh)/(1.0+kappa*thresh*v*(1.0-thresh))
        return np.array([fMax - (fMax-fMin)*response(v) for v in val])

    def _atHost(self,environ):
        '''
        Remove mosquitoes who found a host.

        '''
        #FIXME
        
        

class klinotaxis(mosquitoPopulation):
    '''
    Plume tracking using memory.

    '''
    def __init__(self,initPosx,**kwargs):
        mosquitoPopulation.__init__(initPosx,**kwargs)
        self.previousCO2 = np.zeros(size(initPosx))
        #stub for subclass assignment
        self.previousMotionDir = None 
        # extra params for klinotaxis
        klinParams = {'diffCO2Thresh':(0.01/10)*0.0042/0.0833,'diffCO2Sat':(1.0 - 0.01)/50.,'diffCO2Kappa':0.0,'diffCO2WindowMin':np.pi/36,'diffCO2WindowMax':np.pi} 
        self.mosqParams.update(klinParams)
        self.mosqParams.update(kwargs)

    def _respondInPlume(self,boolarray):
        # calculate mosquito speed
        CO2 = self.currentCO2[boolarray]
        mosqSpeed = self._responseStrength('CO2',CO2)
        # calculate direction influenced by CO2
        diffCO2 = CO2 - self.previousCO2[boolarray]
        mosqCO2Window = self._responseStrength('diffCO2',np.abs(diffCO2))
        mosqCO2dir = self.previousMotionDir + mosqCO2Window*(-1.0 + 2.0*np.random.rand(len(CO2))
        # correct direction if CO2 decreased
        lowerCO2 = diffCO2 < 0.0
        mosqCO2dir[lowerCO2] -= np.pi
        # update CO2 memory
        self.previousCO2[boolarray] = CO2
        # calculate upwind component of motion
        U = self.currentU[boolarray]
        V = self.currentV[boolarray]
        velMag = np.sqrt(U**2 + V**2)
        mosqWindWindow = self._responseStrength('wind',velMag)
        mosqWindDir = np.pi + np.arctan2(V,U) + mosqWindWindow*(-1.0 + 2.0*np.random.rand(len(U)) # the pi term gives upwind flight
        # Advection plus average CO2 and wind responses. 
        dx = self.mosqParams['decisionInterval'] * (U + 0.5 * mosqSpeed * (np.cos(mosqCO2dir) + np.cos(mosqWindDir)))
        dy = self.mosqParams['decisionInterval'] * (V + 0.5 * mosqSpeed * (np.sin(mosqCO2dir) + np.sin(mosqWindDir)))
        # update position memory
        self.previousMotionDir[boolarray] = np.arctan2(dy,dx)
        return dx, dy



class upwind(klinotaxis):
    def __init__(self,simulation,environ,initPosx,**kwargs):
        '''
        simulation is an instance of class numericalSims and environ is
        an instance of class environment

        '''
        klinotaxis.__init__(initPosx,**kwargs)
        self.initPosy = simulation.simsParams['domainLength'] - simulation.simsParams['h']
        self.currentPosy = self.initPosy*np.ones(initPosx.shape)
        self.previousMotionDir = -np.pi/2 * np.ones(initPosx.shape)   
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)

    def _respondWindOnly(self,boolarray):
        # calculate speed
        mosqSpeed = self.mosqParams['spdMax']
        # calculate direction
        U = self.currentU[boolarray]
        V = self.currentV[boolarray]
        velMag = np.sqrt(U**2 + V**2)
        mosqWindWindow = self._responseStrength('wind',velMag)
        mosqWindDir = np.pi + np.arctan2(V,U) + mosqWindWindow*(-1.0 + 2.0*np.random.rand(len(U)) # the pi term gives upwind flight
        # Advection plus wind response. 
        dx = self.mosqParams['decisionInterval'] * (U + mosqSpeed*np.cos(mosqWindDir))
        dy = self.mosqParams['decisionInterval'] * (V + mosqSpeed*np.sin(mosqWindDir))
        # update position memory
        self.previousMotionDir[boolarray] = np.arctan2(dy,dx)
        return dx, dy



class downwind(klinotaxis):
    def __init__(self,simulation,environ,initPosx,**kwargs):
        '''
        simulation is an instance of class numericalSims and environ is
        an instance of class environment

        '''
         klinotaxis.__init__(initPosx,**kwargs)
       self.initPosy = 0.0
        self.currentPosy = np.zeros(initPosx.shape)
        self.previousMotionDir = np.pi/2 * np.ones(initPosx.shape)   
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)

    def _respondWindOnly(self,boolarray):
        # calculate speed
        mosqSpeed = self.mosqParams['spdMax']
        # calculate direction
        U = self.currentU[boolarray]
        V = self.currentV[boolarray]
        velMag = np.sqrt(U**2 + V**2)
        mosqWindWindow = self._responseStrength('wind',velMag)
        mosqWindDir = np.arctan2(V,U) + mosqWindWindow*(-1.0 + 2.0*np.random.rand(len(U)) 
        # Advection plus wind response. 
        dx = self.mosqParams['decisionInterval'] * (U + mosqSpeed*np.cos(mosqWindDir))
        dy = self.mosqParams['decisionInterval'] * (V + mosqSpeed*np.sin(mosqWindDir))
        # update position memory
        self.previousMotionDir[boolarray] = np.arctan2(dy,dx)
        return dx, dy


class crosswind(klinotaxis):
    def __init__(self,simulation,environ,initPosx,**kwargs):
        '''
        simulation is an instance of class numericalSims and environ is
        an instance of class environment

        '''
        klinotaxis.__init__(initPosx,**kwargs)
        self.initPosy = 0.0
        self.currentPosy = np.zeros(initPosx.shape)
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)
        #CW specific parameters
        self.crosswindDuration = 5 + 4.999 * np.random.rand(len(initPosx)).astype(int)
        self.crosswindDirection = np.array([-1.0 if np.random.rand() < 0.5 else 1.0 for _ in range(len(initPosx))])
        self.previousMotionDir = np.pi/2 * self.crosswindDirection  

    def _setDurationDirection(self):
        boolarray = self.crosswindDuration == 0
        self.crosswindDuration[boolarray] = 5 + 4.999 * np.random.rand(len(boolarray)).astype(int)
        self.crosswindDirection[boolarray] = -1.0*self.crosswindDirection[boolarray]


    def _respondWindOnly(self,boolarray):
        # calculate speed
        mosqSpeed = self.mosqParams['spdMax']
        # calculate direction
        U = self.currentU[boolarray]
        V = self.currentV[boolarray]
        velMag = np.sqrt(U**2 + V**2)
        mosqWindWindow = self._responseStrength('wind',velMag)
        mosqWindDir = self.crosswindDirection*np.pi/2 + np.arctan2(V,U) + mosqWindWindow*(-1.0 + 2.0*np.random.rand(len(U)) 
        # Advection plus wind response. 
        dx = self.mosqParams['decisionInterval'] * (U + mosqSpeed*np.cos(mosqWindDir))
        dy = self.mosqParams['decisionInterval'] * (V + mosqSpeed*np.sin(mosqWindDir))
        # update memory
        self.previousMotionDir[boolarray] = np.arctan2(dy,dx)
        self.crosswindDuration[boolarray] = self.crosswindDuration[boolarray] - 1
        self._setDurationDirection()
        return dx, dy



if __name__ == '__main__':
    
    mymosqs = mosquitoPopulation([1.0,2.0,3.0],windKappa=1.0,windThresh=0.01)
    responseStr = 'wind'
    currentVal = np.array([-1.e-12,0.1,1.0+1.e-12])
    print(mymosqs._responseCurve(responseStr,currentVal))
    print(mymosqs.mosqParams)
    
    


