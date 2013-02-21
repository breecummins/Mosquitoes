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
        self.numMosq = len(initPosx)
        # placeholders for subclass assignment
        self.currentPosy = None
        self.currentCO2 = None 
        # construct parameter dictionary
        self.mosqParams = {'startTime':350.0,'decisionInterval':1.0,'hostRadius':5,
                    'CO2Thresh':0.01,'CO2Sat':1.0,'CO2Kappa':0.0,'spdMin':0.4,'spdMax':1.5,
                    'windThresh':0.0,'windSat':0.5,'windKappa':0.0,'windDirMin':np.pi/6,
                    'windDirMax':np.pi/2}
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
        This method determines strength of mosquito response to a signal.
        responseStr is the string that identifies which thresh, sat, and kappa 
        parameter set to use. Examples: 'CO2' or 'wind'.
        currentVal is an array of the numerical values of the wind or CO2 etc. 
        at the location of each mosquito.

        '''
        unscaledSat = self.mosqParams[responseStr+'Sat']
        val = v/unscaledSat
        sat = 1.0
        kappa = self.mosqParams[responseStr+'Kappa']
        thresh = self.mosqParams[responseStr+'ScaledThresh']
        def response(v):
            if v <= thresh:
                return 0.0
            elif v >= sat:
                return 1.0
            else:
                return (1.0+kappa*thresh)*(v - thresh)/(1.0+kappa*thresh*v*(1.0-thresh))
        return map(response,val)




class klinotaxis(mosquitoPopulation):
    '''
    Plume tracking using memory.

    '''
    def __init__(self,initPosx,**kwargs):
        mosquitoPopulation.__init__(initPosx,**kwargs)
        self.previousCO2 = np.zeros(size(initPosx))
        # need dif sat, thresh, and kappa values

    def _respondInPlume(self,boolarray):
        #FIXME - stub
        dx = None
        dy = None
        return dx, dy


class upwind(klinotaxis):
    def __init__(self,simulation,environ,initPosx,**kwargs):
        '''
        simulation is an instance of class numericalSims and environ is
        an instance of class environment

        '''
        mosquitoPopulation.__init__(initPosx,**kwargs)
        self.initPosy = simulation.simsParams['domainLength'] - simulation.simsParams['h']
        self.currentPosy = self.initPosy*np.ones(initPosx.shape)
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)

    def _respondWindOnly(self,boolarray):
        #FIXME - stub
        dx = None
        dy = None
        return dx, dy



class downwind(klinotaxis):
    def __init__(self,simulation,environ,initPosx,**kwargs):
        '''
        simulation is an instance of class numericalSims and environ is
        an instance of class environment

        '''
        mosquitoPopulation.__init__(initPosx,**kwargs)
        self.initPosy = 0.0
        self.currentPosy = np.zeros(initPosx.shape)
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)

    def _respondWindOnly(self,boolarray):
        #FIXME - stub
        dx = None
        dy = None
        return dx, dy


class crosswind(klinotaxis):
    def __init__(self,simulation,environ,initPosx,**kwargs):
        '''
        simulation is an instance of class numericalSims and environ is
        an instance of class environment

        '''
        mosquitoPopulation.__init__(initPosx,**kwargs)
        self.initPosy = 0.0
        self.currentPosy = np.zeros(initPosx.shape)
        self.crosswindDuration = None # placeholder for duration and direction
        self.currentU,self.currentV,self.currentCO2 = environ.getSignal(self.currentPosx,self.currentPosy,simulation)

    def _respondWindOnly(self,boolarray):
        #FIXME - stub
        dx = None
        dy = None
        return dx, dy



if __name__ == '__main__':
    
    mymosqs = mosquitoPopulation([1.0,2.0,3.0],windKappa=1.0,windThresh=0.01)
    responseStr = 'wind'
    currentVal = np.array([-1.e-12,0.1,1.0+1.e-12])
    print(mymosqs._responseCurve(responseStr,currentVal))
    print(mymosqs.mosqParams)
    
    


