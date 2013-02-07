import numpy as np
import environment

class mosquitoPopulation(object):

    '''
    This class represents a group of mosquito agents with the same parameters. 
    It is mandatorily subclassed to choose a plume tracking strategy (subclasses
    tropotaxis and klinotaxis) and subclassed again to choose a plume finding 
    strategy (upwind, downwind, and crosswind).

    '''

    def __init__(self,initPos,**kwargs):
        '''
        initPos is a list of initial x positions of the mosquito population. 
        The length of the list is the number of mosquitoes in the population.
        kwargs are optional arguments that may overwrite any of the default 
        parameter values assigned below.

        '''
        self.initPos = initPos
        self.currentPosx = initPos
        self.currentPosy = None #must be assigned in plume finding subclass
        self.numMosq = len(initPos)
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

    def _updatePosition(self,environment):
        '''
        Move to new position.
        
        '''
        u,v,CO2 = environment.getSignal(self.currentPosx,self.currentPosy)
        self = self._respondToSignal(u,v,CO2)

    def _respondToSignal(self,u,v,CO2):
        '''
        Stub to be defined by subclass

        '''
        return None


class klinotaxis(mosquitoPopulation):
    '''
    Plume tracking using memory.

    '''
    def __init__(self,initPos,**kwargs):
        mosquitoPopulation.__init__(initPos,**kwargs)
        #need extra dict entries here

    def _responseCurve(self,responseStr,currentVal):
        '''
        This method determines strength of mosquito response to a signal.
        responseStr is the string that identifies which thresh, sat, and kappa 
        parameter set to use. Examples: 'CO2' or 'wind'.
        currentVal is a list of the numerical values of the wind, CO2 etc. at 
        the location of each mosquito.

        '''
        unscaledSat = self.mosqParams[responseStr+'Sat']
        val = [v/unscaledSat for v in currentVal]
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

    def _inPlume(self,ind):
        '''
        Put in agent rules code.

        This generator determines which mosquitoes are currently inside an odor
        plume.
        currentCO2 is the CO2 level at the current mosquito locations.

        '''
        if self.currentCO2[ind] >= self.mosqParams['CO2Thresh']:
            return True
        else:
            return False
        


class upwind(klinotaxis):
    pass

class downwind(klinotaxis):
    pass

class crosswind(klinotaxis):
    pass


if __name__ == '__main__':
    
    mymosqs = mosquitoPopulation([1.0,2.0,3.0],windKappa=1.0,windThresh=0.01)
    responseStr = 'wind'
    currentVal = [-1.e-12,0.1,1.0+1.e-12]
    print(mymosqs._responseCurve(responseStr,currentVal))
    print(mymosqs.mosqParams)
    
    


