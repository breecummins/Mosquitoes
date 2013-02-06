import numpy as np

class mosquitoPopulation(object):

    '''
    This class instantiates a group of mosquito agents with the same parameters. 
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

    def _responseCurve(self,responseStr,currentVal):
        '''
        This method determines strength of mosquito response to a signal.
        responseStr is the string that identifies which thresh, sat, and kappa 
        parameter set to use. Examples: 'CO2' or 'wind'.
        currentVal is a list of the numerical values of the wind, CO2 etc. at 
        the location of each mosquito.

        '''
        unscaledThresh = self.mosqParams[responseStr+'Thresh']
        unscaledSat = self.mosqParams[responseStr+'Sat']
        val = [v/unscaledSat for v in currentVal]
        thresh = unscaledThresh/unscaledSat
        sat = 1.0
        kappa = self.mosqParams[responseStr+'Kappa']
        if kappa <= -1.0/thresh:
            raise ValueError(responseStr+'Kappa must be > -1.0 / %0.3f' %thresh)
        def response(v):
            if v <= thresh:
                return 0.0
            elif v >= sat:
                return 1.0
            else:
                return (1.0+kappa*thresh)*(v - thresh)/(1.0+kappa*thresh*v*(1.0-thresh))
        return map(response,val)

    def _outsideDomain(self,bottom,top,left,right):
        '''
        This method determines which mosquitoes are currently outside the 
        CO2 computational domain.

        '''
        outsidex = [1.0 if x > right or x < left else 0.0 for x in self.currentPosx]
        outsidey = [1.0 if y > top or y < bottom else 0.0 for x in self.currentPosy]
        return [(outsidex[k] or outsidey[k]) for k in range(len(outsidex))]
    
    def _plumeFinding(self):
        '''
        Stub to be defined by subclass
        
        '''
        return None

    def _plumeTracking(self):
        '''
        Stub to be defined by subclass
        
        '''
        return None

if __name__ == '__main__':
    print('Working')


