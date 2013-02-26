import numpy as np
from interpFunctions import interpFromGridNumpyArrays 
from interpFunctions import extrapToGrid 

class numericalSims(object):

    def __init__(self,**kwargs):
        '''
        Sets parameters for numerical simulations. Any defaults in simsParams can 
        be overwritten, but not the dimensional parameters.
        Also note that the left and bottom domain edges must both be zero in order
        for the interpolation function interpFromGrid to work properly.

        '''
        self.dimensionalParams = {'mosquitoFlightSpeed (m/s)':1.0,'mosquitoDecisionTime (s)': 0.1,'CO2Sat (units CO2/unit air or 10^6 ppm)':4.e-3}
        self.simsParams = {'domainLength':100.0,'numGridPoints':128}
        self.simsParams.update(kwargs)
        derivedQuantities = {'h':self.simsParams['domainLength']/self.simsParams['numGridPoints']}
        self.simsParams.update(derivedQuantities)
        self._makeGrid(self.simsParams['h'],self.simsParams['domainLength'])


    def _makeGrid(self,h,L):
        self.xg,self.yg = np.mgrid[h/2.0:L:h,h/2.0:L:h]

    def interpFromGrid(self,xy,randVel1,randVel2,CO2):
        ur, vr, c = interpFromGridNumpyArrays(xy,self.simsParams['h'],randVel1,randVel2,CO2)

    def _addSourceTerm(self,environ,CO2):
        cnew = extrapToGrid(environ.hostPositionx,environ.hostPositiony,environ.hostSourceStrength,self.simsParams['h'],CO2.shape)
        return CO2+cnew

    def _upwindScheme(self):
        pass

    def _implicitRK(self):
        pass


