import numpy as np

class numericalSims(object):

    def __init__(self,**kwargs):
        '''
        Sets parameters for numerical simulations. Any defaults in simsParams can 
        be overwritten, but not the dimensional parameters.
        Also note that the left and bottom domain edges must both be zero in order
        for the interpolation function interpFromGrid to work properly.

        '''
        self.dimensionalParams = {'mosquitoFlightSpeed m/s':1.0}
        self.simsParams = {'domainLength':100.0,'numGridPoints':128}
        self.simsParams.update(kwargs)
        derivedQuantities = {'h':self.simsParams['domainLength']/self.simsParams['numGridPoints']}
        self.simsParams.update(derivedQuantities)
        self._makeGrid(self.simsParams['h'],self.simsParams['domainLength'])


    def _makeGrid(self,h,L):
        self.xg,self.yg = np.mgrid[h/2.0:L:h,h/2.0:L:h]

