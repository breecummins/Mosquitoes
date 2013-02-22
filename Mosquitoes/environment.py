import numpy as np

def constantVel(x,y,xmag=0.0,ymag=0.2):
    return xmag*np.ones(x.shape), ymag*np.ones(y.shape)

class environment(object):
    '''
    This class represents the environment in which the mosquitoes fly.

    '''
    def __init__(self,hostPositions,velocityFunctionHandle=constantVel):
        self.hostPositions = hostPositions
        self.velfunc = velocityFunctionHandle
        self.randSeeds = map(int,203863455938475394857*np.random.rand(100000))
        # need host params

    def getSignal(self,x,y,simulation):
        '''
        This function returns three arrays: u,v,CO2 for every (x,y) pair. 
        x and y are arrays of the same length denoting position in 2D.
        simulation is an instance of the numericalSims class.

        '''
        
        # Get bulk flow wind and background CO2
        u,v = self.velfunc(x,y)
        CO2 = np.zeros(x.shape)

        # Get random velocities and CO2 inside domain
        # Assume domain is square with lower left corner at (0,0) and is cell-centered
        L = numericalSims.simsParams['domainLength']
        h = numericalSims.simsParams['h']
        insideDom = [k for k in range(len(x)) if self.currentPosx[k] < (L-h/2.0) and self.currentPosx[k] > h/2.0 and self.currentPosy[k] < (L-h/2.0) and self.currentPosy[k] > h/2.0]
        inx = [self.currentPosx[k] for k in insideDom]
        iny = [self.currentPosy[k] for k in insideDom]
        ur,vr,c = self.interpFromGrid(inx,iny)  

        # Add interpolated values to bulk values
        u[insideDom] = u[insideDom] + ur
        v[insideDom] = v[insideDom] + vr       
        CO2[insideDom] = CO2[insideDom] + c     

        return u,v,CO2


if __name__ == '__main__':
    pass