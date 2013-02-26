import numpy as np
import lib_numericalMethods as nMeth

def constantVel(x,y,xmag=0.0,ymag=0.2):
    '''
    Constant velocity needs to be slower than mosquito flight speed of 1.0 (nondimensional).

    '''
    return xmag*np.ones(x.shape), ymag*np.ones(y.shape)

def constantSourceStrength(dimParams,numHosts):
    '''
    Constant CO2 emission per second from a chicken ((units CO2/unit air)/s) (Brackenbury 1982 notes): ~0.1 ml CO2 / ml air / min for 2 kg chicken.

    '''
    return np.array([(0.1/60)*(dimParams['mosquitoDecisionTime (s)']/dimParams['CO2Sat (units CO2/unit air or 10^6 ppm)'])]*numHosts)


class environment(object):
    '''
    This class represents the environment in which the mosquitoes fly, and 
    the numerical view (grid) into the environment.

    Sets parameters for numerical simulations. Any defaults in simsParams can 
    be overwritten, but not the dimensional parameters.
    Also note that the left and bottom domain edges must both be zero in order
    for the interpolation function interpFromGrid to work properly.

    '''

    def __init__(self,numericalSims,hostPositionx,hostPositiony,hostSourceHandle=constantSourceStrength,velocityFunctionHandle=constantVel,**kwargs):
        # host positions and params (add the other params)
        self.hostPositionx = hostPositionx #numpy array of x positions
        self.hostPositiony = hostPositiony #numpy array of y positions
        self.hostSourceStrength = hostSourceHandle(numericalSims.dimensionalParams,len(hostPositionx)) 
        # velocity parameters
        self.velfunc = velocityFunctionHandle
        self.randSeeds = map(int,203863455938475394857*np.random.rand(100000))
        # dimensional parameters to interpret results (code is nondimensional)
        self.dimensionalParams = {'mosquitoFlightSpeed (m/s)':1.0,'mosquitoDecisionTime (s)': 0.1,'CO2Sat (units CO2/unit air or 10^6 ppm)':4.e-3}
        # numerical parameters for the simulation, may be overwritten with kwargs
        self.simsParams = {'domainLength':100.0,'numGridPoints':128}
        self.simsParams.update(kwargs)
        derivedQuantities = {'h':self.simsParams['domainLength']/self.simsParams['numGridPoints']}
        self.simsParams.update(derivedQuantities)
        self.xg, self.yg = nMeth.makeGrid(self.simsParams['h'],self.simsParams['domainLength'])

    def getSignal(self,x,y,randVel1,randVel2,CO2):
        '''
        This function returns three arrays: u,v,c for every (x,y) pair. 
        x and y are arrays of the same length denoting mosquito position in 2D.
        randVel* and CO2 are the grid values for the random velocity components
        and the CO2.
        
        '''       
        # Get bulk flow wind and background CO2
        u,v = self.velfunc(x,y)
        c = np.zeros(x.shape)
        # Get random velocities and CO2 inside domain
        # Assume domain is square with lower left corner at (0,0) and is cell-centered
        L = self.simsParams['domainLength']
        h = self.simsParams['h']
        insideDom = [k for k in range(len(x)) if x[k] < (L-h/2.0) and x[k] > h/2.0 and y[k] < (L-h/2.0) and y[k] > h/2.0]
        ur,vr,ci = nMeth.interpFromGrid(x[insideDom],y[insideDom],self.simsParams['h'],randVel1,randVel2,CO2)  
        # Add interpolated values to bulk values
        u[insideDom] = u[insideDom] + ur
        v[insideDom] = v[insideDom] + vr       
        c[insideDom] = c[insideDom] + ci     
        return u,v,CO2

    def updateEnvironment(self):
        pass


if __name__ == '__main__':
    pass