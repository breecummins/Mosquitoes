import numpy as np
import lib_numericalMethods as nMeth

def constantVel(x,y,xmag=0.0,ymag=0.2):
    '''
    Constant velocity needs to be slower than mosquito flight speed of 1.0 (nondimensional). Return components in x and y directions.
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

    def __init__(self,hostPositionx,hostPositiony,hostSourceHandle=constantSourceStrength,velocityFunctionHandle=constantVel,**kwargs):
        # host positions and params 
        self.hostPositionx = hostPositionx #numpy array of x positions
        self.hostPositiony = hostPositiony #numpy array of y positions
        self.hostSourceStrength = hostSourceHandle(numericalSims.dimensionalParams,len(hostPositionx)) 
        # velocity parameters
        self.velfunc = velocityFunctionHandle
        self.randSeeds = map(int,203863455938475394857*np.random.rand(100000))
        self.randVelMag = 0.375*0.2 #needs to be smaller than bulk flow
        # dimensional parameters to interpret results (code is nondimensional)
        self.dimensionalParams = {'mosquitoFlightSpeed (m/s)':1.0,'mosquitoDecisionTime (s)': 0.1,'CO2Sat (units CO2/unit air or 10^6 ppm)':4.e-3}
        # numerical parameters for the simulation, may be overwritten with kwargs
        # dt must be 1.0/N where N is an integer, so that the mosquito decisions 
        # occurring every 1.0 happen at a time step boundary.
        # finalTime must be long enough to allow all crosswind mosquitoes to find a host.
        # numGridPoints must be chosen so that grid spacing is on the order of a 
        # single mosquito flight
        self.simsParams = {'domainLength':100.0,'numGridPoints':128,'initialTime':0.0,'finalTime':5000.0,'dt':1.0/10,'randVelSwitch':20.0}
        self.simsParams.update(kwargs)
        derivedQuantities = {'h':self.simsParams['domainLength']/self.simsParams['numGridPoints'],'numRandVelTimeSteps':self.simsParams['randVelSwitch']/self.simsParams['dt']}
        self.simsParams.update(derivedQuantities)
        # grid and grid quantities
        self.xg, self.yg = nMeth.makeGrid(self.simsParams['h'],self.simsParams['domainLength'])
        self.CO2 = np.zeros(self.xg.shape)
        self.randVel1 = np.zeros(self.xg.shape) 
        self.randVel2 = np.zeros(self.xg.shape)
        self.constantU, self.constantV = self.velfunc(self.xg,self.yg) #will have to be changed for time dependent velocity
        self.constantSource = nMeth.extrapToGrid(self.hostPositionx,self.hostPositiony,self.hostSourceStrength,self.simsParams['h'],self.xg.shape) #will have to change for space or time dependent host breathing

    def _setHeavisideRandVel(self,ind):
        np.random.seed(self.randSeeds[ind])
        self.randVel1 = self.randVelMag*np.random.randn(self.simsParams['numGridPoints'],self.simsParams['numGridPoints'])
        self.randVel2 = self.randVelMag*np.random.randn(self.simsParams['numGridPoints'],self.simsParams['numGridPoints'])
        
    def _setContinuousRandomVel(self,ind):
        np.random.seed(self.randSeeds[ind])
        self.randVel1n = self.randVelMag*np.random.randn(self.simsParams['numGridPoints'],self.simsParams['numGridPoints'])
        self.randVel2n = self.randVelMag*np.random.randn(self.simsParams['numGridPoints'],self.simsParams['numGridPoints'])
        np.random.seed(self.randSeeds[ind+1])
        self.randVel1np1 = self.randVelMag*np.random.randn(self.simsParams['numGridPoints'],self.simsParams['numGridPoints'])
        self.randVel2np1 = self.randVelMag*np.random.randn(self.simsParams['numGridPoints'],self.simsParams['numGridPoints'])

    def _continuousRandVel(self,timestep):
        ratio = timestep / self.simsParams['numRandVelTimeSteps']
        self.randVel1 = self.randVel1n + ratio * (self.randVel1np1 - self.randVel1n)
        self.randVel2 = self.randVel2n + ratio * (self.randVel2np1 - self.randVel2n)
        
    def querySignal(self,x,y):
        '''
        This function returns three arrays: u,v,c for every (x,y) pair. 
        x and y are arrays of the same length denoting mosquito position in 2D.
        
        '''       
        # Get bulk flow wind and background CO2
        u,v = self.velfunc(x,y)
        c = np.zeros(x.shape)
        # Get random velocities and CO2 inside domain
        # Assume domain is square with lower left corner at (0,0) and is cell-centered
        L = self.simsParams['domainLength']
        h = self.simsParams['h']
        insideDom = [k for k in range(len(x)) if x[k] < (L-h/2.0) and x[k] > h/2.0 and y[k] < (L-h/2.0) and y[k] > h/2.0]
        ur,vr,ci = nMeth.interpFromGrid(x[insideDom],y[insideDom],self.simsParams['h'],self.randVel1,self.randVel2,self.CO2)  
        # Add interpolated values to bulk values
        u[insideDom] = u[insideDom] + ur
        v[insideDom] = v[insideDom] + vr       
        c[insideDom] = c[insideDom] + ci     
        return u,v,c

    def updateEnvironment(self,currentTime):
        self.CO2 = nMeth.explicitRK4(currentTime,self.CO2,self.simsParams['dt'],self._updateCO2HeavisideRandVel)

    def _updateCO2HeavisideRandVel(self,t,CO2):
        ind,rem = divmod(t,self.simsParams['randVelSwitch'])
        if rem < self.simsParams['dt']/2:
            self._setHeavisideRandVel(ind)
        #FIXME: do upwind scheme
        #add self.constantSource to the velocity term

    def _updateCO2ContinuousRandVel(self,t,CO2):
        ind,rem = divmod(t,self.simsParams['randVelSwitch'])
        if rem < self.simsParams['dt']/2:
            self._setContinuousRandVel(ind)
        timestep,_ = divmod(rem,self.simsParams['dt'])
        self._continuousRandVel(timestep)
        #FIXME: do upwind scheme
        #add self.constantSource to the velocity term

if __name__ == '__main__':
    pass