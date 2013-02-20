import numpy as np
import numericalSims

class environment(object):
    '''
    This class represents the environment in which the mosquitoes fly.

    '''
    def __init__(self,velocityFunctionHandle=lambda x,y: (0,0.2)):
        self.velfunc = velocityFunctionHandle
        self.randSeeds = map(int,203863451938475394857*np.random.rand(100000))

    def getSignal(self,x,y,numericalSims):
        '''
        This function returns three lists: u,v,CO2 for every (x,y) pair. 
        x and y are lists of the same length denoting position in 2D.
        bottom, top, left, and right are the lower y, upper y, lower x, 
        and upper x scalar values of the rectangular domain, respectively.

        '''
        
        # Get bulk flow wind 
        U = map(self.velfunc,x,y)
        #Can I avoid the unpacking? Should I just have lists of tuples all the time?
    
        # Get background CO2
        CO2 = [0]*len(x)

        # Get random velocities and CO2 inside domain
        # Assume domain is square with lower left corner at (0,0) and is cell-centered
        L = numericalSims.simsParams['domainLength']
        h = numericalSims.simsParams['h']
        insideDom = [k for k in range(len(x)) if self.currentPosx[k] < (L-h/2.0) and self.currentPosx[k] > h/2.0 and self.currentPosy[k] < (L-h/2.0) and self.currentPosy[k] > h/2.0]
        ixy = [(self.currentPosx[k],self.currentPosy[k]) for k in insideDom]
        ur,vr,c = self.interpFromGrid(ixy)  

        # Add interpolated values to bulk values
        U = [(U[k][0]+ur[k],U[k][1]+vr[k]) for k in insideDom]       
        CO2 = [CO2[k]+c[k] for k in insideDom]     

        return U,CO2


if __name__ == '__main__':
    myenv = environment()
    mysim = numericalSims()
    u,v,CO2 = myenv.getSignal(23.4,72.8,mysim)
