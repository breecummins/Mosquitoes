import numpy as np

def makeGrid(h,L):
    '''
    Makes a cell-centered square grid of length L and spacing h with lower
    left corner at (0,0). This means the first grid point is located at
    (h/2,h/2). The following interpolation and extrapolation functions
    depend on this grid construction.

    '''
    xg,yg = np.mgrid[h/2.0:L:h,h/2.0:L:h]
    return xg,yg

def _getIndicesNodes(x,y,h):
    '''
    Helper function for interpFromGrid and extrapToGrid. Assumes cell-centered
    grid with lower left corner at (0,0); i.e. lower left grid point of (h/2,h/2).
    (x,y) are 2D positions. h is grid spacing.

    '''
    # Find indices of the closest node to (x,y) TO THE LOWER LEFT.
    # Note that int always rounds toward zero.
    i = (x/h - 0.5).astype('int')
    j = (y/h - 0.5).astype('int') 
    # find remainders to assign proportions to nodes
    rx = x/h - 0.5 - i
    ry = y/h - 0.5 - j        
    # Find the proportion of the value at each node that contributes to the 
    # interpolation at (x,y). nodes = (lowerleft, upperleft, lowerright, upperright)
    nodes = np.array([(1-rx)*(1-ry), (1-rx)*ry, rx*(1-ry), rx*ry])
    return i,j,nodes

def interpFromGrid(x,y,h,randVel1,randVel2,CO2):
    '''
    This function interpolates values located at grid nodes to the 
    locations (x,y). 
    x, y are numpy arrays of positions in the x and y directions. 
    h is (scalar) grid spacing.
    randVel* and CO2 are velocity and CO2 values on the grid (np.array).

    '''
    # get indices and proportional values
    i,j,nodes = _getIndicesNodesNumpyArrays(x,y,h)
    # get the values of the CO2 and random wind at the four closest nodes
    V1 = np.array([randVel1[[i,j]],randVel1[[i,j+1]],randVel1[[i+1,j]],randVel1[[i+1,j+1]]]) 
    V2 = np.array([randVel2[[i,j]],randVel2[[i,j+1]],randVel2[[i+1,j]],randVel2[[i+1,j+1]]]) 
    C = np.array([CO2[[i,j]],CO2[[i,j+1]],CO2[[i+1,j]],CO2[[i+1,j+1]]])
    # perform the interpolation
    ur = np.sum(nodes*V1,0)
    vr = np.sum(nodes*V2,0)
    c = np.sum(nodes*C,0)
    return ur,vr,c

def extrapToGrid(x,y,s,h,size):
    '''
    This function extrapolates exhaled host CO2 to the closest 4 grid nodes. 
    x, y are numpy arrays of positions in the x and y directions. 
    s is a numpy array of the host CO2 values at (x,y).
    h is (scalar) grid spacing.
    size is the shape of the 2D computational grid (tuple).

    '''
    # get indices and proportional values
    i,j,nodes = _getIndicesNodesNumpyArrays(x,y,h)
     # calculate additional CO2 at each node
    sarray = np.zeros(size)
    sarray[[i,i,i+1,i+1],[j,j+1,j,j+1]] = nodes*s 
     return sarray

def implicitRK():
    pass

def explicitRK4(t,y,dt,func):
    # 4th order Runge-Kutta solver for dy/dt = func(y,t)
    # func takes args t (scalar time) and y (array values at time t) 
    # dt is time step
    k1 = dt*func(t,y)
    k2 = dt*func(t+dt/2., y+k1/2.)
    k3 = dt*func(t+dt/2., y+k2/2.)
    k4 = dt*func(t+dt, y+k3)
    return y + (k1 + 2*k2 + 2*k3 + k4)/6.0

def forwardEuler(t,y,dt,func):
    return y + dt*func(t,y)

def upwindScheme(U,V,environ):
    CO2 = environ.CO2
    # add ghost cells to velocity arrays
    uxm = np.vstack([environ.leftedge,U[:-1,:]])
    uxp = np.vstack([U[1:,:],environ.rightedge])
    vxm = np.hstack([environ.bottomedge,V[:,:-1]])
    vxp = np.hstack([V[:,1:],environ.topedge])
    # find values at cell edges
    um = 0.5*(U+uxm)
    up = 0.5*(U+uxp)
    vm = 0.5*(V+vxm)
    vp = 0.5*(V+vxp)
    # find wind direction on cell edges
    bool_upp = (up > 0).astype('int')
    bool_upm = (up <= 0).astype('int')
    bool_ump = (um > 0).astype('int')
    bool_umm = (um <= 0).astype('int')
    bool_vpp = (vp > 0).astype('int')
    bool_vpm = (vp <= 0).astype('int')
    bool_vmp = (vm > 0).astype('int')
    bool_vmm = (vm <= 0).astype('int')            
    # find locations for outflow bcs
    leftind = np.nonzero(um[0,:] < 0)
    rightind = np.nonzero(up[-1,:] > 0)
    bottomind = np.nonzero(vm[:,0] < 0)
    topind = np.nonzero(vp[:,-1] > 0)
    # add ghost cells to CO2 array and adjust for outflow bcs
    z = np.zeros(environ.simsParams['numGridPoints'])
    Cxm = np.vstack([z,CO2[:-1,:]])
    Cxp = np.vstack([CO2[1:,:],z])
    Cym = np.hstack([z[:,np.newaxis],CO2[:,:-1]])
    Cyp = np.hstack([CO2[:,1:],z[:,np.newaxis]])
    Cxm[0,leftind] = 2*CO2[0,leftind] - CO2[1,leftind]
    Cxp[-1,rightind] = 2*CO2[-1,rightind] - CO2[-2,rightind]
    Cym[bottomind,0] = 2*CO2[bottomind,0] - CO2[bottomind,1]
    Cyp[topind,-1] = 2*CO2[topind,-1] - CO2[topind,-2]
    # calculate flux term using upwinding scheme
    Flxx=(CO2*bool_upp + Cxp*bool_upm)*up - (Cxm*bool_ump + CO2*bool_umm)*um
    Flxy=(CO2*bool_vpp + Cyp*bool_vpm)*vp - (Cym*bool_vmp + CO2*bool_vmm)*vm
    return (Flxx+Flxy)/environ.simsParams['h']
 
