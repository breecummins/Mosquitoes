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

def addSourceTerm(environ,CO2):
    '''
    Extrapolate exhaled CO2 to grid and add it to current CO2.

    '''
    cnew = extrapToGrid(environ.hostPositionx,environ.hostPositiony,environ.hostSourceStrength,simsParams['h'],CO2.shape)
    return CO2+cnew

def upwindScheme():
    pass

def implicitRK():
    pass

def explicitRK():
    pass


