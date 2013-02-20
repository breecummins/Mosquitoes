import numpy as np

def interpFromGridWithMaps(xy,h):
    '''
    This function interpolates values located at grid nodes to the 
    locations (x,y). 
    xy is a list of (x,y) tuples. 
    h is (scalar) grid spacing.

    '''
    # Find indices of the closest node to (x,y) TO THE LOWER LEFT assuming 
    # a cell-centered grid with lower left corner of the domain located at 
    # (0,0); i.e. lowest left-most grid point in the domain is (h/2, h/2).
    ij = map(int,np.fix(map(lamba (a,b): (a/h - 0.5, b/h - 0.5), xy)))

    # find the remainders to do the interpolation
    rxy = map(lambda ((a,b),(c,d)): (a/h - 0.5 - c, b/h -0.5 - d), (xy, ij))
    
    # Find the proportion of the value at each node that contributes to the 
    # interpolation at (x,y). nodes = (lowerleft, lowerright, upperleft, upperright)
    nodes = map(lambda (a,b): ( (1-a)*(1-b), a*(1-b), (1-a)*b, a*b ), rxy)

    # get the values of the CO2 and random wind at the four closest nodes
    Vij = map(lambda matinds: map(lambda m: (randVel1[m], randVel2[m], CO2[m]),matinds),ij,[(p,k+1) for (p,k) in ij],[(p+1,k) for (p,k) in ij], [(p+1,k+1) for (p,k) in ij])

    def interpolate(ind):
        return [nodes[k][0]*Vij[0][k][ind] + nodes[k][1]*Vij[1][k][ind] + nodes[k][2]*Vij[2][k][ind] + nodes[k][3]*Vij[3][k][ind] for k in range(len(Vij[0]))]
    
    # perform the interpolation
    ur = interpolate(0)
    vr = interpolate(1)
    c = interpolate(2)

    return ur,vr,c
