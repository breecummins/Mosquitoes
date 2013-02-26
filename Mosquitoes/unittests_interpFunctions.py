import numericalSims as nS
import interpFunctions as iF
import numpy as np

def testaccuracy(xy):
    mysim = nS.numericalSims()
    # bilinear/linear/affine functions should be recovered exactly
    randVel1 = 0.03*mysim.xg + 0.1
    randVel2 = -0.02*mysim.yg
    CO2 = mysim.xg + 2.5*mysim.yg 
    # exact values at (x,y) locations
    ur_exact = [0.03*z[0] + 0.1 for z in xy]
    vr_exact = [-0.02*z[1] for z in xy]
    c_exact = [z[0] + 2.5*z[1] for z in xy] 
    # call interpolation functions and calculate error
    ur,vr,c=iF.interpFromGridWithMaps(xy,mysim.simsParams['h'],randVel1,randVel2,CO2)
    mapserr = [[abs(ur[k]-ur_exact[k]) for k in range(len(xy))],[abs(vr[k]-vr_exact[k]) for k in range(len(xy))],[abs(c[k]-c_exact[k]) for k in range(len(xy))]]
    print('Error in maps method:')
    print(max(max(mapserr)))
    ur,vr,c=iF.interpFromGridSingleForLoop(xy,mysim.simsParams['h'],randVel1,randVel2,CO2)
    forlooperr=[[abs(ur[k]-ur_exact[k]) for k in range(len(xy))],[abs(vr[k]-vr_exact[k]) for k in range(len(xy))],[abs(c[k]-c_exact[k]) for k in range(len(xy))]]
    print('Error in for loop method:')
    print(max(max(forlooperr)))
    # List comprehension is tremendously faster than casting a list of tuples 
    # straight to an array.
    x = np.array([z[0] for z in xy]) 
    y = np.array([z[1] for z in xy])
    ur,vr,c=iF.interpFromGridNumpyArrays(x,y,mysim.simsParams['h'],randVel1,randVel2,CO2)
    nparrayerr=[[abs(ur[k]-ur_exact[k]) for k in range(len(xy))],[abs(vr[k]-vr_exact[k]) for k in range(len(xy))],[abs(c[k]-c_exact[k]) for k in range(len(xy))]]
    print('Error in numpy array method:')
    print(max(max(nparrayerr)))
    ur,vr,c=iF.interpFromGridListComp(xy,mysim.simsParams['h'],randVel1,randVel2,CO2)
    listcomperr=[[abs(ur[k]-ur_exact[k]) for k in range(len(xy))],[abs(vr[k]-vr_exact[k]) for k in range(len(xy))],[abs(c[k]-c_exact[k]) for k in range(len(xy))]]
    print('Error in list comp method:')
    print(max(max(listcomperr)))

def testspeed():
    '''
    To test speed, import cProfile and testspeed from this module into a python 
    interpreter and do
    cProfile.runctx('testspeed()',globals(),locals())

    '''
    for k in range(10):
        x = 1+98*np.random.rand(10000)
        y = 1+98*np.random.rand(10000)
        xy = zip(x,y)
        testaccuracy(xy)

def testextrap(x,y):
    mysim = nS.numericalSims()
    # bilinear/linear/affine functions should be recovered exactly
    s = np.random.rand(len(x))
    # print('s',s)
    # call extrapolation functions and calculate error
    c=iF.extrapToGrid(x,y,s,mysim.simsParams['h'],mysim.xg.shape)
    i,j,nodes = iF.getIndicesNodesNumpyArrays(x,y,mysim.simsParams['h'])
    checksum = c[[i,j]] + c[[i,j+1]] + c[[i+1,j]] + c[[i+1,j+1]]
    # print('checksum',checksum)
    # check that max extrap val occurs at min dist
    match = []
    for k in range(len(x)):
        vals = [c[i[k],j[k]],c[i[k],j[k]+1],c[i[k]+1,j[k]],c[i[k]+1,j[k]+1]]
        valinds = [a[0] for a in sorted(enumerate(vals), key=lambda z:z[1])]
        def mydist(k,inds):
            dist = ( (x[k] - mysim.xg[inds])**2 + (y[k] - mysim.yg[inds])**2 )**0.5
            if dist >= mysim.simsParams['h'] * (2.0**0.5):
                print('Distance out of bounds.')
            return dist
        dists = [mydist(k,inds) for inds in [(i[k],j[k]),(i[k],j[k]+1),(i[k]+1,j[k]),(i[k]+1,j[k]+1)]]
        distinds = [a[0] for a in sorted(enumerate(dists), key=lambda z:z[1], reverse=True)]
        match.append(valinds == distinds)
        if not match[-1]:
            print('vals',vals)
            print('valinds',valinds)
            print('dists',dists)
            print('distinds',distinds)
            print('(x,y)',(x[k],y[k]))
            print('(xg,yg)',(mysim.xg[i[k],j[k]],mysim.yg[i[k],j[k]]),(mysim.xg[i[k],j[k]+1],mysim.yg[i[k],j[k]+1]),(mysim.xg[i[k]+1,j[k]],mysim.yg[i[k]+1,j[k]]),(mysim.xg[i[k]+1,j[k]+1],mysim.yg[i[k]+1,j[k]+1]))
    print('Error in extrapolation:')
    print(np.max(np.max(np.array([np.abs(checksum[k] - s[k]) for k in range(len(x))]))))
    print('Does the max value occur at the closest node? (It should.)')
    if all(match):
        print('yes')
    else:
        chkinds = np.nonzero(checksum < 1.e-10)
        print('No. Number of mismatches is {} and values of extrapolation are {}'.format(len(match)-sum(match)),checksum[chkinds])



if __name__ == '__main__':
    xy = [(48.32,5.02),(16.94,34.43),(69.50,90.98)]
    # testaccuracy(xy)
    # testspeed()
    np.random.seed(3123)
    x = 1+98*np.random.rand(100)
    np.random.seed(84379)
    y = 1+98*np.random.rand(100)
    # x = np.array([x for (x,y) in xy])
    # y = np.array([y for (x,y) in xy])
    testextrap(x,y)