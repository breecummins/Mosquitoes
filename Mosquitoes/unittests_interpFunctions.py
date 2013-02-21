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
    ur,vr,c=iF.interpFromGridNumpyArrays(xy,mysim.simsParams['h'],randVel1,randVel2,CO2)
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


if __name__ == '__main__':
    # xy = [(48.32,5.02),(16.94,34.43),(69.50,90.98)]
    # testaccuracy(xy)
    testspeed()