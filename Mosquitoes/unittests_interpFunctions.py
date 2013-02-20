import numericalSims as nS
import interpFunctions as iF

xy = [(48.32,5.02),(16.94,34.43),(69.50,90.98)]

def testaccuracy():
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

def testspeed():
    pass

if __name__ == '__main__':
    testaccuracy()
    testspeed()