class environment(object):
    '''
    This class represents the environment in which the mosquitoes fly.

    '''
    def __init__(self):
        pass

    def getSignal(self,x,y):
        pass

    def _outsideDomain(self,ind,bottom,top,left,right):
        '''
        Put in environment.getSignal

        This generator determines which mosquitoes are currently outside the 
        CO2 computational domain. 
        bottom, top, left, and right are the lower y, upper y, lower x, and 
        upper x scalar values of the rectangular domain, respectively.

        '''
        valx = self.currentPosx[ind]
        valy = self.currentPosy[ind]
        if valx > right or valx < left or valy > top or valy < bottom:
            return True
        else:
            return False
            





