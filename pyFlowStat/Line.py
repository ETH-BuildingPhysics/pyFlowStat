'''
Line.py
'''


import numpy as np



class Line(object):
    '''
 
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,
                 xyz):
                     
        self.xyz=xyz

    
    # class methods #
    #---------------#
    def component(self,dim):
        '''
        '''
        raise NotImplementedError('Line subclasses should implement component.')

 
    def rawVars(self):
        '''
        '''
        raise NotImplementedError('Line subclasses should implement rawVars.')
    
class LineDict(dict): 
    '''
    A class which has the same behavior than a classic python dict(), but it
    can be filled only with object of type Line.
    
    Usage:
        TODO: change
        >>> lDict = LineDict()
        >>> tsDict['U'] = myTriSurfaceVectorObject
        >>> tsDict['T'] = myTriSurfaceScalarObject
        >>> tsDict['afloat'] = 12,927
        TypeError: item is not of type "Line"
    '''      
    def __setitem__(self, key, item):
        if isinstance(item,Line):
            super(LineDict,self).__setitem__(key, item)
        else:
            raise TypeError('Item is not of type "Line"')  