'''
LineScalar.py
'''

#import re

import numpy as np

import pyFlowStat.Line as Line


class LineScalar(Line.Line):
    '''
    '''
    # constructors #
    #--------------#
    
    def __init__(self,
                 xyz,
                 s):
        '''
        base constructor.
        
        Arguments:
            *xyz*: numpy array of shape (npoints,3).
             position vector
             
            *s*: numpy array of shape (npoints).
             scalar values.
        '''
        super(LineScalar,self).__init__(xyz)

        self.s = np.asarray(s)
        
    # class methods #
    #---------------#
    def __call__(self,dim=0):
        return self.component(dim=dim)
        
    def component(self,dim=0):
        return self.s
        
    def rawVars(self):
        '''
        Return the scalar field defined in the source coordinate system.
        
        Returns:
            *rawData*: numpy array of shape (N,)
        '''
        return self.s