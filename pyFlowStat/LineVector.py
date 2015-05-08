'''
LineVector.py
'''

#import re

import numpy as np

import pyFlowStat.Line as Line


class LineVector(Line.Line):
    '''
    '''
    # constructors #
    #--------------#
    
    def __init__(self,
                 xyz,
                 vx,
                 vy,
                 vz):
        '''
        base constructor.
        
        Arguments:
            *xyz*: numpy array of shape (npoints,3).
             position vector
            *s*: numpy array of shape (npoints).
             scalar values.
        '''
        super(LineVector,self).__init__(xyz)

        self.vx = np.asarray(vx)
        self.vy = np.asarray(vy)
        self.vz = np.asarray(vz)
        
        
    # class methods #
    #---------------#
    def __call__(self,dim=0):
        return self.component(dim=dim)
        
    def component(self,dim):
        if dim==0:
            return self.vx
        if dim==1:
            return self.vy
        if dim==2:
            return self.vz
        
    def rawVars(self):
        '''
        Return the scalar field defined in the source coordinate system.
        
        Returns:
            *rawData*: numpy array of shape (N,)
        '''
        return np.vstack((self.vx,self.vy,self.vz)).T