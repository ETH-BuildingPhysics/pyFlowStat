'''
LineSymmTensor.py
'''

#import re

import numpy as np

import pyFlowStat.Line as Line


class LineSymmTensor(Line.Line):
    '''
    '''
    # constructors #
    #--------------#
    
    def __init__(self,
                 xyz,
                 txx,
                 txy,
                 txz,
                 tyy,
                 tyz,
                 tzz):
        '''
        base constructor.
        
        Arguments:
            *xyz*: numpy array of shape (npoints,3).
             position vector

        '''
        super(LineSymmTensor,self).__init__(xyz)

        self.txx = np.asarray(txx)
        self.txy = np.asarray(txy)
        self.txz = np.asarray(txz)
        self.tyy = np.asarray(tyy)
        self.tyz = np.asarray(tyz)
        self.tzz = np.asarray(tzz)
        
    # class methods #
    #---------------#
    def __call__(self,dim=0):
        return self.component(dim=dim)
        
    def component(self,dim):
        if dim==0:
            return self.txx
        if dim==1:
            return self.txy
        if dim==2:
            return self.txz
        if dim==3:
            return self.tyy
        if dim==4:
            return self.tyz
        if dim==5:
            return self.tzz
        
    def rawVars(self):
        '''
        Return the scalar field defined in the source coordinate system.
        
        Returns:
            *rawData*: numpy array of shape (N,6)
        '''
        return np.vstack((self.txx,self.txy,self.txz,self.tyy,self.tyz,self.tzz)).T