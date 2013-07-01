'''
TriSurface.py

Class TriSruface holds a planar surface, which hols a discret field. The discret field
can by unstructured.

In the foamFiles, the plane in define in base a={a1,a2,a3} (standard normed orthogonal base)
base

'''

import numpy as np
import scipy as sp


class TriSurface(object):
    """TriSurface Class"""
    
    def __init__(self):
        # in base a (default base from input files)
        # not really needed (storage bloat in some extend)
        self.varsa = []
        self.pointsa = []
        self.facesa = []
        self.origb = np.zeros(3)
        # in base b (base of planar surface, all z of spatial corrdinates are zero)
        self.varsb = []
        self.pointsb = []
        self.facesb = []
        self.origb = np.zeros(3)
        # Affine transformation matrix and it's inverse.
        self.aTransMat = np.zeros([4,4])
        self.invaTransMat = np.zeros([4,4])
        # rotation matrix and translation vector which compose the affine transformation
        self.rotMat = np.zeros([3,3])
        self.vecTrans = np.zeros([3,1])

    def readFromFoamFile(self):
        '''read planar surface from OpenFOAM. In foamFile format.'''
        pass
    
    def readFromVTK(self):
        '''read planar surface from OpenFOAM. In vtk format.'''
        pass
    
    def createTriangle(self):
        '''
        create a Delauney triangles from self.points if not exist. Store the
        triangle in self.faces        
        '''
        pass
    
    def getTransMat(self):
        '''
                get affine transformation matrix
        '''
        pass