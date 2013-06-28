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
        # in base a
        self.varsa = []
        self.pointsa = []
        self.facesa = []
        self.origb = np.zeros(3)
        # in base b
        self.varsb = []
        self.pointsb = []
        self.facesb = []
        self.origb = np.zeros(3)
        # in base b, without last dimension
        self.vars = []
        self.points = []
        self.faces = []
        self.orig = np.zeros(2)
        
        self.tranMat = np.zeros([3,3])
        self.invtranMat = np.zeros([3,3])

    def readFromFoamFile(self):
        pass
    
    def readFromVTK(self):
        pass
    
    def createTriang(self):
        pass