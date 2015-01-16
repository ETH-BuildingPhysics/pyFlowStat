'''
TriSurfaceMesh.py
'''

import re

import numpy as np
import matplotlib.tri as tri

import CoordinateTransformation as coorTrans
import TriSurfaceFunctions


class TriSurfaceMesh(object):
    
    # constructors #
    #--------------#
    def __init__(self,
                 x,
                 y,
                 z,
                 triangles=None,
                 mask=None,
                 affTrans=None,
                 linTrans=None):
        '''
        '''                     
        self.triangulation = tri.Triangulation(x, y, triangles=triangles, mask=mask)
    
      
        # "private" member variable. Don't play with them if you are not sure...
        self.__z = z
        
        self.__affTrans = affTrans
        self.__linTrans = linTrans
    
    @classmethod
    def readFromFoamFile(cls,
                         pointsFile,
                         facesFile,
                         viewAnchor,
                         xViewBasis,
                         yViewBasis,
                         srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):

        '''
        Construct from a surface saved  by OpenFOAM in foamFile format.
        '''
        afftrans, lintrans = TriSurfaceFunctions.getTransformation(viewAnchor,
                                                                   xViewBasis,
                                                                   yViewBasis,
                                                                   srcBasisSrc)
                                       
        # get x and y vector (in ptsTgt)
        ptsSrc = TriSurfaceFunctions.parseFoamFile_sampledSurface(pointsFile)
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])

        #get triangles
        triangles = TriSurfaceFunctions.parseFoamFile_sampledSurface(facesFile)[:,1:4]

        # update class member variables
        return cls(x=ptsTgt[:,0],
                   y=ptsTgt[:,1],
                   z=ptsTgt[:,2],
                   triangles=triangles,
                   mask=None,
                   affTrans=afftrans,
                   linTrans=lintrans)

    @classmethod   
    def readFromVTK(cls,
                    vtkFile,
                    viewAnchor,
                    xViewBasis,
                    yViewBasis,
                    srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
        '''
        Construct from a surface saved by OpenFOAM in VTK format.
        '''
        afftrans, lintrans = TriSurfaceFunctions.getTransformation(viewAnchor,
                                                                   xViewBasis,
                                                                   yViewBasis,
                                                                   srcBasisSrc)
       
        # read VTK file
        ptsSrc, triangles, vecsSrc = TriSurfaceFunctions.parseVTK_ugly_sampledSurface(vtkFile)
        
        # Transform the points
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])

        # update class member variables
        return cls(x=ptsTgt[:,0],
                   y=ptsTgt[:,1],
                   z=ptsTgt[:,2],
                   triangles=triangles,
                   mask=None,
                   affTrans=afftrans,
                   linTrans=lintrans)

    # getters #
    #---------#  
    @property
    def x(self):
        '''
        Get x coordinate of the grid points.
        '''
        return self.triangulation.x
        
    @property    
    def y(self):
        '''
        Get y coordinate of the grid points.
        '''
        return self.triangulation.y
        
    @property
    def triangles(self):
        '''
        Get triangles from the grid.
        '''
        return self.triangulation.triangles
        
    @property
    def affTrans(self):
        '''
        Get affine transformation object.
        '''
        return self.__affTrans
    
    @property
    def linTrans(self):
        '''
        Get linear transformation object.
        '''
        return self.__linTrans
        
        
    # class methods #
    #---------------#
        
    def rawPoints(self):
        '''
        Return the grid points in the source coordinate system.
        
        Returns:
            *rawPoints*: numpy array of shape (N,3)
        '''
        surfacePoints = np.vstack((self.x,self.y,self.__z)).T
        rawPoints = np.zeros((surfacePoints.shape[0],surfacePoints.shape[1]))
        for i in range(surfacePoints.shape[0]):
            rawPoints[i,:] = self.affTrans.tgtToSrc(surfacePoints[i,:])
        return rawPoints                   
