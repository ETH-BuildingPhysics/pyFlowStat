'''
TriSurfaceMesh.py
'''

import numpy as np
import matplotlib.tri as tri

import pyFlowStat.CoordinateTransformation as coorTrans
import pyFlowStat.ParserFunctions as ParserFunctions


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
        
        Arguments:
            *pointsFile*: python string.
             Location of the point file. It holds the points (summits) of
             the triangulate grid.
            
            *facesFile*: python string.
             Location of the faces file. It holds the triangles of the grid.
             If facesFile=None, the triangles are created with the Delauney
             method.
             
            *viewAnchor*: python of numpy array. Shape = 3.
             Location of the (0,0) point of your grid. Defined in the source
             basis (meaning: the coordinate system of the OpenFOAM simulation)
             
            *xViewBasis*: python of numpy array. Shape = 3.
            
            *yViewBasis*: python of numpy array. Shape = 3.
            
            *srcBasisSrc*: python of numpy array. Shape = 3,3.
        '''
        afftrans, lintrans = getTransformation(viewAnchor,
                                               xViewBasis,
                                               yViewBasis,
                                               srcBasisSrc)
                                       
        # get x and y vector (in ptsTgt)
        ptsSrc = ParserFunctions.parseFoamFile_sampledSurface(pointsFile)
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])

        #get triangles
        if facesFile==None:
            triangles = None
        else:
            triangles = ParserFunctions.parseFoamFile_sampledSurface(facesFile)[:,1:4]

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
        afftrans, lintrans = getTransformation(viewAnchor,
                                               xViewBasis,
                                               yViewBasis,
                                               srcBasisSrc)
       
        # read VTK file
        ptsSrc, triangles, vecsSrc = ParserFunctions.parseVTK_ugly_sampledSurface(vtkFile)
        
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


def getTransformation(viewAnchor,
                      xViewBasis,
                      yViewBasis,
                      srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
    '''
    Return the affine and linear transfomation object for coordinate
    transformation from a source coordinate system S to a target coordinate 
    system T.
    
    Arguments:
        *viewAnchor*: python list or numpy array of shape (3,).
         Location of the origin of T, defined in S.
         
        *xViewBasis*: python list or numpy array of shape (3,).
         x axis of T, defined in S.
        
        *yViewBasis*: python list or numpy array of shape (3,).
         y axis of T, defined in S.
        
        *srcBasisSrc*: python list or numpy array of shape (3,3).
         x,y,z axis of S, defined in S. For advenced user only.
         Default=[[1,0,0],[0,1,0],[0,0,1]]
          
    Returns:
        *affTrans*: AffineTransformation object.
            
        *linTrans*: LinearTransformation object.
    '''
    # check and convert arguments
    srcBasisSrc = np.array(srcBasisSrc,dtype=float)
    if srcBasisSrc.shape!=(3,3):
        raise ValueError('srcBasis must be a 3x3 matrix')
        
    xViewBasis = np.array(xViewBasis,dtype=float)
    yViewBasis = np.array(yViewBasis,dtype=float)
    if xViewBasis.shape!=(3,) or yViewBasis.shape!=(3,):
        raise ValueError('xViewBasis.shape and yViewBasis. ',
                         'shape must be equal to (3,)')
        
    # get the basis and the transformation object
    tgtBasisSrc = np.zeros((3,3))
    tgtBasisSrc[:,0] = xViewBasis
    tgtBasisSrc[:,1] = yViewBasis
    tgtBasisSrc[:,2] = np.cross(xViewBasis,yViewBasis)
    afftrans = coorTrans.AffineTransformation(srcBasisSrc,tgtBasisSrc,viewAnchor)
    lintrans = coorTrans.LinearTransformation(srcBasisSrc,tgtBasisSrc)
    
    return afftrans, lintrans