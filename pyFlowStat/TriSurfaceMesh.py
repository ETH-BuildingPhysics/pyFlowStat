'''
TriSurfaceMesh.py
'''

#import re

import numpy as np
import matplotlib.tri as tri

import pyFlowStat.ParserFunctions as ParserFunctions
import pyFlowStat.TriSurface as TriSurface

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
                         xViewBasis,
                         yViewBasis=None,
                         viewAnchor=(0,0,0),
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
            
            *yViewBasis*: python of numpy array. Shape = 3. If it's None, it 
             will be estimated. Works for surfaces normal to one coordinate
             direction
            
            *srcBasisSrc*: python of numpy array. Shape = 3,3. Keep default.
        '''
        
                
        ptsSrc = ParserFunctions.parseFoamFile_sampledSurface(pointsFile)
        if not yViewBasis:
            n=getN(ptsSrc)
            yViewBasis=getYBasis(n,xViewBasis)
            
        afftrans, lintrans = TriSurface.getTransformation(viewAnchor,
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
    def readFromHdf5(cls,
                     hdf5Parser,
                     xViewBasis,
                     yViewBasis=None,
                     viewAnchor=(0,0,0),
                     srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
        '''
        Construct from a Hdf5Parser. The source Hdf5File must have at least the
        following structure:
        
        hdf5file
        {
            mesh
            {
                points: list of points. Shape=(N,3)
                triangles: list of triangles (connectivity). Shape=(M,3)
            }
            otherKey1
            otherKey2
            ...
            
        }
        
        Arguments:
            *hdf5Parser*: a h5py parser.
             the hdf5Parser can by created as follow:
             >>> import h5py
             >>> hdf5Parser = h5py.File('sourceHdf5File','r')
             
            *xViewBasis*: python array of shape=3.
             x direction of the triSurface, defined in the source coordinate
             system.
        
            *yViewBasis*: python array of shape=3.
             y direction of the triSurface, defined in the source coordinate
             system. Default=None
        
            *viewAnchor*: python array of shape=3.
             Origin of the surface coordinate system, defined in the source
             coordinate system. Default=(0,0,0)
             
            *srcBasisSrc*: python array of shape=3x3.
             Default=[[1,0,0],[0,1,0],[0,0,1]].
             
        Usage:
             >>> import h5py
             >>> from pyFlowStat.TriSurfaceMesh import TriSurfaceMesh
             >>> hdf5Parser = h5py.File('sourceHdf5File','r')
             >>> tsm = TriSurfaceMesh.readFromHdf5(hdf5Parser,xViewBasis=(1,0,0))
             >>> hdf5Parser.close()
        
        '''
        # get mest data
        gName = 'mesh'
        ptsSrc = hdf5Parser[gName]['points'].value
        triangles = hdf5Parser[gName]['faces'].value
        
        # create the transformation objects
        if not yViewBasis:
            n = getN(ptsSrc)
            yViewBasis = getYBasis(n,xViewBasis)
    
        afftrans, lintrans = TriSurface.getTransformation(viewAnchor=viewAnchor,
                                                          xViewBasis=xViewBasis,
                                                          yViewBasis=yViewBasis,
                                                          srcBasisSrc=srcBasisSrc)
    
        # transform the points from the source basis to the target basis
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
        afftrans, lintrans = TriSurface.getTransformation(viewAnchor,
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


    def area(self):
        '''
        Calculate and return the area of the surface. This algorithm uses
        the Green theorem. From it, one can derived a general formula to
        calculate the area A of any simple polygon of n summits:
        
        $$A = \frac{(x_{0}+x_n)(y_{0}-y_n)}{2} + \sum_{i=0}^{n-1}\frac{(x_{i+1}+x_i)(y_{i+1}-y_i)}{2}$$
        
        Arguments:
            *none*
            
        Returns
            *area*: python float()
             The total area of the mesh.
        '''
        xtri = self.x[self.triangles]
        ytri = self.y[self.triangles]
        
        return np.sum( 0.5*( (xtri[:,1]+xtri[:,0])*(ytri[:,1]-ytri[:,0])
                            +(xtri[:,2]+xtri[:,1])*(ytri[:,2]-ytri[:,1])
                            +(xtri[:,0]+xtri[:,2])*(ytri[:,0]-ytri[:,2])) )

                            
# helper functions #
#------------------#      
                      
def getN(rawPoints):
    '''
    get normal of plane
    
    *rawPoints*: numpy array of shape (N,3)
    
    returns n, numpy array of shape 3
    '''
    norm=0
    while norm==0:
        idx_points=np.random.randint(0,len(rawPoints),size=3)
        p1=rawPoints[idx_points[0]]
        p2=rawPoints[idx_points[1]]
        p3=rawPoints[idx_points[2]]
        v1=p1-p2
        v2=p1-p3
        n=np.cross(v1,v2)
        norm=np.linalg.norm(n)
    return np.round(np.abs(n)/norm)

def getYBasis(n,xViewBasis):
    '''
    get y basis vector from x basis and plane normal
    
    *n*: python of numpy array. Shape = 3.
    
    *xViewBasis*: python of numpy array. Shape = 3.
    
    returns yViewBasis, python of numpy array. Shape = 3.
    '''
    yViewBasis=np.cross(n,xViewBasis)
    norm=np.linalg.norm(yViewBasis)
    if norm==0:
        raise ValueError("n and xViewBasis are collinar")
    return np.round(np.abs(yViewBasis)/norm)

def getBasis(rawPoints,zIsVertical=True):
    '''
    fully automatic basis estimation (currently unused)
    
    *rawPoints*: numpy array of shape (N,3)
    
    returns xViewBasis,yViewBasis
    '''
    n=getN(rawPoints)
    xViewBasis=np.array([1.0,0.0,0.0])
    try:
        yViewBasis=getYBasis(n,xViewBasis)
    except:
        if zIsVertical:
            xViewBasis=np.array([0.0,1.0,0.0])
        else:
            xViewBasis=np.array([0.0,0.0,1.0])
        yViewBasis=getYBasis(n,xViewBasis)

    return xViewBasis,yViewBasis