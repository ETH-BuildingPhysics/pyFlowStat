'''
TriSurfaceVector.py
'''

import re

import numpy as np
import matplotlib.tri as tri

import CoordinateTransformation as coorTrans


class TriSurfaceVector(object):
    '''
    class TriSurfaceVector.
    
    Class which describes a 3D vector field V=(vx,vy,vz) on a 2D (flat)
    triangulated surface S of N points. The triangulation is a grid of M 
    triangles. The x coordinate is the horizontal direction of S, y the
    vertical direction of S and z the off-plan direction. Therefore, vx and vy
    are the in-plane coordinate and vz the off-plane.
    
    Member variables:
        *triangulation* : triangulation object. 
         Triangulation object, which holds the grid. Therefore, it also holds
         the grid points coordinates (x and y), the triangles and some other
         stuffs. See matplotlib.Tri.Trangulation.Triangulation for more
         details.
         
        *vx*: numpy array of shape (N,)
         
        *vy*: numpy array of shape (N,)
        
        *vz*: numpy array of shape (N,)
        
        *vx_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *vy_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *vz_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *data*: python dict
        
        *data_i*: python dict
        
        *__interType*: python string
        
        *__interKind*: python string
        
    Member functions:
        *__init__*: default constructor
        
        *readFromFoamFile*: constructor from a foamFile generate by OpenFOAM
        
        *readFromVTK*: contructor from a VTK file generate by OpenFOAM
        
        *x*
         Returns the x coordinate of the grid points.
        
        *y*
         Returns the y coordinate of the grid points.
        
        *trangles*
        
        *getInterpolator*
        
        *gradient*
        
        *gradient_i*
        
        
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self, x, y, vx, vy, vz, triangles=None, mask=None):
        '''
        base constructor from a list of x, y and z. list of triangles and mask 
        optional.
        
        Arguments:
            *x*: numpy array of shape (npoints).
             x-coordinates of grid points.
    
            *y*: numpy array of shape (npoints).
             y-coordinates of grid points.
             
            *vx*: numpy array of shape (npoints).
             "in-plan" vector componant. Labeled x TriSurfaceVector
             
            *vy*: numpy array of shape (npoints).
             "in-plan" vector componant. Labeled y in TriSurfaceVector
             
             *vz*: numpy array of shape (npoints).
             "off-plan" vector componant. Labeled z in TriSurfaceVector

            *triangles*: integer array of shape (ntri,3).
             For each triangle, the indices of the three points that make
             up the triangle, ordered in an anticlockwise manner. If no
             triangles is passed, a Delauney triangulation is computed. 
             Default=None.
             
            *mask*: optional boolean array of shape (ntri).
             Which triangles are masked out.
             
            *interpoation*: python string. 
             type of interpolation used. Value: "cubic" or "linear".
             Default="cubic". 

            *kind*: python string.
             Definition of the cubic interpolation type. Value: "geom" or
             "min_E". Default="geom". "min_E" is supposed to be the more 
             accurate, but "geom" is way faster.
        '''

        self.triangulation = tri.Triangulation(x, y, triangles=triangles, mask=mask)

        self.vx=np.asarray(vx)
        self.vy=np.asarray(vy)
        self.vz=np.asarray(vz)
        
    
        self.__interType = None
        self.__interKind  = None
        
        self.vx_i = None
        self.vy_i = None
        self.vz_i = None
        
        self.data = dict()
        self.data_i = dict()
        
 
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         pointsFile,
                         facesFile,
                         viewAnchor,
                         xViewBasis,
                         yViewBasis,
                         srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]],
                         mask=None):
        '''
        Construct from a surface saved  by OpenFOAM in foamFile format.
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
        afftrans = coorTrans.AffineTransfomation(srcBasisSrc,tgtBasisSrc,viewAnchor)
        lintrans = coorTrans.LinearTransformation(srcBasisSrc,tgtBasisSrc)

        # get x and y vector (in ptsTgt)
        ptsSrc = parseFoamFile(pointsFile)
        ptsTgt = np.zeros((ptsSrc.shape[0],ptsSrc.shape[1]))
        for i in range(ptsSrc.shape[0]):
            ptsTgt[i,:] = afftrans.srcToTgt(ptsSrc[i,:])

        #get vectors (in vecsTgt)
        vecsSrc = parseFoamFile(varsFile)
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        for i in range(vecsSrc.shape[0]):
            vecsTgt[i,:] = lintrans.srcToTgt(vecsSrc[i,:])
        
        #get triangles
        triangles = parseFoamFile(facesFile)[:,1:4]
        

        # feed x, y and triangles to the base constructor
        return cls(x=ptsTgt[:,0],
                   y=ptsTgt[:,1],
                   vx=vecsTgt[:,0],
                   vy=vecsTgt[:,1],
                   vz=vecsTgt[:,2],
                   triangles=triangles,
                   mask=mask)
 
       
    @classmethod
    def readFromVTK(cls):
        '''
        Construct from a surface saved by OpenFOAM in VTK format.
        '''
        raise NotImplementedError('The method is not implemented')
  
    # getter and setter #
    #-------------------#

    def x(self):
        '''
        Get x coordinate of the grid points.
        '''
        return self.triangulation.x
        
        
    def y(self):
        '''
        Get y coordinate of the grid points.
        '''
        return self.triangulation.y
        
        
    def triangles(self):
        '''
        Get trangles from the grid.
        '''
        return self.triangulation.triangles
        
    
    # class methods #
    #---------------#
    def getInterpolator(self,interpolation='cubic', kind='geom'):
        self.__interType = interpolation
        self.__interKind = kind
        if self.__interType=='cubic':
            self.vx_i = tri.CubicTriInterpolator(self.triangulation, self.vx, kind=self.__interKind)
            self.vy_i = tri.CubicTriInterpolator(self.triangulation, self.vy, kind=self.__interKind)
            self.vz_i = tri.CubicTriInterpolator(self.triangulation, self.vz, kind=self.__interKind)
        elif self.__interType=='linear':
            self.vx_i = tri.LinearTriInterpolator(self.triangulation, self.vx)
            self.vy_i = tri.LinearTriInterpolator(self.triangulation, self.vy)
            self.vz_i = tri.LinearTriInterpolator(self.triangulation, self.vz)
        else:
            raise ValueError('Interpolation must be "cubic" or "linear".')
        
    
    def gradient(self):
        '''
        Calculate and save the gradient at all point of the grid.
        '''
        self.data['dvxdx'],self.data['dvxdy'] = self.vx_i.gradient(self.x(),self.y())
        self.data['dvxdx'],self.data['dvydy'] = self.vy_i.gradient(self.x(),self.y())
 
       
    def gradient_i(self,x,y):
        '''
        Calculate the gradient at the point pt(x,y) and return it.
        
        Arguments:
            *x*: python float.
             x coordinate of the point pt.
            
            *y*: python float.
             y coordinate of the point pt.
             
        Returns:
            *dvxdx, dvxdy, dvydx, dvydy*: python tuple of four float.
             The gradient at point pt(x,y).
        '''
        dvxdx, dvxdy = self.vx_i.gradient(x,y)
        dvydx, dvydy = self.vy_i.gradient(x,y)
        return dvxdx, dvxdy, dvydx, dvydy 


    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionnary "data"
        '''
        return self.data[key]


    def __setitem__(self, key, item):
        '''
        Setter for key "key" on member dictionnary "data"
        '''
        self.data[key] = item

        
    

            


def parseFoamFile(foamFile):
    '''
    Parse a foamFile generated by the OpenFOAM sample tool or sampling library.
    
    Note:
        * It's a primitiv parse, do not add header in your foamFile!
        * Inline comment are allowed only from line start. c++ comment style.
        * It's REALLY a primitive parser!!!
        
    Arguments:
        * foamFile: [str()] Path of the foamFile

    Returns:
        * output: [numpy.array()] data store in foamFile
    '''
    output = []
    catchFirstNb = False
    istream = open(foamFile, 'r')
    for line in istream: 
        # This regex finds all numbers in a given string.
        # It can find floats and integers writen in normal mode (10000) or
        # with power of 10 (10e3).
        match = re.findall('[-+]?\d*\.?\d+e*[-+]?\d*', line)
        if (line.startswith('//')):
            pass
        if (catchFirstNb==False and len(match)==1):
            catchFirstNb = True
        elif (catchFirstNb==True and len(match)>0):
            matchfloat = list()
            for nb in match:                
                matchfloat.append(float(nb))
            output.append(matchfloat)
        else:
            pass
    istream.close()
    return np.array(output)  