'''
TriSurfaceVector.py
'''

import re

import numpy as np
import matplotlib.tri as tri

import CoordinateTransformation as coorTrans
import TriSurfaceMesh as TriSurfaceMesh
import TriSurfaceFunctions


class TriSurfaceVector(object):
    '''
    class TriSurfaceVector.
    
    Class which describes a 3D vector field V=(vx,vy,vz) on a 2D (flat)
    triangulated surface S of N points. The triangulation is a grid of M 
    triangles. 
    The coordinate system of S is the following:
        * x = the horizontal coordinate, pointing East (in-plane)
        * y = the vertical coordinate, pointing North (in-plane)
        * z = the off-plane coordinate, pointing torward you   
    Therefore, vx and vy are the in-plane vector coordinate and vz the
    off-plane one.
    
    To construct a TriSurfaceVector object, use better the @classmethod
    "readFromFoamFile" or "readFromVTK".

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
        
        *__affTrans*:
        
        *__linTrans*:
        
        *__projectedField*: python bool
        
    Member functions:
        *__init__*: base constructor.
        
        *readFromFoamFile*: constructor from a foamFile generate by OpenFOAM.
        
        *readFromVTK*: contructor from a VTK file generate by OpenFOAM.
        
        *x*: Returns the x coordinate of the grid points.
        
        *y*: Returns the y coordinate of the grid points.
        
        *trangles*
        
        *getInterpolator*
        
        *addFields*

        *addFieldFromFoamFile*        
        
        *gradient*
        
        *gradient_i*
        
        
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,
                 vx,
                 vy,
                 vz,
                 time,
                 triSurfaceMesh,
                 projectedField=True,
                 interpolation=None,
                 kind=None):
        '''
        base constructor.
        
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
             
            *projectedField* python bool (default=True)
             Defines if the data fields has to be projected in the basis of the
             surface. Example where such option really matters:
             Let us consider a 3D laminar flow with the follow velocity field: 
             U=(Ux,0,0). Let us consider the surface S1 with the normal n1=(0,ny,0)
             and the surface S2 with the normal n2=(nx,ny,0). On S1, Ux will be
             in the surface plan, therefore S1.gradient() makes sense. On the 
             other hand, Ux in NOT in the plane of S2.
             
            *interpoation*: python string. 
             type of interpolation used. Value: "cubic" or "linear".
             Default="cubic". 

            *kind*: python string.
             Definition of the cubic interpolation type. Value: "geom" or
             "min_E". Default="geom". "min_E" is supposed to be the more 
             accurate, but "geom" is way faster.
        '''
        self.triSurfaceMesh = triSurfaceMesh

        self.vx=np.asarray(vx)
        self.vy=np.asarray(vy)
        self.vz=np.asarray(vz)
        
        self.vx_i = None
        self.vy_i = None
        self.vz_i = None
        
        self.data = dict()
        self.data_i = dict()
        
        self.time = float(time)

        # "private" member variable. Don't play with them if you are not sure...        
        self.__interType = interpolation
        self.__interKind  = kind
        self.__projectedField = projectedField
        
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         triSurfaceMesh,
                         time,
                         projectedField=True):
        '''
        Construct from a surface saved  by OpenFOAM in foamFile format.
        '''

        #get vectors (in vecsTgt)
        vecsSrc = TriSurfaceFunctions.parseFoamFile_sampledSurface(varsFile)
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = lintrans.srcToTgt(vecsSrc[i,:])
        else:
            vecsTgt = vecsSrc

        # update class member variables
        return cls(vx=vecsTgt[:,0],
                   vy=vecsTgt[:,1],
                   vz=vecsTgt[:,2],
                   time=time,
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None)
 
    @classmethod   
    def readFromVTK(cls,
                    vtkFile,
                    triSurfaceMesh,
                    time,
                    projectedField=True):
        '''
        Construct from a surface saved by OpenFOAM in VTK format.
        '''     
        # read VTK file
        ptsSrc, triangles, vecsSrc = TriSurfaceFunctions.parseVTK_ugly_sampledSurface(vtkFile)
            
        # transform the data
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = lintrans.srcToTgt(vecsSrc[i,:])
        else:
            vecsTgt = vecsSrc

        # update class member variables
        return cls(vx=vecsTgt[:,0],
                   vy=vecsTgt[:,1],
                   vz=vecsTgt[:,2],
                   time=time,
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None)
  
    # getters #
    #---------#
    @property
    def x(self):
        return self.triSurfaceMesh.x
        
    @property    
    def y(self):
        return self.triSurfaceMesh.y
        
    @property
    def triangulation(self):
        return self.triSurfaceMesh.triangulation
        
    @property
    def triangles(self):
        return self.triSurfaceMesh.triangles
        
    @property
    def affTrans(self):
        return self.triSurfaceMesh.affTrans
    
    @property
    def linTrans(self):
        return self.triSurfaceMesh.linTrans

    @property
    def rawPoints(self):
        return self.triSurfaceMesh.rawPoints
            
    @property        
    def rawVars(self):
        surfaceData = self.surfaceVars
        rawData = np.zeros((surfaceData.shape[0],surfaceData.shape[1]))
        if self.__projectedField==True:
            for i in range(surfaceData.shape[0]):
                rawData[i,:] = self.linTrans.tgtToSrc(surfaceData[i,:])
        else:
            rawData = surfaceData
        return rawData
        
    @property
    def surfaceVars(self):
        return np.vstack((self.vx,self.vy,self.vz)).T
        


    # setters #
    #---------#


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
    
    
    def addField(self,field,fieldname):
        '''
        Add a field F of dimension d (e.g: d=3 for a vector filed) to the
        current TriSurfaceVector object TSV. The grid of F (N points) must be
        identical, in term of number of points and their location, as the grid
        of TSV. F will be stored in TSV.data['fieldName'] or TSV['fieldName'].
        
        Arguments:
            *field*: numpy array of shape (N,d).
            
            *fieldName*: python string.
        '''
        fieldSrc = field
        fieldShape = fieldSrc.shape
        fieldTgt = np.zeros((fieldShape[0],fieldShape[1]))

        if (self.__projectedField==True and fieldShape[1]>1):
            for i in range(fieldShape[0]):
                fieldTgt[i,:] = self.linTrans.srcToTgt(fieldSrc[i,:])
        else:
            fieldTgt = fieldSrc
        self.data[fieldname] = fieldTgt
            
        
    def addFieldFromFoamFile(self,fieldFile,fieldname):
        '''
        Add a field F of dimension d (e.g: d=3 for a vector filed) to the
        current TriSurfaceVector object TSV. The grid of F (N points) must be
        identical, in term of number of points and their location, as the grid
        of TSV. F will be stored in TSV.data['fieldName'] or TSV['fieldName'].
        
        Arguments:
            *field*: numpy array of shape (N,d).
            
            *fieldName*: python string.
        '''
        #get field
        fieldSrc = parseFoamFile_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
        
    
    def gradient(self):
        '''
        Calculate and save the gradient at all point of the grid.
        '''
        self.data['dvxdx'],self.data['dvxdy'] = self.vx_i.gradient(self.x(),self.y())
        self.data['dvydx'],self.data['dvydy'] = self.vy_i.gradient(self.x(),self.y())
        self.data['dvzdx'],self.data['dvzdy'] = self.vz_i.gradient(self.x(),self.y())
 
       
    def gradient_i(self,x,y):
        '''
        Calculate the gradient at the point pt(x,y) and return it.
        
        Arguments:
            *x*: python float or numpy array.
             x coordinate of the point pt.
            
            *y*: python float or numpy array.
             y coordinate of the point pt.
             
        Returns:
            *dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy *: python tuple of four float.
             The gradient at point pt(x,y).
        '''
        dvxdx, dvxdy = self.vx_i.gradient(x,y)
        dvydx, dvydy = self.vy_i.gradient(x,y)
        dvzdx, dvzdy = self.vz_i.gradient(x,y)
        return dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy


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
