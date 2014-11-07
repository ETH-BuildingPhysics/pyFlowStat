'''
TriSurfaceVector.py
'''

#import re

import numpy as np
import matplotlib.tri as tri

#import CoordinateTransformation as coorTrans
#import TriSurfaceMesh as TriSurfaceMesh
import TriSurfaceFunctions


class TriSurfaceVector(object):
    '''
    class TriSurfaceVector.
    
    Class which describes a 3D vector field V=(vx,vy,vz) on a 2D (flat)
    triangulated surface S of N points. The triangulation is a TriSurfaceMesh
    object of M triangles.
    
    The coordinate system of S is the following:
        * x = the horizontal coordinate, pointing East (in-plane)
        * y = the vertical coordinate, pointing North (in-plane)
        * z = the off-plane coordinate, pointing torward you   
    If "projectedField" is "True" in the constructor, vx and vy are the
    in-plane vector coordinate and vz the off-plane one. If "False", no
    transformation is made on the input vx, vy and vz.
    
    To construct a TriSurfaceVector object, use better the @classmethod
    "readFromFoamFile" or "readFromVTK".

    Member variables:
        *triSurfaceMesh* :  TriSurfaceMesh object.
         
        *vx*: numpy array of shape (N,)
         
        *vy*: numpy array of shape (N,)
        
        *vz*: numpy array of shape (N,)
        
        *vx_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *vy_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *vz_i*: CubicTriInterpolator or LinearTriInterpolation object
        
        *data*: python dict
        
        *data_i*: python dict
  
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
             
            *vx*: numpy array of shape (npoints).
             "in-plan" vector componant. Labeled x TriSurfaceVector
             
            *vy*: numpy array of shape (npoints).
             "in-plan" vector componant. Labeled y in TriSurfaceVector
             
            *vz*: numpy array of shape (npoints).
             "off-plan" vector componant. Labeled z in TriSurfaceVector
             
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=True)
             Defines if the data fields has to be projected in the basis of the
             surface. 
             
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
        
        Arguments:
            *varsFile*: python string.
             Path to the file holding the vector field.
             
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=True)
             Defines if the data fields has to be projected in the basis of the
             surface. 
        '''

        #get vectors (in vecsTgt)
        vecsSrc = TriSurfaceFunctions.parseFoamFile_sampledSurface(varsFile)
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = triSurfaceMesh.linTrans.srcToTgt(vecsSrc[i,:])
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
        
        Arguments:
            *varsFile*: python string.
             Path to the vtk-file.
             
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=True)
             Defines if the data fields has to be projected in the basis of the
             surface. 
        '''     
        # read VTK file
        ptsSrc, triangles, vecsSrc = TriSurfaceFunctions.parseVTK_ugly_sampledSurface(vtkFile)
            
        # transform the data
        vecsTgt = np.zeros((vecsSrc.shape[0],vecsSrc.shape[1]))
        if projectedField==True:
            for i in range(vecsSrc.shape[0]):
                vecsTgt[i,:] = triSurfaceMesh.linTrans.srcToTgt(vecsSrc[i,:])
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


    # setters #
    #---------#


    # class methods #
    #---------------#
    def rawPoints(self):
        '''
        Return the grid points in the source coordinate system.
        
        Returns:
            *rawPoints*: numpy array of shape (N,3)
        '''
        return self.triSurfaceMesh.rawPoints()
            
     
    def rawVars(self):
        '''
        Return the vector field defined the source coordinate system.
        
        Returns:
            *rawData*: numpy array of shape (N,3)
        '''
        surfaceData = self.surfaceVars()
        rawData = np.zeros((surfaceData.shape[0],surfaceData.shape[1]))
        if self.__projectedField==True:
            for i in range(surfaceData.shape[0]):
                rawData[i,:] = self.linTrans.tgtToSrc(surfaceData[i,:])
        else:
            rawData = surfaceData
        return rawData
        

    def surfaceVars(self):
        '''
        Return the vector field as saved in the TriSurfaceVector object.
        
        Returns:
            *surfaceVars*: numpy array of shape (N,3)
        '''
        return np.vstack((self.vx,self.vy,self.vz)).T


    def Umag(self):
        return np.sqrt(np.square(self.vx)+np.square(self.vy)+np.square(self.vz))
      
     
    def VortZ(self):
        if self.vx_i!=None:  
            VortZ = None
            if (self.data.has_key('dvydx')==True and self.data.has_key('dvydx')==True):
                VortZ = self.data['dvydx']-self.data['dvydx']
            else:
                dvydx,dvydy = self.vy_i.gradient(self.x,self.y)
                dvxdx,dvxdy = self.vx_i.gradient(self.x,self.y)
                VortZ = dvydx-dvydx
            
            return VortZ
            
        else:
            raise ValueError('this method needs interpolators. Please run',
                             'method "addInterpolator" first.')
    
    def Q(self):
        '''
        Evaluate the Q criteria on the vector field v. Q make sense only if v
        is a velocity field. Be careful: The full definition of Q needs the
        velocity derivatives in the z direction too (off-plane derivatives).
        As those derivative are not accessible in TriSurfaceVector, Q() returns
        a sort of troncated Q...
        
        Returns:
            *Q*: numpy array of shape (N,).
        '''
        if self.vx_i!=None:
            Q = None
            if (self.data.has_key('dvydx')==True and self.data.has_key('dvydx')==True):
                Q = 0.5*(-2.0*self['dvxdy']*self['dvydx']-self['dvxdx']**2-self['dvydy']**2)
            else:
                dvydx,dvydy = self.vy_i.gradient(self.x,self.y)
                dvxdx,dvxdy = self.vx_i.gradient(self.x,self.y)
                Q = 0.5*(-2.0*dvxdy*dvydx-dvxdx**2-dvydy**2)
            
            return Q
        else:
            raise ValueError('this method needs interpolators. Please run',
                             'method "addInterpolator" first.')


    def gradientxy(self,x,y):
        '''
        Calculate the gradient at the point pt(x,y) and return it.
        
        Arguments:
            *x*: python float or numpy array.
             x coordinate of the point pt.
            
            *y*: python float or numpy array.
             y coordinate of the point pt.
             
        Returns:
            *dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy *: python tuple of six float or numpy array.
             The gradient at point pt(x,y).
        '''
        if self.vx_i!=None: 
            dvxdx, dvxdy = self.vx_i.gradient(x,y)
            dvydx, dvydy = self.vy_i.gradient(x,y)
            dvzdx, dvzdy = self.vz_i.gradient(x,y)
            return dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy
        else:
            raise ValueError('this method needs interpolators. Please run',
                             'method "addInterpolator" first.')
    
    
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
        
    # class methods - adders #
    #------------------------#
        
    def addInterpolator(self,interpolation='cubic', kind='geom'):
        '''
        Add interpolator Object to the vector field.
        '''
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
        fieldTgt = np.zeros(fieldShape)

        if (self.__projectedField==True and len(fieldShape)>1):
            for i in range(fieldShape[0]):
                fieldTgt[i,:] = self.linTrans.srcToTgt(fieldSrc[i,:])
        else:
            fieldTgt = fieldSrc
        self.data[fieldname] = fieldTgt
            
        
    def addFieldFromFoamFile(self,fieldFile,fieldname):
        '''
        Add a field F (shape d) stored in a foamFile to the current
        TriSurfaceVector object. See docstring from self.addField() for more
        information.

        '''
        #get field
        fieldSrc = TriSurfaceFunctions.parseFoamFile_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
        
        
    def addFieldFromVTK(self,fieldFile,fieldname):
        '''
        Add a field F (shape d) stored in a VTK file to the current
        TriSurfaceVector object. See docstring from self.addField() for more
        information.

        '''
        #get field
        points, polygon, fieldSrc = TriSurfaceFunctions.parseVTK_ugly_sampledSurface(fieldFile)
        self.addField(fieldSrc,fieldname)
    
    
    def addGradient(self):
        '''
        Calculate and save the gradient at all point of the grid. As expected,
        the dvidz does not exist.
        '''   
        dvxdx, dvxdy, dvydx, dvydy, dvzdx, dvzdy = self.gradientxy(self.x,self.y)
        self.data['dvxdx'] = dvxdx
        self.data['dvxdy'] = dvxdy
        
        self.data['dvydx'] = dvydx
        self.data['dvydy'] = dvydy
        
        self.data['dvzdx'] = dvzdx
        self.data['dvzdy'] = dvzdy
        

    def addUmag(self):
        '''
        Add velocity magnitude in data  (TriSurfaceVector['Umag']). This method
        makes sense only if vx, vy and vz are velocity componants. 
        '''
        self.data['Umag'] = self.Umag()
        
        
    def addU(self):
        '''
        Add velocity in data  (TriSurfaceVector['U']). This method
        makes sense only if vx, vy and vz are velocity componants.
        '''
        self.data['U'] = self.surfaceVars()
        
        
    def addVortZ(self):
        '''
        Add z componant of the vorticity in data  (TriSurfaceVector['VortZ']). 
        This method makes sense only if vx, vy and vz are velocity componants.
        '''
        self.data['VortZ'] = self.VortZ()
        
    def addQ(self):
        '''
        Add the Q criterion in data  (TriSurfaceVector['Q']). 
        This method makes sense only if vx, vy and vz are velocity componants.
        '''
        self.data['Q'] = self.Q()
        
    
            
    
    
 
