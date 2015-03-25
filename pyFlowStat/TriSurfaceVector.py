'''
TriSurfaceVector.py
'''

#import re

import numpy as np
import matplotlib.tri as tri

import pyFlowStat.TriSurface as TriSurface
import pyFlowStat.ParserFunctions as ParserFunctions


class TriSurfaceVector(TriSurface.TriSurface):
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
                 projectedField=False,
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
        super(TriSurfaceVector,self).__init__(time=time,
                                        triSurfaceMesh=triSurfaceMesh,
                                        projectedField=projectedField,
                                        interpolation=interpolation,
                                        kind=kind)

        self.vx=np.asarray(vx)
        self.vy=np.asarray(vy)
        self.vz=np.asarray(vz)
        
        self.vx_i = None
        self.vy_i = None
        self.vz_i = None

        
    @classmethod
    def readFromFoamFile(cls,
                         varsFile,
                         triSurfaceMesh,
                         time,
                         projectedField=False):
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
             
            *projectedField* python bool (default=False)
             Defines if the data fields has to be projected in the basis of the
             surface. 
        '''

        #get vectors (in vecsTgt)
        vecsSrc = ParserFunctions.parseFoamFile_sampledSurface(varsFile)
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
                    projectedField=False):
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
             
            *projectedField* python bool (default=False)
             Defines if the data fields has to be projected in the basis of the
             surface. 
        '''     
        # read VTK file
        ptsSrc, triangles, vecsSrc = ParserFunctions.parseVTK_ugly_sampledSurface(vtkFile)
            
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


    # setters #
    #---------#


    # class methods #
    #---------------#
    def component(self,dim):
        if dim==0:
            return self.vx
        if dim==1:
            return self.vy
        if dim==2:
            return self.vz
            
    def rawVars(self):
        '''
        Return the vector field defined the source coordinate system.
        
        Returns:
            *rawData*: numpy array of shape (N,3)
        '''
        surfaceData = self.surfaceVars()
        rawData = np.zeros((surfaceData.shape[0],surfaceData.shape[1]))
        if self.projectedField==True:
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

  
    # class methods - adders #
    #------------------------#
        
    def addInterpolator(self,interpolation='cubic', kind='geom'):
        '''
        Add interpolator Object to the vector field.
        '''
        self.interType = interpolation
        self.interKind = kind
        if self.interType=='cubic':
            self.vx_i = tri.CubicTriInterpolator(self.triangulation, self.vx, kind=self.interKind)
            self.vy_i = tri.CubicTriInterpolator(self.triangulation, self.vy, kind=self.interKind)
            self.vz_i = tri.CubicTriInterpolator(self.triangulation, self.vz, kind=self.interKind)
        elif self.interType=='linear':
            self.vx_i = tri.LinearTriInterpolator(self.triangulation, self.vx)
            self.vy_i = tri.LinearTriInterpolator(self.triangulation, self.vy)
            self.vz_i = tri.LinearTriInterpolator(self.triangulation, self.vz)
        else:
            raise ValueError('Interpolation must be "cubic" or "linear".')
            

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

