'''
TriSurfaceSymmTensor.py
'''

#import re

import numpy as np
import matplotlib.tri as tri
import pyFlowStat.TriSurface as TriSurface
import pyFlowStat.ParserFunctions as ParserFunctions


class TriSurfaceSymmTensor(TriSurface.TriSurface):
    '''
    class TriSurfaceSymmTensor.
    
    Class which describes a symmetric tensor field tii on a 2D (flat)
    triangulated surface S of N points. The triangulation is a TriSurfaceMesh
    object of M triangles.
    
    The coordinate system of S is the following:
        * x = the horizontal coordinate, pointing East (in-plane)
        * y = the vertical coordinate, pointing North (in-plane)
        * z = the off-plane coordinate, pointing torward you     
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,
                 txx,
                 txy,
                 txz,
                 tyy,
                 tyz,
                 tzz,
                 time,
                 triSurfaceMesh,
                 projectedField=False,
                 interpolation=None,
                 kind=None):
        '''
        base constructor.
        
        Arguments:          
            *txx*: numpy array of shape (npoints).
             txx values.
             
            *txy*: numpy array of shape (npoints).
             txy values.
            
            *txz*: numpy array of shape (npoints).
             txz values.
             
            *tyy*: numpy array of shape (npoints).
             tyy values.
             
            *tyz*: numpy array of shape (npoints).
             tyz values.
             
            *tzz*: numpy array of shape (npoints).
             tzz values.
         
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=False)
             Unused for the scalar field, but might be used for the fields
             added with the methods "addField", "addFieldFromFoamFile" and
             "addFieldFromVtk"
             
            *interpoation*: python string. 
             type of interpolation used. Value: "cubic" or "linear".
             Default="cubic". 

            *kind*: python string.
             Definition of the cubic interpolation type. Value: "geom" or
             "min_E". Default="geom". "min_E" is supposed to be the more 
             accurate, but "geom" is way faster.
        '''
        super(TriSurfaceSymmTensor,self).__init__(time=time,
                                                  triSurfaceMesh=triSurfaceMesh,
                                                  projectedField=projectedField,
                                                  interpolation=interpolation,
                                                  kind=kind)

        self.txx = np.asarray(txx)
        self.txy = np.asarray(txy)
        self.txz = np.asarray(txz)
        self.tyy = np.asarray(tyy)
        self.tyz = np.asarray(tyz)
        self.tzz = np.asarray(tzz)
        
        self.txx_i = None
        self.txy_i = None 
        self.txz_i = None 
        self.tyy_i = None 
        self.tyz_i = None 
        self.tzz_i = None 

        
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
             Path to the file holding the scalar field.
             
            *time*: python float
             timestep of the surface. If this information does not matter,
             use 0.
             
            *triSurfaceMesh* :  TriSurfaceMesh object.
             TriSurfaceMesh object, which holds the mesh information.
             
            *projectedField* python bool (default=False)
             Unused for the scalar field, but might be used for the fields
             added with the methods "addField", "addFieldFromFoamFile" and
             "addFieldFromVtk"
        '''

        #get scalars
        stensSrc = ParserFunctions.parseFoamFile_sampledSurface(varsFile)
        stensTgt = np.zeros((stensSrc.shape[0],stensSrc.shape[1]))
        if projectedField==True:
            for i in range(stensSrc.shape[0]):
                stenAsMat = TriSurface.mat(stensSrc[i,:])
                stensTgt[i,:] = triSurfaceMesh.linTrans.srcToTgt(stenAsMat)
        else:
            stensTgt = stensSrc

        # update class member variables
        return cls(txx=stensTgt[:,0],
                   txy=stensTgt[:,1],
                   txz=stensTgt[:,2],
                   tyy=stensTgt[:,3],
                   tyz=stensTgt[:,4],
                   tzz=stensTgt[:,5],
                   time=time,
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None)


    @classmethod 
    def readFromHdf5(cls,
                     hdf5Parser,
                     varName,
                     triSurfaceMesh,
                     key,
                     projectedField=False):
        '''
        '''
        time = hdf5Parser[key]['time'].value
        stensSrc = hdf5Parser[key][varName].value
        stensTgt = np.zeros((stensSrc.shape[0],stensSrc.shape[1]))
        if projectedField==True:
            for i in range(stensSrc.shape[0]):
                stenAsMat = TriSurface.mat(stensSrc[i,:])
                stensTgt[i,:] = triSurfaceMesh.linTrans.srcToTgt(stenAsMat)
        else:
            stensTgt = stensSrc
        
        # update class member variables
        return cls(txx=stensTgt[:,0],
                   txy=stensTgt[:,1],
                   txz=stensTgt[:,2],
                   tyy=stensTgt[:,3],
                   tyz=stensTgt[:,4],
                   tzz=stensTgt[:,5],
                   time=float(time),
                   triSurfaceMesh=triSurfaceMesh,
                   projectedField=projectedField,
                   interpolation=None,
                   kind=None) 
                   
                   
#    @classmethod   
#    def readFromVTK(cls,
#                    vtkFile,
#                    triSurfaceMesh,
#                    time,
#                    projectedField=True):
#        '''
#        Construct from a surface saved by OpenFOAM in VTK format.
#        
#        Arguments:
#            *varsFile*: python string.
#             Path to the vtk-file.
#             
#            *time*: python float
#             timestep of the surface. If this information does not matter,
#             use 0.
#             
#            *triSurfaceMesh* :  TriSurfaceMesh object.
#             TriSurfaceMesh object, which holds the mesh information.
#             
#            *projectedField* python bool (default=True)
#             Defines if the data fields has to be projected in the basis of the
#             surface. 
#        '''     
#        # read VTK file
#        ptsSrc, triangles, slrsTgt = TriSurfaceFunctions.parseVTK_ugly_sampledSurface(vtkFile)
#
#
#        # update class member variables
#        return cls(s=slrsTgt,
#                   time=time,
#                   triSurfaceMesh=triSurfaceMesh,
#                   projectedField=projectedField,
#                   interpolation=None,
#                   kind=None)
#  
#    # getters #
#    #---------#      
#
#
#    # setters #
#    #---------#
#
#
#    # class methods #
#    #---------------#

    def __call__(self,dim):
        return self.component(dim)
      
      
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
            
            
    def interpolate(self,x,y,dim):
        try:
            if dim==0:
                return self.txx_i(x,y)
            if dim==1:
                return self.txy_i(x,y)
            if dim==2:
                return self.txz_i(x,y)
            if dim==3:
                return self.tyy_i(x,y)
            if dim==4:
                return self.tyz_i(x,y)
            if dim==5:
                return self.tzz_i(x,y)
        except:
            raise ValueError('this method needs interpolators. Please run',
                                 'method "addInterpolator" first.')

#    def rawVars(self):
#        '''
#        Return the scalar field defined in the source coordinate system.
#        
#        Returns:
#            *rawData*: numpy array of shape (N,)
#        '''
#        return self.s
#        
#
#    def surfaceVars(self):
#        '''
#        Return the scalar field as saved in the TriSurfaceScalar object.
#        
#        Returns:
#            *surfaceVars*: numpy array of shape (N,)
#        '''
#        return self.s
#
#
#    def gradientxy(self,x,y):
#        '''
#        Calculate the gradient at the point pt(x,y) and return it. The method
#        works also for a list of N points. In such case, x and y are arrays.
#        
#        Arguments:
#            *x*: python float or numpy array (shape=N).
#             x coordinate of the point pt.
#            
#            *y*: python float or numpy array (shape=N).
#             y coordinate of the point pt.
#             
#        Returns:
#            *dsdx, dsdy*: python tuple of two float or numpy array.
#             The gradient at point pt(x,y).
#        '''
#        if self.s_i!=None: 
#            dsdx, dsdy = self.s_i.gradient(x,y)
#            return dsdx, dsdy
#        else:
#            raise ValueError('this method needs interpolators. Please run',
#                             'method "addInterpolator" first.')
#
#  
#    # class methods - adders #
#    #------------------------#
#        
    def addInterpolator(self,interpolation='cubic', kind='geom'):
        '''
        Add interpolator Object to the vector field.
        '''
        self.interType = interpolation
        self.interKind = kind
        if self.interType=='cubic':
            self.txx_i = tri.CubicTriInterpolator(self.triangulation, self.txx, kind=self.interKind)
            self.txy_i = tri.CubicTriInterpolator(self.triangulation, self.txy, kind=self.interKind)
            self.txz_i = tri.CubicTriInterpolator(self.triangulation, self.txz, kind=self.interKind)
            self.tyy_i = tri.CubicTriInterpolator(self.triangulation, self.tyy, kind=self.interKind)
            self.tyz_i = tri.CubicTriInterpolator(self.triangulation, self.tyz, kind=self.interKind)
            self.tzz_i = tri.CubicTriInterpolator(self.triangulation, self.tzz, kind=self.interKind)
        elif self.interType=='linear':
            self.txx_i = tri.LinearTriInterpolator(self.triangulation, self.txx)
            self.txy_i = tri.LinearTriInterpolator(self.triangulation, self.txy)
            self.txz_i = tri.LinearTriInterpolator(self.triangulation, self.txz)
            self.tyy_i = tri.LinearTriInterpolator(self.triangulation, self.tyy)
            self.tyz_i = tri.LinearTriInterpolator(self.triangulation, self.tyz)
            self.tzz_i = tri.LinearTriInterpolator(self.triangulation, self.tzz)
        else:
            raise ValueError('Interpolation must be "cubic" or "linear".')
#            
#
#    def addGradient(self):
#        '''
#        Calculate and save the gradient at all point of the grid.
#        '''   
#        dsdx, dsdy = self.gradientxy(self.x,self.y)
#        self.data['dsdx'] = dsdx
#        self.data['dsdy'] = dsdy