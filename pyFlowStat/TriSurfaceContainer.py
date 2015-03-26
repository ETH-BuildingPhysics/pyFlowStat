'''
TriSurfaceContainer.py
'''
import numpy as np
import os
import glob
#import matplotlib.tri as tri
import matplotlib.tri as tri

from pyFlowStat.TriSurfaceVector import TriSurfaceVector
from pyFlowStat.TriSurfaceScalar import TriSurfaceScalar
from pyFlowStat.TriSurfaceVector import TriSurfaceVector
from pyFlowStat.TriSurfaceSymmTensor import TriSurfaceSymmTensor
from pyFlowStat.TriSurfaceMesh import TriSurfaceMesh

import pyFlowStat.ParserFunctions as ParserFunctions


class TriSurfaceContainer(object):
    '''
 
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self,triSurfaceMesh):
        '''
        base constructor.
        '''
        self.triSurfaceMesh = triSurfaceMesh
        self.fields=dict()
        self.data=dict()
        
    @classmethod
    def createFromTriSurface(cls,triSurface,name):
        c=cls(triSurfaceMesh=triSurface.triSurfaceMesh)
        c.fields[name]=triSurface
        return c
    
    @classmethod
    def createFromFoamFolder(cls,pathname,xViewBasis,yViewBasis=None,viewAnchor=(0,0,0),time=0.0,names=[]):
        '''
        This function creates a TriSurfaceContainer for an OpenFOAM sampled surface folder (e.g./postProcessing/surfaces/0/planeXY/) and loads all the fields.
        '''
        pointsFile=os.path.join(pathname,'points')
        facesFile=os.path.join(pathname,'faces')
        tsm=TriSurfaceMesh.readFromFoamFile(pointsFile=pointsFile,
                                           facesFile=facesFile,
                                           viewAnchor=viewAnchor,
                                           xViewBasis=xViewBasis,
                                           yViewBasis=yViewBasis)
        
        c=cls(triSurfaceMesh=tsm)
        c.data['name']=os.path.basename(pathname)
        c.data['pathname']=pathname
        c.addFieldFromFoamFile(names=names,time=time)
        return c
        
    @property
    def triangulation(self):
        return self.triSurfaceMesh.triangulation
        
    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionary "fields"
        '''
        return self.fields[key]
        
    def addTriSurface(self,triSurface,name):
        if triSurface.triSurfaceMesh is not self.triSurfaceMesh:
            raise ValueError("triSurfaceMesh is not identical")
        self.fields[name]=triSurface
        
    def addFoamScalarField(self,fieldpath,name,time,projectedField=False):
        '''
        '''
        tSurf=TriSurfaceScalar.readFromFoamFile(fieldpath,self.triSurfaceMesh,time=time,projectedField=projectedField)
        self.fields[name] = tSurf
    
    def addFoamVectorField(self,fieldpath,name,time,projectedField=False):
        '''
        '''
        tSurf=TriSurfaceVector.readFromFoamFile(fieldpath,self.triSurfaceMesh,time=time,projectedField=projectedField)
        self.fields[name] = tSurf
                               
    def addFoamSymmTensorField(self,fieldpath,name,time,projectedField=False):
        '''
        '''
        tSurf=TriSurfaceSymmTensor.readFromFoamFile(fieldpath,self.triSurfaceMesh,time=time,projectedField=projectedField)           
        self.fields[name] = tSurf
        
    def addFieldFromFoamFolder(self,surfacepath='',names=[],time=0.0):
        '''
        '''
        if surfacepath=='':
            surfacepath=self.data['pathname']
        
        if os.path.exists(surfacepath):
            scalarFields=glob.glob(os.path.join(surfacepath,'scalarField\*'))
            vectorFields=glob.glob(os.path.join(surfacepath,'vectorField\*'))
            symmTensorFields=glob.glob(os.path.join(surfacepath,'symmTensorField\*'))
            if len(names)<1:
                for name in scalarFields:
                    names.append(os.path.basename(name))
                for name in vectorFields:
                    names.append(os.path.basename(name))
                for name in symmTensorFields:
                    names.append(os.path.basename(name))
            
            for field in scalarFields:
                key=os.path.basename(field)
                if key in names:
                    self.addFoamScalarField(field,key,time)
            for field in vectorFields:
                key=os.path.basename(field)
                if key in names:
                    self.addFoamVectorField(field,key,time)
            for field in symmTensorFields:
                key=os.path.basename(field)
                if key in names:
                    self.addFoamSymmTensorField(field,key,time)
        else:
            raise IOError("Folder does not exist")
        