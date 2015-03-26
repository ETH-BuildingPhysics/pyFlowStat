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
from pyFlowStat.TriSurfaceSymmTensor import TriSurfaceSymmTensor
from pyFlowStat.TriSurfaceMesh import TriSurfaceMesh


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
        c.addFieldFromFoamFolder(names=names,time=time)
        return c

    @classmethod
    def createFromHdf5(cls,
                        hdf5Parser,
                        xViewBasis,
                        yViewBasis=None,
                        viewAnchor=(0,0,0),
                        srcBasisSrc=[[1,0,0],[0,1,0],[0,0,1]]):
        '''
        '''
        tsm = TriSurfaceMesh.readFromHdf5(hdf5Parser=hdf5Parser,
                                          xViewBasis=xViewBasis,
                                          yViewBasis=yViewBasis,
                                          viewAnchor=viewAnchor,
                                          srcBasisSrc=srcBasisSrc)
        return cls(triSurfaceMesh=tsm)


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

        
    # class methods #
    #---------------#
    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionary "fields"
        '''
        return self.fields[key]

        
    def addTriSurface(self,triSurface,name):
        '''
        Add an existing TriSurface<type> (for example: TriSurfaceScalar or 
        TriSurfaceVector) to TriSurfaceContainer.
        
        Arguments:
            *triSurface*: TriSurface<type> object.
             The TrisSurface<type> to add.
            
            *name*: string.
             Name of the added field.
        '''
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
            
    def addFieldFromHdf5(self,hdf5Parser,names=[],time=0,projectedField=False):
        '''
        '''
        gTime = str(time)
        if len(names)==0:
            names = hdf5Parser[gTime].keys()
            names.pop(names.index('time'))
        
        for name in names:
            try:
                dataShape = hdf5Parser[gTime][name].value.shape
                if len(dataShape)==1:  #data is a scalar
                    tss = TriSurfaceScalar.readFromHdf5(hdf5Parser=hdf5Parser,
                                                        varName=name,
                                                        triSurfaceMesh=self.triSurfaceMesh,
                                                        time=time,
                                                        projectedField=projectedField)
                    self.fields[name] = tss
                elif len(dataShape)==2 and dataShape[1]==3:  #data is a vector
                    tsv = TriSurfaceVector.readFromHdf5(hdf5Parser=hdf5Parser,
                                                        varName=name,
                                                        triSurfaceMesh=self.triSurfaceMesh,
                                                        time=time,
                                                        projectedField=projectedField)
                    self.fields[name] = tsv
                elif len(dataShape)==2 and dataShape[1]==6:  #data is a symmtensor
                    tsst = TriSurfaceSymmTensor.readFromHdf5(hdf5Parser=hdf5Parser,
                                                             varName=name,
                                                             triSurfaceMesh=self.triSurfaceMesh,
                                                             time=time,
                                                             projectedField=projectedField)
                    self.fields[name] = tsst
                else:
                    raise IOError('variable of name "'+name+'" is not a'
                                  'scalar, not a vector, not a symmTensor.')
            except:
                raise IOError('time "'+gTime+'" and/or name "'+name+'" does not '
                              'exist as key in the HDF5 parser.')
        