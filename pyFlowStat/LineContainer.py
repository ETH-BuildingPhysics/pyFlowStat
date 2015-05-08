'''
LineContainer.py
'''
import numpy as np
import os
import glob

from pyFlowStat.Line import LineDict
from pyFlowStat.LineScalar import LineScalar
from pyFlowStat.LineVector import LineVector
from pyFlowStat.LineSymmTensor import LineSymmTensor



class LineContainer(object):
    '''
 
    '''
    
    # constructors #
    #--------------#
    
    def __init__(self):
        '''
        base constructor.
        '''
        self.lines=LineDict()
        self.data=dict()
        self.time=0.0
        
    @classmethod
    def createFromFoamFolder(cls,pathname,time=0.0,names=[]):
        c=cls()
        c.data['name']=os.path.basename(pathname)
        c.data['pathname']=pathname
        c.time=time
        #TODO: find names here?
        c.addLinesFromFoamFolder(names=names,time=time)
        return c
        
    # class methods #
    #---------------#
    def __getitem__(self, key):
        '''
        Getter for key "key" on member dictionary "fields"
        '''
        return self.lines[key]
        
    def __setitem__(self, key, item):
        '''
        Add an existing TriSurface<type> (for example: TriSurfaceScalar or 
        TriSurfaceVector) to TriSurfaceContainer.
        '''
        self.addLine(item,key)

        
    def addLine(self,line,name):
        '''
        Add an existing Line<type> (for example: LineScalar or 
        LineVector) to LineContainer.
        
        Arguments:
            *line*: Line<type> object.
             The Line<type> to add.
            
            *name*: string.
             Name of the added field.
        '''

        self.lines[name]=line
        

    def addFoamScalarLine(self,linepath,name):
        '''
        '''
        #TODO: change
        tLine=LineScalar.readFromFoamFile(linepath)
        self.lines[name] = tLine
        
    
    def addFoamVectorLine(self,linepath,name):
        '''
        '''
        #TODO: change
        tLine=LineVector.readFromFoamFile(linepath)
        self.lines[name] = tLine

                               
    def addFoamSymmTensorLine(self,linepath,name):
        '''
        '''
        #TODO: change
        tLine=LineSymmTensor.readFromFoamFile(linepath)           
        self.lines[name] = tLine
        
    def addLinesFromFoamFolder(self,linepath='',names=[]):
        '''
        '''
        if linepath=='':
            linepath=self.data['pathname']

        if os.path.exists(linepath):
            #TODO: change
            '''
            for line in scalarLines:
                key=os.path.basename(field)
                if key in names:
                    self.addFoamScalarField(field,key)
            for field in vectorLines:
                key=os.path.basename(field)
                if key in names:
                    self.addFoamVectorField(field,key)
            for field in symmTensorLines:
                key=os.path.basename(field)
                if key in names:
                    self.addFoamSymmTensorField(field,key)
            '''
            pass
        else:
            raise IOError("Folder does not exist")