import numpy as np
import modred

class POD(object):
    def __init__(self,scalarFieldList):
        '''
        Arguments:
            *scalarFieldList*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example one component of a SurfaceList).
        '''
        self.surfaces=scalarFieldList
        self.surfShape=self.surfaces.shape
        self.result=dict()
        pass
    
    def decompose(self,nMode,PODmethod='snap'):
        '''
        Compute the Porper Orthogonal Decomposition (POD) from a list of N snapshots. 
        The snapshots are PIV surfaces of size (surfX,surfY) of field F. F can be any
        scalar field (Ux, ux, T, vorticity, R11,...)
        
        Arguments:

            
            *nMode*: python integer.
             Number of modes of the POD.
            
            *PODmethod*: python string. Default='snap'
             Type of POD algorithm used. For the snapshot method, use PODmethod='snap' and
             for the direct method, use PODmewthod='direct'. The default value is 'snap'.
             
             How to choose between 'snap' and 'direct': if surfX*surfY > N**2, use the snap
             method, it should increase the coputational speed for only a little loss in 
             precision.
             
        Returns:
            None
        
        '''

        nSnap = self.surfShape[0]
        
        # reshape the input
        surfacesPOD = np.nan_to_num(self.surfaces.reshape(self.surfShape[0],self.surfShape[1]*self.surfShape[2]).T)
        
        if PODmethod=='snap':
            modesPOD, eigVals = modred.compute_POD_matrices_snaps_method(surfacesPOD, range(nMode))
        elif PODmethod=='direct':
            modesPOD, eigVals = modred.compute_POD_matrices_direct_method(surfacesPOD, range(nMode))
        else:
            print('error: argument '+str(PODmethod)+' is not valid. Use \'snap\' or \'direct\'')
        
        # convert modesPOD in an numpy array
        modesPOD = np.asarray(modesPOD)
        
        # reshape modes: from 1D array to 2D array (an image)
        modes = np.asarray(modesPOD).T.reshape(nMode,self.surfShape[1],self.surfShape[2])
        
        # compute the time dependent coefficients ai
        ai = np.zeros((nSnap,nMode))
        i = 0
        for i in range(0,nMode):
            ai[:,i] = np.inner(surfacesPOD.T,modesPOD[:,i]).T

        self.result['modes']=modes
        self.result['eigVals']=eigVals
        self.result['ai']=ai
        
    def plotEnegry(self,ax,start,end,**kwargs):
        ax.plot(self.result['eigVals'][start:end]/np.sum(self.result['eigVals'][start:]),**kwargs)
        
    def plotCumulativeEnegry(self,ax,start,end,**kwargs):
        ax.plot(np.cumsum(self.result['eigVals'][start:end])/np.sum(self.result['eigVals'][start:]),**kwargs)
        
    def reconstructFrame(self,frame,modeIdxList):
        res=np.zeros(shape=self.result['modes'].shape[1:])
        for idx in modeIdxList:
            res=res+self.result['ai'][frame,idx]*self.result['modes'][idx]
        return res
            
    def reconstructFrames(self,frameList,modeIdxList):
        res=np.zeros(shape=(len(frameList),self.result['modes'].shape[1],self.result['modes'].shape[2]))
        for frame in frameList:
            for idx in modeIdxList:
                res[frame,:]=res[frame,:]+self.result['ai'][frame,idx]*self.result['modes'][idx]
        return res
                
    def clearInputData(self):
        self.surfaces=None
        
        
class DMD(object):
    def __init__(self):
        pass