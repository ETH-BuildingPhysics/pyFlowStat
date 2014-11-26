import numpy as np
import modred

class PODPiv(object):
    def __init__(self,surfaceList):
        '''
        Arguments:
            *scalarFieldList*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example one component of a SurfaceList).
        '''
        FourD=np.array([[s.data['Ux'],s.data['Uy'],s.data['Uz']] for s in surfaceList])
        inputShape=FourD.shape
        self.inputShape=inputShape
        self.vecs=FourD.reshape((inputShape[0],inputShape[1]*inputShape[2]*inputShape[3])).T

        self.result=dict()
        pass

    def decompose(self,nMode,method='snap',subtractMean=False):
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

        nSnap = self.vecs.shape[1]
        
        # reshape the input
        #vecs = self.surfaces.reshape(self.surfShape[0],self.surfShape[1]*self.surfShape[2]).T

        if np.any(self.vecs):
            self.vecs = np.nan_to_num(self.vecs)
        
        if subtractMean:
            self.vecs=self.vecs-np.mean(self.vecs,axis=1,keepdims=True)
            
        if method=='snap':
            modesPOD, eigVals = modred.compute_POD_matrices_snaps_method(self.vecs, range(nMode))
        elif method=='direct':
            modesPOD, eigVals = modred.compute_POD_matrices_direct_method(self.vecs, range(nMode))
        else:
            print('error: argument '+str(method)+' is not valid. Use \'snap\' or \'direct\'')
        
        # convert modesPOD in an numpy array
        modesPOD = np.asarray(modesPOD)
        
        # reshape modes: from 1D array to 2D array (an image)
        #modes = np.asarray(modesPOD).T.reshape(nMode,self.surfShape[1],self.surfShape[2])
        modes = np.asarray(modesPOD).T.reshape(nMode,self.inputShape[1],self.inputShape[2],self.inputShape[3])
        
        # compute the time dependent coefficients ai
        ai = np.zeros((nSnap,nMode))
        i = 0
        for i in range(0,nMode):
            ai[:,i] = np.inner(self.vecs.T,modesPOD[:,i]).T

        self.result['modes']=modes
        self.result['eigVals']=eigVals
        self.result['ai']=ai
        
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
    
    def decompose(self,nMode,method='snap',subtractMean=False):
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
        vecs = self.surfaces.reshape(self.surfShape[0],self.surfShape[1]*self.surfShape[2]).T

        if np.any(vecs):
            vecs = np.nan_to_num(vecs)
        
        if subtractMean:
            vecs=vecs-np.mean(vecs,axis=1,keepdims=True)
            
        if method=='snap':
            modesPOD, eigVals = modred.compute_POD_matrices_snaps_method(vecs, range(nMode))
        elif method=='direct':
            modesPOD, eigVals = modred.compute_POD_matrices_direct_method(vecs, range(nMode))
        else:
            print('error: argument '+str(method)+' is not valid. Use \'snap\' or \'direct\'')
        
        # convert modesPOD in an numpy array
        modesPOD = np.asarray(modesPOD)
        
        # reshape modes: from 1D array to 2D array (an image)
        modes = np.asarray(modesPOD).T.reshape(nMode,self.surfShape[1],self.surfShape[2])
        
        # compute the time dependent coefficients ai
        ai = np.zeros((nSnap,nMode))
        i = 0
        for i in range(0,nMode):
            ai[:,i] = np.inner(vecs.T,modesPOD[:,i]).T

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
    def __init__(self,scalarFieldList,dt):
        '''
        Arguments:
            *scalarFieldList*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example one component of a SurfaceList).
        '''
        self.surfaces=scalarFieldList
        self.surfShape=self.surfaces.shape
        self.result=dict()
        self.dt=dt
    
    def eig(self):
        return self.result['ritz_vals']
        
    def eigReal(self):
        return np.real(self.result['ritz_vals'])
        
    def eigImag(self):
        return np.imag(self.result['ritz_vals'])    
        
    def eigAbs(self):
        return np.absolute(self.result['ritz_vals'])
        
    def eigAngle(self):
        return np.angle(self.result['ritz_vals'])
        
    def decompose(self,nMode=0,method='snap',subtractMean=False):
        '''
        Compute the Dynamic Mode Decomposition/Koopman Mode Decomposition (DMD) from 
        a list of N snapshots. 
        The snapshots are PIV surfaces of size (surfX,surfY) of field F. F can be any
        scalar field (Ux, ux, T, vorticity, R11,...)
        
        Arguments:
            *surfaces*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example).
            
            *nMode*: python integer.
             Number of modes of theDMD.
            
            *DMDmethod*: python string. Default='snap'
             Type of DMD algorithm used. For the snapshot method, use DMDmethod='snap' and
             for the direct method, use DMDmethod='direct'. The default value is 'snap'.
             
             How to choose between 'snap' and 'direct': if surfX*surfY > N**2, use the snap
             method, it should increase the coputational speed for only a little loss in 
             precision.
             
            *returnDMD*: python bool. Default=False.
             If True, the 1D DMD vectors are also returned
             
        Returns:
            *tons of cool stuff...*
        
        '''
        if nMode>self.surfShape[0]-1 or nMode==0:
            nMode=self.surfShape[0]-1
        # reshape the input
        surfaces = self.surfaces.reshape(self.surfShape[0],self.surfShape[1]*self.surfShape[2]).T

        if np.any(surfaces):
            surfaces = np.nan_to_num(surfaces)
            
        if subtractMean:
            surfaces=surfaces-np.mean(surfaces,axis=1,keepdims=True)
            
        if method=='snap':
            modes, ritz_vals, mode_norms = modred.compute_DMD_matrices_snaps_method(surfaces, range(nMode))
        elif method=='direct':
            modes, ritz_vals, mode_norms = modred.compute_DMD_matrices_snaps_method(surfaces, range(nMode))
        else:
            print('error: argument '+str(method)+' is not valid. Use \'snap\' or \'direct\'')

        self.result['modes']=modes
        self.result['ritz_vals']=ritz_vals
        self.result['mode_norms']=mode_norms
        
    def getMode(self,nMode):
        return self.getModeLst()[nMode]
        
    def getModeLst(self):
        return np.asarray(self.result['modes']).T.reshape(self.result['modes'].shape[1],self.surfShape[1],self.surfShape[2])
    
    def reconstructDMD(self,k,idx=None,verbose=False):
        tmp=np.zeros(shape=self.result['modes'].shape[0])
        if idx==None:
            idx=range(len(self.result['ritz_vals']))
        for i in idx:
            t=self.result['ritz_vals'][i]**k*self.result['modes'][:,i]
            #print t,dmd_ritz_vals[i]
            tmp=tmp+np.array(t)[:,0]
        if verbose:
            print 'Frq Band:',np.imag(self.result['ritz_vals'][idx[0]])*600/(2*np.pi),'Hz -',np.imag(self.result['ritz_vals'][idx[-1]])*600/(2*np.pi),'Hz'
        tmp_real=np.asarray(np.real(tmp)).T.reshape(self.surfShape[1],self.surfShape[2])
        tmp_imag=np.asarray(np.imag(tmp)).T.reshape(self.surfShape[1],self.surfShape[2])
        return tmp, tmp_real, tmp_imag
    
    def getFrqSortedIdx(self,verbose=False):
        '''
        Helper Function
        returns frq_idx_pos,frq_idx_neg
        '''
        frq_idx=np.argsort(np.angle(self.result['ritz_vals']))
        frq_idx_pos=[i for i in frq_idx if np.angle(self.result['ritz_vals'][i])>0]
        frq_idx_neg=[i for i in frq_idx if np.angle(self.result['ritz_vals'][i])<0][::-1]
        if len(frq_idx_pos)!=len(frq_idx_neg):
            if verbose:
                print 'getFrqSortedIdx: Not same length! Pos:', len(frq_idx_pos),'Neg:',len(frq_idx_neg)
            l=min(len(frq_idx_pos),len(frq_idx_neg))
            frq_idx_pos=frq_idx_pos[:l+1]
            frq_idx_neg=frq_idx_neg[:l+1]
        return frq_idx_pos,frq_idx_neg,self.getFrqList(frq_idx_pos)
        
    def getNormSortedIdx(self):
        '''
        Helper Function
        returns frq_idx_pos,frq_idx_neg
        '''
        norm_idx=np.argsort(self.result['mode_norms'])
        return norm_idx
        
    def getIdxforFrq(self,f,verbose=False):
        '''
        Get the dmd index corresponding to a frequency
        '''
        array=np.angle(self.result['ritz_vals'])/self.dt/(2.0*np.pi)
        idx = (np.abs(array-f)).argmin()
        if verbose:
            print 'Target f:',f,'Hz, Closest Match:',array[idx],'Hz'
        return idx

    def __getIdxforFrq__(self,frqs,f,verbose=False):
        '''
        Get the dmd index corresponding to a frequency
        '''
        #array=np.angle(self.result['ritz_vals'])/self.dt/(2.0*np.pi)
        array=np.array(frqs)
        idx = (np.abs(array-f)).argmin()
        if verbose:
            print 'Target f:',f,'Hz, Closest Match:',array[idx],'Hz'
        return idx
        
    def getFrq(self,idx):
        return np.angle(self.result['ritz_vals'][idx])/self.dt/(2.0*np.pi)
        
    def getFrqStr(self,idx,digits=2):
        return str(np.round(self.getFrq(idx),digits))+' Hz'
        
    def getFrqList(self,idxLst=None):
        if idxLst==None:
            return np.angle(self.result['ritz_vals'])/self.dt/(2.0*np.pi)
        else:
            return [self.getFrq(i) for i in idxLst] 
        
    def getIdxFrqBand(self,f1,f2,verbose=False):
        '''
        Get the dmd indieces corresponding to a frequency band
        '''
        if verbose:
            print '**getIdxFrqBand**'
        frq_idx_pos,frq_idx_neg,f=self.getFrqSortedIdx()

        idx_f1=self.__getIdxforFrq__(f,f1,verbose=verbose)
        idx_f2=self.__getIdxforFrq__(f,f2,verbose=verbose)
    
        res=np.sort(np.hstack([frq_idx_pos[idx_f1:idx_f2+1],frq_idx_neg[idx_f1:idx_f2+1]]))
        if verbose:
            print 'SubIdx:',idx_f1,'up to',idx_f2,',length',len(range(idx_f1,idx_f2+1)),'/',len(frq_idx_pos)
            print 'Total Points:',len(res),'\n'
        return res
