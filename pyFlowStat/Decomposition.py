import numpy as np
import modred
import scipy.signal as sp
import h5py

class POD(object):
    def __init__(self,vecs):
        '''
        Arguments:
            *vecs*: numpy array of shape (elements,snapshots).
        '''
        self.vecs = vecs
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
        
        self.result['nMode']=nMode
        self.result['nSnap']=nSnap
        
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
        self.result['raw_modes']=modesPOD
        self.result['eigVals']=eigVals
        
        ai=self.projectOnMode(self.vecs)
        self.result['ai']=ai
    
    def projectOnMode(self,vec):
        nSnap = self.result['nSnap']
        nSnap=vec.shape[1]
        nMode = self.result['nMode']
        
        # compute the time dependent coefficients ai
        ai = np.zeros((nSnap,nMode))
        i = 0
        for i in range(0,nMode):
            ai[:,i] = np.inner(vec.T,self.result['raw_modes'][:,i]).T
        return ai
        
        
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
        
    def saveResult(self,filename):
        '''
        experimental: works with DMDPiv class
        '''
        
        self.result['inputShape']=self.inputShape
        #self.result['dt']=self.dt
        
        keys = self.result.keys()
        fwm = h5py.File(filename, 'w')
        gDict = fwm.create_group('dict')
        for k in keys:
            gDict.create_dataset(k,data=self.result[k])
        
        fwm.close()
        
    def loadResult(self,filename):
        '''
        experimental: works with DMDPiv class
        '''

        fwm = h5py.File(filename, 'r')
        
        keys = fwm['dict'].keys()
        
        for k in keys:
            self.result[k]=fwm['dict'][k].value
            
        fwm.close()
        #self.result['modes']=np.matrix(self.result['modes'])
        #self.dt=self.result['dt']
        self.inputShape=tuple(self.result['inputShape'])

class PODPiv(POD):
    def __init__(self,surfaceList,subslice=None):
        '''
        Arguments:
            *surfaceList*: numpy array of Surfaces (N).
            *dt*: float, timestep
            *subslice*: numpy index tuples, created wit np.s_
        '''
        if subslice is None:
            FourD=np.array([[s.data['Ux'],s.data['Uy'],s.data['Uz']] for s in surfaceList])
        else:
            FourD=np.array([[s.data['Ux'][subslice],s.data['Uy'][subslice],s.data['Uz'][subslice]] for s in surfaceList])
        inputShape=FourD.shape
        vecs = FourD.reshape((inputShape[0],np.prod(inputShape[1:]))).T
        
        super(PODPiv,self).__init__(vecs)
        self.inputShape=inputShape
   
    def decompose(self,nMode,method='snap',subtractMean=False):
        super(PODPiv,self).decompose(nMode=nMode,method=method,subtractMean=subtractMean)
        # reshape modes: from 1D array to 2D array (an image)
        modes = np.asarray(self.result['raw_modes']).T.reshape(nMode,self.inputShape[1],self.inputShape[2],self.inputShape[3])
        self.result['modes']=modes 
        
class PODscalarFieldList(POD):
    def __init__(self,scalarFieldList):
        '''
        Arguments:
            *scalarFieldList*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example one component of a SurfaceList).
        '''
        surfShape=scalarFieldList.shape
        vecs = scalarFieldList.reshape(surfShape[0],surfShape[1]*surfShape[2]).T
        super(PODscalarFieldList,self).__init__(vecs)
        self.surfaces=scalarFieldList
        self.surfShape=self.surfaces.shape
        
    def decompose(self,nMode,method='snap',subtractMean=False):
        
        super(PODscalarFieldList,self).decompose(nMode=nMode,method=method,subtractMean=subtractMean)
        # reshape modes: from 1D array to 2D array (an image)
        modes = np.asarray(self.result['raw_modes']).T.reshape(nMode,self.surfShape[1],self.surfShape[2])
        self.result['modes']=modes 
        
class DMD(object):
    def __init__(self,vecs,dt):
        '''
        Arguments:
            *scalarFieldList*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example one component of a SurfaceList).
        '''
        self.vecs=vecs
        self.result=dict()
        self.dt=dt
        
    def decompose(self,nMode=0,method='snap',subtractMean=False):
        '''
        Compute the Dynamic Mode Decomposition/Koopman Mode Decomposition (DMD)
        
        Arguments:
            *nMode*: python integer.
             Number of modes of theDMD.
            
            *method*: python string. Default='snap'
             Type of DMD algorithm used. For the snapshot method, use DMDmethod='snap' and
             for the direct method, use DMDmethod='direct'. The default value is 'snap'.
             
             How to choose between 'snap' and 'direct': if surfX*surfY > N**2, use the snap
             method, it should increase the coputational speed for only a little loss in 
             precision.
             
            *subtractMean*: python bool. Default=False.
             If True, the DMD reduces to a DFT, see:
             Chen, K., Tu, J., & Rowley, C. (2012). Variants of dynamic mode 
             decomposition: boundary condition, Koopman, and Fourier analyses.
             Journal of Nonlinear Science.
            
        Remarks:
        use the key 'mode_norms' to plot the power spectral density, since:
        myDMD_piv.result['mode_norms'] =  np.linalg.norm(myDMD_piv.result['modes'],axis=0)**2
        
        '''

            
        nSnap = self.vecs.shape[1]
        nVecs=self.vecs.shape[0]
        
        if nMode>nSnap-1 or nMode==0:
            nMode=nSnap-1

        if nMode>nVecs:
            nMode=nVecs
        
        if np.any(self.vecs):
            self.vecs = np.nan_to_num(self.vecs)
        
        if subtractMean:
            self.vecs=self.vecs-np.mean(self.vecs,axis=1,keepdims=True)
        #print nMode
        if method=='snap':
            modes, ritz_vals, mode_norms = modred.compute_DMD_matrices_snaps_method(self.vecs, range(nMode))
        elif method=='direct':
            modes, ritz_vals, mode_norms = modred.compute_DMD_matrices_snaps_method(self.vecs, range(nMode))
        else:
            print('error: argument '+str(method)+' is not valid. Use \'snap\' or \'direct\'')
            
        self.result['modes']=modes
        self.result['ritz_vals']=ritz_vals
        self.result['mode_norms']=mode_norms
        self.result['ai']=np.array([self.result['ritz_vals']**t for t in range(self.vecs.shape[1])])
        self.result['m']=self.result['modes'].shape[1]
        
    def saveResult(self,filename):
        '''
        experimental: works with DMDPiv class
        '''
        
        self.result['inputShape']=self.inputShape
        self.result['dt']=self.dt
        
        keys = self.result.keys()
        fwm = h5py.File(filename, 'w')
        gDict = fwm.create_group('dict')
        for k in keys:
            gDict.create_dataset(k,data=self.result[k])
        
        fwm.close()
        
    def loadResult(self,filename):
        '''
        experimental: works with DMDPiv class
        '''

        fwm = h5py.File(filename, 'r')
        
        keys = fwm['dict'].keys()
        
        for k in keys:
            self.result[k]=fwm['dict'][k].value
            
        fwm.close()
        self.result['modes']=np.matrix(self.result['modes'])
        self.dt=self.result['dt']
        self.inputShape=tuple(self.result['inputShape'])
        
    def getKey(self,key,idx=None):
        if idx==None:
            return self.result[key]
        else:
            return self.result[key][idx]    
            
    def getEig(self,idx=None):
        return self.getKey('ritz_vals',idx=idx)
            
    def getEigReal(self,idx=None):
        return np.real(self.getEig(idx=idx))

    def getEigImag(self,idx=None):
        return np.imag(self.getEig(idx=idx))    
        
    def getEigAbs(self,idx=None):
        return np.absolute(self.getEig(idx=idx))
        
    def getEigAngle(self,idx=None):
        return np.angle(self.getEig(idx=idx))
        
    def getMode(self,nMode):
        return self.getModeLst()[nMode]
        
    def getModeLst(self):
        return np.asarray(self.result['modes']).T.reshape((self.result['modes'].shape[1],)+self.inputShape[1:])

    def getResidualVec(self):
        k=self.result['m']
        xm=self.vecs[:,k]
        xm_rec=DMD.reconstructDMD(self,k)
        return xm-xm_rec
        
    def getResidual(self):
        return np.linalg.norm(self.getResidualVec())
        
    def getResidualRelative(self):
        k=self.result['m']
        xm=self.vecs[:,k]
        xm_rec=DMD.reconstructDMD(self,k)
        return np.linalg.norm(xm-xm_rec)/np.linalg.norm(xm)

    def reconstructDMD(self,k,idx=None,verbose=False):
        '''
        TODO: correct 
        '''
        tmp=np.zeros(shape=self.result['modes'].shape[0])
        if idx==None:
            idx=range(self.result['modes'].shape[1])
        for i in idx:
            t=self.result['ritz_vals'][i]**k*self.result['modes'][:,i]
            tmp=tmp+np.array(t)[:,0]

        return tmp
    
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
            
    def getNorm(self,idx=None):
        return self.getKey('mode_norms',idx=idx)
        
    def getNormSortedIdx(self,idx=None):
        '''
        returns norm_idx
        '''
        norm_idx=np.argsort(self.getNorm())
        return norm_idx

    def getSortedIdx(self,function,idx=None):
        '''
        returns norm_idx
        '''
        if idx==None:
            sorted_idx=np.argsort(function())
            return sorted_idx
        else:
            tmp=function()
            #print tmp
            sorted_idx=idx[np.argsort(tmp[idx])]
            return sorted_idx
        
    def getGrowthSortedIdx(self):
        '''
        returns growth_idx
        '''
        growth_idx=np.argsort(self.getEigAbs())
        return growth_idx

    def getStabSortedIdx(self):
        '''
        returns indices of modes sorted by stability
        with 0 being the closest to the unity circle
        '''
        stab_idx=np.argsort(np.abs(np.log(self.getEigAbs())))
        return stab_idx
        
    def getFilteredIndex(self,min_norm=1.0,min_abs=0.001,positiveOnly=False):
        '''
        returns indices that according to filter
        '''
        idx_norm=np.where(self.result['mode_norms']>=min_norm)
        idx_stab=np.where(np.abs(1.0-self.getEigAbs())<=min_abs)
        idx_intersect=np.intersect1d(idx_norm[0],idx_stab[0])
        if positiveOnly:
            idx_intersect=idx_intersect[np.where(self.getEigImag(idx_intersect)>=0)]
        return self.getSortedIdx(self.getNorm,idx=idx_intersect)
    
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
        
class DMDPiv(DMD):
    def __init__(self,surfaceList,dt,subslice=None,filter_kernel_size=0):
        '''
        Arguments:
            *surfaceList*: numpy array of Surfaces (N).
            *dt*: float, timestep
            *subslice*: numpy index tuples, created wit np.s_
        '''
        if subslice is None:
            FourD=np.array([[s.data['Ux'],s.data['Uy'],s.data['Uz']] for s in surfaceList])
        else:
            FourD=np.array([[s.data['Ux'][subslice],s.data['Uy'][subslice],s.data['Uz'][subslice]] for s in surfaceList])
            
        if filter_kernel_size>0:
            print 'start smoothing'
            for s in FourD:
                s[0]=sp.medfilt(s[0],kernel_size=filter_kernel_size)
                s[1]=sp.medfilt(s[1],kernel_size=filter_kernel_size)
                s[2]=sp.medfilt(s[2],kernel_size=filter_kernel_size)
                '''
                '''
            print 'stop'
        inputShape=FourD.shape
        vecs = FourD.reshape((inputShape[0],np.prod(inputShape[1:]))).T
        
        super(DMDPiv,self).__init__(vecs,dt)
        self.inputShape=inputShape
        
    def reconstructDMD(self,k,idx=None,verbose=False):
        tmp = DMD.reconstructDMD(self,k=k,idx=idx,verbose=verbose)
        #tmp = super(DMDscalarFieldList,self).reconstructDMD(k=k,idx=idx,verbose=verbose)
        return tmp.T.reshape(self.inputShape[1:])
        
    def getMode(self,nMode,component=0):
        return DMD.getMode(self,nMode)[component]

class DMDscalarFieldList(DMD):
    def __init__(self,scalarFieldList,dt,subslice=None):
        '''
        Arguments:
            *scalarFieldList*: numpy array of shape (N,surfX,surfY).
             Surfaces from a PIV measurment (for example one component of a SurfaceList).
        '''
        inputShape=scalarFieldList[subslice].shape
        vecs = scalarFieldList[subslice].reshape(inputShape[0],np.prod(inputShape[1:])).T
        super(DMDscalarFieldList,self).__init__(vecs,dt)
        #self.surfaces=scalarFieldList
        self.inputShape=inputShape
        
    def reconstructDMD(self,k,idx=None,verbose=False):
        tmp = DMD.reconstructDMD(self,k=k,idx=idx,verbose=verbose)
        #tmp = super(DMDscalarFieldList,self).reconstructDMD(k=k,idx=idx,verbose=verbose)
        return tmp.T.reshape(self.inputShape[1:])
        
