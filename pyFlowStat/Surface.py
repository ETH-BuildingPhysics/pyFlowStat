#===========================================================================#
# load modules
#===========================================================================#
#standard modules
#import sys

#scientific modules
import numpy as np
#import scipy as sp
import os
import matplotlib.tri as tri
from pyFlowStat.TriSurfaceMesh import TriSurfaceMesh
from pyFlowStat.TriSurfaceVector import TriSurfaceVector
from pyFlowStat.TriSurfaceScalar import TriSurfaceScalar
from pyFlowStat.TriSurfaceSymmTensor import TriSurfaceSymmTensor
#from pyFlowStat.TriSurface import parseFoamFile

# special modules
from ctypes import *

TypeName = ["Image", "2D-PIV-Vector (header, 4x(Vx,Vy))",
            "2D-Vector (Vx,Vy)", "2D-PIV+p.ratio (header, 4x(Vx,Vy), peakratio)",
          "3D-Vector (Vx,Vy,Vz)", "3D-Vector+p.ratio (header, 4x(Vx,Vy), peakratio)"]
          
WORD=c_ushort
 
class AttributeList(Structure):
    '''
    ctypes wrapper Davis VC7 struct AttributeList 
    
    #typedef struct AttributeList
    {
       char*          name;
       char*          value;
       AttributeList* next;
    } AttributeList;
    
    '''
    def __getattr__(self, key):
        if key=='pairs':
            self.get_pairs()
            return self.pairs
        if key=='dict':
            return self.as_dict()
        else:
            raise AttributeError(u"Does not have %s atribute" % key)
            
    def get_pairs(self):
        att = self
        self.pairs = []
        while att!=0:
            try:
                self.pairs.append((att.name, att.value))
                att = att.next[0]
            except ValueError:
                break
    
    def as_dict(self):
        self.get_pairs()
        return dict(self.pairs)
        
    def delete(self):
        del_attributelist(self)
        
AttributeList._fields_=[("name",c_char_p),("value",c_char_p),("next",POINTER(AttributeList))]


class _bufarray(Union):
    '''
    ctypes wrapper Davis VC7 union
    union
    {
          float*   floatArray;
          Word*    wordArray;
    };
    '''
    _fields_=[("floatArray",POINTER(c_float)),("wordArray",POINTER(WORD))]


class BufferScaleType(Structure):
    '''
    ctypes wrapper for Davis VC7 struct BufferScaleType
    
    typedef struct
    {
    	float	factor;
    	float offset;
    	char	description[16];
    	char	unit[16];
    } BufferScaleType;
    
    '''
    _fields_=[("factor",c_float),("offset",c_float),("description",c_char*16),
              ("unit",c_char*16)]


class BufferType(Structure):
    '''
    ctypes wrapper for Davis VC7 class BufferType
    
    typedef struct
    {
      int         isFloat;
      int         nx,ny,nz,nf;
      int         totalLines;
    	int			vectorGrid;			// 0 for images
    	int			image_sub_type;	// BufferFormat_t
      union
    	{
          float*   floatArray;
          Word*    wordArray;
      };
    	BufferScaleType	scaleX;		// x-scale
    	BufferScaleType	scaleY;		// y-scale
    	BufferScaleType	scaleI;		// intensity scale
    	bool*			bMaskArray;			// mask array, NULL if no mask exists
    } BufferType;
    
    '''
    _anonymous_ = ("bufarray",)
    _fields_=[("isFloat",c_int),("ny",c_int),("nx",c_int),("nz",c_int),("nf",c_int),
             ("totalLines",c_int),("vectorGrid",c_int),("image_sub_type",c_int),
             ("bufarray",_bufarray),("scaleX",BufferScaleType),("scaleY",BufferScaleType),
             ("scaleI",BufferScaleType),("bMaskArray",POINTER(c_bool))]

def getMode(buf,theX_,theY_,width_,frameOffset):
    '''
    helper method to get mode from Davis VC7 buffer
    '''
    mode = int(buf.floatArray[theX_ + theY_*width_ + frameOffset])
    if mode<0:
        return -1
    elif mode>4:
        #// interpolated or filled vector
        mode = 4
    mode=mode-1
    return mode

class Surface(object):
    '''
    Holds 2D data on a equidistant,cartesian grid
    
    Attributes:
      * vx,vy,vz (numpy ndarray): velocity data, saved on cell center.
      * dx,dy (float): cell size, in mm.
      * minX,maxX,minY,maxY (float): min/max position of cell centers in mm.
      * extent (list of floats): [minX-dx/2,maxX+dx/2,minY-dy/2,maxY+dy/2] in mm.
      * data (dict): dictionary to hold processed data, created by createDataDict().
        
    Note: Units for distances have to be in mm (dx,dy,minX,maxX,minY,maxY and extent)
    in order for gradients to be calculated correctly.
    '''
    def __init__(self):
        self.vx=[]
        self.vy=[]
        self.vz=[]
        
        self.dx = float()
        self.dy = float()        

        self.minX = float()
        self.maxX = float()
        self.minY = float()
        self.maxY = float()
        self.extent = []

        self.data=dict()
        return

    def createDataDict(self):
        '''
        Creates the "data" dictionnary from member variables vx, vy, vz, dx, dy

        Member variable data (python ditionary) is created.

        By default, the following keys are included in data:
            Ux:  [numpy.array.shape=(ny,nx)] Velocity Ux
            Uy:  [numpy.array.shape=(ny,nx)] Velocity Uy
            Uz:  [numpy.array.shape=(ny,nx)] Velocity Uz
            dx:    [float] spacing in x dirction
            dy:    [float] spacing in y dirction
        '''

        self.data = dict()
        self.data['Ux'] = self.vx
        self.data['Uy'] = self.vy
        self.data['Uz'] = self.vz
        self.data['dx'] = self.dx
        self.data['dy'] = self.dy

    def emptyCopy(self):
        s=Surface()
        s.vx=np.zeros(self.data['Ux'].shape)
        s.vy=np.zeros(self.data['Uy'].shape)
        s.vz=np.zeros(self.data['Uz'].shape)
        s.dx=self.dx
        s.dy=self.dy
        s.minX=self.minX
        s.maxY=self.maxY
        s.maxX=self.maxX
        s.minY=self.minY
        s.extent=self.extent
        return s

    def generateUmag(self):
        Umag = np.zeros(self.data['Ux'].shape)
        Umag = np.sqrt(self.data['Ux']**2+self.data['Uy']**2+self.data['Uz']**2)
        self.data['Umag']=Umag
        
    def generateUmagFluct(self):
        try:
            Umag = np.zeros(self.data['ux'].shape)
            Umag = np.sqrt(self.data['ux']**2+self.data['uy']**2+self.data['uz']**2)
            self.data['umag']=Umag
        except KeyError as err:
            print err.message
            print 'add fluctuating field using addReynoldsDecomposition()'
        
    def generateUmag2D(self):
        Umag2D = np.zeros(self.data['Ux'].shape)
        Umag2D = np.sqrt(self.data['Ux']**2+self.data['Uy']**2)
        self.data['Umag2D']=Umag2D
        
    def generateFields(self):
        '''
        Generates additional dictionary entries.
        '''
        self.generateUmag()
        self.generateUmag2D()

        self.computeGradients()
        self.computeVorticity()
        self.computeQ()
        self.computeSignedQ()
        self.computeOWQ()
        self.computeLambda2()
        
        dudx=self.data['dudx']
        dvdy=self.data['dvdy']
        self.data['Div2D']=dudx+dvdy
        
        self.data['KE']=0.5*(self.vx**2+self.vy**2+self.vz**2)
        
        #self.computeGradients(method='r')
        #self.computeGradients(method='ls')
        
        
#        tensorS= np.empty(self.data['Ux'].shape)
#        tensorW= np.empty(self.data['Ux'].shape)
#        tensorS= 0.5*[[dudx+dudx,dudy+dvdx],[dvdx+dudy,dvdy+dvdy]]
#        tensor2= 0.5*[[0.0,dudy-dvdx],[dvdx-dudy,0.0]]
        
    def computeGradients(self,method='numpy'):
        if method=='numpy':
            dudy,dudx=np.gradient(self.vx,-self.dy/1000,self.dx/1000)
            dvdy,dvdx=np.gradient(self.vy,-self.dy/1000,self.dx/1000)
            self.data['dudy']=dudy
            self.data['dudx']=dudx
            self.data['dvdy']=dvdy
            self.data['dvdx']=dvdx
        elif method=='r':
            dudx_r = np.zeros(self.data['Ux'].shape)
            for i in range(2,self.data['Ux'].shape[0]-2):
                for j in range(2,self.data['Ux'].shape[1]-2):
                    dudx_r[i,j]=(self.data['Ux'][i,j-2]-8.0*self.data['Ux'][i,j-1]+8.0*self.data['Ux'][i,j+1]-self.data['Ux'][i,j+2])/(12.0*self.dx/1000.0)
            self.data['dudx_r']=dudx_r
        elif method=='ls':
            dudx_ls = np.zeros(self.data['Ux'].shape)
            for i in range(2,self.data['Ux'].shape[0]-2):
                for j in range(2,self.data['Ux'].shape[1]-2):
                    dudx_ls[i,j]=(2.0*self.data['Ux'][i,j+2]+self.data['Ux'][i,j+1]-self.data['Ux'][i,j-1]-2.0*self.data['Ux'][i,j-2])/(10.0*self.dx/1000.0)
            self.data['dudx_ls']=dudx_ls
            
    def removeGradients(self):
        for k in self.data.keys():
            if k.startswith('dudx'):
                self.data.pop(k)
            if k.startswith('dudy'):
                self.data.pop(k)
            if k.startswith('dvdx'):
                self.data.pop(k)
            if k.startswith('dvdy'):
                self.data.pop(k)
        
    def computeQ(self):
        
        dudy=self.data['dudy']
        dudx=self.data['dudx']
        dvdy=self.data['dvdy']
        dvdx=self.data['dvdx']
        #self.data['Q']=np.zeros(self.data['Ux'].shape)
        self.data['Q']=0.5*(-2.0*dudy*dvdx-dudx**2-dvdy**2)
        
    def computeSignedQ(self):
        Q_sign=self.data['Q'].copy()
        Q_sign[Q_sign<0]=0.0
        Q_sign[self.data['VortZ']<0]=Q_sign[self.data['VortZ']<0]*-1.0
        self.data['Q_sign']=Q_sign
        
    def computeSwirlingStrength(self):
        
        dudy=self.data['dudy']
        dudx=self.data['dudx']
        dvdy=self.data['dvdy']
        dvdx=self.data['dvdx']
        #self.data['Q']=np.zeros(self.data['Ux'].shape)
        self.data['SwirlingStrength^2']=(1.0/(4.0*dudx))**2+(1.0/(4.0*dvdy))**2-0.5*dudx*dvdy+dvdx*dudy
        
    def computeOWQ(self):
        '''
        Okubo-Weiss
        '''
        dudy=self.data['dudy']
        dudx=self.data['dudx']
        dvdy=self.data['dvdy']
        dvdx=self.data['dvdx']
        #self.data['Q']=np.zeros(self.data['Ux'].shape)
        self.data['OW-Q']=(dudx-dvdy)**2+(dudy+dvdx)**2-(dvdx-dudy)**2
        
    def computeLambda2(self):
        dudy=self.data['dudy']
        dudx=self.data['dudx']
        dvdy=self.data['dvdy']
        dvdx=self.data['dvdx']
        self.data['lambda2'] = self.getLambda2(dudx,dudy,dvdx,dvdy)
        
    def computeVorticity(self):
        
        dudy=self.data['dudy']
        dvdx=self.data['dvdx']
        vort_z=dvdx-dudy
        self.data['VortZ']=vort_z
        
    def getLambda2(self,dudx,dudy,dvdx,dvdy):
        S11 = dudx
        S12 = 0.5*(dudy+dvdx)
        S21 = 0.5*(dvdx+dudy)
        S22 = dvdy
        S13=S23=S33=S31=S32 = np.zeros(dudx.shape)
        W13=W23=W33=W31=W32 = np.zeros(dudx.shape)
        W11 = np.zeros(dudx.shape)
        W12 = 0.5*(dudy-dvdx)
        W21 = 0.5*(dvdx-dudy)
        W22 = np.zeros(dudx.shape)

        P11=S11*S11+S12*S12+S13*S13-W12*W12-W13*W13
        P12=S12*(S11+S22)+S13*S23-W13*W23
        P13=S13*(S11+S33)+S12*S23+W12*W23
        P22=S12*S12+S22*S22+S23*S23-W12*W12-W23*W23
        P23=S23*(S22+S33)+S12*S13-W12*W13
        P33=S13*S13+S23*S23+S33*S33-W13*W13-W23*W23

        a=-1.0
        b=P11+P22+P33
        c=P12*P12+P13*P13+P23*P23-P11*P22-P11*P33-P22*P33
        d=P11*P22*P33+2.0*P12*P13*P23-P12*P12*P33-P13*P13*P22-P23*P23*P11

        x=((3.0*c/a)-b*b/(a*a))/3.0
        y=(2.0*b*b*b/(a*a*a)-9.0*b*c/(a*a)+27.0*d/a)/27.0
        z=y*y/4.0+x*x*x/27.0

        i=np.sqrt(y*y/4.0-z)
        j=-pow(i,1.0/3.0)
        k=np.arccos(-(y/(2.0*i)))
        m=np.cos(k/3.0)
        n=np.sqrt(3.0)*np.sin(k/3.0)
        p=b/(3.0*a)

        lam1=2.0*j*m+p;
        lam2=-j*(m+n)+p;
        lam3=-j*(m-n)+p;
        lam=np.zeros(lam1.shape)
        row,col=lam1.shape
        for arow in range(0,row):
            for acol in range(0,col):
                l=[]
                l.append(lam1[arow,acol])
                l.append(lam2[arow,acol])
                l.append(lam3[arow,acol])
                l.sort()
                lam[arow,acol]=l[1]
        return lam*-1.0

    def addReynoldsDecomposition(self,MeanFlowSurface,addReStresses=True):
        '''
        Generate fluctuations by subtracting the mean flow (surface of same size)
        Adds fluctuation fields ux,uy,uz and correleations uu,vv,ww,uv,uw and TKE
        '''
        self.data['ux']=self.data['Ux']-MeanFlowSurface.data['Ux']
        self.data['uy']=self.data['Uy']-MeanFlowSurface.data['Uy']
        self.data['uz']=self.data['Uz']-MeanFlowSurface.data['Uz']
        if addReStresses:
            self.data['uu']=self.data['ux']**2
            self.data['vv']=self.data['uy']**2
            self.data['ww']=self.data['uz']**2
            self.data['uv']=self.data['ux']*self.data['uy']
            self.data['uw']=self.data['ux']*self.data['uz']
            self.data['vw']=self.data['uy']*self.data['uz']
            self.data['TKE']=0.5*(self.data['uu']+self.data['vv']+self.data['ww'])

    def readFromVC7(self,filename,v=False):
        '''
        reads PIV vector data in tha Davis format, using the 64bit windows DLL
        '''
        dllpath = os.path.dirname(os.path.realpath(__file__))
        ReadIMX64 = cdll.LoadLibrary(dllpath+"\ReadIMX64.dll")

        tmpBuffer = BufferType()
        attributeLst = AttributeList()
        self.vx=[]
        self.vy=[]
        self.vz=[]

        res = ReadIMX64.ReadIM7(filename, byref(tmpBuffer), byref(attributeLst))
        
        #print len(self.attributeLst)
        #attv=att.as_dict()['_SCALE_X']
        #print attv
        #print res
        if res>0:
            print "Error reading image"
            return

        if v:
            print "Size (ny, nx)"
            print tmpBuffer.ny
            print tmpBuffer.nx
            print TypeName[tmpBuffer.image_sub_type]

        theFrame=0
        frameOffset = theFrame * tmpBuffer.nx * tmpBuffer.ny * tmpBuffer.nz;
        width = tmpBuffer.nx;
        height = tmpBuffer.ny;
        componentOffset = width * height;


        self.vx=np.empty((width,height), dtype=float)
        self.vx[:] = np.NAN
        self.vy=np.empty((width,height), dtype=float)
        self.vy[:] = np.NAN
        self.vz=np.empty((width,height), dtype=float)
        self.vz[:] = np.NAN

        if tmpBuffer.image_sub_type == 3 or tmpBuffer.image_sub_type == 1:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(tmpBuffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX] = np.sign(tmpBuffer.scaleX.factor)*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+1)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vy[theY,theX] = np.sign(tmpBuffer.scaleY.factor)*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+2)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    else:
                        pass
        if tmpBuffer.image_sub_type == 4:
            for theY in range(0,width):
                for theX in range(0,height):
                    self.vx[theY,theX] = np.sign(tmpBuffer.scaleX.factor)*(tmpBuffer.floatArray[theX + theY*height + frameOffset]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    self.vy[theY,theX] = np.sign(tmpBuffer.scaleY.factor)*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    self.vz[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*2]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
        if tmpBuffer.image_sub_type == 5:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(tmpBuffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX]=np.sign(tmpBuffer.scaleX.factor)*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+1)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vy[theY,theX]=np.sign(tmpBuffer.scaleY.factor)*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+2)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vz[theY,theX]=(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+3)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    else:
                        pass
#        print tmpBuffer.scaleX.factor
#        print tmpBuffer.scaleX.offset
#        print tmpBuffer.scaleY.factor
#        print tmpBuffer.scaleY.offset
        #self.maxX=
        #self.maxY=
        if np.isnan(self.vz).all():
            self.vz.fill(0)

        self.dx=abs(tmpBuffer.scaleX.factor*tmpBuffer.vectorGrid)
        self.dy=abs(tmpBuffer.scaleY.factor*tmpBuffer.vectorGrid)
        self.minX=tmpBuffer.scaleX.factor*tmpBuffer.vectorGrid*(0.5)+tmpBuffer.scaleX.offset
        self.maxY=tmpBuffer.scaleY.factor*tmpBuffer.vectorGrid*(0.5)+tmpBuffer.scaleY.offset
        self.maxX=tmpBuffer.scaleX.factor*tmpBuffer.vectorGrid*(tmpBuffer.ny-0.5)+tmpBuffer.scaleX.offset
        self.minY=tmpBuffer.scaleY.factor*tmpBuffer.vectorGrid*(tmpBuffer.nx-0.5)+tmpBuffer.scaleY.offset
        self.extent=[self.minX-(self.dx/2),self.maxX+(self.dx/2),self.minY-(self.dy/2),self.maxY+(self.dy/2)]
        ReadIMX64.DestroyBuffer(tmpBuffer)
        self.createDataDict()
        #plot(vx)
        #plot(vy)
        #plot(vz)

#        			mode = (int) theBuffer.floatArray[ theX + theY*width + frameOffset ];
#			if (mode<=0)
#			{	// disabled vector
#				return true;
#			}
#			if (mode>4)
#			{	// interpolated or filled vector
#				mode = 4;
#			}
#			mode--;
#			vx = theBuffer.floatArray[ theX + theY*width + frameOffset + componentOffset*(mode*2+1) ];
#			vy = theBuffer.floatArray[ theX + theY*width + frameOffset + componentOffset*(mode*2+2) ];

    def readFromIM7(self,filename,key,frame=0,v=False,scale=1.0):
        '''
        reads PIV image data in tha Davis format, using the 64bit windows DLL
        '''

        dllpath = os.path.dirname(os.path.realpath(__file__))
        ReadIMX64 = cdll.LoadLibrary(dllpath+"\ReadIMX64.dll")

        tmpBuffer = BufferType()
        self.attributeLst = AttributeList()
        res = ReadIMX64.ReadIM7(filename, byref(tmpBuffer), byref(self.attributeLst))
        s=None
        
        if tmpBuffer.image_sub_type < 0:
            if v:
                print "Size (ny, nx)"
                print tmpBuffer.ny
                print tmpBuffer.nx
                print tmpBuffer.nf
                print 'type:',tmpBuffer.image_sub_type

            theFrame=frame
            width = tmpBuffer.ny;
            height = tmpBuffer.nx;
            if tmpBuffer.isFloat:
                s=np.empty((height,width), dtype=float)
                s[:] = np.NAN
            else:
                s=np.empty((height,width), dtype=int)
                s[:] = np.NAN
            if(tmpBuffer.isFloat):
                for y in range(0,tmpBuffer.ny):
                    for x in range(0,tmpBuffer.nx): 
                        s[x,y] = (tmpBuffer.floatArray[theFrame*tmpBuffer.nx*tmpBuffer.ny+x*tmpBuffer.ny+y]*tmpBuffer.scaleI.factor*scale+tmpBuffer.scaleI.offset)
            else:
                for y in range(0,tmpBuffer.ny):
                    for x in range(0,tmpBuffer.nx):
                        s[x,y] = (tmpBuffer.wordArray[theFrame*tmpBuffer.nx*tmpBuffer.ny+x*tmpBuffer.ny+y]*tmpBuffer.scaleI.factor*scale+tmpBuffer.scaleI.offset)
        
        ReadIMX64.DestroyBuffer(tmpBuffer)
        self.data[key]=s

    def interpolateField(self,values,grid_x,grid_y,triangulation,method='cubic',kind='min_E'):
        '''
        helper function
        method=linear,cubic (default)
        kind = geom, min_E (default)
        '''
        if method=='cubic':
            itp=tri.CubicTriInterpolator(triangulation,values,kind=kind)
        elif method=='linear':
            itp=tri.LinearTriInterpolator(triangulation,values)
        else:
            itp=tri.CubicTriInterpolator(triangulation,values,kind=kind)
        zi_ma = itp(grid_x, grid_y)
        zi=zi_ma.filled(np.nan)

        return zi

    def readFromFoamFile(self,
                         pointsFile,
                         facesFile,
                         velFile,
                         scalarFileList=[],
                         symTensorFileList=[],
                         viewAnchor=(0,0,0),
                         xViewBasis=(1,0,0),
                         yViewBasis=(0,1,0),
                         dx=None,
                         dy=None,
                         interpolationMethod='cubic',
                         kind='min_E'):
        '''
        Read an OpenFOAM surface (triangulated grid) in the current Surface
        object (cartesian grid). As the "grid" change (tri to cartesian), the
        value must be interpolated.
        
        
        Arguments:
            *pointFile*: python string.
             Point file  generate by OpenFOAM. This is the grid point
             coordinates.
            
            *facesFile*: python string.
             Face file generate by OpenFOAM. It is a list of triangles, which
             compose the grid.
            
            *velFile*: python string.
             Vector file generate by OpenFOAM. This is the data associated with
             each grid point.
            
            *scalarFileList*: python list.
            
            *symTensorFileList*: python list.
            
            *dx*: python float.
             Physical size of a pixel in the Surface class (x discretisation).
             Must be given in mm.
            
            *dy*: python float.
             Physical size of a pixel in the Surface class (y discretisation).
             Must be given in mm.
            
            *interpolationMethod*: python string. 
             Interpolation method used to interpolate from the triangulated
             grid to the cartesian grid. "cubic" or "linear". Default="cubic"
             
            *kind*: python string.
             Defines the algorithm used for the cubic interpolation. Choices:
             "min_E" or "geom". "min_E" should be the more accurate, but it is 
             also the most time time consuming.
             
        Returns:
            none
        '''

        print 'Reading Velocity'

        tsm = TriSurfaceMesh.readFromFoamFile(pointsFile=pointsFile,
                                              facesFile=facesFile,
                                              viewAnchor=viewAnchor,
                                              xViewBasis=xViewBasis,
                                              yViewBasis=yViewBasis)
                                              
        tsv = TriSurfaceVector.readFromFoamFile(varsFile=velFile,
                                                triSurfaceMesh=tsm,
                                                time=0,
                                                projectedField=False)                  

        points = np.vstack((tsv.x,tsv.y)).T
        
        print 'Creating Grid and Interpolator'
        if dx==None:
            dxlist=[a for a in np.abs(np.diff(points[:,0])) if a>0]
            dx=np.min(dxlist)
        if dy==None:
            dylist=[a for a in np.abs(np.diff(points[:,1])) if a>0]
            dy=np.min(dylist)

        MaxX=np.max(points[:,0])
        MinX=np.min(points[:,0])
        MaxY=np.max(points[:,1])
        MinY=np.min(points[:,1])
        extent=[MinX-dx/2,MaxX+dx/2,MinY-dy/2,MaxY+dy/2]

        cellsX=int((MaxX-MinX)/dx)+1
        cellsY=int((MaxY-MinY)/dy)+1

        grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
        triang = tsv.triangulation

        print 'Interpolating Velocity'
        vx_i=self.interpolateField(tsv.vx,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
        vy_i=self.interpolateField(tsv.vy,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
        vz_i=self.interpolateField(tsv.vz,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
        self.vx=np.flipud(vx_i)
        self.vy=np.flipud(vy_i)
        self.vz=np.flipud(vz_i)

        self.dx=dx
        self.dy=dy
        self.minX=MinX
        self.maxX=MaxX
        self.minY=MinY
        self.maxY=MaxY
        self.extent=extent
        self.createDataDict()

        for scalarFile in scalarFileList:
            varName=os.path.basename(scalarFile)
            print 'Reading Scalar',varName
            tsv.addFieldFromFoamFile(fieldFile=scalarFile,fieldname=varName)
            scalar_i=self.interpolateField(tsv[varName],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            self.data[varName]=np.flipud(scalar_i)

        for symTensorFile in symTensorFileList:
            varName=os.path.basename(symTensorFile)
            print 'Reading Tenstor',varName
            tsv.addFieldFromFoamFile(fieldFile=symTensorFile,fieldname=varName)
            tensor_11=self.interpolateField(tsv[varName][:,0],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            tensor_12=self.interpolateField(tsv[varName][:,1],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            tensor_13=self.interpolateField(tsv[varName][:,2],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            tensor_22=self.interpolateField(tsv[varName][:,3],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            tensor_23=self.interpolateField(tsv[varName][:,4],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            tensor_33=self.interpolateField(tsv[varName][:,5],grid_x, grid_y, triang, method=interpolationMethod, kind=kind)

            tensor_11=np.flipud(tensor_11)
            tensor_12=np.flipud(tensor_12)
            tensor_13=np.flipud(tensor_13)
            tensor_22=np.flipud(tensor_22)
            tensor_23=np.flipud(tensor_23)
            tensor_33=np.flipud(tensor_33)

            if varName=='UPrime2Mean':
                print 'Adding UPrime2Mean'
                self.data['uu_bar']=tensor_11
                self.data['uv_bar']=tensor_12
                self.data['uw_bar']=tensor_13
                self.data['vv_bar']=tensor_22
                self.data['vw_bar']=tensor_23
                self.data['ww_bar']=tensor_33
                self.data['TKE_bar']=0.5*(self.data['uu_bar']+self.data['vv_bar']+self.data['ww_bar'])
            else:
                print 'Adding symTensor',varName
                self.data[varName+'_ii']=[tensor_11,tensor_12,tensor_13,tensor_22,tensor_23,tensor_33]



    def readVelFromFoamFile(self,
                            varsFile,
                            pointsFile,
                            facesFile,
                            viewAnchor=(0,0,0),
                            xViewBasis=(1,0,0),
                            yViewBasis=(0,1,0),
                            dx=None,
                            dy=None,
                            interpolationMethod='cubic',
                            kind='min_E'):
        '''
        '''

        tsm = TriSurfaceMesh.readFromFoamFile(pointsFile=pointsFile,
                                              facesFile=facesFile,
                                              viewAnchor=viewAnchor,
                                              xViewBasis=xViewBasis,
                                              yViewBasis=yViewBasis)
                                              
        tsv = TriSurfaceVector.readFromFoamFile(varsFile=varsFile,
                                                triSurfaceMesh=tsm,
                                                time=0,
                                                projectedField=False)                  

        points = np.vstack((tsv.x,tsv.y)).T
        
        print 'Creating Grid and Interpolator'
        if dx==None:
            dxlist=[a for a in np.abs(np.diff(points[:,0])) if a>0]
            dx=np.min(dxlist)
        if dy==None:
            dylist=[a for a in np.abs(np.diff(points[:,1])) if a>0]
            dy=np.min(dylist)

        MaxX=np.max(points[:,0])
        MinX=np.min(points[:,0])
        MaxY=np.max(points[:,1])
        MinY=np.min(points[:,1])
        extent=[MinX-dx/2,MaxX+dx/2,MinY-dy/2,MaxY+dy/2]

        cellsX=int((MaxX-MinX)/dx)+1
        cellsY=int((MaxY-MinY)/dy)+1

        grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
        triang = tsv.triangulation

        vx_i=self.interpolateField(tsv.vx,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
        vy_i=self.interpolateField(tsv.vx,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
        vz_i=self.interpolateField(tsv.vx,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)

        self.vx=np.flipud(vx_i)
        self.vy=np.flipud(vy_i)
        self.vz=np.flipud(vz_i)
        
        self.dx=dx
        self.dy=dy
        self.minX=MinX
        self.maxX=MaxX
        self.minY=MinY
        self.maxY=MaxY
        self.createDataDict()
        self.extent=extent

    def readScalarFromFoamFile(self,
                               varsFile,
                               pointsFile,
                               facesFile,
                               viewAnchor=(0,0,0),
                               xViewBasis=(1,0,0),
                               yViewBasis=(0,1,0),
                               dx=None,
                               dy=None,
                               interpolationMethod='cubic',
                               kind='min_E'):
        '''
        '''
        varName=os.path.basename(varsFile)        
        tsm = TriSurfaceMesh.readFromFoamFile(pointsFile=pointsFile,
                                              facesFile=facesFile,
                                              viewAnchor=viewAnchor,
                                              xViewBasis=xViewBasis,
                                              yViewBasis=yViewBasis)
                                              
        tss = TriSurfaceScalar.readFromFoamFile(varsFile=varsFile,
                                                triSurfaceMesh=tsm,
                                                time=0,
                                                projectedField=False)                  

        points = np.vstack((tss.x,tss.y)).T

        #if not hasattr(self,'data'):
            #print 'dict does not exists'
        if not self.data.has_key('dx') or self.data.has_key('dy'):
            print 'keys dx and dy does not exist'
            if dx==None:
                dxlist=[a for a in np.abs(np.diff(points[:,0])) if a>0]
                dx=np.min(dxlist)
            if dy==None:
                dylist=[a for a in np.abs(np.diff(points[:,1])) if a>0]
                dy=np.min(dylist)

            MaxX=np.max(points[:,0])
            MinX=np.min(points[:,0])
            MaxY=np.max(points[:,1])
            MinY=np.min(points[:,1])
            extent=[MinX,MaxX,MinY,MaxY]
            #print MinX,MaxX,MinY,MaxY

            cellsX=int((MaxX-MinX)/dx)
            cellsY=int((MaxY-MinY)/dy)
            #print cellsX,cellsY
            grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
            triang = tss.triangulation
#            scalar_i=doInterp(triang,tss.s,grid_x, grid_y)
            scalar_i=self.interpolateField(tss.s,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            vx_i=np.empty(scalar_i.shape)
            vy_i=np.empty(scalar_i.shape)
            vz_i=np.empty(scalar_i.shape)
            vx_i[:]=np.NAN
            vy_i[:]=np.NAN
            vz_i[:]=np.NAN

            self.vx=np.flipud(vx_i)
            self.vy=np.flipud(vy_i)
            self.vz=np.flipud(vz_i)
            self.extent=extent
            self.minX=MinX
            self.maxX=MaxX
            self.minY=MinY
            self.maxY=MaxY
            self.dx=dx
            self.dy=dy
            self.createDataDict()

            self.data[varName]=np.flipud(scalar_i)
        else:
            print 'dict exists'
            MaxX=self.extent[1]
            MinX=self.extent[0]
            MaxY=self.extent[3]
            MinY=self.extent[2]

            cellsX=int((MaxX-MinX)/self.dx)
            cellsY=int((MaxY-MinY)/self.dy)
            #print cellsX,cellsY
            grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
            triang = tss.triangulation
            scalar_i=self.interpolateField(tss.s,grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            print 'adding scalar',varName
            self.data[varName]=np.flipud(scalar_i)


    def readReStressFromFoamFile(self,
                                 varsFile,
                                 pointsFile,
                                 facesFile,
                                 viewAnchor=(0,0,0),
                                 xViewBasis=(1,0,0),
                                 yViewBasis=(0,1,0),
                                 dx=None,
                                 dy=None,
                                 interpolationMethod='cubic',
                                 kind='min_E'):
        '''
        '''
        tsm = TriSurfaceMesh.readFromFoamFile(pointsFile=pointsFile,
                                              facesFile=facesFile,
                                              viewAnchor=viewAnchor,
                                              xViewBasis=xViewBasis,
                                              yViewBasis=yViewBasis)
                                              
        tsst = TriSurfaceSymmTensor.readFromFoamFile(varsFile=varsFile,
                                                     triSurfaceMesh=tsm,
                                                     time=0,
                                                     projectedField=False)                  

        points = np.vstack((tsst.x,tsst.y)).T

        #if not hasattr(self,'data'):
        if not self.data.has_key('dx') or self.data.has_key('dy'):
            print 'keys dx and dy does not exist'
            if dx==None:
                dxlist=[a for a in np.abs(np.diff(points[:,0])) if a>0]
                dx=np.min(dxlist)
            if dy==None:
                dylist=[a for a in np.abs(np.diff(points[:,1])) if a>0]
                dy=np.min(dylist)

            MaxX=np.max(points[:,0])
            MinX=np.min(points[:,0])
            MaxY=np.max(points[:,1])
            MinY=np.min(points[:,1])
            extent=[MinX,MaxX,MinY,MaxY]
            #print MinX,MaxX,MinY,MaxY



            cellsX=int((MaxX-MinX)/dx)
            cellsY=int((MaxY-MinY)/dy)
            #print cellsX,cellsY
            grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
            triang = tsst.triangulation            
            uu_bar=self.interpolateField(tsst.txx, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            uv_bar=self.interpolateField(tsst.txy, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            uw_bar=self.interpolateField(tsst.tyy, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            vv_bar=self.interpolateField(tsst.tyy, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            vw_bar=self.interpolateField(tsst.tyz, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            ww_bar=self.interpolateField(tsst.tzz, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            vx_i=np.empty(uu_bar.shape)
            vy_i=np.empty(uu_bar.shape)
            vz_i=np.empty(uu_bar.shape)
            vx_i[:]=np.NAN
            vy_i[:]=np.NAN
            vz_i[:]=np.NAN

            self.vx=np.flipud(vx_i)
            self.vy=np.flipud(vy_i)
            self.vz=np.flipud(vz_i)
            self.extent=extent
            self.dx=dx
            self.dy=dy
            self.minX=MinX
            self.maxX=MaxX
            self.minY=MinY
            self.maxY=MaxY
            self.createDataDict()

            print 'adding Tensor'
            self.data['uu_bar']=np.flipud(uu_bar)
            self.data['uv_bar']=np.flipud(uv_bar)
            self.data['uw_bar']=np.flipud(uw_bar)
            self.data['vv_bar']=np.flipud(vv_bar)
            self.data['vw_bar']=np.flipud(vw_bar)
            self.data['ww_bar']=np.flipud(ww_bar)
            self.data['TKE_bar']=0.5*(self.data['uu_bar']+self.data['vv_bar']+self.data['ww_bar'])

        else:
            print 'dict exists'
            MaxX=self.maxX
            MinX=self.minX
            MaxY=self.maxY
            MinY=self.minY

            cellsX=int((MaxX-MinX)/self.dx)
            cellsY=int((MaxY-MinY)/self.dy)
            #print cellsX,cellsY
            grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
            triang = tsst.triangulation
            uu_bar=self.interpolateField(tsst.txx, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            uv_bar=self.interpolateField(tsst.txy, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            uw_bar=self.interpolateField(tsst.tyy, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            vv_bar=self.interpolateField(tsst.tyy, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            vw_bar=self.interpolateField(tsst.tyz, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            ww_bar=self.interpolateField(tsst.tzz, grid_x, grid_y, triang, method=interpolationMethod, kind=kind)
            print 'adding Tensor'
            self.data['uu_bar']=np.flipud(uu_bar)
            self.data['uv_bar']=np.flipud(uv_bar)
            self.data['uw_bar']=np.flipud(uw_bar)
            self.data['vv_bar']=np.flipud(vv_bar)
            self.data['vw_bar']=np.flipud(vw_bar)
            self.data['ww_bar']=np.flipud(ww_bar)
            self.data['TKE_bar']=0.5*(self.data['uu_bar']+self.data['vv_bar']+self.data['ww_bar'])


def getVC7SurfaceList(directory,nr=0,step=1):
    '''
    Get a list of Surfaces read from PIV data
    '''
    filelist=getVC7filelist(directory,nr,step)

    surfaces=[]
    surfaces=[Surface()]*len(filelist)
    #os.chdir(directory)
    if nr==0:
        nr=len(filelist)
    for i in range(0,min(len(filelist),nr)):
        print("reading " + filelist[i])
        surfaces[i]=Surface()
        surfaces[i].readFromVC7(os.path.join(directory,filelist[i]))
    return surfaces

def getVC7filelist(directory,nr=0,step=1):
    '''
    Get a list of filenames of PIV vetor data files
    '''
    filelist=[]
    if os.path.exists(directory):
        for files in os.listdir(directory):
            if files.endswith(".vc7"):
                filelist.append(files)
        filelist.sort()
        filelist=filelist[0::step]
        if nr==0:
            nr=len(filelist)
        filelist=filelist[0:min(len(filelist),nr)]

    return filelist
    
def getIM7filelist(directory,nr=0,step=1):
    '''
    Get a list of filenames of PIV vetor data files
    '''
    filelist=[]
    if os.path.exists(directory):
        for files in os.listdir(directory):
            if files.endswith(".im7"):
                filelist.append(files)
        filelist.sort()
        filelist=filelist[0::step]
        if nr==0:
            nr=len(filelist)
        filelist=filelist[0:min(len(filelist),nr)]

    return filelist
    
def getIM7SurfaceList(directory,nr=0,step=1):
    '''
    Get a list of Surfaces read from PIV data
    '''
    filelist=getIM7filelist(directory,nr,step)

    surfaces=[]
    surfaces=[Surface()]*len(filelist)
    #os.chdir(directory)
    if nr==0:
        nr=len(filelist)
    for i in range(0,min(len(filelist),nr)):
        print("reading " + filelist[i])
        surfaces[i]=Surface()
        surfaces[i].readFromIM7(os.path.join(directory,filelist[i]))
    return surfaces

class rect(object):
    '''
    Defines a rectangle using the cooridnate of two points
    '''
    def __init__(self,x0,x1,y0,y1,name=''):
        self.x0=x0
        self.x1=x1
        self.y0=y0
        self.y1=y1
        self.name=name

    def width(self):
        return np.abs(self.x1-self.x0)

    def height(self):
        return np.abs(self.y1-self.y0)

    def p1(self):
        '''
        returns lower left point
        '''
        xmin=np.min([self.x0,self.x1])
        ymin=np.min([self.y0,self.y1])

        return (xmin,ymin)

class IM7(object):
    def __init__(self,filename):
        '''
        reads PIV image data in tha Davis format, using the 64bit windows DLL
        '''

        dllpath = os.path.dirname(os.path.realpath(__file__))
        self.ReadIMX64 = cdll.LoadLibrary(dllpath+"\ReadIMX64.dll")

        self.myBuffer = BufferType()
        self.attributeLst = AttributeList()
        self.res = self.ReadIMX64.ReadIM7(filename, byref(self.myBuffer), byref(self.attributeLst))
        
    def nx(self):
        return self.myBuffer.nx
        
    def ny(self):
        return self.myBuffer.ny
        
    def nf(self):
        return self.myBuffer.nf
        
    def imageSubType(self):
        BufferFormat_t=dict()
        BufferFormat_t["-2"]='BUFFER_FORMAT_MEMPACKWORD'
        BufferFormat_t["-3"]='BUFFER_FORMAT_FLOAT'
        BufferFormat_t["-4"]='BUFFER_FORMAT_WORD'
        return BufferFormat_t[str(self.myBuffer.image_sub_type)]
        
    def getData(self,frame,v=False):
        s=None
        if self.myBuffer.image_sub_type < 0:
            if v:
                print "Size (ny, nx)"
                print self.myBuffer.ny
                print self.myBuffer.nx
                print self.myBuffer.nf
                print 'type:',self.myBuffer.image_sub_type

            theFrame=frame
            width = self.myBuffer.ny;
            height = self.myBuffer.nx;
            if self.myBuffer.isFloat:
                s=np.empty((height,width), dtype=float)
                s[:] = np.NAN
            else:
                s=np.empty((height,width), dtype=int)
                s[:] = np.NAN
            if(self.myBuffer.isFloat):
                for y in range(0,self.myBuffer.ny):
                    for x in range(0,self.myBuffer.nx): 
                        s[x,y] = (self.myBuffer.floatArray[theFrame*self.myBuffer.nx*self.myBuffer.ny+x*self.myBuffer.ny+y])
            else:
                for y in range(0,self.myBuffer.ny):
                    for x in range(0,self.myBuffer.nx):
                        s[x,y] = (self.myBuffer.wordArray[theFrame*self.myBuffer.nx*self.myBuffer.ny+x*self.myBuffer.ny+y])
        return s
        
    def __del__(self):        
        self.ReadIMX64.DestroyBuffer(self.myBuffer)
    

class SurfaceTimeSeries(object):
    '''
    Transfroms a list of Surfaces into a 3D array
    '''
    def __init__(self):
        self.vx=[]
        self.vy=[]
        self.vz=[]
        self.t=[]
        
        self.dx = float()
        self.dy = float()

        self.minX = float()
        self.maxX = float()
        self.minY = float()
        self.maxY = float()
        self.extent = []

        self.data=dict()
        return
        
    def loadFromSurfaceList(self,slist,frq):
        self.vx=np.array([s.data['Ux'] for s in slist])
        self.vy=np.array([s.data['Uy'] for s in slist])
        self.vz=np.array([s.data['Uz'] for s in slist])
        self.dx=slist[0].dx
        self.dy=slist[0].dy        
        self.minX=slist[0].minX
        self.maxX=slist[0].maxX
        self.minY=slist[0].minY
        self.maxY=slist[0].maxY
        self.extent=slist[0].extent
        
        self.data['frq']=frq
        self.data['dt']=1.0/frq       
        self.t=np.linspace(0,(self.vx.shape[0]-1)/self.data['frq'],self.vx.shape[0])
        self.data['t'] = self.t