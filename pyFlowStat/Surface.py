#===========================================================================#
# load modules
#===========================================================================#
#standard modules
import sys

#scientific modules
import numpy as np
import scipy as sp
import os
import matplotlib.tri as tri
from pyFlowStat.TriSurface import TriSurface
from pyFlowStat.TriSurface import parseFoamFile

# special modules
from ctypes import *

TypeName = ["Image", "2D-PIV-Vector (header, 4x(Vx,Vy))",
            "2D-Vector (Vx,Vy)", "2D-PIV+p.ratio (header, 4x(Vx,Vy), peakratio)",
          "3D-Vector (Vx,Vy,Vz)", "3D-Vector+p.ratio (header, 4x(Vx,Vy), peakratio)"]

WORD=c_ushort

#typedef struct AttributeList
#{
#   char*          name;
#   char*          value;
#   AttributeList* next;
#} AttributeList;
class AttributeList(Structure):
    pass
AttributeList._fields_=[("name",c_char_p),("value",c_char_p),("next",POINTER(AttributeList))]

#   union
#	{
#      float*   floatArray;
#      Word*    wordArray;
#   };
class _bufarray(Union):
    _fields_=[("floatArray",POINTER(c_float)),("wordArray",POINTER(WORD))]



#typedef struct
#{
#   int         isFloat;
#   int         nx,ny,nz,nf;
#   int         totalLines;
#	int			vectorGrid;			// 0 for images
#	int			image_sub_type;	// BufferFormat_t
#   union
#	{
#      float*   floatArray;
#      Word*    wordArray;
#   };
#	BufferScaleType	scaleX;		// x-scale
#	BufferScaleType	scaleY;		// y-scale
#	BufferScaleType	scaleI;		// intensity scale
#	bool*			bMaskArray;			// mask array, NULL if no mask exists
#} BufferType;
class BufferScaleType(Structure):
    _fields_=[("factor",c_float),("offset",c_float),("description",c_char*16),
              ("unit",c_char*16)]

#typedef struct
#{
#	float	factor;
#	float offset;
#	char	description[16];
#	char	unit[16];
#} BufferScaleType;
class BufferType(Structure):
    _anonymous_ = ("bufarray",)
    _fields_=[("isFloat",c_int),("ny",c_int),("nx",c_int),("nz",c_int),("nf",c_int),
             ("totalLines",c_int),("vectorGrid",c_int),("image_sub_type",c_int),
             ("bufarray",_bufarray),("scaleX",BufferScaleType),("scaleY",BufferScaleType),
             ("scaleI",BufferScaleType),("bMaskArray",POINTER(c_bool))]

def getMode(buf,theX_,theY_,width_,frameOffset):
    mode = int(buf.floatArray[theX_ + theY_*width_ + frameOffset])
    if mode<0:
        return -1
    elif mode>4:
        #// interpolated or filled vector
        mode = 4
    mode=mode-1
    return mode

#Read file of type IMG/IMX/VEC, returns error code ImReadError_t
#extern "C" int EXPORT ReadIMX ( const char* theFileName, BufferType* myBuffer, AttributeList** myList );
class Surface(object):
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


    def generateFields(self):
        '''
        Generates additional dictionary entries.
        '''
        Umag = np.zeros(self.data['Ux'].shape)
        Umag2D = np.zeros(self.data['Ux'].shape)
        Umag = np.sqrt(self.data['Ux']**2+self.data['Uy']**2+self.data['Uz']**2)
        Umag2D = np.sqrt(self.data['Ux']**2+self.data['Uy']**2)
        self.data['Umag']=Umag
        self.data['Umag2D']=Umag2D

        dudy,dudx=np.gradient(self.vx,-self.dy/1000,self.dx/1000)
        dvdy,dvdx=np.gradient(self.vy,-self.dy/1000,self.dx/1000)
        vort_z=dvdx-dudy
        self.data['dudy']=dudy
        self.data['dudx']=dudx
        self.data['dvdy']=dvdy
        self.data['dvdx']=dvdx
        self.data['VortZ']=vort_z
        self.data['KE']=0.5*(self.vx**2+self.vy**2+self.vz**2)
        self.data['Div2D']=dudx+dvdy

        #self.data['SwirlingStrength^2']=np.zeros(self.data['Ux'].shape)
        self.data['Q']=np.zeros(self.data['Ux'].shape)
        #self.data['SwirlingStrength^2']=(1.0/(4.0*dudx))**2+(1.0/(4.0*dvdy))**2-0.5*dudx*dvdy+dvdx*dudy
        self.data['Q']=0.5*(-2.0*dudy*dvdx-dudx**2-dvdy**2)
#        tensorS= np.empty(self.data['Ux'].shape)
#        tensorW= np.empty(self.data['Ux'].shape)
#        tensorS= 0.5*[[dudx+dudx,dudy+dvdx],[dvdx+dudy,dvdy+dvdy]]
#        tensor2= 0.5*[[0.0,dudy-dvdx],[dvdx-dudy,0.0]]
        self.data['lambda2'] = self.getLambda2(dudx,dudy,dvdx,dvdy)

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

    def addReynoldsDecomposition(self,MeanFlowSurface):
        '''
        Generate fluctuations by subtracting the mean flow (surface of same size)
        Adds fluctuation fields ux,uy,uz and correleations uu,vv,ww,uv,uw and TKE
        '''
        self.data['ux']=self.data['Ux']-MeanFlowSurface.data['Ux']
        self.data['uy']=self.data['Uy']-MeanFlowSurface.data['Uy']
        self.data['uz']=self.data['Uz']-MeanFlowSurface.data['Uz']
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
                        self.vx[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+1)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vy[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*2+2)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    else:
                        pass
        if tmpBuffer.image_sub_type == 4:
            for theY in range(0,width):
                for theX in range(0,height):
                    self.vx[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    self.vy[theY,theX] = -1*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                    self.vz[theY,theX] = (tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*2]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
        if tmpBuffer.image_sub_type == 5:
            for theY in range(0,width):
                for theX in range(0,height):
                    mode = getMode(tmpBuffer,theX,theY,height,frameOffset)
                    if mode >= 0:
                        self.vx[theY,theX]=(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+1)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
                        self.vy[theY,theX]=-1*(tmpBuffer.floatArray[theX + theY*height + frameOffset + componentOffset*(mode*3+2)]*tmpBuffer.scaleI.factor+tmpBuffer.scaleI.offset)
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

    def interpolateField(self,values,grid_x,grid_y,triangulation,method='cubic'):
        '''
        helper function
        methode=linear,cubic (default)
        '''
        if method=='cubic':
            itp=tri.CubicTriInterpolator(triangulation,values)
        elif method=='linear':
            itp=tri.LinearTriInterpolator(triangulation,values)
        else:
            itp=tri.CubicTriInterpolator(triangulation,values)
        zi_ma = itp(grid_x, grid_y)
        zi=zi_ma.filled(np.nan)

        return zi

    def readFromFoamFile(self,pointsFile,facesFile,velFile,scalarFileList=[],symTensorFileList=[],viewAnchor=(0,0,0),xViewBasis=(1,0,0),yViewBasis=(0,1,0),dx=None,dy=None,interpolationMethod='cubic'):

        print 'Reading Velocity'

        s=TriSurface()
        #s.storeMesh=False
        s.readFromFoamFile(varsFile=velFile,pointsFile=pointsFile,facesFile=facesFile,viewAnchor=viewAnchor,xViewBasis=xViewBasis,yViewBasis=yViewBasis)

        points=s.xys
        faces=s.faces
        #points=parseFoamFile(pointsFile)
        #faces = parseFoamFile(facesFile)[:,1:4]

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
        #extent=[MinX,MaxX,MinY,MaxY]

        cellsX=int((MaxX-MinX)/dx)+1
        cellsY=int((MaxY-MinY)/dy)+1

        grid_y, grid_x = np.mgrid[MinY:MaxY:np.complex(0,cellsY),MinX:MaxX:np.complex(0,cellsX)]
        triang = tri.Triangulation(points[:,0], points[:,1], faces)

        print 'Interpolating Velocity'
        vx_i=self.interpolateField(s.vars[:,0],grid_x, grid_y, triang,method=interpolationMethod)
        vy_i=self.interpolateField(s.vars[:,1],grid_x, grid_y, triang,method=interpolationMethod)
        vz_i=self.interpolateField(s.vars[:,2],grid_x, grid_y, triang,method=interpolationMethod)
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
            s.vars=parseFoamFile(scalarFile)
            scalar_i=self.interpolateField(s.vars[:,0],grid_x, grid_y, triang,method=interpolationMethod)
            self.data[varName]=np.flipud(scalar_i)

        for symTensorFile in symTensorFileList:
            varName=os.path.basename(symTensorFile)
            print 'Reading Tenstor',varName
            s.vars=parseFoamFile(symTensorFile)
            tensor_11=self.interpolateField(s.vars[:,0],grid_x, grid_y, triang,method=interpolationMethod)
            tensor_12=self.interpolateField(s.vars[:,1],grid_x, grid_y, triang,method=interpolationMethod)
            tensor_13=self.interpolateField(s.vars[:,2],grid_x, grid_y, triang,method=interpolationMethod)
            tensor_22=self.interpolateField(s.vars[:,3],grid_x, grid_y, triang,method=interpolationMethod)
            tensor_23=self.interpolateField(s.vars[:,4],grid_x, grid_y, triang,method=interpolationMethod)
            tensor_33=self.interpolateField(s.vars[:,5],grid_x, grid_y, triang,method=interpolationMethod)

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



    def readVelFromFoamFile(self,varsFile,pointsFile,facesFile,viewAnchor=(0,0,0),xViewBasis=(1,0,0),yViewBasis=(0,1,0),dx=None,dy=None):


        s=TriSurface()
        s.readFromFoamFile(varsFile,pointsFile,facesFile,viewAnchor,xViewBasis,yViewBasis)

        points=s.xys
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
        triang = tri.Triangulation(points[:,0], points[:,1], s.faces)

        vx_i=self.interpolateField(s.vars[:,0],grid_x, grid_y,triang)
        vy_i=self.interpolateField(s.vars[:,1],grid_x, grid_y,triang)
        vz_i=self.interpolateField(s.vars[:,2],grid_x, grid_y,triang)

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

    def readScalarFromFoamFile(self,varsFile,pointsFile,facesFile,viewAnchor=(0,0,0),xViewBasis=(1,0,0),yViewBasis=(0,1,0),dx=None,dy=None):


        varName=os.path.basename(varsFile)
        s=TriSurface()
        s.readFromFoamFile(varsFile,pointsFile,facesFile,viewAnchor,xViewBasis,yViewBasis)
        points=s.xys

        #if not hasattr(self,'data'):
            #print 'dict does not exists'
        if not data.has_key('dx') or data.has_key('dy'):
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
            triang = tri.Triangulation(points[:,0], points[:,1], s.faces)
            scalar_i=doInterp(triang,s.vars[:,0],grid_x, grid_y)
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
            triang = tri.Triangulation(points[:,0], points[:,1], s.faces)
            scalar_i=doInterp(triang,s.vars[:,0],grid_x, grid_y)
            print 'adding scalar',varName
            self.data[varName]=np.flipud(scalar_i)


    def readReStressFromFoamFile(self,varsFile,pointsFile,facesFile,viewAnchor=(0,0,0),xViewBasis=(1,0,0),yViewBasis=(0,1,0),dx=None,dy=None):


        s=TriSurface()
        s.readFromFoamFile(varsFile,pointsFile,facesFile,viewAnchor,xViewBasis,yViewBasis)
        points=s.xys

        #if not hasattr(self,'data'):
        if not data.has_key('dx') or data.has_key('dy'):
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
            triang = tri.Triangulation(points[:,0], points[:,1], s.faces)
            uu_bar=doInterp(triang,s.vars[:,0],grid_x, grid_y)
            uv_bar=doInterp(triang,s.vars[:,1],grid_x, grid_y)
            uw_bar=doInterp(triang,s.vars[:,2],grid_x, grid_y)
            vv_bar=doInterp(triang,s.vars[:,3],grid_x, grid_y)
            vw_bar=doInterp(triang,s.vars[:,4],grid_x, grid_y)
            ww_bar=doInterp(triang,s.vars[:,5],grid_x, grid_y)
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
            triang = tri.Triangulation(points[:,0], points[:,1], s.faces)
            uu_bar=doInterp(triang,s.vars[:,0],grid_x, grid_y)
            uv_bar=doInterp(triang,s.vars[:,1],grid_x, grid_y)
            uw_bar=doInterp(triang,s.vars[:,2],grid_x, grid_y)
            vv_bar=doInterp(triang,s.vars[:,3],grid_x, grid_y)
            vw_bar=doInterp(triang,s.vars[:,4],grid_x, grid_y)
            ww_bar=doInterp(triang,s.vars[:,5],grid_x, grid_y)
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

class rect(object):
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
        xmin=np.min([self.x0,self.x1])
        ymin=np.min([self.y0,self.y1])

        return (xmin,ymin)

class SurfaceTimeSeries(object):
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
        self.data['t']=np.linspace(0,(self.vx.shape[0]-1)/self.data['frq'],self.vx.shape[0])