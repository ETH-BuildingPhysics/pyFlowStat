import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyFlowStat.Surface as Surface

def PlotField(ax,surface,field,vmin,vmax,offset=[0,0],interpolation='nearest',modifier=None,**kwargs):
    extent=[0,0,0,0]
    extent[0]=surface.extent[0]-offset[0]
    extent[1]=surface.extent[1]-offset[0]
    extent[2]=surface.extent[2]-offset[1]
    extent[3]=surface.extent[3]-offset[1]
    if not modifier:
        im=ax.imshow(surface.data[field],vmin=vmin,vmax=vmax,interpolation=interpolation,extent=extent,**kwargs)
    else:
        im=ax.imshow(modifier(surface.data[field]),vmin=vmin,vmax=vmax,interpolation=interpolation,extent=extent,**kwargs)
    return im
    
def PlotContour(ax,surface,field,vmin,vmax,offset=[0,0], contourlevels=21, contourlabels=11,**kwargs):

    X,Y = surface.getMeshgrid(offset=offset)
    contour_levels = np.linspace(vmin, vmax, contourlevels)
    contour_levels_label = np.linspace(vmin, vmax, contourlabels)
    print field
    cts=ax.contourf(X,Y,surface.data[field],contour_levels,alpha=.75,**kwargs)
#    ax.colorbar()

    C = ax.contour(X,Y,surface.data[field], contour_levels_label, colors='black', linewidth=.5)
    ax.clabel(C, inline=1, fontsize=10)
    ax.set_aspect('equal')
    return cts

#    xfill=[-100,0,0,100,100,200,200,-100]
#    yfill=[0,0,-100,-100,0,0,-120,-120]
#    fill(xfill,yfill,'k')
#    
#    xlim([-30,130])
#    ylim([-105,130])
def PlotStreamLine(ax,surface,vmin,vmax,density=10,offset=[0,0],**kwargs):
    X,Y = surface.getMeshgrid(offset=offset)

    cnorm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)

    return ax.streamplot(X,Y,surface.data['Ux'],surface.data['Uy'],density=density,norm=cnorm,color='k',**kwargs)
#    return ax.streamplot(X,Y,surface.data['Ux'],surface.data['Uy'],density=density,norm=cnorm,color=u)

def PlotColoredStreamLine(ax,surface,vmin,vmax,density=10,offset=[0,0]):
    X,Y = surface.getMeshgrid(offset=offset)
    
    u=np.nan_to_num(np.sqrt(surface.data['Ux']**2+surface.data['Uy']**2))

    cnorm=mpl.colors.Normalize(vmin=0,vmax=np.max(u))


    return ax.streamplot(X,Y,surface.data['Ux'],surface.data['Uy'],density=density,norm=cnorm,color=u)
#    return ax.streamplot(X,Y,surface.data['Ux'],surface.data['Uy'],density=density,norm=cnorm,color=u)

def PlotVelocityVectors(ax,surface,scale=1,offset=[0,0],spacing=1,**kwargs):
 
    if 'width' not in kwargs:
        kwargs['width']=0.1
        
    X,Y = surface.getMeshgrid(offset=offset)
    return plt.quiver(X[::spacing,::spacing],Y[::spacing,::spacing],surface.data['Ux'][::spacing,::spacing],surface.data['Uy'][::spacing,::spacing],scale=scale,angles='uv',units='xy',**kwargs)