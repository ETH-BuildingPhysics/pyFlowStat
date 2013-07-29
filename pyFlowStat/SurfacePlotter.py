import matplotlib.pyplot as plt
import numpy as np
import pyFlowStat.Surface as Surface

def PlotField(ax,surface,field,vmin,vmax):
    im=plt.imshow(surface.data[field],vmin=vmin,vmax=vmax,interpolation='nearest',extent=surface.extent)
    return im