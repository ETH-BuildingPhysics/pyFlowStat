#! /usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyFlowStat.PointProbe import PointProbe

path=r'Data/10/U'
pt = PointProbe()
pt.readFromOpenFoam([0.3,-0.3,0.05],path)

plt.figure()
ax=plt.plot(pt.data['t'],pt.data['U'])
plt.title('Time resolved velocity')
plt.grid(True)
plt.xlabel('Time $[s]$')
plt.ylabel('Velocity $[m/s]$')
plt.legend(['$U_x$','$U_y$','$U_z$'])
plt.show()