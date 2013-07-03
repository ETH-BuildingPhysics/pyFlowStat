#! /usr/bin/env python

from pylab import *
from pyFlowStat.PointProbe import PointProbe

probeFile1 = '10/U'
probeFile2 = '19.82/U'
point = [0.5,-0.3,0.1]   #must be included in the file defined in "path"

#Create object pt of class PointProbe
#pt = PointProbe(point,path)
pt = PointProbe(point,probeFile2)

#access to data:
# velocity
print(pt.data['U'])
# time
print(pt.data['t'])
# time serie frequency
print(pt.data['frq'])
# get a list of all available datas
print(pt.data.keys())

# generate aditional statistics
pt.generateStatistics()

# auto-correlation coefficients of fluctuation in direction 1 (u1) and related lag
print(pt.data['r11'])
print(pt.data['taur11'])

#plot auto-correlation coefficient
figure()
plot(pt.data['taur11'],pt.data['r11'],color='red',label='auto-corr coeff u1')
title('My title')
xlabel('my x label')
ylabel('my y label')
grid(True)
legend()
show()

#plot velocity field in the three directions and their average
figure()
plot(pt.data['t'],pt.data['U'][:,0],color='blue',label='U1')
plot(pt.data['t'],pt.data['U'][:,1],color='red',label='U2')
plot(pt.data['t'],pt.data['U'][:,2],color='green',label='U3')
plot(pt.data['t'],pt.data['Uoo'][:,0],c='darkblue',lw=2,label='Uave1')
plot(pt.data['t'],pt.data['Uoo'][:,1],c='darkred',lw=2,label='Uave2')
plot(pt.data['t'],pt.data['Uoo'][:,2],c='darkgreen',lw=2,label='Uave3')
grid(True)
legend()
show()




#Create object pt of class PointProbe
pt = PointProbe()
pt.readFromOpenFoam(point,probeFile1)

# generate aditional statistics
pt.generateStatistics()

#plot velocity field in the three directions and their average
figure()
plot(pt.data['t'],pt.data['U'][:,0],color='blue',label='U1')
plot(pt.data['t'],pt.data['U'][:,1],color='red',label='U2')
plot(pt.data['t'],pt.data['U'][:,2],color='green',label='U3')
plot(pt.data['t'],pt.data['Uoo'][:,0],c='darkblue',lw=2,label='Uave1')
plot(pt.data['t'],pt.data['Uoo'][:,1],c='darkred',lw=2,label='Uave2')
plot(pt.data['t'],pt.data['Uoo'][:,2],c='darkgreen',lw=2,label='Uave3')
xlabel('time [s]')
ylabel('velocity [m/s]')
title('U and Ubar at location '+str(point))
grid(True)
legend()
show()
