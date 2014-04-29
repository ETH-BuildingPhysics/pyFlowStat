#! /usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyFlowStat.PointProbe import PointProbe

path=r'Data/10/U'
pt = PointProbe()
pt.readFromOpenFoam([0.3,-0.3,0.05],path)

plt.figure('Time resolved velocity')
ax=plt.plot(pt.data['t'],pt.data['U'])
plt.title('Time resolved velocity')
plt.grid(True)
plt.xlabel('Time $[s]$')
plt.ylabel('Velocity $[m/s]$')
plt.legend(['$U_x$','$U_y$','$U_z$'])

# In the "pt" object of class "PointProbe", the data are stored in python ditionnary.
# Therefore, all the methods of a python dictonnary can be used. All single data
# (like U, Umag or the time) are stored in a numpy.array. The following lines
# show some examples.

# get a list of all the keys
print('data in pt:')
print(pt.data.keys())

#print Ux like this
print(pt.data['U'][:,0])
# or like this
print(pt['U'][:,0])

#print Uy
print(pt['U'][:,1])

# print the time
print(pt['t'])

# generate statistics (like cross corelation, fft, Reynolds stress tensor):
pt.generateStatistics()

# get a list of all the keys (lots of new keys)
print('new data in pt:')
print(pt.data.keys())

# plot velocity U and the mean
plt.figure('Velocities and means')
plt.plot(pt['t'],pt['U'][:,0],c='blue',label='$U_x$')
plt.plot(pt['t'],pt['U'][:,1],c='red',label='$U_y$')
plt.plot(pt['t'],pt['U'][:,2],c='green',label='$U_z$')
plt.plot(pt['t'],pt['Uoo'][:,0],c='blue',lw=2,label=r'$\langle U_x \rangle$')
plt.plot(pt['t'],pt['Uoo'][:,1],c='red',lw=2,label=r'$\langle U_y \rangle$')
plt.plot(pt['t'],pt['Uoo'][:,2],c='green',lw=2,label=r'$\langle U_z \rangle$')
plt.title('Velocities and means')
plt.grid(True)
plt.xlabel('Time $[s]$')
plt.ylabel('Velocity $[m/s]$')
plt.legend()

# plot auto-correlation
plt.figure('auto-correlation')
plt.plot(pt['taur11'],pt['r11'],label='$r_{11}$')
plt.plot(pt['taur22'],pt['r22'],label='$r_{22}$')
plt.plot(pt['taur33'],pt['r33'],label='$r_{33}$')
plt.xlim([0,1000])
plt.xlabel('step shift')
plt.ylabel('$r_{ii}$')
plt.title('auto-correlation')
plt.grid(True)
plt.legend()

# As the data are store in a dictonnary in the object pt, a new entry can be
# added easily

# My new entry is 2*Umean
newData = 2*pt['Uoo']
pt['2timeUoo'] = newData
print('list of data in pt:')
print(pt.data.keys())
print('Here is my new entry "2timeUoo":')
print(pt['2timeUoo'])








plt.show()

