#! /usr/bin/env python
from pylab import *
from pyFlowStat.PointProbe import PointProbe

a = PointProbe()

path='19.82/U'
lst=PointProbe.readOfRuntime([0.5,-0.3,0.1],path)
print lst
