#-------------------------------------------------------------------------------
# Name:        Velocity Estimation R&D chain
# Purpose:
#
# Author:      vfomin
#
# Created:     16/06/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Numpy, Matplotplib
import numpy as np
import matplotlib.pyplot as plt

# Project
import Tools.GlobalConsts as GC
import Sensors.GlobalVariablesCSK_HI as Sensor
import SARSim.SARSim_proto2 as SARSim

xc = 0.0
yc = Sensor.H * np.tan(55.0 * np.pi / 180.0)
s = SARSim.SARSim(Sensor,verbose=False)
# define SAR acquisition parameters:
params = [[-4200, 4200], [54.995, 55.005], [xc, yc]]
s.setParams(params)
# set targets
# a) moving targets
targets = [
            [xc-122.0, yc-100.0],
            [xc-100.0, yc-80.0, 50.0, 0.0],
            [xc-80.0, yc-60.0, 0.0, 50.0],
            [xc-60.0, yc-40.0, 25.0, 27.0],
            [xc-40.0, yc-20.0, -25.0, 27.0],
            [xc-20.0, yc, -25.0, -27.0],
            [xc, yc+20.0, 25.0, -27.0],
            [xc+20.0, yc+40.0, -45.0, 0.0],
            [xc+40.0, yc+60.0, 0.0, -45.0]
        ]
targetXMin = -350
targetXMax = 200
targetYMin = yc-150
targetYMax = yc+100

s.setTargets(targets)

print "Compute raw data ..."
Srx, xArray, tArray = s.computeRawData()

# data compression
print "Compress raw data ..."
S_final, xArray, tArray = s.compressRawData()
##    SARSim.showResponse(S_final, [tArray[0],tArray[-1],xArray[0],xArray[-1]], 'Compressed image', False)
# linear interpolation from Range time to OY coordinates
S_final, yArray = SARSim.warpFromRangeTimeToRegularOYGrid(S_final, tArray, Sensor.H)


xIndx = np.argwhere((xArray < targetXMax) & (xArray > targetXMin)).ravel()
yIndx = np.argwhere((yArray < targetYMax) & (yArray > targetYMin)).ravel()
win = [xIndx[0], xIndx[-1], yIndx[0], yIndx[-1]]
print "Window :", win

##SARSim.showGrayResponse(S_final, [yArray[0]-yc, yArray[-1]-yc,xArray[0],xArray[-1]], 'Compressed image - Ground range')
SARSim.showResponse(S_final, None, 'Compressed image - Ground range')






