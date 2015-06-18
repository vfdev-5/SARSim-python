#-------------------------------------------------------------------------------
# Name:        CSK Satellite global variables
# Purpose:      Defines global CSK satellite (HImage) & sensor parameters
#               File used : CSKS1_SCS_U_HI_24_VV_RD_SF_20130509205142_20130509205150.h5.aux.xml
# Author:      vfomin
#
# Created:     18/05/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Numpy
import numpy as np

# Project
import Tools.GlobalConsts as GC




# GLOBAL VARIABLES :
Freq = 9600000000 # in Hz, "Radar_Frequency"
Lambda = GC.C/Freq # in meters,  Carrier wavelength
H = 625222.52152799 # height of the satellite in m, "Satellite_Height"
Vsat = 7632.8601880878696  # in m/s, mean of values from "ECEF_Satellite_Velocity"
K = -1631469726562.5 # Hz/s, "S01_Range_Chirp_Rate"
Tp = 4e-005 # in seconds, "S01_Range_Chirp_Length"
La = 5.6 # in meters, "Antenna_Length"
PRF = 3063.72549019608 # Hz, "S01_PRF"

Name = "CSK HImage"


def _computeMeanVelocity(vString, delimiter=" "):
    """
    Method to compute mean velocity using data from xml file
    vString = "-325.6660555 -1830.2709704 -7403.3921748 -386.3346157 -1774.236534 -7414.0778361 ..."
    """
    vArray = np.array(vString.split(delimiter)).astype(np.float)
    vArray = vArray.reshape((len(vArray)/3, 3))

    vel = np.sqrt(np.sum(vArray * vArray, 1))
    return vel.mean()



def _readAnalogSignalReconstructionLevels(vString, delimiter=","):
    """
    Method to read string value into array and replace Nan with zeros
    vString = "NaN,NaN,NaN,NaN,NaN,-1226.9649658203125,-766.2490234375 ..."
    """
    vArray = np.array(vString.split(delimiter)).astype(np.float)
    vArray[np.isnan(vArray)] = 0.0
    return vArray



def _dataWriter(dataRe, dataIm, lut, outputFilename):
    """
    Write I-Q (8 bits) data into a npy (np.complex array) file
    """
    data = np.zeros(dataRe.shape,dtype=np.complex)
    for i in range(data.shape[0]):
        re = lut[dataRe[i,:]]
        im = lut[dataIm[i,:]]
        data[i,:] = re + 1j * im

    np.save(outputFilename, data)
