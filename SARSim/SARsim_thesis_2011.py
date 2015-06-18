#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      vfomin
#
# Created:     18/05/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# SOURCE :
# Shaharyar Khwaja. Fast Raw Data Generation of Realistic Environments for a SAR System
# Simulator. Signal and Image processing. Universite Rennes 1, 2008. English.
# <tel-00371992>
# https://tel.archives-ouvertes.fr/tel-00371992


# Python
import sys

# Numpy
import numpy as np
import matplotlib.pyplot as plt

# Opencv
import cv2

# Project
import Tools.GlobalConsts as GC
import Sensors.GlobalVariablesTSX as TSX

# GLOBAL VARIABLES :
Y0 = TSX.R0 * np.sin(TSX.Theta) # in meters
T_p = TSX.Tp
Y_d = Y0*0.001 # in meters, target is illuminated for a total azimuth distance



def d_a_0(y):
    """
    Radar-target distance of a stationary target
    d_a_0(y) = sqrt( R0**2 + (y-Y0)**2 )
    """
    return np.sqrt( TSX.R0**2 + (y - Y0)**2 )

def chirp(t):
    """
    SAR emitting pulse
    """
    return GC.rect(t / T_p) * np.exp(1j*np.pi*( 2.0* TSX.Freq * t + TSX.K * t**2 ) )

def sarResp(t,y,d_a,E_y):
    """
    The demodulated received data for a single point at slant-range and azimuth position of
    r_a and y_a, respectively is given as

    s_a(t,y) = E_t(t - 2*d_a(y)/c) * E_y(y - Y0) * exp( -4*pi*1j * f_c * d_a(y) / c + 1j * pi * K_r * (t - 2*d_a(y) / c)**2 )

    E_t(t) = rect(t/T_p)
    T_p - chirp duration time
    E_y - antenna beam pattern

    """
    tt = t - 2.0*d_a(y) / GC.C
    return GC.rect(tt / T_p) * E_y(y - Y0) * np.exp(-4.0*1j*np.pi * TSX.Freq * d_a(y) / GC.C + 1j*np.pi * TSX.K * (tt)**2 )




if __name__ == '__main__':

    yArray=Y0 + np.arange(-Y_d,Y_d,Y_d*0.001)
    tArray=np.arange(-T_p*0.6,T_p*0.8,T_p*0.001)


    # Chirp signal
##    emitSig=chirp(tArray)
##    plt.subplot(121)
##    plt.title("Chirp module")
##    plt.plot(tArray,np.absolute(emitSig))
##    plt.subplot(122)
##    plt.title("Chirp phase")
##    plt.plot(tArray,np.angle(emitSig))
##    plt.show()

    # Radar-Target distance plot
##    dArray = d_a_0(yArray)
##    plt.plot(yArray,dArray,'-b')
##    plt.show()


    # SAR raw data
    respSig = np.zeros((len(tArray), len(yArray)), dtype=complex)
    E_y = lambda y : GC.rect(y/Y_d)
    for index, t in enumerate(tArray):
        respSig[index,:] = sarResp(t, yArray[:],d_a_0,E_y)

    plt.subplot(121)
    plt.imshow(np.absolute(respSig))
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(np.angle(respSig))
    plt.colorbar()
    plt.show()




















