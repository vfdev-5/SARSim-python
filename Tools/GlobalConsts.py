#-------------------------------------------------------------------------------
# Name:        global variables
# Purpose:    some universal constants and math functions
#
# Author:      vfomin
#
# Created:     18/05/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# DEFINE GLOBAL CONSTS :
C = 299792458.0 # in m/s

# Numpy
import numpy as np


def hs(x):
    """
    Heavyside function :
        hv(x) = {1, if x>=0; 0, otherwise}
    """
    return 0.5 * (np.sign(x) + 1.0)

def ihs(x):
    """
    Inverse Heavyside function:
        ihv(x) = {0, if x>=0; 1, otherwise}
    """

    return 0.5 * (1.0 - np.sign(x))

def rect(x):
    """
    Rectangle function:
        rect(x) = {1, if |x|<= 0.5; 0, otherwise}
    """
    return hs(x + 0.5) * ihs(x - 0.5)


# -------------------------FFT / IFFT-------------------------------------------

def fft(s,n=None):
    """
    Improved fft that gives the correct phase
    Input signal should be such that the element s[middle] corresponds to the time = 0
    """
    return np.fft.fftshift(np.fft.fft(np.fft.fftshift(s), n))

def ifft(s, n=None):
    """
    Improved ifft that gives the correct phase
    Input signal should be such that the element s[middle] corresponds to the frequency = 0
    """
    return np.fft.ifftshift(np.fft.ifft(np.fft.ifftshift(s), n))


def freq(length, bandwidth):
    return np.fft.fftshift(np.fft.fftfreq( length, 1.0 / bandwidth))


def fft2(img):
    """
    Improved 2D fft
    """
    out = np.zeros(img.shape, dtype=complex)
    for i in range(out.shape[1]):
        # get range fixed column
        out[:,i] = fft(img[:,i])
    for j in range(out.shape[0]):
        out[j,:] = fft(out[j,:])
    return out


def ifft2(img):
    """
    Improved 2D ifft
    """
    out = np.zeros(img.shape, dtype=complex)
    for i in range(out.shape[1]):
        out[:,i] = ifft(img[:,i])
    for j in range(out.shape[0]):
        out[j,:] = ifft(out[j,:])
    return out


# -----------------------------------------------------------------------------

def chirpSignal(t,tau_p,K_r):
    """
    Create a chirp signal :
        S_{tx}(t) = rect(t/tau_p) * exp(1j*pi*K_r*t^2)
    """
    return rect(t / tau_p) * np.exp(1j*np.pi*K_r * t**2)


def pulseSignal(t, tau, T, delayTime=0):
    """
    Create a pulse signal

             _______
      ______|       |_______

    -T/2 -tau/2+d tau/2+d T/2

    time = -T*0.5 + step*index

    """
    sArray = GC.rect((t-delayTime)/tau)
    return sArray, tArray

def pulseSignal2(tau, T, N, delayTime=0):
    """
    Create a pulse signal

             _______
      ______|       |_______

    -T/2 -tau/2+d tau/2+d T/2

    time = -T*0.5 + step*index

    """
    step = T*1.0/N
    tArray = np.arange(-T*0.5,T*0.5,step)
    sArray = GC.rect((tArray-delayTime)/tau)
    return sArray, tArray
