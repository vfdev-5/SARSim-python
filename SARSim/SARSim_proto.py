#-------------------------------------------------------------------------------
# Name:        SARSim
# Purpose:     SAR Image simulator
#
# Author:      vfomin
#
# Created:     19/05/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

"""
SAR Imaging tutorial

1) Emitted signal :

    S_{tx}(t) = rect(t/tau_p) * exp(2*pi*1j* f_c * t + 1j*pi*K_r*t^2),

    where
        tau_p - pulse duration (s)
        f_c - carrier frequency (Hz)
        K_r - chirp rate (Hz/s)

        t - is range-time (quick/fast)

    rect(x) = {1, if |x|<= 0.5; 0, otherwise }


    Pulse emission representation

       Emission
        |   Reception
        |       |
    ---||||----------||||----------||||----------> time
       <-->          <------------->
       tau_p              1/PRF


    PRF - pulse repetition frequency


    Frequency time evolution of the signal :
        f(t) = f_c + K_r*t
    and bandwidth B = K_r*tau_p


2) Received signal :

    Let a raw received signal be S_r(t, x) at x (azimuth) position
    of the satellite, from a point stationary target
    S_r(t, x) = rect( (t-2*d(x)/c) / tau_p) * Wa(x) * sigma
                    * exp( -4*pi*1j* f_c * d(x)/c + 1j*pi*K_r*(t-2*d(x)/c)^2 )

    where
        d(x) - distance until the target when the satellite is at x (azimuth) position
            d(x)^2 = H^2 + y_c^2 + (x-x_c)^2
            H - height of the satellite orbit
            x_c - target zone OX (azimuth) position
            y_c - target zone OY (range) position
        Wa(x) - antenna beam pattern, can be approximated
            a) ideal : Wa(x) = rect( (x - x_c)/deltaX )
                deltaX - azimuth flight distance covered by satellite during elimination and reception
                deltaX = deltaT * V_{sat}, V_{sat} - satellite velocity
            b) approximation of a uniform aperture (see SARAS: A synthetic Aperture Radar (SAR) Raw Signal Simulator, G.Franceschetti)
                Wa(x) = sinc^2 ( (x - x_c)/deltaX )

3) Range compression :

    Matching filter is applied on received signal S_r(t)
    S_{out}(t) = Integral_{-inf}^{inf} S_r(u) * H(t-u) du,
    where H(t) - duplicated complex-conjugated signal of the original signal
    H(t) = rect(t/tau_p) * exp(1j*pi*K_r*t^2)

    S_{out}(t,x) = iFFT_t( FFT_t(S_r(t,x)) * H_fft)
    where H_fft is fft of emitted signal


4) Range Cell Migration Correction :



5) Azimuth compression
    Matching filter is applied on received signal S_r(t,x)
    where H(x) = rect((x-xc)/deltaX) * exp(1j*pi*K_a*(x-xc)^2)

"""


# Numpy, Matplotplib
import numpy as np
import matplotlib.pyplot as plt


# Project
import Tools.GlobalConsts as GC
##import GlobalVariablesTSX as Sat
import Sensors.GlobalVariablesAerial as Sat


def emittedSignal(t,tau_p, f_c, K_r):
    """
    Emitted signal :
        S_{tx}(t) = rect(t/tau_p) * exp(2*pi*1j* f_c * t + 1j*pi*K_r*t^2)
    """
    return GC.rect(t / tau_p) * np.exp(1j*np.pi*( 2.0* f_c * t + K_r * t**2 ) )

def chirpSignal(t,tau_p,K_r):
    """
    Emitted signal :
        S_{tx}(t) = rect(t/tau_p) * exp(1j*pi*K_r*t^2)
    """
    return GC.rect(t / tau_p) * np.exp(1j*np.pi*K_r * t**2)


def received2DSignal(t, x, dist, tau_p, f_c, K_r, Wa, sigma=1,) :
    """
    SAR received signal at range time t and azimuth position x
    S_{out}(t,x) = rect( (t-2*d(x)/c) / tau_p) * Wa(x) * sigma
                    * exp( -4*pi*1j* f_c * d(x)/c + 1j*pi*K_r*(t-2*d(x)/c)^2 )


    """
    tt = t - 2.0*dist(x) / GC.C
    phase = -4.0*f_c*dist(x) / GC.C + K_r*(tt**2)
    return GC.rect( tt/tau_p ) * Wa(x) * sigma * np.exp(1j*np.pi*phase)


def distSST(x,H,xc,yc):
    """
    Distance function for a stationary ground target zone
    d(x)^2 = H^2 + yc^2 + (x-xc)^2
    xc, yc is the center of the aimed zone
    """
    return np.sqrt(H**2 + yc**2 + (x-xc)**2)


def deltaR(x,xc,dist):
    """
    delta R for a stationary ground target zone
    deltaR = R(x) - Rc,
    where R(x) is the distance when sensor is at azimuth point x
    Rc=(xc,yc) is the smallest distance to the center of the zone

    - x is the numpy array of azimuth positions
    - xc is the azimuth position of the targeted zone
    - dist is the distance function

    """
    return dist(x) - dist(xc)



def idealWa(x,x_c,deltaX):
    """
    Antenna Beam pattern
    ideal : Wa(x) = rect( (x - x_c)/deltaX ) where
        deltaX - azimuth flight distance covered by satellite during elimination and reception
        deltaX = deltaT * V_{sat}, V_{sat} - satellite velocity
    """
    return GC.rect( (x-x_c) / deltaX )

def approxWa(x,x_c,deltaX):
    """
    Antenna Beam pattern
    Approximation of a uniform aperture (see SARAS: A synthetic Aperture Radar (SAR) Raw Signal Simulator, G.Franceschetti)
    Wa(x) = sinc^2 ( (x - x_c)/deltaX ) where

        deltaX - azimuth flight distance covered by satellite during elimination and reception
        deltaX = deltaT * V_{sat}, V_{sat} - satellite velocity
    """
    wa = np.sinc( (x-x_c) / deltaX )
    return wa**2


def azimuthCompression(Srx, h):
    """
    Method to apply matched filter h on azimuth dimension
    """
    # pad received and emitted signals on common time stemp
    # to obtain common fft
    n1 = len(h)
    n2 = Srx.shape[0]

    Nsa = n1 + n2 - 1
##    Nsa2 = int(np.power(2.0, np.ceil(np.log2( Nsa ))))
    Nsa2 = Nsa if Nsa % 2 == 0 else Nsa+1

    SrxMF=np.zeros((Nsa2, Srx.shape[1]),dtype=complex)
    hMF = np.concatenate((np.zeros(np.floor((Nsa2-n1)*0.5)), h, np.zeros(np.ceil((Nsa2-n1)*0.5))))
    HMF = GC.fft(hMF)

    for i in range(Srx.shape[1]):
        col = Srx[:,i]
        colMF = np.concatenate((np.zeros(np.floor((Nsa2-len(col))*0.5)), col, np.zeros(np.ceil((Nsa2-len(col))*0.5))))
        ColMF = GC.fft(colMF)
        ColMF = ColMF * HMF
        SrxMF[:,i] = GC.ifft(ColMF)

    return SrxMF


def rangeCompression(Srx, h):
    """

    Method to apply matched filter h on range dimension

    """

    # pad received and emitted signals on common time stemp
    # to obtain common fft
    n1 = len(h)
    n2 = Srx.shape[1]

    Nsa = n1 + n2 - 1
##    Nsa2 = int(np.power(2.0, np.ceil(np.log2( Nsa ))))
    Nsa2 = Nsa if Nsa % 2 == 0 else Nsa+1

    SrxMF=np.zeros((Srx.shape[0], Nsa2),dtype=complex)
    hMF = np.concatenate((np.zeros(np.floor((Nsa2-n1)*0.5)), h, np.zeros(np.ceil((Nsa2-n1)*0.5))))
    HMF = GC.fft(hMF)

    for i in range(Srx.shape[0]):
        row = Srx[i,:]
        rowMF = np.concatenate((np.zeros(np.floor((Nsa2-len(row))*0.5)), row, np.zeros(np.ceil((Nsa2-len(row))*0.5))))
        RowMF = GC.fft(rowMF)
        RowMF = RowMF * HMF
        SrxMF[i,:] = GC.ifft(RowMF)

    return SrxMF

def rangeMigration(Srx, dRI):
    """
    RCMC algorithm
    """
    Srm = np.zeros(Srx.shape, dtype=complex)
    rcm_max = dRI.max()

    # Loop on rows
    ind1=np.arange(int(Srm.shape[1] - rcm_max),dtype=int)
    for i in range(Srm.shape[0]):
        # Loop on row elements:
        ind2 = ind1+dRI[i]
        Srm[i,ind1[:]] = Srx[i,ind2]

    return Srm


def showEmittedSignal(tArray, Stx, signalName=""):
##    Stx = emittedSignal(tArray,Sat.Tp, Sat.Freq, Sat.K)
    StxFFT = GC.fft(Stx)
    fArray = Sat.Freq + Sat.K * tArray
    plt.subplot(221)
    plt.title(signalName + " signal, Re")
    plt.plot(tArray, np.real(Stx))
    plt.subplot(222)
    plt.title(signalName + " signal, Im")
    plt.plot(tArray, np.imag(Stx))
    plt.subplot(223)
    plt.title(signalName + " signal FFT, module")
    plt.plot(fArray, np.absolute(StxFFT))
    plt.subplot(224)
    plt.title(signalName + " signal FFT, phase")
    plt.plot(fArray, np.sign(np.absolute(StxFFT)) * np.angle(StxFFT))
    plt.show()


def singleTarget2(xArray, x_c, y_c, deltaX, dt):
    """
    Function to modelize SAR response from a zone at (x_c,y_x)
    """
    dist = lambda x : distSST(x, Sat.H, x_c, y_c)
    tmin = np.min(2.0*dist(xArray)) / GC.C - 1.6*Sat.Tp
    tmax = np.max(2.0*dist(xArray)) / GC.C + 1.6*Sat.Tp
    tArray=np.arange(tmin,tmax,dt)

    if (len(tArray) % 2 == 1):
        tArray = np.resize(tArray, (len(tArray)-1))

    Srx = np.zeros((len(xArray), len(tArray)),dtype=complex)
    Wa = lambda x : approxWa(x,x_c,deltaX)
##    Wa = lambda x : idealWa(x,x_c,deltaX)
    for index, t in enumerate(tArray):
        Srx[:, index] = received2DSignal(t,xArray[:], dist, Sat.Tp, Sat.Freq, Sat.K, Wa)

    return Srx, tArray

##def singleTarget(tArray, xArray, x_c, y_c, deltaX):
##
##    dist = lambda x : distSST(x, Sat.H, x_c, y_c)
##    Srx = np.zeros((len(xArray), len(tArray)),dtype=complex)
##    Wa = lambda x : approxWa(x,x_c,deltaX)
##    for index, t in enumerate(tArray):
##        Srx[:, index] = received2DSignal(t,xArray[:], dist, Sat.Tp, Sat.Freq, Sat.K, Wa)
##
##    return Srx


def showResponse(Srx):
    plt.subplot(221)
    plt.title("SAR response, module")
    plt.xlabel("Range bins")
    plt.ylabel("Azimuth bins")
    plt.imshow(np.absolute(Srx),aspect='auto')
    plt.colorbar()
    plt.subplot(222)
    plt.title("SAR response, phase")
    plt.imshow(np.angle(Srx),aspect='auto')
    plt.colorbar()
    plt.subplot(223)
    plt.title("SAR response, Re")
    plt.imshow(np.real(Srx),aspect='auto')
    plt.colorbar()
    plt.subplot(224)
    plt.title("SAR response, Im")
    plt.imshow(np.imag(Srx),aspect='auto')
    plt.colorbar()

    plt.show()



if __name__ == '__main__':

    # Range time
    B = Sat.Tp*Sat.K
    Nsa = max(2*B*Sat.Tp, 512)
    dt = Sat.Tp/Nsa

    # 1) Emitted signal
##    tArray=np.arange(-Sat.Tp*0.6,Sat.Tp*0.8,dt)
##    Stx = emittedSignal(tArray,Sat.Tp, Sat.Freq, Sat.K)
##    showEmittedSignal(tArray,Stx)


    # 2) 2D received signal from single target
    deltaT = 1.0; # acquisition time in seconds
    Xb = Sat.Vsat * deltaT * 0.5
    dx=Sat.Vsat / Sat.PRF
    # Azimuth positions
    xArray=np.arange(-Xb*2.5,Xb*2.5,dx)
    x_c = Xb*0.11
    y_c = Sat.H * np.tan(75.0 * np.pi / 180.0)
    Srx, tArrayR = singleTarget2(xArray,x_c,y_c,2.0*Xb,dt)

    showResponse(Srx)

    # 3) Range compression :
    tArray0=np.arange(-Sat.Tp*0.5,Sat.Tp*0.5,dt)
    h=chirpSignal(tArray0, Sat.Tp, -Sat.K)
    Srx_RC = rangeCompression(Srx,h)

##    Nsa = SrxMF.shape[1]
##    tmin = tArrayR[0] - len(h)*0.5*dt
##    tmax = tmin + dt*Nsa
##    tArrayMF = np.arange(tmin,tmax,dt)
##    tArrayMF = np.resize(tArrayMF,Nsa)
    showResponse(Srx_RC)


    # 4) Azimuth FFT
    SrxAF_RC = np.zeros(Srx_RC.shape, dtype=complex)
    for i in range(SrxAF_RC.shape[1]):
        # get range fixed column
        SrxAF_RC[:,i] = GC.fft(Srx_RC[:,i])
    showResponse(SrxAF_RC)

    exit()

    # 5) Range cell migration Compensation (RCMC)
    # from range time to range distance :
    dist = lambda x : distSST(x, Sat.H, x_c, y_c)
    dR = deltaR(xArray, x_c, dist)
    dRI = np.round(dR).astype(np.int)
    Srx_RC_RCMC = rangeMigration(Srx_RC, dRI)




    # 6) Azimuth compression :
    # azimuth reference signal
    Rc = dist(x_c)
    Ka = 2.0 / Sat.Lambda * 1.0 / Rc
    h = GC.rect((xArray-x_c)/(2.0*Xb)) * np.exp(1j * np.pi * Ka * (xArray - x_c)**2)
    Srx_RC_RCMC_AC = azimuthCompression(Srx_RC_RCMC, h)
    showResponse(Srx_RC_RCMC_AC)

##    showResponse(Srx_final)
##    Srx_final = np.zeros(SrxAF_RCMC.shape, dtype=complex)
##    for i in range(SrxAF_RCMC.shape[1]):
##        # get range fixed column
##        Srx_final[:,i] = GC.ifft(SrxAF_RCMC[:,i])
##
##    showResponse(Srx_final)






