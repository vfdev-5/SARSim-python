#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      vfomin
#
# Created:     20/05/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Numpy & Matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Project
import GlobalConsts as GC

def createPulse(tau, T, N, delayTime=0):
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


def compareFFTs(s0,tArray0,s1,tArray1):

    B0 = N/T0
    B1 = N/T1
    fArray0 = np.arange(-B0*0.5,B0*0.5,B0/N)
    fArray1 = np.arange(-B1*0.5,B1*0.5,B1/N)
    sfft0 = GC.fft(s0)
    sfft1 = GC.fft(s1)


    plt.subplot(221)
    plt.title("Pulse of duration " + str(tau) + "s during time T=" + str(T0) + "s")
    plt.xlabel("Time, s")
    plt.plot(tArray0, s0)
    plt.subplot(222)
    plt.title("Pulse of duration " + str(tau) + "s during time T=" + str(T1) + "s")
    plt.xlabel("Time, s")
    plt.plot(tArray1, s1)
    plt.subplot(223)
    plt.title("FFT of Pulse of duration " + str(tau) + "s during time T=" + str(T0) + "s")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray0, np.absolute(sfft0))
    plt.subplot(224)
    plt.title("FFT of Pulse of duration " + str(tau) + "s during time T=" + str(T1) + "s")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray1, np.absolute(sfft1))
    plt.show()

def showSignal(tArray, s):

    B = len(tArray)*1.0/(tArray.max()-tArray.min())
    fArray = np.arange(-B*0.5,B*0.5,B/len(tArray))
    sfft = GC.fft(s)

    plt.subplot(131)
    plt.title("Input Signal")
    plt.xlabel("Time, s")
    plt.plot(tArray, s)
    plt.subplot(132)
    plt.title("Abs FFT of Input Signal")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray, np.absolute(sfft))
    plt.subplot(133)
    plt.title("Phase FFT of Input Signal")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray, np.angle(sfft))

    plt.show()

def plotReIm(x,y):
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(x, np.real(y), 'b', label='R()')
    ax.plot(x, np.imag(y), 'r:', label='I()')
    ax.plot(x, np.abs(y), 'k--', label='abs()')
    ax.legend()


tau = 0.5
T0 = 2.0
T1 = 3.5
N = 512

s0, tArray0 = createPulse(tau,T0,N,tau*0.51)
s1, tArray1 = createPulse(tau,T1,N,tau*0.51)

##showSignal(tArray1, s1)
##compareFFTs(s0,tArray0,s1,tArray1)

delay = tau*1.7
s11, tArray1 =  createPulse(tau,T1,N,delay)
s11 = s11 + 0.07 * np.random.randn(len(s11))

##plt.plot(tArray1, s11)
##plt.show()
##showSignal(tArray1, s11)

# Matched filter :
n1 = len(s1)
n2 = len(s11)
Nsa = n1 + n2 - 1
Nsa2 = int(np.power(2.0, np.ceil(np.log2( Nsa ) ) ) )
##Nsa2 = Nsa

dt = tArray1[1] - tArray1[0]
tmin = tArray1[0] - (Nsa2-n1)*0.5*dt
tmax = tmin + dt*Nsa2
tArray2 = np.arange(tmin,tmax,dt)

# a) FFT
##Nsa2 = Nsa
h = np.concatenate((np.zeros((Nsa2-n1)/2), s1[::-1].conj(), np.zeros((Nsa2-n1)/2)))
s11e = np.concatenate((np.zeros((Nsa2-n2)/2), s11, np.zeros((Nsa2-n2)/2)))
hfft = GC.fft(h)
s11fft = GC.fft(s11e)

fArray = np.fft.fftshift(np.fft.fftfreq(Nsa2, dt))
plotReIm(fArray, hfft)
plotReIm(fArray, s11fft)
plt.show()

Smf = np.absolute( GC.ifft(s11fft * hfft) )

##showSignal(tArray2, Smf)


# b) convolution
hfilter=s1[::-1].conj()
##hfilter = np.zeros((len(s1)))
##for i in range(len(s1)):
##    hfilter[i] = np.conjugate( s1[len(s1) - i - 1] )
Smf2 = np.convolve(s11,hfilter)
tmin=tArray1[0] - n1*0.5*dt
tArray3 = np.arange(tmin,tmin + dt*Nsa,dt)
##showSignal(tArray3,Smf2)

#c) FFT with H(f)
Nsa2 = max(512, n2*2 - 1)
fArray = np.fft.fftshift( np.fft.fftfreq(Nsa2, dt) )
HMF = np.sqrt(2.0*np.pi)/dt *  T1 * np.sinc(T1*fArray)
s11e = np.concatenate((np.zeros(np.floor((Nsa2-n2)*0.5)), s11, np.zeros(np.ceil((Nsa2-n2)*0.5))))
s11fft = GC.fft(s11e)
Smf3 = np.absolute( GC.ifft(s11fft * HMF) )

plotReIm(fArray, HMF)
plotReIm(fArray, s11fft)
plt.show()

##tmin = tArray1[0] - (Nsa2-n1)*0.5*dt
##tArray4 = np.arange(tmin,tmin + dt*Nsa2,dt)
##showSignal(np.arange(len(Smf3)),Smf3)



plt.title("Received & Matched filtered signal")
plt.xlabel("Time, s")
plt.plot(tArray1, s11)
plt.plot(tArray2, Smf)
plt.plot(tArray3, Smf2)
##
plt.show()

