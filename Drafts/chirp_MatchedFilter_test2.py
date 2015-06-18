#-------------------------------------------------------------------------------
# Name:        module2
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


##def showSignal(tArray, s):
##
##    B = len(tArray)*1.0/(tArray.max()-tArray.min())
##    fArray = np.arange(-B*0.5,B*0.5,B/len(tArray))
##    sfft = GC.fft(s)
##
##    plt.subplot(131)
##    plt.title("Input Signal")
##    plt.xlabel("Time, s")
##    plt.plot(tArray, s)
##    plt.subplot(132)
##    plt.title("Abs FFT of Input Signal")
##    plt.xlabel("Freq, Hz")
##    plt.plot(fArray, np.absolute(sfft))
##    plt.subplot(133)
##    plt.title("Phase FFT of Input Signal")
##    plt.xlabel("Freq, Hz")
##    plt.plot(fArray, np.angle(sfft))
##
##    plt.show()

def plotReIm(x,y):
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(x, np.real(y), 'b', label='R()')
    ax.plot(x, np.imag(y), 'r:', label='I()')
    ax.plot(x, np.abs(y), 'k--', label='abs()')
    ax.legend()


# Initial chirp centered at t=0
K=2.5
Tp=4.0
Fs = 250.0
dt=1.0/Fs
tArray0 = np.arange(-15.0, 15.0, dt)
s0 = GC.chirpSignal(tArray0,Tp,K)

fArray0 = np.fft.fftshift( np.fft.fftfreq(len(tArray0),dt) )
S0 = GC.fft(s0)

##plotReIm(tArray0, s0)
##plotReIm(fArray0, S0)
##plt.show()

# Response chirp + noise
delay = 6.0
s1 = GC.chirpSignal(tArray0-delay,Tp,K)
s1 += 0.02 * np.random.randn(len(s1))
##plotReIm(tArray0, s1)
##plt.show()

# Matched filter :
h = s0[::-1].conj()
n1 = len(h)
n2 = len(s1)
Nsa = n1 + n2 - 1
tmin=tArray0[0] - n1*0.5*dt
tArray1 = np.arange(tmin,tmin + dt*Nsa,dt)

# a) convolution
sMF1 = np.convolve(s1,h)

##plt.plot(tArray0, np.abs(s0)/np.abs(s0).max(), 'g-')
##plt.plot(tArray0, np.abs(s1)/np.abs(s1).max(), 'r-')
##plt.plot(tArray1, np.abs(sMF1)/np.abs(sMF1).max(), 'y-')
##plt.show()

# b) FFT
h1 = np.concatenate((np.zeros(np.floor((Nsa-n1)*0.5)), h, np.zeros(np.ceil((Nsa-n1)*0.5))))
s1e = np.concatenate((np.zeros(np.floor((Nsa-n2)*0.5)), s1, np.zeros(np.ceil((Nsa-n2)*0.5))))

H1 = GC.fft(h1)
S1 = GC.fft(s1e) * H1
sMF2 = GC.ifft(S1)

plt.plot(tArray0, np.abs(s0)/np.abs(s0).max(), 'g-')
plt.plot(tArray0, np.abs(s1)/np.abs(s1).max(), 'r-')
plt.plot(tArray1, np.abs(sMF2)/np.abs(sMF2).max(), 'y-')
plt.show()


