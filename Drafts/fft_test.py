#-------------------------------------------------------------------------------
# Name:        FFT test

# Purpose:
#
# Author:      vfomin
#
# Created:     04/06/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


"""

Questions :

1) How FFT and analytic Fourier transform expression are related precisely ?
2)


"""




# Numpy & Matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Project
import GlobalConsts as GC

def showReSignal(tArray, s):
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


def showCSignal(tArray, s):
    B = len(tArray)*1.0/(tArray.max()-tArray.min())
    sfft = GC.fft(s)
    fArray = np.arange(-B*0.5,B*0.5,B/len(tArray))[0:len(tArray)]

    plt.subplot(221)
    plt.title("Input Signal, Re")
    plt.xlabel("Time, s")
    plt.plot(tArray, np.real(s))
    plt.subplot(222)
    plt.title("Input Signal, Im")
    plt.xlabel("Time, s")
    plt.plot(tArray, np.imag(s))

    plt.subplot(223)
    plt.title("Abs FFT of Input Signal")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray, np.absolute(sfft))
    plt.subplot(224)
    plt.title("Phase FFT of Input Signal")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray, np.angle(sfft))

    plt.show()

def showCSignal2(tArray, s):
    B = len(tArray)*1.0/(tArray.max()-tArray.min())
    sfft = GC.fft(s)
    fArray = np.arange(-B*0.5,B*0.5,B/len(tArray))[0:len(tArray)]


    plt.subplot(221)
    plt.title("Input Signal, Re")
    plt.xlabel("Time, s")
    plt.plot(tArray, np.real(s))
    plt.subplot(222)
    plt.title("Input Signal, Im")
    plt.xlabel("Time, s")
    plt.plot(tArray, np.imag(s))

    plt.subplot(223)
    plt.title("FFT of Input Signal : Re")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray, np.real(sfft))
    plt.subplot(224)
    plt.title("FFT of Input Signal : Im")
    plt.xlabel("Freq, Hz")
    plt.plot(fArray, np.imag(sfft))

    plt.show()


def test1():
    # Create complex shirp signal
    T = 10
    N = 2048
    Tp = 2
    K = 25.0

    step = T*1.0/N
    tArray = np.arange(-T*0.5,T*0.5,step)
    delay = T*0.15
    s0 = GC.rect((tArray-delay)/Tp) * np.exp(1j * np.pi * K * (tArray-delay)**2)
    ##s0 = GC.rect((tArray-delay)/Tp) * np.real(np.exp(1j * np.pi * K * (tArray-delay)**2))


    showCSignal(tArray,s0)
    ##showCSignal2(tArray,s0)

    # Compare ffts: ifft(fft(s)) - s
    ##s1 = np.fft.ifft(np.fft.fft(s0))
    ##diff1 = s1 - s0
    ##showCSignal2(tArray,diff1)

    # Compare ffts: nfft(.) = fftshift( fft( fftshift(.) ) ) and nifft(.) = fftshift( ifft( fftshift(.) ) )
    s2 = GC.ifft(GC.fft(s0))
    showCSignal(tArray,s2)
    diff2 = s2 - s0

    plt.plot(tArray,np.real(diff2))
    plt.plot(tArray,np.imag(diff2))
    plt.plot(tArray,np.abs(diff2))
    plt.show()

def plotReIm(x,y):
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(x, np.real(y), 'b', label='R()')
    ax.plot(x, np.imag(y), 'r:', label='I()')
    ax.plot(x, np.abs(y), 'k--', label='abs()')
    ax.legend()


def test2():
    dt = 0.02
    T = 30
    t = np.arange(-T*0.5,T*0.5,dt)
    k = np.fft.fftshift(np.fft.fftfreq(len(t), dt))

    a=3.0
    b=5.0
    f = np.cos(2*np.pi*a*(t-b)) * np.exp(-np.pi*(t-b)**2)
    plotReIm(t,f)


    F = dt/np.sqrt(2*np.pi) *  GC.fft(f)

    w = 2*np.pi*k
    Ft = 1.0 / np.sqrt(2.0*np.pi) * np.cosh(a*w) * np.exp(-np.pi * a**2 - w**2/(4.0*np.pi)) * np.exp(-1j*b*w)

    diff = np.sum(np.abs(F-Ft)**2)
    print diff

    plotReIm(k,F)
    plotReIm(k,Ft)
    plt.show()


def test3():
    a = 0.35
    b = 0.0

    dx = 0.015
    x = np.arange(-12, 12, dx)

    y = np.exp(-a*(x-b)**2)

    k = np.fft.fftshift(np.fft.fftfreq(len(x), dx))
    Y_k = dx/np.sqrt(2*np.pi) *  GC.fft( y )

    w = 2*np.pi*k
    Yt_k = 1.0 / np.sqrt(2.0*a) * np.exp(-w**2 / (4.0*a)) * np.exp(-1j * b * w)

    diff = np.sum(np.abs(Y_k-Yt_k)**2)
    print diff

##    plotReIm(x,y)
##    plotReIm(k,Y_k)
##    plotReIm(k,Yt_k)
##    plt.show()






##test1()
##test2()
##test3()

dt = 0.01
Tp = 5.0
t = np.arange(-15.0, 15, dt)
y = GC.rect(t/Tp)

k = np.fft.fftshift(np.fft.fftfreq(len(t), dt))
Y_k = dt/np.sqrt(2*np.pi)*GC.fft( y )

w = 2*np.pi*k
Yt_k = Tp/np.sqrt(2*np.pi) * np.sinc(0.5*Tp*w)

diff = np.sum(np.abs(Y_k-Yt_k)**2)
print diff

plotReIm(t,y)
plotReIm(k,Y_k)
plotReIm(k,Yt_k)
plt.show()
