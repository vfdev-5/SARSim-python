#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      vfomin
#
# Created:     27/04/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Python
import sys

# Numpy
import numpy as np
import matplotlib.pyplot as plt

# Opencv
import cv2

# Project
#import Common.ImageTools as ImageTools


# GLOBAL VARIABLES :
c = 299792458.0 # in m/s
Freq = 5.331e9 # Hz (http://fr.wikipedia.org/wiki/ENVISAT)
Lambda = c/Freq # in meters
H = 785000.0 # height of the satellite in m (http://fr.wikipedia.org/wiki/ENVISAT)
Theta = 35.0 * np.pi / 180.0 # Elevation angle in radians
R0 = H / np.cos(Theta) # m
Period = 100.6 * 60.0 # orbit period time in seconds (http://fr.wikipedia.org/wiki/ENVISAT)
Vsat = 7126.14  # in m/s (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)
K = 4.17788e+11 # Hz/s (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)
La = 10 # in m (taken from https://earth.esa.int/handbooks/asar/CNTR3-1-1.html)
Y0 = R0 * np.sin(Theta) # in meters
Vrg = 12.0 # in m/s
Vaz = 5.0 # in m/s
Arg = 2.0 # in m/s^2

FS = 18.9624680e6 # Hz (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)
PRF = 1679.9020 # Hz (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)

print "---- Satellite acquisition parameters: "
print "- Elevation angle, Theta (rad) :", Theta
print "- Range distance, R0 (m) :", R0
print "- Wavelength, Lambda (m) :", Lambda
print "- Satellite velocity, Vsat (m/s) :", Vsat
print "- Chirp slope, K :", K

print "---- Ground moving target parameters: "
print "- Ground range coordinate, Y0 (m) :", Y0
print "- Groung velocity, Vrg, Vaz (m/s) :", Vrg, ",",  Vaz
print "- Groung range acceleration, Arg (m/s^2) :", Arg


def R(s):
    """
    Slant range in function of the azimuth time s in (seconds)

    R(t) = R0 + y0*Vrg/R0 * t + (1/(2*R0)) * ( (Vaz - Vsat)^2 + Vrg^2 * (1 - y0^2/R0^2 + y0 * a_rg) ) * t^2

    y0 : Ground range coordinate of the moving target, y0 = R0 * sin(theta), theta is local elevation angle
    Vrg, Vaz : Groung velocity of the moving target
    a_rg : Ground range acceleration of the moving target
    Vsat : Satellite velocity
    """
    return R0 + Y0*Vrg/R0 * s + 0.5 / R0 * ( (Vaz - Vsat)**2 + Vrg**2 * (1.0 - (Y0/R0)**2 + Y0 *Arg)) * s**2


def sarResponse(s, t, R):
    """
     s : Azimuth time at zero Doppler
     t : Range time (fast time)

     The received signal in base-band has the following expression:
     v(s,t) = 0.5 * exp(-j * 4 * pi * R(s) / lambda ) * exp( j * pi * K * (t - 2*(R(s)-R0)/c)^2 )

     c : Light velocity
     K : Chirp slope
     R0 : Slant range at the zero Doppler
     lambda : Wavelenght of the transmitted signal

    """
    return 0.5*np.exp(-4.0*np.pi*1j*R(s)/Lambda ) * np.exp(np.pi*1j*K * (t - 2*(R(s) - R0)/c)**2)


def W(s):
    """
    Antenna pattern
    Wa(s) = sinc( (La / Lambda * arctan(s / R0) )^2
    La : Physical length of the antenna
    """
    return np.power(np.sinc( La/Lambda * np.arctan2(s,R0)), 2.0)



if __name__ == '__main__':

    Rvect = np.vectorize(R)

    # Sample number in range direction
    Nchirp = 1024
    # Sample number of the synthetic apeture length
    Nsa = 1024
    Ncenter = Nsa/2

    # range/quick time
    t = np.arange(-Nchirp*0.5,Nchirp*0.5,1)/FS
    tfull = np.arange(-Nchirp,Nchirp,1)/FS
    # azimuth/slow time
    s = (np.arange(0,Nsa,1) - Ncenter) / PRF
    sfull = (np.arange(0,2*Nsa,1) - 2*Ncenter) / PRF

    plt.plot(W(s),'-b')
    plt.show()
    print W(s).min(), W(s).max()

    #plt.figure()
    #plt.subplot(121)
    #plt.plot(np.absolute(W(s)*sarResponse(s,0,Rvect)),'-b')
    #plt.subplot(122)
    #plt.plot(np.angle(W(s)*sarResponse(s,0,Rvect)),'-b')
    #plt.show()
    sys.exit()


    data = np.zeros((len(sfull),len(tfull)),dtype=complex)
    target = np.zeros((len(s),len(t)),dtype=complex)
    for index, value in enumerate(t):
        target[:,index] = W(s[:]) * sarResponse(s[:],value,Rvect)


    plt.subplot(121)
    #plt.plot(np.absolute(target[:,0]),'-b')
    plt.imshow(np.absolute(target))
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(np.angle(target))
    plt.colorbar()
    plt.show()
    print np.absolute(target).min(), np.absolute(target).max()
    sys.exit()

    data[400:400+len(s),500:500+len(t)]=target

    plt.subplot(121)
    plt.imshow(np.angle(data))
    plt.subplot(122)
    plt.imshow(np.absolute(data))
    plt.show()
    sys.exit()

    # Range & Azimuth compression

    rc = 0.5*np.exp(1j*np.pi*K*tfull*tfull)
    rcfft = np.fft.fft(rc)

##    plt.plot(np.angle(rcfft))
##    plt.show()

    dataRC = np.zeros(data.shape, dtype=np.complex64)
    for s in range(dataRC.shape[0]):
        t0 = data[s,:]
        t1 = np.fft.fft(t0)
        t1 = t1 / rcfft
        dataRC[s,:] = np.fft.ifft(t1)



    f0 = plt.figure(0)
    f0.add_subplot(121)
    plt.imshow(np.angle(dataRC))
    f0.add_subplot(122)
    plt.imshow(np.absolute(dataRC))
    plt.show()

    #------
    #cv2.destroyAllWindows()