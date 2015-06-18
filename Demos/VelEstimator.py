#-------------------------------------------------------------------------------
# Name:     GMTI - Velocity estimation
# Purpose:
#
# Author:      vfomin
#
# Created:     15/06/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Numpy, Matplotplib
import numpy as np
import matplotlib.pyplot as plt

# Project
import Tools.GlobalConsts as GC
##import GlobalVariablesCSK_HI as Sensor
import Sensors.GlobalVariablesAerial2 as Sensor
import SARSim.SARSim_proto2 as SARSim


def azimuthFFT(data):
    procDataAF = np.zeros(data.shape, dtype=complex)
    for i in range(procDataAF.shape[1]):
        # get range fixed column
        procDataAF[:,i] = GC.fft(data[:,i])
    return procDataAF

def azimuthIFFT(dataAF):
    procData = np.zeros(dataAF.shape, dtype=complex)
    for i in range(procData.shape[1]):
        # get range fixed column
        procData[:,i] = GC.ifft(dataAF[:,i])
    return procData

def azimuthTransform(dataAF, transformType, fArray, yArray, Lambda, Vsat, vAzimuth = 0.0):
    """
    Function to perform azimuth image compression/decompression
    - dataAF is azimuth FFT of the compressed data
    - transformType is 'compress' or 'decompress'
    - fArray is a np.array of frequencies corresponding to azimuth times
    - Lambda is sensor wavelength
    - yArray is a np.array of ground range positions
    - Vsat is satellite velocity
    - vAzimuth is azimuth velocity of the transformed target
    """
    sign = +1 if transformType is 'compress' else -1 if transformType is 'decompress' else 0
    assert sign == +1 or sign == -1, logPrint("Parameter transformType should be 'compress' or 'decompress'")

    Dfreq = np.sqrt(1.0 - (0.5 * Lambda * fArray / (Vsat - vAzimuth))**2)
##    Dfreq = - 0.25 * (Lambda * fArray / (Vsat - vAzimuth))**2

    procDataAF = np.zeros(dataAF.shape, dtype=complex)
    for i in range(dataAF.shape[1]):
        HMF = np.exp( sign * 4*np.pi*1j * yArray[i] / Lambda * Dfreq)
        procDataAF[:,i] = dataAF[:,i] * HMF
    return procDataAF


if __name__ == '__main__':

##    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\SARSim_CSK_HI_v1_0_0__v2_10_0.npy"
##    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\SARSim_CSK_HI_v1_50_0.npy"
##    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\SARSim_CSK_HI_vs_10_15_30_50.npy"
    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\SARSim_CSK_HI_v1_0_0__v2_0_5.npy"

    img_set = np.load(filename)

    data = img_set[0] # compressed data
    xArray = img_set[1] # azimuth positions
    yArray = img_set[2] # ground range positions
    dx = (xArray[-1] - xArray[0])/(len(xArray)-1)
    Vsat = Sensor.Vsat
    Lambda = Sensor.Lambda

##    # TEST ON THE WHOLE IMAGE
##    dataAF = azimuthFFT(data)
##    fArray = GC.freq(dataAF.shape[0], Vsat/dx)
##    dataAF = azimuthTransform(dataAF, 'decompress', fArray, yArray, Lambda, Vsat)
##    dataAF2 = azimuthTransform(dataAF, 'compress', fArray, yArray, Lambda, Vsat, 50.0)
##    outData = azimuthIFFT(dataAF2)
##    SARSim.showResponse(outData, None, 'Compressed image - Ground range')
##    exit()


##    win = [1842, 1500, 150, data.shape[0]-1500]

##    win = [1742, 0, 250, data.shape[0]]
##    index = 134

##    data = data[win[1]:win[1]+win[3], win[0]:win[0]+win[2]]
##    yArray = yArray[win[0]:win[0]+win[2]]

    SARSim.showResponse(data, None, 'Compressed image - Ground range')


    # Decompress image at azimuth
    dataAF = azimuthFFT(data)
    fArray = GC.freq(dataAF.shape[0], Vsat/dx)
    dataAF = azimuthTransform(dataAF, 'decompress', fArray, yArray, Lambda, Vsat)

    # Compress at different velocities
    vset = np.arange(0.0, 10.0, .2)
    stack3D = np.zeros((data.shape[0], data.shape[1], len(vset)), dtype=np.float32)
    for i, v in enumerate(vset):
        print "Compress at velocity:", v
        dataAF2 = azimuthTransform(dataAF, 'compress', fArray, yArray, Lambda, Vsat, v)
        outData = azimuthIFFT(dataAF2)
        stack3D[:,:,i] = np.abs(outData[:,:])


    # User interaction :
    data0 =  stack3D[:,:,0]
    fig, ax = plt.subplots(1,2)
    fig.suptitle("Choose a target to estimate azimuth velocity")
    plt.title("Compressed data at v=0")
    ax[0].imshow(data0, picker=1)
    ax[1].imshow(stack3D[:,0,:])

    indices=[]
    def on_pick(event):
        artist = event.artist
        xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
        print 'Artist picked:', event.artist
        print 'x, y of mouse: {:.2f},{:.2f}'.format(xmouse, ymouse)
        print 'data value : ', data0[int(ymouse), int(xmouse)]

        indices.append(int(xmouse))
##        ax[1].imshow(stack3D[:,index,:], extent=[vset[0],vset[-1], 0, stack3D[:,index,:].shape[0]], interpolation='none')
##        sliceData = stack3D[:,index,:]
##        maxindex = np.unravel_index( sliceData.argmax(), sliceData.shape )
##        print "Estimated velocity : ", vset[maxindex[1]]


    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()

    for index in indices:
        plt.imshow(stack3D[:,index,:], extent=[vset[0],vset[-1], 0, stack3D[:,index,:].shape[0]], interpolation='none', aspect='auto')
        plt.show()

        sliceData = stack3D[:,index,:]
        maxindex = np.unravel_index( sliceData.argmax(), sliceData.shape )
        print "Estimated velocity : ", vset[maxindex[1]]




##    plt.subplot(241)
##    plt.imshow(stack3D[:,:,0], interpolation='none')
##    plt.subplot(242)
##    plt.imshow(stack3D[:,:,1], interpolation='none')
##    plt.subplot(243)
##    plt.imshow(stack3D[:,:,2], interpolation='none')
##    plt.subplot(244)
##    plt.imshow(stack3D[:,:,3], interpolation='none')
##    plt.subplot(245)
##    plt.imshow(stack3D[:,:,4], interpolation='none')
##    plt.subplot(246)
##    plt.imshow(stack3D[:,:,5], interpolation='none')
##    plt.subplot(247)
##    plt.imshow(stack3D[:,:,6], interpolation='none')
##    plt.show()


